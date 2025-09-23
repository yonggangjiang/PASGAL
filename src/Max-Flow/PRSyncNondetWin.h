#include <queue>
#include <limits>
#include <algorithm>

#include "graph.h"
#include "hashbag.h"
#include "parlay/delayed_sequence.h"
#include "parlay/sequence.h"

using namespace std;
using namespace parlay;

template <class Graph>
class PRSyncNondetWin {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;
  using Edge = typename Graph::Edge;
  using FlowTy = typename Edge::FlowTy;

  // cannot use const Graph& because e.w is modified
  int n;
  NodeId source, sink;
  Graph& G;
  hashbag<NodeId> bag;
  sequence<bool> in_frontier;
  sequence<NodeId> frontier;
  sequence<int> heights;
  sequence<FlowTy> excess;
  
  // Optimized fields based on maxFlow.C
  sequence<FlowTy> excess_new;  // For atomic updates
  sequence<bool> in_working_set;
  sequence<pair<NodeId, FlowTy>> pushes;  // Push operations buffer
  sequence<NodeId> working_set;
  size_t working_set_size;
  size_t work_counter;
  size_t work_threshold;
  size_t global_relabel_counter;
  bool first_iteration;
  
  // Constants for work computation (from maxFlow.C)
  static constexpr size_t ALPHA = 6;
  static constexpr size_t BETA = 12;

 public:
  PRSyncNondetWin() = delete;
  PRSyncNondetWin(Graph& _G) : G(_G), bag(G.n) {
    n = G.n;
    in_frontier = sequence<bool>(G.n);
    frontier = sequence<NodeId>::uninitialized(G.n);
    heights = sequence<int>::uninitialized(G.n);
    excess = sequence<FlowTy>::uninitialized(G.n);
    
    // Initialize optimized fields
    excess_new = sequence<FlowTy>::uninitialized(G.n);
    in_working_set = sequence<bool>(G.n);
    pushes = sequence<pair<NodeId, FlowTy>>(G.m * 2 + G.n); // Buffer for push operations
    working_set = sequence<NodeId>::uninitialized(G.n);
    working_set_size = 0;
    work_threshold = ALPHA * n + G.m / 2; // From maxFlow.C
    work_counter = 0;
    global_relabel_counter = 0;
    first_iteration = true;
  }
  void add_to_bag(NodeId u) {
    if (u != source && u != sink &&
        compare_and_swap(&in_frontier[u], false, true)) {
      bag.insert(u);
    }
  }

  // Global Relabel using reverse BFS from sink - optimized version
  void global_relabel() {
    global_relabel_counter++;
    
    // Reset all heights to n (unreachable)
    parallel_for(0, G.n, [&](size_t i) {
      if (i != sink) {
        heights[i] = n;
      }
    });
    
    // Sequential BFS for correctness
    sequence<NodeId> current_level(G.n);
    sequence<NodeId> next_level(G.n);
    
    heights[sink] = 0;
    current_level[0] = sink;
    size_t current_size = 1;
    
    int distance = 0;
    while (current_size > 0) {
      size_t next_size = 0;
      
      // Process current level
      for (size_t i = 0; i < current_size; i++) {
        NodeId u = current_level[i];
        
        // Look at incoming edges (reverse edges in residual graph)
        for (size_t j = G.offsets[u]; j < G.offsets[u + 1]; j++) {
          auto& e = G.edges[j];
          NodeId v = e.v;
          
          // Check if reverse edge has capacity and v is not yet labeled
          auto& rev_e = G.edges[e.rev];
          if (rev_e.w > 0 && heights[v] == n && v != source) {
            heights[v] = distance + 1;
            next_level[next_size++] = v;
          }
        }
      }
      
      // Swap levels
      for (size_t i = 0; i < next_size; i++) {
        current_level[i] = next_level[i];
      }
      current_size = next_size;
      distance++;
    }
    
    // Reset work counter
    work_counter = 0;
  }

  void init(NodeId s) {
    parallel_for(0, G.n, [&](size_t i) {
      heights[i] = 0;
      excess[i] = 0;
      excess_new[i] = 0;
      in_working_set[i] = false;
    });
    
    heights[s] = n;
    working_set_size = 0;
    first_iteration = true;
    
    // Saturate source-adjacent edges - optimized version
    sequence<NodeId> temp_working_set(G.offsets[s + 1] - G.offsets[s]);
    size_t temp_size = 0;
    
    for (size_t i = G.offsets[s]; i < G.offsets[s + 1]; i++) {
      auto& e = G.edges[i];
      if (e.w > 0) {
        G.edges[e.rev].w += e.w;
        excess[e.v] += e.w;
        e.w = 0;
        
        // Add to temporary working set if not source or sink
        if (e.v != source && e.v != sink && excess[e.v] > 0) {
          temp_working_set[temp_size++] = e.v;
        }
      }
    }
    
    // Remove duplicates and set up initial working set
    sort(temp_working_set.begin(), temp_working_set.begin() + temp_size);
    temp_size = unique(temp_working_set.begin(), temp_working_set.begin() + temp_size) - temp_working_set.begin();
    
    for (size_t i = 0; i < temp_size; i++) {
      working_set[working_set_size++] = temp_working_set[i];
    }
  }

  void relabel(NodeId u) {
    int new_height = std::numeric_limits<int>::max();
    for (size_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      auto& e = G.edges[i];
      if (e.w > 0) {
        new_height = min(new_height, heights[e.v] + 1);
      }
    }
    heights[u] = new_height;
  }

  // Process a single node - optimized based on maxFlow.C
  pair<size_t, size_t> process_node(NodeId u, pair<NodeId, FlowTy>* push_buffer) {
    if (u == source || u == sink || excess[u] <= 0 || heights[u] >= n) {
      return {0, 0};
    }
    
    size_t work_done = 0;
    size_t push_count = 0;
    FlowTy original_excess = excess[u];
    FlowTy local_excess = original_excess;
    
    // Add self to push buffer first (for excess delta tracking)
    push_buffer[push_count++] = {u, 0};
    
    // Main discharge loop with relabeling - inspired by maxFlow.C RELABEL_TACTIC==3
    while (local_excess > 0 && heights[u] < n) {
      int min_height = n;
      bool found_admissible = false;
      // Scan all edges for pushing and compute new height
      for (size_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
        auto& e = G.edges[i];
        NodeId v = e.v;
        
        if (local_excess <= 0) break;
        
        // Update min height for potential relabeling
        if (e.w > 0) {
          min_height = min(min_height, heights[v] + 1);
        }
        
        // Check if edge is admissible for pushing
        if (e.w > 0 && heights[u] == heights[v] + 1) {
          found_admissible = true;
          FlowTy delta = min(e.w, local_excess);
          
          // Perform push operation
          e.w -= delta;
          G.edges[e.rev].w += delta;
          local_excess -= delta;
          
          // Add to push buffer  
          push_buffer[push_count++] = {v, delta};
          
          work_done += BETA; // Push work
        }
      }
      
      work_done += G.offsets[u + 1] - G.offsets[u]; // Edge scan work
      
      // Relabel if no admissible edges found and still have excess
      if (!found_admissible && local_excess > 0 && min_height < n) {
        heights[u] = min_height;
        work_done += BETA + (G.offsets[u + 1] - G.offsets[u]); // Relabel work
      } else {
        break; // No more work to do or reached unreachable height
      }
      
      // Break if we relabeled to unreachable height
      if (heights[u] >= n) {
        break;
      }
    }
    
    // Store excess delta in first buffer position
    push_buffer[0].second = local_excess - original_excess;
    
    return {work_done, push_count};
  }

  FlowTy max_flow(NodeId s, NodeId t) {
    source = s;
    sink = t;
    if (source == sink) {
      return 0;
    }
    
    init(source);
    
    // Initial global relabel
    if (first_iteration) {
      global_relabel();
      first_iteration = false;
    }
    
    size_t iteration = 0;
    
    // Main algorithm loop - optimized based on maxFlow.C
    while (working_set_size > 0) {
      iteration++;
      
      // Check for global relabel based on work threshold
      const float GLOBAL_UPDATE_FREQ = 3.0;
      if (work_counter * GLOBAL_UPDATE_FREQ > work_threshold) {
        global_relabel();
        
        // Rebuild working set after global relabel
        working_set_size = 0;
        for (NodeId v = 0; v < G.n; v++) {
          if (v != source && v != sink && excess[v] > 0 && heights[v] < n) {
            working_set[working_set_size++] = v;
          }
        }
        
        if (working_set_size == 0) break;
      }
      
      // Adaptive thresholds based on working set size (inspired by maxFlow.C)
      bool use_parallel = working_set_size > 50;
      size_t parallel_threshold = 300;
      
      // Initialize excess_new for this iteration
      parallel_for(0, G.n, [&](size_t i) {
        excess_new[i] = 0;
        in_working_set[i] = false;
      });
      
      // Allocate push buffers with better size estimation
      sequence<size_t> buffer_sizes(working_set_size);
      parallel_for(0, working_set_size, [&](size_t i) {
        NodeId u = working_set[i];
        size_t degree = G.offsets[u + 1] - G.offsets[u];
        buffer_sizes[i] = degree + 2; // +2 for self and safety
      });
      
      // Compute prefix sums for buffer allocation
      sequence<size_t> buffer_offsets = sequence<size_t>::uninitialized(working_set_size + 1);
      buffer_offsets[0] = 0;
      for (size_t i = 0; i < working_set_size; i++) {
        buffer_offsets[i + 1] = buffer_offsets[i] + buffer_sizes[i];
      }
      
      size_t total_buffer_size = buffer_offsets[working_set_size];
      sequence<pair<NodeId, FlowTy>> all_pushes(total_buffer_size);
      
      // Process all nodes in working set 
      sequence<pair<size_t, size_t>> work_and_pushes(working_set_size);
      
      if (use_parallel && working_set_size > parallel_threshold) {
        parallel_for(0, working_set_size, [&](size_t i) {
          NodeId u = working_set[i];
          pair<NodeId, FlowTy>* push_buffer = &all_pushes[buffer_offsets[i]];
          work_and_pushes[i] = process_node(u, push_buffer);
          in_working_set[u] = true;
        });
      } else {
        for (size_t i = 0; i < working_set_size; i++) {
          NodeId u = working_set[i];
          pair<NodeId, FlowTy>* push_buffer = &all_pushes[buffer_offsets[i]];
          work_and_pushes[i] = process_node(u, push_buffer);
          in_working_set[u] = true;
        }
      }
      
      // Accumulate work counter
      for (size_t i = 0; i < working_set_size; i++) {
        work_counter += work_and_pushes[i].first;
      }
      
      // Apply pushes using race-mode approach (inspired by maxFlow.C MODE_RACE)
      for (size_t i = 0; i < working_set_size; i++) {
        NodeId u = working_set[i];
        size_t push_count = work_and_pushes[i].second;
        pair<NodeId, FlowTy>* push_buffer = &all_pushes[buffer_offsets[i]];
        
        // Apply excess change to source node
        excess[u] += push_buffer[0].second;
        
        // Apply pushes to target nodes  
        for (size_t j = 1; j < push_count; j++) {
          NodeId v = push_buffer[j].first;
          FlowTy delta = push_buffer[j].second;
          
          if (v == sink) {
            write_add(&excess[v], delta);
          } else {
            write_add(&excess_new[v], delta);
          }
        }
      }
      
      // Apply accumulated excess updates
      parallel_for(0, G.n, [&](size_t i) {
        if (excess_new[i] > 0) {
          excess[i] += excess_new[i];
        }
      });
      
      // Build new working set - collect candidate nodes efficiently
      sequence<NodeId> candidate_set(total_buffer_size);
      size_t candidate_size = 0;
      
      // Collect all nodes that received pushes
      for (size_t i = 0; i < working_set_size; i++) {
        size_t push_count = work_and_pushes[i].second;
        pair<NodeId, FlowTy>* push_buffer = &all_pushes[buffer_offsets[i]];
        
        for (size_t j = 0; j < push_count; j++) {
          NodeId v = push_buffer[j].first;
          if (v != source && v != sink && heights[v] < n) {
            candidate_set[candidate_size++] = v;
          }
        }
      }
      
      // Filter candidates efficiently (inspired by maxFlow.C filtering)
      if (candidate_size > 1000) {
        // Use parallel filtering for large sets - manual implementation since sequence::filter needs template args
        working_set_size = 0;
        parallel_for(0, candidate_size, [&](size_t i) {
          if (excess[candidate_set[i]] > 0) {
            // Use atomic increment for thread safety
            size_t pos = __sync_fetch_and_add(&working_set_size, 1);
            if (pos < G.n) {
              working_set[pos] = candidate_set[i];
            }
          }
        });
      } else {
        // Use sequential filtering for smaller sets
        sort(candidate_set.begin(), candidate_set.begin() + candidate_size);
        candidate_size = unique(candidate_set.begin(), candidate_set.begin() + candidate_size) 
                         - candidate_set.begin();
        
        working_set_size = 0;
        for (size_t i = 0; i < candidate_size; i++) {
          NodeId v = candidate_set[i];
          if (excess[v] > 0) {
            working_set[working_set_size++] = v;
          }
        }
      }
      
      // Reset working set flags
      for (size_t i = 0; i < working_set_size; i++) {
        in_working_set[working_set[i]] = false;
      }
      
      // Safety check to avoid infinite loops
      if (iteration > 10000) {
        break;
      }
    }
    
    return excess[sink];
  }
};

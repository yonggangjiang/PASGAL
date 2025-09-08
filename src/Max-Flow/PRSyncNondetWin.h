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
    FlowTy local_excess = excess[u];
    
    // Add self to push buffer first
    push_buffer[push_count++] = {u, 0};
    
    while (local_excess > 0 && heights[u] < n) {
      int min_height = n;
      bool found_admissible = false;
      
      // Scan all edges for pushing and relabeling
      for (size_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
        auto& e = G.edges[i];
        NodeId v = e.v;
        
        if (local_excess <= 0) break;
        
        // Update min height for relabeling
        if (e.w > 0) {
          min_height = min(min_height, heights[v] + 1);
        }
        
        // Check if edge is admissible for pushing
        if (e.w > 0 && heights[u] == heights[v] + 1) {
          found_admissible = true;
          FlowTy delta = min(e.w, local_excess);
          
          // Push flow
          e.w -= delta;
          G.edges[e.rev].w += delta;
          local_excess -= delta;
          
          // Add to push buffer
          push_buffer[push_count++] = {v, delta};
          
          work_done += BETA;
        }
      }
      
      work_done += G.offsets[u + 1] - G.offsets[u]; // Edge scan work
      
      if (!found_admissible && local_excess > 0 && min_height < n) {
        // Relabel
        heights[u] = min_height;
        work_done += BETA + (G.offsets[u + 1] - G.offsets[u]);
      } else {
        break;
      }
    }
    
    // Update remaining excess
    push_buffer[0].second = local_excess - excess[u]; // Delta from original
    
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
    
    // Allocate push buffers for all threads
    sequence<pair<NodeId, FlowTy>> all_pushes(G.m * 2 + G.n);
    sequence<size_t> push_offsets(max(working_set_size + 1, size_t(1)));
    
    size_t iteration = 0;
    while (working_set_size > 0) {
      iteration++;
      
      // Check if we need global relabel based on work threshold
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
        
        // Reallocate push buffers if needed
        if (push_offsets.size() < working_set_size + 1) {
          push_offsets = sequence<size_t>(working_set_size + 1);
        }
      }
      
      // Initialize excess_new for this iteration
      parallel_for(0, G.n, [&](size_t i) {
        excess_new[i] = 0;
        in_working_set[i] = false;
      });
      
      // Estimate buffer sizes needed per node (conservative estimate)
      parallel_for(0, working_set_size, [&](size_t i) {
        NodeId u = working_set[i];
        size_t degree = G.offsets[u + 1] - G.offsets[u];
        push_offsets[i] = degree + 1; // +1 for self
      });
      
      // Compute prefix sums for push buffer offsets
      size_t total_pushes = 0;
      for (size_t i = 0; i < working_set_size; i++) {
        size_t current = push_offsets[i];
        push_offsets[i] = total_pushes;
        total_pushes += current;
      }
      push_offsets[working_set_size] = total_pushes;
      
      // Ensure we have enough space
      if (all_pushes.size() < total_pushes) {
        all_pushes = sequence<pair<NodeId, FlowTy>>(total_pushes * 2); // Extra space
      }
      
      // Process all nodes in working set in parallel
      sequence<pair<size_t, size_t>> work_and_pushes(working_set_size);
      parallel_for(0, working_set_size, [&](size_t i) {
        NodeId u = working_set[i];
        pair<NodeId, FlowTy>* push_buffer = &all_pushes[push_offsets[i]];
        work_and_pushes[i] = process_node(u, push_buffer);
        in_working_set[u] = true;
      });
      
      // Sum up work done
      for (size_t i = 0; i < working_set_size; i++) {
        work_counter += work_and_pushes[i].first;
      }
      
      // Apply all pushes
      for (size_t i = 0; i < working_set_size; i++) {
        NodeId u = working_set[i];
        size_t push_count = work_and_pushes[i].second;
        pair<NodeId, FlowTy>* push_buffer = &all_pushes[push_offsets[i]];
        
        // Apply self excess update first
        excess[u] += push_buffer[0].second;
        
        // Apply pushes to other nodes
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
      
      // Apply excess updates
      parallel_for(0, G.n, [&](size_t i) {
        if (excess_new[i] > 0) {
          excess[i] += excess_new[i];
        }
      });
      
      // Build new working set - optimized approach inspired by maxFlow.C
      // Collect all vertices that received pushes
      sequence<NodeId> candidate_set(G.n);
      size_t candidate_size = 0;
      
      // Add vertices from push operations
      for (size_t i = 0; i < working_set_size; i++) {
        size_t push_count = work_and_pushes[i].second;
        pair<NodeId, FlowTy>* push_buffer = &all_pushes[push_offsets[i]];
        
        for (size_t j = 0; j < push_count; j++) {
          NodeId v = push_buffer[j].first;
          if (v != source && v != sink && heights[v] < n) {
            candidate_set[candidate_size++] = v;
          }
        }
      }
      
      // Sort and remove duplicates
      if (candidate_size > 0) {
        sort(candidate_set.begin(), candidate_set.begin() + candidate_size);
        candidate_size = unique(candidate_set.begin(), candidate_set.begin() + candidate_size) - candidate_set.begin();
      }
      
      // Filter to only include vertices with positive excess
      working_set_size = 0;
      for (size_t i = 0; i < candidate_size; i++) {
        NodeId v = candidate_set[i];
        if (excess[v] > 0) {
          working_set[working_set_size++] = v;
        }
      }
      
      // Safety check to avoid infinite loops
      if (iteration > 10000) {
        break;
      }
    }
    
    return excess[sink];
  }
};

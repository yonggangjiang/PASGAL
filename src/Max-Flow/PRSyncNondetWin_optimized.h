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
  
  // Critical optimizations based on maxFlow.C exact implementation
  sequence<FlowTy> enew;  // Atomic excess updates - CRITICAL!
  sequence<bool> in_working_set;
  sequence<pair<NodeId, FlowTy>> push_buffer;  // Global push buffer
  sequence<NodeId> working_set;
  sequence<NodeId> working_set_tmp;
  size_t working_set_size;
  size_t work_counter;
  size_t work_threshold;
  size_t global_relabel_counter;
  bool first_iteration;
  
  // Constants for work computation (from maxFlow.C)
  static constexpr size_t ALPHA = 6;
  static constexpr size_t BETA = 12;
  static constexpr float GLOBAL_UPDATE_FREQ = 3.0;

 public:
  PRSyncNondetWin() = delete;
  PRSyncNondetWin(Graph& _G) : G(_G), bag(G.n) {
    n = G.n;
    in_frontier = sequence<bool>(G.n);
    frontier = sequence<NodeId>::uninitialized(G.n);
    heights = sequence<int>::uninitialized(G.n);
    excess = sequence<FlowTy>::uninitialized(G.n);
    
    // Initialize critical fields exactly like maxFlow.C
    enew = sequence<FlowTy>::uninitialized(G.n);
    in_working_set = sequence<bool>(G.n);
    push_buffer = sequence<pair<NodeId, FlowTy>>(G.m + G.n); // Buffer for push operations
    working_set = sequence<NodeId>::uninitialized(G.n);
    working_set_tmp = sequence<NodeId>::uninitialized(G.n);
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

  // Global Relabel using reverse BFS from sink - exactly like maxFlow.C
  void global_relabel() {
    global_relabel_counter++;
    work_counter = 0;
    
    // Reset all heights to n (unreachable) except sink
    parallel_for(0, G.n, [&](size_t i) {
      if (i != sink) {
        heights[i] = n;
      }
    });
    
    // BFS from sink using frontier approach like maxFlow.C
    sequence<NodeId> current_frontier(G.n);
    current_frontier[0] = sink;
    heights[sink] = 0;
    size_t frontier_size = 1;
    
    int distance = 0;
    while (frontier_size > 0) {
      sequence<NodeId> next_frontier(G.n);
      size_t next_size = 0;
      
      // Process current frontier level
      for (size_t i = 0; i < frontier_size; i++) {
        NodeId u = current_frontier[i];
        
        // Look at all edges from u
        for (size_t j = G.offsets[u]; j < G.offsets[u + 1]; j++) {
          auto& e = G.edges[j];
          NodeId v = e.v;
          
          // Check if reverse edge has capacity (can reach u from v)
          auto& rev_e = G.edges[e.rev];
          if (rev_e.w > 0 && heights[v] == n && v != source) {
            heights[v] = distance + 1;
            next_frontier[next_size++] = v;
          }
        }
      }
      
      // Move to next level
      for (size_t i = 0; i < next_size; i++) {
        current_frontier[i] = next_frontier[i];
      }
      frontier_size = next_size;
      distance++;
    }
  }

  void init(NodeId s) {
    parallel_for(0, G.n, [&](size_t i) {
      heights[i] = 0;
      excess[i] = 0;
      enew[i] = 0;  // Critical initialization
      in_working_set[i] = false;
    });
    
    heights[s] = n;
    working_set_size = 0;
    first_iteration = true;
    
    // Saturate source-adjacent edges - following maxFlow.C exactly
    sequence<NodeId> temp_working_set(G.offsets[s + 1] - G.offsets[s]);
    size_t temp_size = 0;
    
    for (size_t i = G.offsets[s]; i < G.offsets[s + 1]; i++) {
      auto& e = G.edges[i];
      if (e.w > 0) {
        FlowTy delta = e.w;
        // Push flow from source to neighbors
        e.w = 0;
        G.edges[e.rev].w += delta;
        excess[e.v] += delta;
        
        // Add to temporary working set if eligible
        if (e.v != source && e.v != sink && excess[e.v] > 0) {
          temp_working_set[temp_size++] = e.v;
        }
      }
    }
    
    // Remove duplicates and build initial working set
    if (temp_size > 0) {
      sort(temp_working_set.begin(), temp_working_set.begin() + temp_size);
      temp_size = unique(temp_working_set.begin(), temp_working_set.begin() + temp_size) 
                  - temp_working_set.begin();
      
      for (size_t i = 0; i < temp_size; i++) {
        working_set[working_set_size++] = temp_working_set[i];
      }
    }
  }

  // Process a single node - following maxFlow.C processNode exactly for RELABEL_TACTIC != 4
  pair<size_t, size_t> process_node(NodeId u, pair<NodeId, FlowTy>* local_push_buffer) {
    if (u == source || u == sink) {
      return {0, 0};
    }
    
    if (excess[u] <= 0 || heights[u] >= n) {
      local_push_buffer[0] = {u, 0};
      return {0, 1};
    }
    
    size_t work_done = 0;
    size_t push_count = 0;
    FlowTy original_excess = excess[u];
    
    // Mark in working set and add self to buffer (like maxFlow.C)
    if (!in_working_set[u]) {
      in_working_set[u] = true;
      local_push_buffer[push_count++] = {u, 0};
    }
    
    // Main discharge loop with multiple relabels per iteration (RELABEL_TACTIC == 3)
    while (excess[u] > 0 && heights[u] < n) {
      int new_height = n;
      bool skipped = false;
      
      // Scan all edges for pushes and relabel computation
      for (size_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
        auto& e = G.edges[i];
        NodeId v = e.v;
        
        if (excess[u] <= 0) break;
        
        // Update minimum height for relabeling
        if (e.w > 0) {
          new_height = min(new_height, heights[v] + 1);
        }
        
        // Check admissibility and push if possible
        if (e.w > 0 && heights[u] == heights[v] + 1) {
          FlowTy delta = min(e.w, excess[u]);
          
          if (delta > 0) {
            // Perform push
            excess[u] -= delta;
            e.w -= delta;
            G.edges[e.rev].w += delta;
            
            // Add to push buffer for later application via enew
            local_push_buffer[push_count++] = {v, delta};
            
            work_done += BETA;
          }
        }
      }
      
      work_done += G.offsets[u + 1] - G.offsets[u]; // Edge scan work
      
      if (excess[u] <= 0 || skipped) break;
      
      // Relabel if still have excess and found admissible edges
      if (new_height > heights[u] && new_height < n) {
        heights[u] = new_height;
        work_done += BETA + (G.offsets[u + 1] - G.offsets[u]);
      } else {
        break;
      }
      
      if (heights[u] >= n) break;
    }
    
    // Add excess delta to buffer (like maxFlow.C)
    local_push_buffer[0].second = excess[u] - original_excess;
    
    return {work_done, push_count};
  }

  bool need_global_update() {
    return work_counter * GLOBAL_UPDATE_FREQ > work_threshold;
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
    
    // Main algorithm loop - following maxFlow.C run() exactly
    while (working_set_size > 0) {
      iteration++;
      
      if (iteration % 1000 == 0) {
        cout << "Iteration " << iteration << ", working set size: " << working_set_size << endl;
      }
      
      // Check for global relabel
      if (need_global_update()) {
        cout << "Global relabel at iteration " << iteration << endl;
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
      
      // Reset enew array (critical for race mode)
      parallel_for(0, G.n, [&](size_t i) {
        enew[i] = 0;
        in_working_set[i] = false;
      });
      
      // Process all nodes and collect push operations
      sequence<size_t> buffer_sizes(working_set_size);
      sequence<pair<size_t, size_t>> work_and_pushes(working_set_size);
      
      // Compute buffer sizes needed
      parallel_for(0, working_set_size, [&](size_t i) {
        NodeId u = working_set[i];
        buffer_sizes[i] = (G.offsets[u + 1] - G.offsets[u]) + 2; // +2 for safety
      });
      
      // Compute offsets for push buffers
      sequence<size_t> buffer_offsets = sequence<size_t>::uninitialized(working_set_size + 1);
      buffer_offsets[0] = 0;
      for (size_t i = 0; i < working_set_size; i++) {
        buffer_offsets[i + 1] = buffer_offsets[i] + buffer_sizes[i];
      }
      
      size_t total_buffer_size = buffer_offsets[working_set_size];
      sequence<pair<NodeId, FlowTy>> all_pushes(total_buffer_size);
      
      // Process all nodes in parallel (like maxFlow.C)
      bool use_parallel = working_set_size > 50;
      
      if (use_parallel) {
        parallel_for(0, working_set_size, [&](size_t i) {
          NodeId u = working_set[i];
          pair<NodeId, FlowTy>* local_buffer = &all_pushes[buffer_offsets[i]];
          work_and_pushes[i] = process_node(u, local_buffer);
        });
      } else {
        for (size_t i = 0; i < working_set_size; i++) {
          NodeId u = working_set[i];
          pair<NodeId, FlowTy>* local_buffer = &all_pushes[buffer_offsets[i]];
          work_and_pushes[i] = process_node(u, local_buffer);
        }
      }
      
      // Accumulate work
      for (size_t i = 0; i < working_set_size; i++) {
        work_counter += work_and_pushes[i].first;
      }
      
      // Apply pushes using RACE mode (enew mechanism) - critical optimization!
      size_t total_working_set_size = 0;
      
      for (size_t i = 0; i < working_set_size; i++) {
        NodeId u = working_set[i];
        size_t push_count = work_and_pushes[i].second;
        pair<NodeId, FlowTy>* local_buffer = &all_pushes[buffer_offsets[i]];
        
        // Apply excess delta to source node
        if (push_count > 0) {
          excess[u] += local_buffer[0].second;
        }
        
        // Apply all pushes via enew array (race mode)
        for (size_t j = 1; j < push_count; j++) {
          NodeId v = local_buffer[j].first;
          FlowTy delta = local_buffer[j].second;
          
          if (v == sink) {
            write_add(&excess[v], delta);
          } else {
            write_add(&enew[v], delta);
            // Add to next working set if first push to this node
            if (!in_working_set[v] && compare_and_swap(&in_working_set[v], false, true)) {
              size_t pos = __sync_fetch_and_add(&total_working_set_size, 1);
              if (pos < G.n) {
                working_set_tmp[pos] = v;
              }
            }
          }
        }
      }
      
      // Apply enew updates and filter working set (like maxFlow.C race mode)
      parallel_for(0, total_working_set_size, [&](size_t i) {
        if (i < total_working_set_size) {
          NodeId v = working_set_tmp[i];
          excess[v] += enew[v];
          enew[v] = 0;
          in_working_set[v] = false;
        }
      });
      
      // Filter working set for next iteration
      if (total_working_set_size > 1000) {
        working_set_size = 0;
        parallel_for(0, total_working_set_size, [&](size_t i) {
          NodeId v = working_set_tmp[i];
          if (excess[v] > 0 && heights[v] < n && v != sink) {
            size_t pos = __sync_fetch_and_add(&working_set_size, 1);
            if (pos < G.n) {
              working_set[pos] = v;
            }
          }
        });
      } else {
        working_set_size = 0;
        for (size_t i = 0; i < total_working_set_size; i++) {
          NodeId v = working_set_tmp[i];
          if (excess[v] > 0 && heights[v] < n && v != sink) {
            working_set[working_set_size++] = v;
          }
        }
      }
      
      // Safety check
      if (iteration > 100000) {
        break;
      }
    }
    
    return excess[sink];
  }
};

#include <queue>

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
  
  // New fields for PRSyncNondetWin
  sequence<FlowTy> added_excess;
  sequence<bool> is_discovered;
  sequence<vector<NodeId>> discovered_vertices;
  sequence<int> new_heights;
  size_t work_counter;
  size_t work_threshold;
  bool first_pulse;
  sequence<NodeId> working_set;
  size_t working_set_size;

 public:
  PRSyncNondetWin() = delete;
  PRSyncNondetWin(Graph& _G) : G(_G), bag(G.n) {
    n = G.n;
    in_frontier = sequence<bool>(G.n);
    frontier = sequence<NodeId>::uninitialized(G.n);
    heights = sequence<int>::uninitialized(G.n);
    excess = sequence<FlowTy>::uninitialized(G.n);
    
    // Initialize new fields
    added_excess = sequence<FlowTy>::uninitialized(G.n);
    is_discovered = sequence<bool>(G.n);
    discovered_vertices = sequence<vector<NodeId>>(G.n);
    new_heights = sequence<int>::uninitialized(G.n);
    work_threshold = 4 * n; // Typical threshold
    work_counter = 0;
    first_pulse = true;
    working_set = sequence<NodeId>::uninitialized(G.n);
    working_set_size = 0;
  }
  void add_to_bag(NodeId u) {
    if (u != source && u != sink &&
        compare_and_swap(&in_frontier[u], false, true)) {
      bag.insert(u);
    }
  }

  // Global Relabel using reverse BFS from sink
  void global_relabel() {
    // Reset all heights to n (unreachable)
    parallel_for(0, G.n, [&](size_t i) {
      if (i != sink) {
        heights[i] = n;
      }
    });
    
    // Start BFS from sink
    sequence<NodeId> current_level(G.n);
    sequence<NodeId> next_level(G.n);
    size_t current_size = 0;
    size_t next_size = 0;
    
    heights[sink] = 0;
    current_level[current_size++] = sink;
    
    int distance = 0;
    while (current_size > 0) {
      next_size = 0;
      
      // Process all nodes at current distance
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
    
    // Reset work counter and first pulse flag
    work_counter = 0;
    first_pulse = false;
  }

  void init(NodeId s) {
    parallel_for(0, G.n, [&](size_t i) {
      heights[i] = 0;
      excess[i] = 0;
      added_excess[i] = 0;
      is_discovered[i] = false;
      discovered_vertices[i].clear();
    });
    
    heights[s] = n;
    working_set_size = 0;
    
    // Saturate source-adjacent edges
    for (size_t i = G.offsets[s]; i < G.offsets[s + 1]; i++) {
      auto& e = G.edges[i];
      if (e.w > 0) {
        G.edges[e.rev].w += e.w;
        excess[e.v] += e.w;
        e.w = 0;
        
        // Add to working set if not source or sink
        if (e.v != source && e.v != sink && excess[e.v] > 0) {
          working_set[working_set_size++] = e.v;
        }
      }
    }
  }

  void relabel(NodeId u) {
    heights[u] = std::numeric_limits<int>::max();
    for (size_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      auto& e = G.edges[i];
      if (e.w > 0) {
        heights[u] = min(heights[u], heights[e.v]);
      }
    }
    heights[u]++;
  }

  FlowTy max_flow(NodeId s, NodeId t) {
    source = s;
    sink = t;
    if (source == sink) {
      return 0;
    }
    
    init(source);
    first_pulse = true;
    
    while (true) {
      // Check conditions for global relabel
      if (first_pulse || work_counter >= work_threshold || working_set_size == 0) {
        global_relabel();
        
        // Rebuild working set with all vertices that have excess > 0 and height < n
        working_set_size = 0;
        for (NodeId v = 0; v < G.n; v++) {
          if (v != source && v != sink && excess[v] > 0 && heights[v] < n) {
            working_set[working_set_size++] = v;
          }
        }
      }
      
      if (working_set_size == 0) {
        break;
      }
      
      // Initialize per-vertex data structures
      parallel_for(0, working_set_size, [&](size_t i) {
        NodeId v = working_set[i];
        discovered_vertices[v].clear();
        new_heights[v] = heights[v];
        added_excess[v] = 0;
        is_discovered[v] = false;
      });
      
      // Main discharge loop for each vertex in working set
      parallel_for(0, working_set_size, [&](size_t i) {
        NodeId v = working_set[i];
        FlowTy e = excess[v]; // local copy of excess
        
        while (e > 0) {
          int new_label = n;
          bool skipped = false;
          
          // Scan all residual edges
          for (size_t j = G.offsets[v]; j < G.offsets[v + 1]; j++) {
            auto& edge = G.edges[j];
            NodeId w = edge.v;
            
            if (e == 0) break; // vertex is already discharged completely
            
            bool admissible = (new_heights[v] == heights[w] + 1);
            
            // Check if edge is shared between two active vertices
            if (excess[w] > 0) {
              bool win = (heights[v] == heights[w] + 1) || 
                        (heights[v] < heights[w] - 1) || 
                        (heights[v] == heights[w] && v < w);
              if (admissible && !win) {
                skipped = true;
                continue; // skip to next residual edge
              }
            }
            
            // Push flow if edge is admissible and has capacity
            if (admissible && edge.w > 0) {
              FlowTy delta = min(edge.w, e);
              edge.w -= delta;
              G.edges[edge.rev].w += delta;
              e -= delta;
              
              // For sink, update excess immediately; for others use added_excess
              if (w == sink) {
                write_add(&excess[w], delta);
              } else {
                write_add(&added_excess[w], delta);
              }
              
              // Add w to discovered vertices if not sink and not yet discovered
              if (w != sink && compare_and_swap(&is_discovered[w], false, true)) {
                discovered_vertices[v].push_back(w);
              }
            }
            
            // Update new label for relabeling
            if (edge.w > 0 && heights[w] >= new_heights[v]) {
              new_label = min(new_label, heights[w] + 1);
            }
          }
          
          if (e == 0 || skipped) {
            break;
          }
          
          // Relabel
          new_heights[v] = new_label;
          if (new_heights[v] == n) {
            break;
          }
        }
        
        // Store remaining excess
        added_excess[v] = e - excess[v];
        
        // Add self to discovered if has excess
        if (e > 0 && compare_and_swap(&is_discovered[v], false, true)) {
          discovered_vertices[v].push_back(v);
        }
      });
      
      // Update heights and excess values
      parallel_for(0, working_set_size, [&](size_t i) {
        NodeId v = working_set[i];
        heights[v] = new_heights[v];
        excess[v] += added_excess[v];
        added_excess[v] = 0;
        is_discovered[v] = false;
      });
      
      // Collect all discovered vertices for next working set (line 54: Concat)
      sequence<NodeId> new_working_set(G.n);
      size_t new_working_set_size = 0;
      
      for (size_t i = 0; i < working_set_size; i++) {
        NodeId v = working_set[i];
        for (NodeId w : discovered_vertices[v]) {
          new_working_set[new_working_set_size++] = w;
        }
      }
      
      // Filter working set (line 55: workingSet = [v | v ← workingSet, d(v) < n])
      size_t filtered_size = 0;
      for (size_t i = 0; i < new_working_set_size; i++) {
        NodeId v = new_working_set[i];
        if (heights[v] < n) {
          new_working_set[filtered_size++] = v;
        }
      }
      new_working_set_size = filtered_size;
      
      // Update working set
      for (size_t i = 0; i < new_working_set_size; i++) {
        working_set[i] = new_working_set[i];
      }
      working_set_size = new_working_set_size;
      
      // Lines 57-60: Update excess for nodes in new working set
      parallel_for(0, working_set_size, [&](size_t i) {
        NodeId v = working_set[i];
        excess[v] += added_excess[v];
        added_excess[v] = 0;
        is_discovered[v] = false;
      });
      
      // Update work counter
      work_counter += working_set_size;
    }
    
    return excess[sink];
  }
};

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
  }

  void add_to_bag(NodeId u) {
    if (u != source && u != sink &&
        compare_and_swap(&in_frontier[u], false, true)) {
      bag.insert(u);
    }
  }

  // Global relabeling using BFS from sink
  void global_relabel() {
    parallel_for(0, G.n, [&](size_t i) {
      heights[i] = n;
    });
    
    queue<NodeId> q;
    heights[sink] = 0;
    q.push(sink);
    
    while (!q.empty()) {
      NodeId u = q.front();
      q.pop();
      
      for (size_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
        auto& e = G.edges[i];
        NodeId v = e.v;
        Edge& rev_e = G.edges[e.rev];
        
        // If reverse edge has capacity and v is not labeled
        if (rev_e.w > 0 && heights[v] == n) {
          heights[v] = heights[u] + 1;
          q.push(v);
        }
      }
    }
    
    work_counter = 0; // Reset work counter after global relabel
  }

  // Test and set for discovered vertices
  bool test_and_set_discovered(NodeId u) {
    return compare_and_swap(&is_discovered[u], false, true);
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
    parallel_for(G.offsets[s], G.offsets[s + 1], [&](size_t i) {
      auto& e = G.edges[i];
      if (e.w) {
        G.edges[e.rev].w += e.w;
        // Assume no duplicate edges. Otherwise it needs atomic addition
        excess[e.v] += e.w;
        e.w = 0;
        add_to_bag(e.v);
      }
    });
  }

  // Winning condition for conflict resolution
  bool wins(NodeId v, NodeId w) {
    return (heights[v] == heights[w] + 1) || 
           (heights[v] < heights[w] - 1) || 
           (heights[v] == heights[w] && v < w);
  }

  // Discharge function implementing the pseudocode logic
  void discharge_vertex(NodeId v) {
    discovered_vertices[v].clear();
    new_heights[v] = heights[v];
    FlowTy e = excess[v]; // local copy
    
    while (e > 0) {
      int new_label = n;
      bool skipped = false;
      
      for (size_t i = G.offsets[v]; i < G.offsets[v + 1]; i++) {
        auto& fwd_e = G.edges[i];
        NodeId w = fwd_e.v;
        
        if (e == 0) break; // vertex is already discharged completely
        
        bool admissible = (new_heights[v] == heights[w] + 1);
        
        if (excess[w] > 0) { // is the edge shared between two active vertices?
          bool win = wins(v, w);
          if (admissible && !win) {
            skipped = true;
            continue; // skip to next residual edge
          }
        }
        
        if (admissible && fwd_e.w > 0) { // edge is admissible
          FlowTy delta = min(fwd_e.w, e);
          fwd_e.w -= delta;
          G.edges[fwd_e.rev].w += delta;
          e -= delta;
          write_add(&added_excess[w], delta); // atomic fetch-and-add
          
          if (w != sink && test_and_set_discovered(w)) {
            discovered_vertices[v].push_back(w);
          }
        }
        
        if (fwd_e.w > 0 && heights[w] >= new_heights[v]) {
          new_label = min(new_label, heights[w] + 1);
        }
        
        if (e == 0 || skipped) {
          break;
        }
      }
      
      new_heights[v] = new_label;
      if (new_heights[v] == n) {
        break;
      }
    }
    
    added_excess[v] = e - excess[v];
    if (excess[v] > 0 && test_and_set_discovered(v)) {
      discovered_vertices[v].push_back(v);
    }
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
      // Check conditions for global relabel (lines 5-8)
      size_t frontier_size = bag.pack_into(frontier);
      if (first_pulse || work_counter > work_threshold || frontier_size == 0) {
        global_relabel();
        first_pulse = false;
        
        // Filter working set: keep only vertices with d(v) < n (line 9)
        sequence<NodeId> filtered_frontier = 
          parlay::filter(parlay::make_slice(frontier.data(), frontier.data() + frontier_size),
                        [&](NodeId v) { return heights[v] < n; });
        
        parallel_for(0, frontier_size, [&](size_t i) { 
          in_frontier[frontier[i]] = false; 
        });
        
        // Update frontier with filtered vertices
        frontier_size = filtered_frontier.size();
        parallel_for(0, frontier_size, [&](size_t i) {
          frontier[i] = filtered_frontier[i];
          in_frontier[frontier[i]] = true;
        });
      }
      
      // Check if working set is empty (lines 11-12)
      if (frontier_size == 0) {
        break;
      }
      
      // Phase 1: Discharge all vertices in parallel (lines 14-46)
      parallel_for(0, frontier_size, [&](size_t i) {
        NodeId v = frontier[i];
        discharge_vertex(v);
      });
      
      // Phase 2: Update heights and excess (lines 48-52)
      parallel_for(0, frontier_size, [&](size_t i) {
        NodeId v = frontier[i];
        heights[v] = new_heights[v];
        excess[v] += added_excess[v];
        added_excess[v] = 0;
        is_discovered[v] = false;
      });
      
      // Collect discovered vertices (line 54)
      vector<NodeId> new_working_set_vec;
      for (size_t i = 0; i < frontier_size; i++) {
        NodeId v = frontier[i];
        for (NodeId w : discovered_vertices[v]) {
          new_working_set_vec.push_back(w);
        }
      }
      
      // Filter new working set: keep only vertices with d(v) < n (line 55)
      sequence<NodeId> new_working_set = 
        parlay::filter(parlay::to_sequence(new_working_set_vec), 
                      [&](NodeId v) { return heights[v] < n; });
      
      // Reset frontier flags
      parallel_for(0, frontier_size, [&](size_t i) { 
        in_frontier[frontier[i]] = false; 
      });
      
      // Phase 3: Final update for discovered vertices (lines 57-60)
      parallel_for(0, new_working_set.size(), [&](size_t i) {
        NodeId v = new_working_set[i];
        excess[v] += added_excess[v];
        added_excess[v] = 0;
        is_discovered[v] = false;
      });
      
      // Update working set
      frontier_size = new_working_set.size();
      parallel_for(0, frontier_size, [&](size_t i) {
        frontier[i] = new_working_set[i];
        in_frontier[frontier[i]] = true;
      });
      
      work_counter += frontier_size; // Track work done
    }
    return excess[sink];
  }
};

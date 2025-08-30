#include <queue>

#include "graph.h"
#include "hashbag.h"
#include "parlay/delayed_sequence.h"
#include "parlay/sequence.h"

using namespace std;
using namespace parlay;

namespace Basic {

template <class Graph>
class PushRelabel {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;
  using Edge = typename Graph::Edge;
  using FlowTy = typename Edge::FlowTy;

  // cannot use const Graph& because e.w is modified
  size_t n;
  NodeId source, sink;
  Graph& G;
  hashbag<NodeId> bag;
  sequence<bool> in_frontier;
  sequence<NodeId> frontier;
  sequence<int> heights;
  sequence<FlowTy> excess;

 public:
  PushRelabel() = delete;
  PushRelabel(Graph& _G) : G(_G), bag(G.n) {
    n = G.n;
    in_frontier = sequence<bool>(G.n);
    frontier = sequence<NodeId>::uninitialized(G.n);
    heights = sequence<int>::uninitialized(G.n);
    excess = sequence<FlowTy>::uninitialized(G.n);
  }

  void add_to_bag(NodeId u) {
    if (u != source && u != sink &&
        compare_and_swap(&in_frontier[u], false, true)) {
      bag.insert(u);
    }
  }

  void init(NodeId s) {
    parallel_for(0, G.n, [&](size_t i) {
      heights[i] = 0;
      excess[i] = 0;
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

  void relabel(NodeId u) {
    heights[u] = std::numeric_limits<int>::max();
    for (size_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      auto& e = G.edges[i];
      if (e.w > 0) {
        heights[u] = min(heights[u], heights[e.v]);
      }
    }
    heights[u]++;
    // printf("node %u relabeled to %d\n", u, heights[u]);
  }

  void push(NodeId u, Edge& fwd_e, Edge& rev_e) {
    NodeId v = fwd_e.v;
    if (heights[u] != heights[v] + 1) {
      return;
    }
    FlowTy excess_u = excess[u];
    FlowTy send = min(fwd_e.w, excess_u);
    if (send) {
      write_add(&excess[u], -send);
      assert(excess[u] >= 0);
      write_add(&excess[v], send);
      fwd_e.w -= send;
      rev_e.w += send;
      add_to_bag(v);
    }
  }

  void discharge(NodeId u) {
    bool pushed = false;
    for (size_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      auto& e = G.edges[i];
      if (e.w > 0 && heights[u] == heights[e.v] + 1) {
        pushed = true;
        push(u, e, G.edges[e.rev]);
        if (!excess[u]) {
          break;
        }
      }
    }
    if (!pushed) {
      // printf("relabel %d\n", u);
      relabel(u);
      add_to_bag(u);
    }
    if (excess[u] > 0) {
      // printf("node %d has excess %d remaining\n", u, excess[u]);
      add_to_bag(u);
    }
  }

  FlowTy max_flow(NodeId s, NodeId t) {
    source = s;
    sink = t;
    if (source == sink) {
      return 0;
    }
    init(source);
    while (true) {
      size_t frontier_size = bag.pack_into(frontier);
      if (frontier_size == 0) {
        break;
      }
      parallel_for(0, frontier_size,
                   [&](size_t i) { in_frontier[frontier[i]] = false; });
      parallel_for(0, frontier_size, [&](size_t i) {
        NodeId u = frontier[i];
        discharge(u);
      });
    }
    return excess[sink];
  }
};

} // namespace Basic

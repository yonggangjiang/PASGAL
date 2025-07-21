#include <queue>

#include "graph.h"
#include "parlay/sequence.h"

using namespace std;
using namespace parlay;

template <class Graph>
class Dinic {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;
  using Edge = typename Graph::Edge;
  using FlowTy = typename Edge::FlowTy;

  // cannot use const Graph& because e.w is modified
  Graph& G;
  sequence<int> level;
  // 32 should change to 64 if one vertex has more than 2^32 neighbors
  sequence<uint32_t> iter;

 public:
  Dinic() = delete;
  Dinic(Graph& _G) : G(_G) {
    level = sequence<int>::uninitialized(G.n);
    iter = sequence<uint32_t>::uninitialized(G.n);
  }

  bool bfs(NodeId s, NodeId t) {
    fill(level.begin(), level.end(), -1);
    level[s] = 0;
    queue<NodeId> q;
    q.push(s);
    while (!q.empty()) {
      NodeId v = q.front();
      q.pop();
      for (EdgeId i = G.offsets[v]; i < G.offsets[v + 1]; i++) {
        auto& e = G.edges[i];
        if (e.w > 0 && level[e.v] < 0) {
          level[e.v] = level[v] + 1;
          q.push(e.v);
        }
      }
    }
    return level[t] != -1;
  }

  FlowTy dfs(NodeId u, NodeId t, FlowTy upTo) {
    if (u == t) {
      return upTo;
    }
    FlowTy flow = 0;
    for (auto& edge_id = iter[u]; edge_id < G.offsets[u + 1] - G.offsets[u];
         edge_id++) {
      auto& e = G.edges[G.offsets[u] + edge_id];
      if (e.w > 0 && level[e.v] == level[u] + 1) {
        FlowTy d = dfs(e.v, t, min(upTo - flow, e.w));
        if (d > 0) {
          e.w -= d;
          G.edges[e.rev].w += d;
          flow += d;
          if (flow == upTo) {
            break;
          }
        }
      }
    }
    if (flow == 0) {
      level[u] = -1;
    }
    return flow;
  }

  FlowTy max_flow(NodeId s, NodeId t) {
    if (s == t) {
      return 0;
    }
    FlowTy flow = 0;
    while (bfs(s, t)) {
      fill(iter.begin(), iter.end(), 0);
      FlowTy f;
      while ((f = dfs(s, t, numeric_limits<FlowTy>::max())) > 0) {
        flow += f;
      }
    }
    return flow;
  }
};

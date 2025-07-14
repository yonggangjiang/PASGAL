#pragma once
#include <vector>
#include <queue>
#include <algorithm>
#include <limits>

class MaxFlow {
public:
    using NodeId = int;
    using Capacity = long long;

    struct Edge {
        NodeId to, rev;
        Capacity cap;
    };

    MaxFlow(int n) : N(n), G(n) {}

    void add_edge(NodeId u, NodeId v, Capacity cap) {
        G[u].push_back({v, (int)G[v].size(), cap});
        G[v].push_back({u, (int)G[u].size() - 1, 0}); // reverse edge
    }

    // Ford-Fulkerson using BFS (Edmonds-Karp)
    Capacity max_flow(NodeId s, NodeId t) {
        Capacity flow = 0;
        std::vector<int> level(N);
        std::vector<Edge*> parent(N);
        while (true) {
            // BFS to find augmenting path
            std::fill(level.begin(), level.end(), -1);
            std::queue<NodeId> q;
            q.push(s);
            level[s] = 0;
            parent[s] = nullptr;
            while (!q.empty()) {
                NodeId u = q.front(); q.pop();
                for (auto &e : G[u]) {
                    if (level[e.to] < 0 && e.cap > 0) {
                        level[e.to] = level[u] + 1;
                        parent[e.to] = &e;
                        q.push(e.to);
                        if (e.to == t) break;
                    }
                }
            }
            if (level[t] < 0) break; // no more augmenting path
            // Find bottleneck
            Capacity bottleneck = std::numeric_limits<Capacity>::max();
            for (NodeId v = t; v != s;) {
                Edge* e = parent[v];
                NodeId u = -1;
                for (NodeId i = 0; i < N; ++i) {
                    for (auto &ed : G[i]) {
                        if (&ed == e) { u = i; break; }
                    }
                    if (u != -1) break;
                }
                bottleneck = std::min(bottleneck, e->cap);
                v = u;
            }
            // Augment flow
            for (NodeId v = t; v != s;) {
                Edge* e = parent[v];
                NodeId u = -1;
                for (NodeId i = 0; i < N; ++i) {
                    for (auto &ed : G[i]) {
                        if (&ed == e) { u = i; break; }
                    }
                    if (u != -1) break;
                }
                e->cap -= bottleneck;
                G[e->to][e->rev].cap += bottleneck;
                v = u;
            }
            flow += bottleneck;
        }
        return flow;
    }

private:
    int N;
    std::vector<std::vector<Edge>> G;
};

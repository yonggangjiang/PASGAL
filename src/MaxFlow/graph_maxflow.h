#pragma once
#include <vector>
#include <cstdint>
#include <string>
#include <fstream>
#include <iostream>
#include <cassert>

class MaxFlowGraph {
public:
    using NodeId = uint32_t;
    using EdgeId = uint64_t;
    struct Edge {
        NodeId to;
        int64_t cap;
        EdgeId rev;
    };
    std::vector<std::vector<Edge>> adj;
    size_t n, m;

    MaxFlowGraph() : n(0), m(0) {}

    void read_binary_format(const char* filename) {
        std::ifstream in(filename, std::ios::binary);
        if (!in) {
            std::cerr << "Cannot open file: " << filename << std::endl;
            std::abort();
        }
        uint64_t n64, m64, sizes;
        in.read(reinterpret_cast<char*>(&n64), 8);
        in.read(reinterpret_cast<char*>(&m64), 8);
        in.read(reinterpret_cast<char*>(&sizes), 8);
        n = n64;
        m = m64;
        assert(sizes == (n + 1) * 8 + m * 4 + 3 * 8);
        std::vector<uint64_t> offsets(n + 1);
        std::vector<uint32_t> edges(m);
        in.read(reinterpret_cast<char*>(offsets.data()), (n + 1) * 8);
        in.read(reinterpret_cast<char*>(edges.data()), m * 4);
        adj.assign(n, {});
        for (NodeId u = 0; u < n; ++u) {
            for (uint64_t i = offsets[u]; i < offsets[u + 1]; ++i) {
                NodeId v = edges[i];
                add_edge(u, v, 1); // default capacity 1, can be changed
            }
        }
    }

    void add_edge(NodeId u, NodeId v, int64_t cap) {
        adj[u].push_back({v, cap, (EdgeId)adj[v].size()});
        adj[v].push_back({u, 0, (EdgeId)adj[u].size() - 1});
    }
};

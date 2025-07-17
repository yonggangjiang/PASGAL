#include <iostream>
#include <chrono>
#include <fstream>
#include <queue>
#include <vector>
#include <algorithm>
#include "graph.h"

using namespace std;

// Optimized push-relabel implementation based on the reference algorithm
uint64_t push_relabel_maxflow(Graph<uint32_t, uint64_t, uint64_t>& G, uint32_t source, uint32_t sink) {
    const uint32_t n = G.n;
    const uint32_t m = 2 * G.m;
    
    // Arc structure for residual graph
    struct Arc {
        uint32_t to;
        uint32_t rev;
        uint64_t cap;
    };
    
    vector<Arc> arcs(m);
    vector<uint32_t> first(n + 1);
    
    // Build residual graph with correct arc indexing
    vector<uint32_t> deg(n, 0);
    for (uint32_t u = 0; u < n; u++) {
        for (uint64_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
            deg[u]++;
            deg[G.edges[i].v]++;
        }
    }
    
    first[0] = 0;
    for (uint32_t u = 0; u < n; u++) {
        first[u + 1] = first[u] + deg[u];
    }
    
    vector<uint32_t> pos = first;
    for (uint32_t u = 0; u < n; u++) {
        for (uint64_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
            uint32_t v = G.edges[i].v;
            uint64_t cap = G.edges[i].w;
            
            uint32_t fwd_idx = pos[u]++;
            uint32_t rev_idx = pos[v]++;
            
            arcs[fwd_idx] = {v, rev_idx, cap};
            arcs[rev_idx] = {u, fwd_idx, 0};
        }
    }
    
    // Node data structures
    vector<uint64_t> excess(n, 0);
    vector<uint32_t> height(n, 0);
    vector<uint32_t> current(n);
    vector<uint32_t> count(2 * n + 1, 0);
    
    for (uint32_t u = 0; u < n; u++) {
        current[u] = first[u];
    }
    
    // Global relabeling using reverse BFS
    auto global_relabel = [&]() {
        fill(height.begin(), height.end(), n);
        fill(count.begin(), count.end(), 0);
        
        height[sink] = 0;
        count[0] = 1;
        count[n] = n - 1;
        
        queue<uint32_t> q;
        q.push(sink);
        
        while (!q.empty()) {
            uint32_t u = q.front();
            q.pop();
            
            for (uint32_t i = first[u]; i < first[u + 1]; i++) {
                uint32_t v = arcs[i].to;
                if (height[v] == n && arcs[arcs[i].rev].cap > 0) {
                    count[n]--;
                    height[v] = height[u] + 1;
                    count[height[v]]++;
                    q.push(v);
                }
            }
        }
        
        // Reset current pointers
        for (uint32_t u = 0; u < n; u++) {
            current[u] = first[u];
        }
    };
    
    // Initial global relabel
    global_relabel();
    
    if (height[source] == n) {
        count[n]--;
    }
    height[source] = n;
    count[n]++;
    
    // Working set for active nodes
    vector<uint32_t> work_set, next_work_set;
    work_set.reserve(n);
    next_work_set.reserve(n);
    
    // Saturate all arcs from source
    for (uint32_t i = first[source]; i < first[source + 1]; i++) {
        Arc& arc = arcs[i];
        if (arc.cap > 0) {
            uint32_t v = arc.to;
            uint64_t flow = arc.cap;
            
            excess[v] += flow;
            arc.cap = 0;
            arcs[arc.rev].cap += flow;
            
            if (v != sink) {
                work_set.push_back(v);
            }
        }
    }
    
    // Remove duplicates and sort work set
    sort(work_set.begin(), work_set.end());
    work_set.erase(unique(work_set.begin(), work_set.end()), work_set.end());
    
    uint64_t work_since_update = 0;
    const uint64_t global_update_freq = 6 * n + m / 2;
    
    // Push operation
    auto push = [&](uint32_t u, uint32_t arc_idx) -> bool {
        Arc& arc = arcs[arc_idx];
        uint32_t v = arc.to;
        uint64_t flow = min(excess[u], arc.cap);
        
        if (flow == 0) return false;
        
        excess[u] -= flow;
        excess[v] += flow;
        arc.cap -= flow;
        arcs[arc.rev].cap += flow;
        
        if (v != sink && v != source && excess[v] == flow) {
            next_work_set.push_back(v);
        }
        
        return arc.cap == 0;
    };
    
    // Relabel operation
    auto relabel = [&](uint32_t u) {
        count[height[u]]--;
        
        // Gap heuristic
        if (count[height[u]] == 0 && height[u] < n) {
            for (uint32_t v = 0; v < n; v++) {
                if (v != source && height[v] > height[u] && height[v] < n) {
                    count[height[v]]--;
                    height[v] = n;
                    count[n]++;
                }
            }
        }
        
        uint32_t new_height = 2 * n;
        for (uint32_t i = first[u]; i < first[u + 1]; i++) {
            if (arcs[i].cap > 0) {
                new_height = min(new_height, height[arcs[i].to] + 1);
            }
        }
        
        height[u] = new_height;
        count[height[u]]++;
        current[u] = first[u];
        
        work_since_update += 12 + (first[u + 1] - first[u]);
    };
    
    // Discharge operation
    auto discharge = [&](uint32_t u) {
        while (excess[u] > 0) {
            if (current[u] >= first[u + 1]) {
                relabel(u);
            } else {
                Arc& arc = arcs[current[u]];
                if (arc.cap > 0 && height[u] == height[arc.to] + 1) {
                    push(u, current[u]);
                } else {
                    current[u]++;
                }
            }
        }
        work_since_update += 6 + (first[u + 1] - first[u]);
    };
    
    // Main algorithm loop
    while (!work_set.empty()) {
        next_work_set.clear();
        
        // Process all nodes in current work set
        for (uint32_t u : work_set) {
            if (excess[u] > 0) {
                discharge(u);
            }
        }
        
        // Sort and remove duplicates from next work set
        if (!next_work_set.empty()) {
            sort(next_work_set.begin(), next_work_set.end());
            next_work_set.erase(unique(next_work_set.begin(), next_work_set.end()), 
                               next_work_set.end());
        }
        
        work_set = std::move(next_work_set);
        
        // Periodic global relabeling
        if (work_since_update > global_update_freq) {
            work_since_update = 0;
            global_relabel();
            
            // Rebuild work set with only active nodes
            work_set.clear();
            for (uint32_t u = 0; u < n; u++) {
                if (u != source && u != sink && excess[u] > 0 && height[u] < n) {
                    work_set.push_back(u);
                }
            }
        }
    }
    
    return excess[sink];
}

int main(int argc, char* argv[]) {

    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }

    const char* input_file = argv[1];
    Graph<uint32_t, uint64_t, uint64_t> G;
    G.read_graph(input_file);

    // Use source and sink directly from the graph structure
    uint32_t source = G.source;
    uint32_t sink = G.sink;

    auto start = std::chrono::high_resolution_clock::now();
    uint64_t flow = push_relabel_maxflow(G, source, sink);
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();

    cout << "Max flow value: " << flow << endl;
    cout << "Elapsed time (s): " << elapsed << endl;

    return 0;
}

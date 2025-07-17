#include <iostream>
#include <chrono>
#include <fstream>
#include <queue>
#include <vector>
#include <algorithm>
#include <atomic>
#include <mutex>
#include "graph.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "utils.h"

using namespace std;
using namespace parlay;

// Number of threads (can be set via command line)
int num_threads = 1;

void set_num_threads(int n) {
    num_threads = n;
    // Note: parlay threading is controlled via environment variables or different API
    // For now, we'll rely on the parallel_for implementation to use available cores
}

// Corrected parallel push-relabel implementation
uint64_t push_relabel_maxflow(Graph<uint32_t, uint64_t, uint64_t>& G, uint32_t source, uint32_t sink) {
    const uint32_t n = G.n;
    const uint32_t m = 2 * G.m;
    
    // Debug output
    cout << "Graph: n=" << n << ", m=" << G.m << ", source=" << source << ", sink=" << sink << endl;
    
    // Arc structure for residual graph 
    struct Arc {
        uint32_t to;
        uint32_t rev;
        uint64_t cap;
    };
    
    // Node structure for parallel processing
    struct Node {
        atomic<uint64_t> excess;
        atomic<uint32_t> height;
        atomic<bool> in_work_set;
        uint32_t current_arc;
        
        Node() : excess(0), height(0), in_work_set(false), current_arc(0) {}
    };
    
    vector<Arc> arcs(m);
    vector<uint32_t> first(n + 1);
    vector<Node> nodes(n);
    
    // Build residual graph correctly
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
    
    // Per-arc mutexes for thread safety
    vector<mutex> arc_mutexes(m);
    
    // Node data structures
    vector<uint32_t> count(2 * n + 1, 0);
    atomic<uint64_t> work_since_update(0);
    
    // Sequential global relabeling using reverse BFS (more reliable)
    auto global_relabel_old = [&]() {
        // Reset all heights to n (unreachable)
        for (uint32_t u = 0; u < n; u++) {
            nodes[u].height.store(n);
        }
        
        fill(count.begin(), count.end(), 0);
        nodes[sink].height.store(0);
        count[0] = 1;
        count[n] = n - 1;
        
        queue<uint32_t> q;
        q.push(sink);
        
        while (!q.empty()) {
            uint32_t u = q.front();
            q.pop();
            
            for (uint32_t i = first[u]; i < first[u + 1]; i++) {
                uint32_t v = arcs[i].to;
                uint32_t rev_arc = arcs[i].rev;
                
                if (arcs[rev_arc].cap > 0 && nodes[v].height.load() == n) {
                    count[n]--;
                    uint32_t new_height = nodes[u].height.load() + 1;
                    nodes[v].height.store(new_height);
                    count[new_height]++;
                    q.push(v);
                }
            }
        }
        
        // Reset current arc pointers
        for (uint32_t u = 0; u < n; u++) {
            nodes[u].current_arc = first[u];
        }
    };
    
    // Initial global relabel
    global_relabel_old();
    nodes[source].height.store(n);
    if (count[n] > 0) count[n]--;
    count[n]++;
    
    // Working set for active nodes
    vector<uint32_t> work_set;
    work_set.reserve(n);
    
    // Saturate all arcs from source
    for (uint32_t i = first[source]; i < first[source + 1]; i++) {
        Arc& arc = arcs[i];
        if (arc.cap > 0) {
            uint32_t v = arc.to;
            uint64_t cap = arc.cap;
            
            nodes[v].excess.fetch_add(cap);
            arc.cap = 0;
            arcs[arc.rev].cap += cap;
            
            if (v != sink && !nodes[v].in_work_set.exchange(true)) {
                work_set.push_back(v);
            }
        }
    }
    
    const uint64_t global_update_freq = 6 * n + m / 2;
    
    // Global relabeling using BFS from sink (re-enabled for efficiency)
    auto global_relabel = [&]() {
        // Reset all heights to n (unreachable)
        for (uint32_t u = 0; u < n; u++) {
            nodes[u].height.store(n);
        }
        
        fill(count.begin(), count.end(), 0);
        nodes[sink].height.store(0);
        count[0] = 1;
        count[n] = n - 1;
        
        queue<uint32_t> q;
        q.push(sink);
        
        while (!q.empty()) {
            uint32_t u = q.front();
            q.pop();
            
            for (uint32_t i = first[u]; i < first[u + 1]; i++) {
                uint32_t v = arcs[i].to;
                uint32_t rev_arc = arcs[i].rev;
                
                if (arcs[rev_arc].cap > 0 && nodes[v].height.load() == n) {
                    count[n]--;
                    uint32_t new_height = nodes[u].height.load() + 1;
                    nodes[v].height.store(new_height);
                    count[new_height]++;
                    q.push(v);
                }
            }
        }
        
        // Reset current arc pointers
        for (uint32_t u = 0; u < n; u++) {
            nodes[u].current_arc = first[u];
        }
        
        cout << "Global relabel completed" << endl;
    };
    
    // Thread-safe push operation
    auto safe_push = [&](uint32_t u, uint32_t arc_idx) -> bool {
        Arc& arc = arcs[arc_idx];
        uint32_t v = arc.to;
        
        lock_guard<mutex> lock(arc_mutexes[arc_idx]);
        
        uint64_t u_excess = nodes[u].excess.load();
        uint64_t flow = min(u_excess, arc.cap);
        
        if (flow == 0) return false;
        
        nodes[u].excess.fetch_sub(flow);
        nodes[v].excess.fetch_add(flow);
        arc.cap -= flow;
        arcs[arc.rev].cap += flow;
        
        return flow > 0;
    };
    
    // Relabel operation
    auto relabel = [&](uint32_t u) -> bool {
        uint32_t old_height = nodes[u].height.load();
        uint32_t new_height = 2 * n;
        
        for (uint32_t i = first[u]; i < first[u + 1]; i++) {
            if (arcs[i].cap > 0) {
                new_height = min(new_height, nodes[arcs[i].to].height.load() + 1);
            }
        }
        
        if (new_height > old_height) {
            nodes[u].height.store(new_height);
            nodes[u].current_arc = first[u];
            work_since_update.fetch_add(12 + (first[u + 1] - first[u]));
            return true;
        }
        return false;
    };
    
    // Process a single node (discharge operation)
    auto process_node = [&](uint32_t u) -> vector<uint32_t> {
        vector<uint32_t> new_active;
        
        while (nodes[u].excess.load() > 0) {
            uint32_t &cur = nodes[u].current_arc;
            bool found_push = false;
            
            // First scan for admissible edges from current position
            for (uint32_t i = cur; i < first[u + 1]; i++) {
                Arc& arc = arcs[i];
                uint32_t v = arc.to;
                
                if (arc.cap > 0 && 
                    nodes[u].height.load() == nodes[v].height.load() + 1) {
                    
                    if (safe_push(u, i)) {
                        if (v != sink && v != source && 
                            !nodes[v].in_work_set.exchange(true)) {
                            new_active.push_back(v);
                        }
                        found_push = true;
                        cur = i; // Advance current arc
                        break; // Continue the while loop to try more pushes
                    }
                }
            }
            
            if (!found_push) {
                // No admissible arc found from current position
                // Reset to beginning and try relabel
                cur = first[u];
                
                // Try scanning from beginning once more before relabel
                for (uint32_t i = first[u]; i < first[u + 1]; i++) {
                    Arc& arc = arcs[i];
                    uint32_t v = arc.to;
                    
                    if (arc.cap > 0 && 
                        nodes[u].height.load() == nodes[v].height.load() + 1) {
                        
                        if (safe_push(u, i)) {
                            if (v != sink && v != source && 
                                !nodes[v].in_work_set.exchange(true)) {
                                new_active.push_back(v);
                            }
                            found_push = true;
                            cur = i;
                            break;
                        }
                    }
                }
                
                if (!found_push) {
                    // Still no push possible, need to relabel
                    if (!relabel(u)) {
                        break; // Can't relabel further, done with this node
                    }
                    // After relabel, reset current arc and continue
                    cur = first[u];
                }
            }
        }
        
        work_since_update.fetch_add(6 + (first[u + 1] - first[u]));
        return new_active;
    };
    
    // Main algorithm loop
    uint32_t iterations = 0;
    while (!work_set.empty()) {
        iterations++;
        if (iterations % 100 == 0) {
            cout << "Iteration " << iterations << ", active nodes: " << work_set.size() << endl;
        }
        
        vector<vector<uint32_t>> new_active_lists(work_set.size());
        
        // Process all active nodes 
        if (num_threads > 1 && work_set.size() > 10) {
            parallel_for(0, work_set.size(), [&](uint32_t i) {
                uint32_t u = work_set[i];
                nodes[u].in_work_set.store(false);
                
                if (nodes[u].excess.load() > 0) {
                    new_active_lists[i] = process_node(u);
                }
            });
        } else {
            // Sequential processing for small work sets
            for (uint32_t i = 0; i < work_set.size(); i++) {
                uint32_t u = work_set[i];
                nodes[u].in_work_set.store(false);
                
                if (nodes[u].excess.load() > 0) {
                    new_active_lists[i] = process_node(u);
                }
            }
        }
        
        // Collect new active nodes
        vector<uint32_t> new_work_set;
        for (const auto& list : new_active_lists) {
            for (uint32_t v : list) {
                new_work_set.push_back(v);
            }
        }
        
        // Remove duplicates
        if (!new_work_set.empty()) {
            sort(new_work_set.begin(), new_work_set.end());
            new_work_set.erase(unique(new_work_set.begin(), new_work_set.end()), 
                              new_work_set.end());
        }
        
        work_set = std::move(new_work_set);
        
        // Periodic global relabeling with proper frequency
        if (work_since_update.load() > global_update_freq) {
            work_since_update.store(0);
            global_relabel();
            
            // Rebuild work set with only active nodes
            work_set.clear();
            for (uint32_t u = 0; u < n; u++) {
                nodes[u].in_work_set.store(false);
                if (u != source && u != sink && 
                    nodes[u].excess.load() > 0 && 
                    nodes[u].height.load() < n) {
                    
                    nodes[u].in_work_set.store(true);
                    work_set.push_back(u);
                }
            }
            
            // if (debug) {
            //     cout << "Global relabel: work set size = " << work_set.size() << endl;
            // }
        }
    }
    
    cout << "Completed in " << iterations << " iterations" << endl;
    return nodes[sink].excess.load();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <input_file> [num_threads]" << endl;
        cout << "  num_threads: number of parallel threads to use (default: 1)" << endl;
        return 1;
    }

    const char* input_file = argv[1];
    
    // Set number of threads if specified
    if (argc >= 3) {
        num_threads = atoi(argv[2]);
        if (num_threads < 1) num_threads = 1;
        set_num_threads(num_threads);
        cout << "Using " << num_threads << " threads" << endl;
    } else {
        cout << "Using 1 thread (sequential)" << endl;
    }
    
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

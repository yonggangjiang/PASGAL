#include <iostream>
#include <queue>
#include <vector>
#include <algorithm>
#include <climits>
#include "src/graph.h"

using namespace std;

// BFS to compute shortest path distance
int bfs_distance(const Graph<>& g, int source, int target) {
    vector<int> dist(g.n, -1);
    queue<int> q;
    
    dist[source] = 0;
    q.push(source);
    
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        
        if (u == target) {
            return dist[target];
        }
        
        for (size_t i = g.offsets[u]; i < g.offsets[u + 1]; i++) {
            int v = g.edges[i].v;
            if (dist[v] == -1) {
                dist[v] = dist[u] + 1;
                q.push(v);
            }
        }
    }
    
    return -1; // No path found
}

// Compute undirected diameter using BFS from a single sampled node
int compute_undirected_diameter(const Graph<>& g) {
    if (g.n == 0) return 0;
    
    // For fast approximation, just sample one node (middle node)
    int start_node = g.n / 2;
    
    // Create undirected adjacency list
    vector<vector<int>> adj(g.n);
    for (size_t u = 0; u < g.n; u++) {
        for (size_t i = g.offsets[u]; i < g.offsets[u + 1]; i++) {
            int v = g.edges[i].v;
            if (v != u) { // Avoid self-loops
                adj[u].push_back(v);
                adj[v].push_back(u); // Add reverse edge for undirected
            }
        }
    }
    
    // Remove duplicates
    for (size_t u = 0; u < g.n; u++) {
        sort(adj[u].begin(), adj[u].end());
        adj[u].erase(unique(adj[u].begin(), adj[u].end()), adj[u].end());
    }
    
    // BFS from single starting node to find maximum distance
    vector<int> dist(g.n, -1);
    queue<int> q;
    
    dist[start_node] = 0;
    q.push(start_node);
    
    int max_distance = 0;
    
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        
        for (int v : adj[u]) {
            if (dist[v] == -1) {
                dist[v] = dist[u] + 1;
                q.push(v);
                max_distance = max(max_distance, dist[v]);
            }
        }
    }
    
    return max_distance;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <graph_file>" << endl;
        return 1;
    }
    
    string filename = argv[1];
    Graph<> g;
    
    try {
        g.read_graph(filename.c_str());
        
        cout << "Nodes: " << g.n << endl;
        cout << "Edges: " << g.m << endl;
        cout << "Density: " << (double)g.m / g.n << endl;
        
        // Compute source to sink distance (assuming source=0, sink=1)
        int source_sink_dist = bfs_distance(g, g.source, g.sink);
        cout << "Source-to-sink distance: " << source_sink_dist << endl;
        
        // Compute undirected diameter (sample-based for large graphs)
        int diameter = compute_undirected_diameter(g);
        cout << "Undirected diameter: " << diameter << endl;
        
    } catch (const exception& e) {
        cerr << "Error reading graph: " << e.what() << endl;
        return 1;
    } catch (...) {
        cerr << "Unknown error reading graph" << endl;
        return 1;
    }
    
    return 0;
}

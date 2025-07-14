#include "graph_maxflow.h"
#include <iostream>
#include <random>
#include <chrono>
#include <fstream>
#include <cstdlib>

// Edmonds-Karp/Push-Relabel or any maxflow implementation can be used here
#include "maxflow.h"

long long run_maxflow(MaxFlowGraph &G, int s, int t, double &cpp_time, long long &flow_out) {
    MaxFlow mf(G.n);
    // Copy edges from G to mf
    for (size_t u = 0; u < G.n; ++u) {
        for (const auto &e : G.adj[u]) {
            if (e.cap > 0) // Only forward edges
                mf.add_edge(u, e.to, e.cap);
        }
    }
    auto start = std::chrono::high_resolution_clock::now();
    long long flow = mf.max_flow(s, t);
    auto end = std::chrono::high_resolution_clock::now();
    cpp_time = std::chrono::duration<double>(end - start).count();
    flow_out = flow;
    return flow;
}

void write_edge_list(MaxFlowGraph &G, const std::string &filename) {
    std::ofstream ofs(filename);
    ofs << G.n << " " << G.m << "\n";
    for (size_t u = 0; u < G.n; ++u) {
        for (const auto &e : G.adj[u]) {
            if (e.cap > 0) // Only forward edges
                ofs << u << " " << e.to << " " << e.cap << "\n";
        }
    }
    ofs.close();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file> [source] [sink]\n";
        return 1;
    }
    MaxFlowGraph G;
    G.read_binary_format(argv[1]);
    std::cout << "Graph size: |V|=" << G.n << ", |E|=" << G.m << std::endl;
    long long flow = 0;
    int s, t;
    if (argc < 4) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, G.n-1);
        s = dis(gen);
        t = dis(gen);
        while (t == s) t = dis(gen);
        std::cout << "Random source: " << s << ", sink: " << t << std::endl;
    } else {
        s = std::stoi(argv[2]);
        t = std::stoi(argv[3]);
    }
    double cpp_time = 0.0;
    run_maxflow(G, s, t, cpp_time, flow);
    std::cout << "Max flow from " << s << " to " << t << ": " << flow << std::endl;
    std::cout << "C++ running time: " << cpp_time << " seconds" << std::endl;

    // Write edge list for Python verification
    std::string edgelist_file = "maxflow_edgelist.txt";
    write_edge_list(G, edgelist_file);
    // Call Python script for verification
    std::string cmd = "python3 /home/yonggang/PASGAL/src/MaxFlow/verify_maxflow.py " + edgelist_file + " " + std::to_string(s) + " " + std::to_string(t) + " " + std::to_string(flow);
    std::cout << "Running: " << cmd << std::endl;
    int ret = std::system(cmd.c_str());
    if (ret != 0) {
        std::cerr << "Python verification failed!\n";
    }
    return 0;
}

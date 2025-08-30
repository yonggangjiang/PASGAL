#include "push-relabel.h"
#include "PRSyncNondetWin.h"

#include <queue>
#include <type_traits>
#include <string>

#include "graph.h"

using namespace std;
using namespace parlay;
using namespace Basic;

typedef uint32_t NodeId;
typedef uint64_t EdgeId;
#ifdef FLOAT
typedef float FlowTy;
#else
typedef int32_t FlowTy;
#endif

constexpr int NUM_SRC = 5;
constexpr int NUM_ROUND = 5;
constexpr int LOG2_WEIGHT = 6;
constexpr int WEIGHT_RANGE = 1 << LOG2_WEIGHT;

template <class _FlowTy = int32_t>
class FlowEdge : public BaseEdge<NodeId> {
 public:
  using FlowTy = _FlowTy;
  FlowTy w;  // capacity
  EdgeId rev;
  FlowEdge() = default;
  FlowEdge(NodeId _v, FlowTy _w) : BaseEdge<NodeId>(_v), w(_w) {}
  bool operator<(const FlowEdge &rhs) const {
    return this->v < rhs.v || (this->v == rhs.v && w < rhs.w);
  }
  bool operator==(const FlowEdge &rhs) const {
    return this->v == rhs.v && w == rhs.w;
  }
};

template <class Graph>
Graph generate_reverse_edges(const Graph &G) {
  using Edge = typename Graph::Edge;
  parlay::sequence<std::tuple<NodeId, NodeId, EdgeId>> edgelist(G.m * 2);
  parlay::parallel_for(0, G.n, [&](NodeId u) {
    parlay::parallel_for(G.offsets[u], G.offsets[u + 1], [&](EdgeId i) {
      // forward edge
      std::get<0>(edgelist[i * 2]) = u;
      std::get<1>(edgelist[i * 2]) = G.edges[i].v;
      std::get<2>(edgelist[i * 2]) = i * 2;
      // reverse edge
      std::get<0>(edgelist[i * 2 + 1]) = G.edges[i].v;
      std::get<1>(edgelist[i * 2 + 1]) = u;
      std::get<2>(edgelist[i * 2 + 1]) = i * 2 + 1;
    });
  });
  parlay::sort_inplace(parlay::make_slice(edgelist), [&](auto &a, auto &b) {
    return std::get<0>(a) < std::get<0>(b) ||
           (std::get<0>(a) == std::get<0>(b) &&
            std::get<1>(a) < std::get<1>(b));
  });

  Graph G_rev;
  G_rev.n = G.n;
  G_rev.m = G.m * 2;
  G_rev.offsets = parlay::sequence<EdgeId>(G.n + 1, G.m * 2);
  G_rev.edges = parlay::sequence<Edge>(G.m * 2);
  parlay::parallel_for(0, G.m * 2, [&](EdgeId i) {
    EdgeId id = std::get<2>(edgelist[i]);
    G_rev.edges[i].v = std::get<1>(edgelist[i]);
    G_rev.edges[i].w = (id & 1) ? 0 : G.edges[id / 2].w;
    if (i == 0 || std::get<0>(edgelist[i]) != std::get<0>(edgelist[i - 1])) {
      G_rev.offsets[std::get<0>(edgelist[i])] = i;
    }
  });
  parlay::scan_inclusive_inplace(
      parlay::make_slice(G_rev.offsets.rbegin(), G_rev.offsets.rend()),
      parlay::minm<EdgeId>());
  parlay::sequence<tuple<uint32_t, uint32_t>> rev_id(G.m);
  parlay::parallel_for(0, G_rev.n, [&](NodeId u) {
    parlay::parallel_for(G_rev.offsets[u], G_rev.offsets[u + 1], [&](EdgeId i) {
      EdgeId id = std::get<2>(edgelist[i]);
      if (id & 1) {
        std::get<1>(rev_id[id / 2]) = i;
      } else {
        std::get<0>(rev_id[id / 2]) = i;
      }
    });
  });

  parlay::parallel_for(0, G.m * 2, [&](EdgeId i) {
    EdgeId id = std::get<2>(edgelist[i]);
    if (id & 1) {
      G_rev.edges[i].rev = std::get<0>(rev_id[id / 2]);
    } else {
      G_rev.edges[i].rev = std::get<1>(rev_id[id / 2]);
    }
  });

  // validate
  parlay::parallel_for(0, G_rev.n, [&](NodeId u) {
    parlay::parallel_for(G_rev.offsets[u], G_rev.offsets[u + 1], [&](EdgeId i) {
      FlowTy w = G_rev.edges[i].w;
      EdgeId rev = G_rev.edges[i].rev;
      assert(G_rev.edges[rev].v == u);
      if (w) {
        assert(G_rev.edges[rev].w == 0);
      }
    });
  });
  return G_rev;
}
template <class Algo, class Graph, class NodeId = typename Graph::NodeId>
void run(Algo &algo, [[maybe_unused]] Graph &G, NodeId s, NodeId t) {
  cout << "source " << s << ", target " << t << endl;
  double total_time = 0;
  FlowTy max_flow = 0;
  parlay::sequence<FlowTy> flows = parlay::tabulate(G.m, [&](EdgeId i) {
    return G.edges[i].w;
  });
  for (int i = 0; i <= NUM_ROUND; i++) {
    parlay::parallel_for(0, G.m, [&](EdgeId j) {
      G.edges[j].w = flows[j];
    });
    internal::timer tm;
    max_flow = algo.max_flow(s, t);
    tm.stop();
    if (i == 0) {
      cout << "Warmup Round: " << tm.total_time() << endl;
    } else {
      cout << "Round " << i << ": " << tm.total_time() << endl;
      total_time += tm.total_time();
    }
    cout << "Max flow: " << max_flow << endl;
  }
  double average_time = total_time / NUM_ROUND;
  cout << "Average time: " << average_time << endl;

  ofstream ofs("push-relabel.tsv", ios_base::app);
  ofs << s << '\t' << t << '\t' << max_flow << '\t' << average_time << '\n';
  ofs.close();
}

template <class Algo, class Graph>
void run(Algo &algo, Graph &G) {
  using NodeId = typename Graph::NodeId;
  for (int v = 0; v < NUM_SRC; v++) {
    NodeId s = hash32(v) % G.n;
    NodeId t = hash32(v + NUM_SRC) % G.n;
    run(algo, G, s, t);
  }
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    cout << "Usage: " << argv[0] << " [-i input_file] [-s] [-a algorithm]\n"
         << "Options:\n"
         << "\t-i,\tinput file path\n"
         << "\t-s,\tsymmetrized input graph\n"
         << "\t-a,\talgorithm (default: basic, options: basic, nondet)\n";
    return 0;
  }
  char c;
  char const *input_path = nullptr;
  bool symmetrized = false;
  uint32_t source = UINT_MAX;
  uint32_t target = UINT_MAX;
  string algorithm = "basic";
  while ((c = getopt(argc, argv, "i:sr:t:a:")) != -1) {
    switch (c) {
      case 'i':
        input_path = optarg;
        break;
      case 's':
        symmetrized = true;
        break;
      case 'r':
        source = atol(optarg);
        break;
      case 't':
        target = atol(optarg);
        break;
      case 'a':
        algorithm = optarg;
        break;
      default:
        std::cerr << "Error: Unknown option " << optopt << std::endl;
        abort();
    }
  }

  cout << "Reading graph..." << endl;
  Graph<NodeId, EdgeId, FlowEdge<FlowTy>> G;
  G.read_graph(input_path);
  G.symmetrized = symmetrized;
  if (!G.weighted) {
    cout << "Generating edge weights..." << endl;
    G.generate_random_weight(1, WEIGHT_RANGE);
  }
  // for (NodeId i = 0; i < G.n; i++) {
  //   for (EdgeId j = G.offsets[i]; j < G.offsets[i + 1]; j++) {
  //     printf("(%u,%u,%d)\n", i, G.edges[j].v, G.edges[j].w);
  //   }
  // }
  // printf("\n");

  G = generate_reverse_edges(G);

  cout << "Running on " << input_path << ": |V|=" << G.n << ", |E|=" << G.m
       << ", num_src=" << NUM_SRC << ", num_round=" << NUM_ROUND << endl;

  cout << "Using algorithm: " << algorithm << endl;
  
  using GraphType = Graph<NodeId, EdgeId, FlowEdge<FlowTy>>;
  
  if (algorithm == "nondet") {
    auto solver = PRSyncNondetWin<GraphType>(G);
    if (source == UINT_MAX || target == UINT_MAX) {
      run(solver, G);
    } else {
      run(solver, G, source, target);
    }
  } else if (algorithm == "basic") {
    auto solver = PushRelabel<GraphType>(G);
    if (source == UINT_MAX || target == UINT_MAX) {
      run(solver, G);
    } else {
      run(solver, G, source, target);
    }
  } else {
    std::cerr << "Error: Unknown algorithm '" << algorithm << "'. Available: basic, nondet" << std::endl;
    return 1;
  }
  return 0;
}

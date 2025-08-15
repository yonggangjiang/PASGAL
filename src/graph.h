#ifndef GRAPH_H
#define GRAPH_H

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <cassert>
#include <fstream>
#include <type_traits>
#include <vector>

#include "parlay/io.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/utilities.h"
#include "utils.h"

class Empty {};

template <class NodeId, class EdgeTy>
class WEdge {
 public:
  NodeId v;
  [[no_unique_address]] EdgeTy w;
  WEdge() {}
  WEdge(NodeId _v) : v(_v) {}
  WEdge(NodeId _v, EdgeTy _w) : v(_v), w(_w) {}

  bool operator<(const WEdge &rhs) const {
    if constexpr (std::is_same_v<EdgeTy, Empty>) {
      return v < rhs.v;
    } else {
      if (v != rhs.v) {
        return v < rhs.v;
      }
      return w < rhs.w;
    }
  }

  bool operator==(const WEdge &rhs) const {
    if constexpr (std::is_same_v<EdgeTy, Empty>) {
      return v == rhs.v;
    } else {
      return v == rhs.v && w == rhs.w;
    }
  }

  bool operator!=(const WEdge &rhs) const {
    return !(*this == rhs);
  }
};

template <class _NodeId = uint32_t, class _EdgeId = uint64_t,
          class _EdgeTy = Empty>
class Graph {
 public:
  using NodeId = _NodeId;
  using EdgeId = _EdgeId;
  using EdgeTy = _EdgeTy;
  using Edge = WEdge<NodeId, EdgeTy>;

  size_t n;
  size_t m;
  bool symmetrized;
  bool weighted;
  parlay::sequence<EdgeId> offsets;
  parlay::sequence<Edge> edges;
  parlay::sequence<EdgeId> in_offsets;
  parlay::sequence<Edge> in_edges;
  // Add source and sink for flow problems
  _NodeId source;
  _NodeId sink;

  Graph() {
    n = m = 0;
    symmetrized = weighted = false;
    source = sink = 0;
  }

  auto in_neighors(NodeId u) const {
    if (symmetrized) {
      return edges.cut(offsets[u], offsets[u + 1]);
    } else {
      return in_edges.cut(in_offsets[u], in_offsets[u + 1]);
    }
  }

  void make_inverse() {
    parlay::sequence<std::pair<NodeId, Edge>> edgelist(m);
    parlay::parallel_for(0, n, [&](NodeId u) {
      parlay::parallel_for(offsets[u], offsets[u + 1], [&](EdgeId i) {
        edgelist[i] = std::make_pair(edges[i].v, Edge(u, edges[i].w));
      });
    });
    parlay::sort_inplace(
        parlay::make_slice(edgelist),
        [](const std::pair<NodeId, Edge> &a, const std::pair<NodeId, Edge> &b) {
          if (a.first != b.first) {
            return a.first < b.first;
          }
          return a.second.v < b.second.v;
        });
    in_offsets = parlay::sequence<EdgeId>(n + 1, m);
    in_edges = parlay::sequence<Edge>(m);
    parlay::parallel_for(0, m, [&](size_t i) {
      in_edges[i] = edgelist[i].second;
      if (i == 0 || edgelist[i].first != edgelist[i - 1].first) {
        in_offsets[edgelist[i].first] = i;
      }
    });
    parlay::scan_inclusive_inplace(
        parlay::make_slice(in_offsets.rbegin(), in_offsets.rend()),
        parlay::minm<EdgeId>());
  }

  void read_pbbs_format(char const *filename) {
    auto chars = parlay::chars_from_file(std::string(filename));
    auto tokens_seq = tokens(chars);
    auto header = tokens_seq[0];
    n = chars_to_ulong_long(tokens_seq[1]);
    m = chars_to_ulong_long(tokens_seq[2]);
    bool weighted_input;
    if (header == parlay::to_chars("WeightedAdjacencyGraph")) {
      weighted_input = true;
      assert(tokens_seq.size() == n + m + m + 3);
    } else if (header == parlay::to_chars("AdjacencyGraph")) {
      weighted_input = false;
      assert(tokens_seq.size() == n + m + 3);
    } else {
      std::cerr << "Unrecognized header" << std::endl;
      abort();
    }
    offsets = parlay::sequence<EdgeId>(n + 1);
    edges = parlay::sequence<Edge>(m);
    parlay::parallel_for(0, n, [&](size_t i) {
      offsets[i] = parlay::internal::chars_to_int_t<NodeId>(
          make_slice(tokens_seq[i + 3]));
    });
    offsets[n] = m;
    parlay::parallel_for(0, m, [&](size_t i) {
      edges[i].v = parlay::internal::chars_to_int_t<NodeId>(
          make_slice(tokens_seq[i + n + 3]));
    });
    if (weighted_input) {
      if constexpr (std::is_same_v<EdgeTy, Empty>) {
        weighted = false;
        std::cout << "Warning: skipping edge weights in file" << std::endl;
      } else {
        weighted = true;
        parlay::parallel_for(0, m, [&](size_t i) {
          if constexpr (std::is_integral_v<EdgeTy>) {
            edges[i].w = parlay::internal::chars_to_int_t<NodeId>(
                make_slice(tokens_seq[i + n + m + 3]));
          } else if constexpr (std::is_floating_point_v<EdgeTy>) {
            edges[i].w = parlay::chars_to_double(tokens_seq[i + n + m + 3]);
          } else {
            std::cerr << "Error: EdgeTy is not arithmetic" << std::endl;
            abort();
          }
        });
      }
    } else {
      weighted = false;
      // Add random weights if EdgeTy is not Empty
      if constexpr (!std::is_same_v<EdgeTy, Empty>) {
        generate_random_weight(1, n);
        std::cout << "Added random weights in range [1, " << n << "]" << std::endl;
      }
    }
    // Set default source and sink for non-flow graphs
    source = 0;
    sink = (n > 1) ? 1 : 0;
  }

  void read_binary_format(char const *filename) {
    std::cout << "Opening binary file..." << std::flush;
    // Uses mmap to accelerate reading
    struct stat sb;
    int fd = open(filename, O_RDONLY);
    if (fd == -1) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      abort();
    }
    if (fstat(fd, &sb) == -1) {
      std::cerr << "Error: Unable to acquire file stat" << std::endl;
      abort();
    }
    std::cout << " mapping..." << std::flush;
    char *data =
        static_cast<char *>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
    size_t len = sb.st_size;
    n = reinterpret_cast<uint64_t *>(data)[0];
    m = reinterpret_cast<uint64_t *>(data)[1];
    size_t sizes = reinterpret_cast<uint64_t *>(data)[2];
    assert(sizes == (n + 1) * 8 + m * 4 + 3 * 8);
    std::cout << " allocating..." << std::flush;
    offsets = parlay::sequence<EdgeId>::uninitialized(n + 1);
    edges = parlay::sequence<Edge>::uninitialized(m);
    std::cout << " reading offsets..." << std::flush;
    parlay::parallel_for(0, n + 1, [&](size_t i) {
      offsets[i] = reinterpret_cast<uint64_t *>(data + 3 * 8)[i];
    });
    std::cout << " reading edges..." << std::flush;
    parlay::parallel_for(0, m, [&](size_t i) {
      edges[i].v = reinterpret_cast<uint32_t *>(data + 3 * 8 + (n + 1) * 8)[i];
    });
    std::cout << " unmapping..." << std::flush;
    if (data) {
      const void *b = data;
      munmap(const_cast<void *>(b), len);
    }
    std::cout << " done" << std::endl;
    
    // Binary format doesn't include weights, so add random weights if EdgeTy is not Empty
    weighted = false;
    if constexpr (!std::is_same_v<EdgeTy, Empty>) {
      generate_random_weight(1, n);
      std::cout << "Added random weights in range [1, " << n << "]" << std::endl;
    }
    
    // Set default source and sink for non-flow graphs
    source = 0;
    sink = (n > 1) ? 1 : 0;
  }

  void read_hyperlink2012(const char *filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      abort();
    }
    size_t sizes;
    ifs.read(reinterpret_cast<char *>(&n), sizeof(size_t));
    ifs.read(reinterpret_cast<char *>(&m), sizeof(size_t));
    ifs.read(reinterpret_cast<char *>(&sizes), sizeof(size_t));
    assert(sizes == (n + 1) * 8 + m * 4 + 3 * 8);
    offsets = parlay::sequence<EdgeId>::uninitialized(n + 1);
    edges = parlay::sequence<Edge>::uninitialized(m);
    ifs.read(reinterpret_cast<char *>(offsets.begin()), (n + 1) * 8);
    ifs.read(reinterpret_cast<char *>(edges.begin()), m * 4);

    ifs.read(reinterpret_cast<char *>(&n), sizeof(size_t));
    ifs.read(reinterpret_cast<char *>(&m), sizeof(size_t));
    ifs.read(reinterpret_cast<char *>(&sizes), sizeof(size_t));
    in_offsets = parlay::sequence<EdgeId>::uninitialized(n + 1);
    in_edges = parlay::sequence<Edge>::uninitialized(m);
    ifs.read(reinterpret_cast<char *>(in_offsets.begin()), (n + 1) * 8);
    ifs.read(reinterpret_cast<char *>(in_edges.begin()), m * 4);
    if (ifs.peek() != EOF) {
      std::cerr << "Error: Bad input graph" << std::endl;
      abort();
    }
    ifs.close();
    
    // Binary format doesn't include weights, so add random weights if EdgeTy is not Empty
    weighted = false;
    if constexpr (!std::is_same_v<EdgeTy, Empty>) {
      generate_random_weight(1, n);
      std::cout << "Added random weights in range [1, " << n << "]" << std::endl;
    }
    
    // Set default source and sink for non-flow graphs
    source = 0;
    sink = (n > 1) ? 1 : 0;
  }

  void read_flowflow_format(const char *filename) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs.is_open()) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      abort();
    }
    
    // Read and verify header
    char header[9];
    ifs.read(header, 8);
    header[8] = '\0';
    if (std::string(header) != "FLOWFLOW") {
      std::cerr << "Error: Invalid FLOWFLOW file header" << std::endl;
      abort();
    }
    
    // Read basic graph info
    uint64_t n_val, m_val, source, sink;
    ifs.read(reinterpret_cast<char*>(&n_val), sizeof(uint64_t));
    ifs.read(reinterpret_cast<char*>(&m_val), sizeof(uint64_t));
    ifs.read(reinterpret_cast<char*>(&source), sizeof(uint64_t));
    ifs.read(reinterpret_cast<char*>(&sink), sizeof(uint64_t));
    
    n = n_val;
    m = m_val;
    this->source = static_cast<NodeId>(source);
    this->sink = static_cast<NodeId>(sink);
    
    // Read vertex offsets
    parlay::sequence<uint64_t> vertex_offsets(n);
    for (size_t i = 0; i < n; ++i) {
      uint64_t offset;
      ifs.read(reinterpret_cast<char*>(&offset), sizeof(uint64_t));
      vertex_offsets[i] = offset;
    }
    
    // Read edge data (neighbor, capacity pairs)
    parlay::sequence<std::pair<uint64_t, uint64_t>> edge_data(m);
    for (size_t i = 0; i < m; ++i) {
      uint64_t neighbor, capacity;
      ifs.read(reinterpret_cast<char*>(&neighbor), sizeof(uint64_t));
      ifs.read(reinterpret_cast<char*>(&capacity), sizeof(uint64_t));
      edge_data[i] = {neighbor, capacity};
    }
    
    ifs.close();
    
    // Convert to our graph format
    offsets = parlay::sequence<EdgeId>(n + 1);
    edges = parlay::sequence<Edge>(m);
    
    // Set up offsets array
    offsets[0] = 0;
    for (size_t i = 0; i < n; ++i) {
      if (i + 1 < n) {
        offsets[i + 1] = vertex_offsets[i + 1];
      } else {
        offsets[i + 1] = m;
      }
    }
    
    // Set up edges array
    parlay::parallel_for(0, m, [&](size_t i) {
      if constexpr (std::is_same_v<EdgeTy, Empty>) {
        edges[i] = Edge(static_cast<NodeId>(edge_data[i].first));
      } else {
        edges[i] = Edge(static_cast<NodeId>(edge_data[i].first), 
                       static_cast<EdgeTy>(edge_data[i].second));
      }
    });
    
    weighted = !std::is_same_v<EdgeTy, Empty>;
  }

  void read_graph(const char *filename) {
    std::string str_filename(filename);
    if (str_filename.find("hyperlink2012.bin") != std::string::npos) {
      // hack for hyperlink2012
      read_hyperlink2012(filename);
      return;
    }
    
    // Check for FLOWFLOW format by reading header
    std::ifstream test_file(filename, std::ios::binary);
    if (test_file.is_open()) {
      char header[9];
      test_file.read(header, 8);
      header[8] = '\0';
      test_file.close();
      if (std::string(header) == "FLOWFLOW") {
        read_flowflow_format(filename);
        return;
      }
    }
    
    size_t idx = str_filename.find_last_of('.');
    if (idx == std::string::npos) {
      std::cerr << "Error: No graph extension provided" << std::endl;
      abort();
    }
    std::string subfix = str_filename.substr(idx + 1);
    if (subfix == "adj") {
      read_pbbs_format(filename);
    } else if (subfix == "bin") {
      read_binary_format(filename);
    } else {
      std::cerr << "Error: Invalid graph extension" << std::endl;
      abort();
    }
  }

  void write_pbbs_format(char const *filename) {
    parlay::chars chars;
    if constexpr (std::is_same_v<EdgeTy, Empty>) {
      chars.insert(chars.end(), parlay::to_chars("AdjacencyGraph\n"));
    } else {
      chars.insert(chars.end(), parlay::to_chars("WeightedAdjacencyGraph\n"));
    }
    chars.insert(chars.end(), parlay::to_chars(std::to_string(n) + "\n"));
    chars.insert(chars.end(), parlay::to_chars(std::to_string(m) + "\n"));
    chars.insert(chars.end(),
                 parlay::flatten(parlay::tabulate(n * 2, [&](size_t i) {
                   if (i % 2 == 0) {
                     return parlay::to_chars(offsets[i / 2]);
                   } else {
                     return parlay::to_chars('\n');
                   }
                 })));
    chars.insert(chars.end(),
                 parlay::flatten(parlay::tabulate(m * 2, [&](size_t i) {
                   if (i % 2 == 0) {
                     return parlay::to_chars(edges[i / 2].v);
                   } else {
                     return parlay::to_chars('\n');
                   }
                 })));
    if constexpr (!std::is_same_v<EdgeTy, Empty>) {
      chars.insert(chars.end(),
                   parlay::flatten(parlay::tabulate(m * 2, [&](size_t i) {
                     if (i % 2 == 0) {
                       return parlay::to_chars(std::to_string(edges[i / 2].w));
                     } else {
                       return parlay::to_chars('\n');
                     }
                   })));
    }
    chars_to_file(chars, std::string(filename));
  }

  void write_flowflow_format(const char *filename, NodeId source = 0, NodeId sink = 1) {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs.is_open()) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      abort();
    }
    
    // Write header
    ofs.write("FLOWFLOW", 8);
    
    // Write basic graph info
    uint64_t n_val = n, m_val = m;
    uint64_t source_val = source, sink_val = sink;
    ofs.write(reinterpret_cast<const char*>(&n_val), sizeof(uint64_t));
    ofs.write(reinterpret_cast<const char*>(&m_val), sizeof(uint64_t));
    ofs.write(reinterpret_cast<const char*>(&source_val), sizeof(uint64_t));
    ofs.write(reinterpret_cast<const char*>(&sink_val), sizeof(uint64_t));
    
    // Write vertex offsets
    for (size_t i = 0; i < n; ++i) {
      uint64_t offset = offsets[i];
      ofs.write(reinterpret_cast<const char*>(&offset), sizeof(uint64_t));
    }
    
    // Write edge data (neighbor, capacity pairs)
    for (size_t i = 0; i < m; ++i) {
      uint64_t neighbor = static_cast<uint64_t>(edges[i].v);
      uint64_t capacity;
      if constexpr (std::is_same_v<EdgeTy, Empty>) {
        capacity = 1; // Default capacity for unweighted graphs
      } else {
        capacity = static_cast<uint64_t>(edges[i].w);
      }
      ofs.write(reinterpret_cast<const char*>(&neighbor), sizeof(uint64_t));
      ofs.write(reinterpret_cast<const char*>(&capacity), sizeof(uint64_t));
    }
    
    ofs.close();
  }

  void write_binary_format(char const *filename) {
    assert(sizeof(EdgeId) == sizeof(uint64_t));
    assert(sizeof(NodeId) == sizeof(uint32_t));
    size_t sizes = (n + 1) * 8 + m * 4 + 3 * 8;
    auto tmp_edges =
        parlay::tabulate<NodeId>(m, [&](size_t i) { return edges[i].v; });

    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      abort();
    }
    ofs.write(reinterpret_cast<char *>(&n), sizeof(size_t));
    ofs.write(reinterpret_cast<char *>(&m), sizeof(size_t));
    ofs.write(reinterpret_cast<char *>(&sizes), sizeof(size_t));
    ofs.write(reinterpret_cast<char *>(offsets.begin()),
              sizeof(EdgeId) * (n + 1));
    ofs.write(reinterpret_cast<char *>(tmp_edges.begin()), sizeof(NodeId) * m);
    ofs.close();
  }

  // Generates integral edge weights in range [l, r) based on edge endpoints
  // This ensures the same graph structure always produces the same weights
  void generate_random_weight(uint32_t l, uint32_t r) {
    if constexpr (std::is_same_v<EdgeTy, Empty>) {
      std::cerr << "Error: Graph instance does not have a edge weight field"
                << std::endl;
      abort();
    }
    if (weighted) {
      std::cerr << "Warning: Overwriting original edge weights\n" << std::endl;
    } else {
      weighted = true;
    }
    
    std::cout << "Generating deterministic weights for " << m << " edges..." << std::flush;
    uint32_t range = r - l + 1;
    
    // Generate weights based on edge endpoints for deterministic results
    parlay::parallel_for(0, n, [&](NodeId u) {
      parlay::parallel_for(offsets[u], offsets[u + 1], [&](EdgeId i) {
        NodeId v = edges[i].v;
        // Use hash of both endpoints to generate deterministic weight
        uint32_t hash_val = (u * 2654435761U) ^ (v * 2654435761U) ^ ((u + v) * 2654435761U);
        edges[i].w = (hash_val % range) + l;
      });
    });
    
    std::cout << " done" << std::endl;
  }

  void count_self_loop_and_parallel_edges() {
    size_t self_loops = 0, parallel_edges = 0;
    parlay::parallel_for(0, n, [&](size_t i) {
      size_t pre = n + 1;
      for (size_t j = offsets[i]; j < offsets[i + 1]; j++) {
        NodeId v = edges[j].v;
        if (i == v) {
          write_add(&self_loops, 1);
        }
        if (v == pre) {
          write_add(&parallel_edges, 1);
        }
        pre = v;
      }
    });
    printf("num of self-loops: %zu\n", self_loops);
    printf("num of parallel edges: %zu\n", parallel_edges);
  }

  // sanity check of input graph
  void validate() {
    auto check_offsets = [&](size_t i) {
      if (i == 0) {
        assert(offsets[i] == 0);
      } else if (i == n) {
        assert(offsets[i] == m);
      } else {
        assert(offsets[i] >= offsets[i - 1]);
      }
    };

    auto check_edges = [&](size_t i) {
      NodeId u = edges[i].v;
      assert(u >= 0 && u < n);
    };

    auto check_sorted = [&](NodeId u) {
      bool sorted = true;
      for (size_t i = offsets[u] + 1; i < offsets[u + 1]; ++i) {
        if (edges[i].v < edges[i - 1].v) {
          sorted = false;
          break;
        }
      }
      if (!sorted) {
        printf("Warning: Edges are not sorted\n");
      }
    };

    parlay::parallel_for(0, n + 1, check_offsets);
    parlay::parallel_for(0, m, check_edges);
    parlay::parallel_for(0, n, check_sorted);

    count_self_loop_and_parallel_edges();

    if (symmetrized) {
      auto check_symmetry = [&](NodeId u) {
        for (size_t i = offsets[u]; i < offsets[u + 1]; ++i) {
          NodeId v = edges[i].v;
          auto e = WEdge(u, edges[i].w);
          if (*std::lower_bound(edges.begin() + offsets[v],
                                edges.begin() + offsets[v + 1], e) != e) {
            symmetrized = false;
            break;
          }
        }
      };

      parlay::parallel_for(0, n, check_symmetry);

      if (!symmetrized) {
        std::cerr << "Error: edges are not symmetrized" << std::endl;
        abort();
      }
    }

    printf("Passed!\n");
  }

  void generate_random_graph(size_t _n = 100000000, size_t _m = 8000000000) {
    n = _n;
    m = _m;
    auto edgelist =
        parlay::sequence<std::pair<NodeId, NodeId>>::uninitialized(m);
    parlay::parallel_for(0, m / 2, [&](size_t i) {
      NodeId u = parlay::hash32(i * 2) % n;
      NodeId v = parlay::hash32(i * 2 + 1) % n;
      if (symmetrized) {
        edgelist[i * 2] = std::make_pair(u, v);
        edgelist[i * 2 + 1] = std::make_pair(v, u);
      } else {
        edgelist[i * 2] = std::make_pair(u, v);
        edgelist[i * 2 + 1] =
            std::make_pair(parlay::hash32(u) % n, parlay::hash32(v) % n);
      }
    });
    parlay::sort_inplace(parlay::make_slice(edgelist),
                         [](const std::pair<NodeId, NodeId> &a,
                            const std::pair<NodeId, NodeId> &b) {
                           if (a.first != b.first) {
                             return a.first < b.first;
                           }
                           return a.second < b.second;
                         });
    offsets = parlay::sequence<EdgeId>(n + 1, m);
    edges = parlay::sequence<Edge>(m);
    parlay::parallel_for(0, m, [&](size_t i) {
      edges[i].v = edgelist[i].second;
      if (i == 0 || edgelist[i].first != edgelist[i - 1].first) {
        offsets[edgelist[i].first] = i;
      }
    });
    parlay::scan_inclusive_inplace(
        parlay::make_slice(offsets.rbegin(), offsets.rend()),
        parlay::minm<EdgeId>());
  }
};

template <class NodeId = uint32_t>
class Forest {
 public:
  size_t num_trees;
  Graph<NodeId, NodeId> G;
  parlay::sequence<NodeId> vertex;
  parlay::sequence<NodeId> offsets;
};

template <class NodeId, class EdgeId, class EdgeTy>
Graph<NodeId, EdgeId, EdgeTy> edgelist2graph(
    parlay::sequence<std::pair<NodeId, WEdge<NodeId, EdgeTy>>> &edgelist,
    size_t n, size_t m) {
  using Edge = WEdge<NodeId, EdgeTy>;
  Graph<NodeId, EdgeId, EdgeTy> G;
  G.n = n;
  G.m = m;
  parlay::sort_inplace(
      parlay::make_slice(edgelist),
      [](const std::pair<NodeId, Edge> &a, const std::pair<NodeId, Edge> &b) {
        if (a.first != b.first) {
          return a.first < b.first;
        }
        return a.second.v < b.second.v;
      });
  G.offsets = parlay::sequence<EdgeId>(n + 1, m);
  G.edges = parlay::sequence<Edge>(m);
  parlay::parallel_for(0, m, [&](size_t i) {
    G.edges[i] = edgelist[i].second;
    if (i == 0 || edgelist[i].first != edgelist[i - 1].first) {
      G.offsets[edgelist[i].first] = i;
    }
  });
  parlay::scan_inclusive_inplace(
      parlay::make_slice(G.offsets.rbegin(), G.offsets.rend()),
      parlay::minm<EdgeId>());
  return G;
}

template <class Graph>
Graph make_symmetrized(const Graph &G) {
  size_t n = G.n;
  size_t m = G.m;
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;
  using EdgeTy = typename Graph::EdgeTy;
  using Edge = typename Graph::Edge;
  parlay::sequence<std::pair<NodeId, Edge>> edgelist(m * 2);
  parlay::parallel_for(0, n, [&](NodeId u) {
    parlay::parallel_for(G.offsets[u], G.offsets[u + 1], [&](EdgeId i) {
      NodeId v = G.edges[i].v;
      EdgeTy w = G.edges[i].w;
      edgelist[i * 2] = std::make_pair(u, Edge(v, w));
      edgelist[i * 2 + 1] = std::make_pair(v, Edge(u, w));
    });
  });
  sort_inplace(make_slice(edgelist));
  auto pred = parlay::delayed_seq<bool>(m * 2, [&](size_t i) {
    if (i != 0 && edgelist[i].first == edgelist[i - 1].first &&
        edgelist[i].second.v == edgelist[i - 1].second.v) {
      return false;
    } else if (edgelist[i].first == edgelist[i].second.v) {
      return false;
    }
    return true;
  });
  edgelist = parlay::pack(make_slice(edgelist), pred);
  return edgelist2graph<NodeId, EdgeId, EdgeTy>(edgelist, n, edgelist.size());
}

template <class Graph>
Graph Transpose(const Graph &G) {
  size_t n = G.n;
  size_t m = G.m;
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;
  using EdgeTy = typename Graph::EdgeTy;
  using Edge = typename Graph::Edge;
  parlay::sequence<std::pair<NodeId, Edge>> edgelist(m);
  parlay::parallel_for(0, n, [&](NodeId u) {
    parlay::parallel_for(G.offsets[u], G.offsets[u + 1], [&](EdgeId i) {
      edgelist[i] = std::make_pair(G.edges[i].v, Edge(u, G.edges[i].w));
    });
  });
  return edgelist2graph<NodeId, EdgeId, EdgeTy>(edgelist, n, m);
}

#endif  // GRAPH_H

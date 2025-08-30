// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef _BENCH_GRAPH_IO
#define _BENCH_GRAPH_IO

#include <iostream>
#include <stdint.h>
#include <cstring>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "parallel.h"
#include "IO.h"

using namespace benchIO;

static uint32_t generate_random_weight(uint64_t edge_id) {
  // Use hash of edge ID to generate deterministic weight
  uint32_t hash_val = edge_id * 2654435761U;
  hash_val ^= hash_val >> 16;
  hash_val *= 2654435761U;
  hash_val ^= hash_val >> 16;
  return (hash_val % 10000) + 1;
}
template <class intT>
int xToStringLen(edge<intT> a) { 
  return xToStringLen(a.u) + xToStringLen(a.v) + 1;
}

template <class intT>
void xToString(char* s, edge<intT> a) { 
  int l = xToStringLen(a.u);
  xToString(s, a.u);
  s[l] = ' ';
  xToString(s+l+1, a.v);
}

template <class intT>
int xToStringLen(wghEdge<intT> a) { 
  return xToStringLen(a.u) + xToStringLen(a.v) + xToStringLen(a.weight) + 2;
}

template <class intT>
void xToString(char* s, wghEdge<intT> a) { 
  int lu = xToStringLen(a.u);
  int lv = xToStringLen(a.v);
  xToString(s, a.u);
  s[lu] = ' ';
  xToString(s+lu+1, a.v);
  s[lu+lv+1] = ' ';
  xToString(s+lu+lv+2, a.weight);
}

namespace benchIO {
  using namespace std;

  string AdjGraphHeader = "AdjacencyGraph";
  string EdgeArrayHeader = "EdgeArray";
  string WghEdgeArrayHeader = "WeightedEdgeArray";
  string WghAdjGraphHeader = "WeightedAdjacencyGraph";

  template <class intT>
  int writeGraphToFile(graph<intT> G, char* fname) {
    intT m = G.m;
    intT n = G.n;
    intT totalLen = 2 + n + m;
    intT *Out = newA(intT, totalLen);
    Out[0] = n;
    Out[1] = m;
    parallel_for (intT i=0; i < n; i++) {
      Out[i+2] = G.V[i].degree;
    }
    intT total = sequence::scan(Out+2,Out+2,n,utils::addF<intT>(),(intT)0);
    for (intT i=0; i < n; i++) {
      intT *O = Out + (2 + n + Out[i+2]);
      vertex<intT> v = G.V[i];
      for (intT j = 0; j < v.degree; j++) 
	O[j] = v.Neighbors[j];
    }
    int r = writeArrayToFile(AdjGraphHeader, Out, totalLen, fname);
    free(Out);
    return r;
  }

  template <class intT>
  int writeWghGraphToFile(wghGraph<intT> G, char* fname) {
    intT m = G.m;
    intT n = G.n;
    intT totalLen = 2 + n + m*2;
    intT *Out = newA(intT, totalLen);
    Out[0] = n;
    Out[1] = m;
    parallel_for (intT i=0; i < n; i++) {
      Out[i+2] = G.V[i].degree;
    }
    intT total = sequence::scan(Out+2,Out+2,n,utils::addF<intT>(),(intT)0);
    for (intT i=0; i < n; i++) {
      intT *O = Out + (2 + n + Out[i+2]);
      wghVertex<intT> v = G.V[i];
      for (intT j = 0; j < v.degree; j++) {
	O[j] = v.Neighbors[j];
	O[j+m] = v.nghWeights[j];
      }
    }
    int r = writeArrayToFile(WghAdjGraphHeader, Out, totalLen, fname);
    free(Out);
    return r;
  }


  template <class intT>
  int writeEdgeArrayToFile(edgeArray<intT> EA, char* fname) {
    intT m = EA.nonZeros;
    int r = writeArrayToFile(EdgeArrayHeader, EA.E, m, fname);
    return r;
  }

  template <class intT>
  int writeWghEdgeArrayToFile(wghEdgeArray<intT> EA, char* fname) {
    uintT m = EA.m;
    int r = writeArrayToFile(WghEdgeArrayHeader, EA.E, m, fname);
    return r;
  }

  template <class intT>
  edgeArray<intT> readEdgeArrayFromFile(char* fname) {
    _seq<char> S = readStringFromFile(fname);
    words W = stringToWords(S.A, S.n);
    if (W.Strings[0] != EdgeArrayHeader) {
      cout << "Bad input file" << endl;
      abort();
    }
    long n = (W.m-1)/2;
    edge<intT> *E = newA(edge<intT>,n);
    {parallel_for(long i=0; i < n; i++)
      E[i] = edge<intT>(atol(W.Strings[2*i + 1]), 
		  atol(W.Strings[2*i + 2]));}
    //W.del(); // to deal with performance bug in malloc

    intT maxR = 0;
    intT maxC = 0;
    for (long i=0; i < n; i++) {
      maxR = max<intT>(maxR, E[i].u);
      maxC = max<intT>(maxC, E[i].v);
    }
    intT maxrc = max<intT>(maxR,maxC) + 1;
    return edgeArray<intT>(E, maxrc, maxrc, n);
  }

  template <class intT>
  wghEdgeArray<intT> readWghEdgeArrayFromFile(char* fname) {
    _seq<char> S = readStringFromFile(fname);
    words W = stringToWords(S.A, S.n);
    if (W.Strings[0] != WghEdgeArrayHeader) {
      cout << "Bad input file" << endl;
      abort();
    }
    long n = (W.m-1)/3;
    wghEdge<intT> *E = newA(wghEdge<intT>,n);
    {parallel_for(long i=0; i < n; i++)
      E[i] = wghEdge<intT>(atol(W.Strings[3*i + 1]), 
			   atol(W.Strings[3*i + 2]),
			   atof(W.Strings[3*i + 3]));}
    //W.del(); // to deal with performance bug in malloc

    intT maxR = 0;
    intT maxC = 0;
    for (long i=0; i < n; i++) {
      maxR = max<intT>(maxR, E[i].u);
      maxC = max<intT>(maxC, E[i].v);
    }
    return wghEdgeArray<intT>(E, max<intT>(maxR,maxC)+1, n);
  }

  template <class intT>
  graph<intT> readGraphFromFile(char* fname) {
    _seq<char> S = readStringFromFile(fname);
    words W = stringToWords(S.A, S.n);
    if (W.Strings[0] != AdjGraphHeader) {
      cout << "Bad input file: missing header: " << AdjGraphHeader << endl;
      abort();
    }

    long len = W.m -1;
    uintT * In = newA(uintT, len);
    {parallel_for(long i=0; i < len; i++) In[i] = atol(W.Strings[i + 1]);}
    //W.del(); // to deal with performance bug in malloc
    
    long n = In[0];
    long m = In[1];

    if (len != n + m + 2) {
      cout << "Bad input file: length = "<<len<< " n+m+2 = " << n+m+2 << endl;
      abort();
    }
    vertex<intT> *v = newA(vertex<intT>,n);
    uintT* offsets = In+2;
    uintT* edges = In+2+n;

    parallel_for (uintT i=0; i < n; i++) {
      uintT o = offsets[i];
      uintT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
      v[i].degree = l;
      v[i].Neighbors = (intT*)(edges+o);
    }
    return graph<intT>(v,(intT)n,(uintT)m,(intT*)In);
  }

  template <class intT, class intE>
    graphC<intT, intE> readGraphCFromFile(char* fname) {
    _seq<char> S = readStringFromFile(fname);
    words W = stringToWords(S.A, S.n);
    if (W.Strings[0] != AdjGraphHeader) {
      cout << "Bad input file: missing header: " << AdjGraphHeader << endl;
      abort();
    }

    long len = W.m -1;
    long n = atol(W.Strings[1]);
    long m = atol(W.Strings[2]);

    if (len != n + m + 2) {
      cout << "Bad input file: length = "<<len<< " n+m+2 = " << n+m+2 << endl;
      abort();
    }

    intT * offsets = newA(intT, n+1);
    intE * edges = newA(intE,m);

    {parallel_for(long i=0; i < n; i++) offsets[i] = atol(W.Strings[i+3]);}
    offsets[n] = m;
    {parallel_for(long i=0; i < m; i++) edges[i] = atol(W.Strings[n+i+3]);}

    //W.del(); // to deal with performance bug in malloc
    
    return graphC<intT,intE>(offsets,edges,n,m);
  }

  template <class intT>
    wghGraph<intT> readWghGraphFromFile(char* fname) {
    _seq<char> S = readStringFromFile(fname);
    words W = stringToWords(S.A, S.n);
    if (W.Strings[0] != WghAdjGraphHeader) {
      cout << "Bad input file" << endl;
      abort();
    }
    
    long len = W.m -1;
    intT * In = newA(intT, len);
    {parallel_for(long i=0; i < len; i++) In[i] = atol(W.Strings[i + 1]);}
    //W.del(); // to deal with performance bug in malloc
    
    long n = In[0];
    long m = In[1];

    if (len != n + 2*m + 2) {
      cout << "Bad input file" << endl;
      abort();
    }
    wghVertex<intT> *v = newA(wghVertex<intT>,n);
    uintT* offsets = (uintT*)In+2;
    uintT* edges = (uintT*)In+2+n;
    intT* weights = In+2+n+m;
    parallel_for (uintT i=0; i < n; i++) {
      uintT o = offsets[i];
      uintT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
      v[i].degree = l;
      v[i].Neighbors = (intT*)(edges+o);
      v[i].nghWeights = (weights+o);
    }
    return wghGraph<intT>(v,(intT)n,(uintT)m,(intT*)In,weights);
  }

  void errorOut(const char* s) {
    cerr << s << endl;
    throw s;
  }
  void packInt64(int64_t x, uint8_t buf[8]) {
    uint64_t xu = x;
    for (int i = 0; i < 8; ++i)
      buf[i] = (xu >> (8 * i)) & 0xff;
  }
  int64_t unpackInt64(const uint8_t buf[8]) {
    uint64_t xu = 0;
    for (int i = 0; i < 8; ++i)
      xu |= ((uint64_t)buf[i]) << (i * 8);
    return (int64_t)xu;
  }

  void writeInt(ostream& out, char buf[8], int64_t x) {
    packInt64(x, (uint8_t*)buf);
    out.write(buf, 8);
  }
  int64_t readInt(istream& in, char buf[8]) {
    in.read(buf, 8);
    return unpackInt64((uint8_t*)buf);
  }
  template<typename intT>
  void writeFlowGraph(ostream& out, FlowGraph<intT> g) {
    char buf[8];
    out.write("FLOWFLOW", 8);
    writeInt(out, buf, g.g.n);
    writeInt(out, buf, g.g.m);
    writeInt(out, buf, g.source);
    writeInt(out, buf, g.sink);
    intT o = 0;
    for (intT i = 0; i < g.g.n; ++i) {
      writeInt(out, buf, o);
      o += g.g.V[i].degree;
    }
    for (intT i = 0; i < g.g.n; ++i) {
      wghVertex<intT>& v = g.g.V[i];
      for (intT j = 0; j < v.degree; ++j) {
        writeInt(out, buf, v.Neighbors[j]);
        writeInt(out, buf, v.nghWeights[j]);
      }
    }
  }
  template<typename intT>
  FlowGraph<intT> readFlowGraph(istream& in) {
    cout << "Starting readFlowGraph..." << endl;
    char buf[10];
    in.read(buf, 8);
    buf[8] = 0;
    cout << "Read header: " << buf << endl;
    
    // Check if it's FLOWFLOW format
    if (strcmp(buf, "FLOWFLOW") == 0) {
      // Original FLOWFLOW format reading
      intT n = readInt(in, buf);
      intT m = readInt(in, buf);
      intT S = readInt(in, buf);
      intT T = readInt(in, buf);
      intT *offset = newA(intT, n);
      intT* adj = newA(intT, m);
      intT* weights = newA(intT, m);
      wghVertex<intT>* v = newA(wghVertex<intT>, n);
      for (intT i = 0; i < n; ++i) {
        offset[i] = readInt(in, buf);
        v[i].Neighbors = adj + offset[i];
        v[i].nghWeights = weights + offset[i];
        if (i > 0)
          v[i - 1].degree = offset[i] - offset[i - 1];
      }
      v[n - 1].degree = m - offset[n - 1];
      free(offset);
      for (intT i = 0; i < m; ++i) {
        adj[i] = readInt(in, buf);
        weights[i] = readInt(in, buf);
      }
      return FlowGraph<intT>(wghGraph<intT>(v, n, m, adj, weights), S, T);
    } else {
      cout << "Binary format detected, reading header..." << endl;
      // Assume binary format - rewind and read as binary
      in.seekg(0, ios::beg);
      
      // Read binary format header
      uint64_t n, m, sizes;
      in.read(reinterpret_cast<char*>(&n), sizeof(uint64_t));
      in.read(reinterpret_cast<char*>(&m), sizeof(uint64_t));
      in.read(reinterpret_cast<char*>(&sizes), sizeof(uint64_t));
      
      cout << "Binary graph: n=" << n << ", m=" << m << ", sizes=" << sizes << endl;
      
      // Verify expected size
      if (sizes != (n + 1) * 8 + m * 4 + 3 * 8) {
        errorOut("Invalid binary graph file format");
      }
      
      // Read offsets
      cout << "Reading offsets array..." << endl;
      uint64_t* offsets_raw = newA(uint64_t, n + 1);
      in.read(reinterpret_cast<char*>(offsets_raw), (n + 1) * sizeof(uint64_t));
      
      // Read edge neighbors
      cout << "Reading edge neighbors..." << endl;
      uint32_t* neighbors_raw = newA(uint32_t, m);
      in.read(reinterpret_cast<char*>(neighbors_raw), m * sizeof(uint32_t));
      
      // Convert to flow graph format
      intT* adj = newA(intT, m);
      intT* weights = newA(intT, m);
      wghVertex<intT>* v = newA(wghVertex<intT>, n);
      
      for (intT i = 0; i < n; ++i) {
        v[i].Neighbors = adj + offsets_raw[i];
        v[i].nghWeights = weights + offsets_raw[i];
        v[i].degree = offsets_raw[i + 1] - offsets_raw[i];
      }
      
      for (intT i = 0; i < m; ++i) {
        adj[i] = static_cast<intT>(neighbors_raw[i]);
        // Generate deterministic weight based on edge ID
        weights[i] = generate_random_weight(i);
      }
      
      free(offsets_raw);
      free(neighbors_raw);
      
      // Default source and sink
      intT S = 0;
      intT T = (n > 1) ? 1 : 0;
      
      cout << "Read binary graph with " << n << " nodes, " << m << " edges" << endl;
      cout << "Added random weights in range [1, 10000]" << endl;
      cout << "Using default source=" << S << ", sink=" << T << endl;
      
      return FlowGraph<intT>(wghGraph<intT>(v, n, m, adj, weights), S, T);
    }
  }

  const char nl = '\n';
  template <typename intT>
  FlowGraph<intT> writeFlowGraphDimacs(ostream& out, FlowGraph<intT> g) {
    out << "c DIMACS flow network description" << nl;
    out << "c (problem-id, nodes, arcs)" << nl;
    out << "p max " << g.g.n << " " << g.g.m << nl;

    out << "c source" << nl;
    out << "n " << g.source + 1 << " s" << nl;
    out << "c sink" << nl;
    out << "n " << g.sink + 1 << " t" << nl;

    out << "c arc description (from, to, capacity)" << nl;

    for (intT i = 0; i < g.g.n; ++i) {
      wghVertex<intT>& v = g.g.V[i];
      for (intT j = 0; j < v.degree; ++j) {
        out << "a " << i + 1 << " " << v.Neighbors[j] + 1 << " "
            << v.nghWeights[j] << nl;
      }
    }
  }

  template<typename intT>
  struct intWghEdge {
    intT from, to, w;
  };
  int readDimacsLinePref(istream& in, const char* expected) {
    char type;
    while (in >> type) {
      if (type == 'c') {
        while (in.peek() != EOF && in.peek() != '\n')
          in.ignore();
        in >> ws;
        continue;
      } else if (!strchr(expected, type)) {
        errorOut((string("Unexpected DIMACS line (expected 'c' or one of '")
                      + expected + "')").c_str());
      }
      return type;
    }
    return EOF;
  }

  template <typename intT>
  FlowGraph<intT> readFlowGraphDimacs(istream& in) {
    string tmp;
    intT n, m;
    int type = readDimacsLinePref(in, "p");
    if (type == EOF)
      errorOut("Unexpected EOF while reading DIMACS file");
    in >> tmp >> n >> m;
    intWghEdge<intT>* edges = newA(intWghEdge<intT>, m);
    intT edgei = 0;
    intT* pos = newA(intT, n + 1);
    intT S = -1, T = -1;
    while (EOF != (type = readDimacsLinePref(in, "an"))) {
      if (type == 'n') {
        intT x;
        char st;
        in >> x >> st;
        x--;
        if (st == 's') S = x;
        else T = x;
      } else { // type == 'a'
        intT from, to, cap;
        in >> from >> to >> cap;
        from--; to--;
        edges[edgei] = (intWghEdge<intT>) { from, to, cap };
        edgei++;
        pos[from + 1]++;
      }
    }
    if (S < 0)
      errorOut("No source was specified in DIMACS input file");
    if (T < 0)
      errorOut("No sink was specified in DIMACS input file");
    if (m != edgei)
      errorOut("Inconsistent edge count in DIMACS input file");
    intT* adj = newA(intT, m);
    intT* weights = newA(intT, m);
    wghVertex<intT>* v = newA(wghVertex<intT>, n);
    for (intT i = 0; i < n; ++i) {
      pos[i + 1] += pos[i];
      v[i].Neighbors = adj + pos[i];
      v[i].nghWeights = weights + pos[i];
      v[i].degree = pos[i + 1] - pos[i];
    }
    for (intT i = 0; i < m; ++i) {
      intT& p = pos[edges[i].from];
      adj[p] = edges[i].to;
      weights[p] = edges[i].w;
      p++;
    }
    free(edges);
    free(pos);
    return FlowGraph<intT>(wghGraph<intT>(v, n, m, adj, weights), S, T);
  }

  template<typename intT>
  FlowGraph<intT> readFlowGraphBinary(const char* filename) {
    cout << "Reading binary format with memory mapping..." << endl;
    
    // Use memory mapping for fast reading
    struct stat sb;
    int fd = open(filename, O_RDONLY);
    if (fd == -1) {
      errorOut("Cannot open file for memory mapping");
    }
    if (fstat(fd, &sb) == -1) {
      errorOut("Unable to acquire file stat");
    }
    
    char *data = static_cast<char *>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
    if (data == MAP_FAILED) {
      errorOut("Memory mapping failed");
    }
    
    // Read binary format header from memory mapped data
    uint64_t n = reinterpret_cast<uint64_t *>(data)[0];
    uint64_t m = reinterpret_cast<uint64_t *>(data)[1];
    uint64_t sizes = reinterpret_cast<uint64_t *>(data)[2];
    
    cout << "Memory mapped binary graph: n=" << n << ", m=" << m << ", sizes=" << sizes << endl;
    
    // Verify expected size
    if (sizes != (n + 1) * 8 + m * 4 + 3 * 8) {
      errorOut("Invalid binary graph file format");
    }
    
    // Allocate arrays for the graph
    intT* adj = newA(intT, m);
    intT* weights = newA(intT, m);
    wghVertex<intT>* v = newA(wghVertex<intT>, n);
    
    // Get pointers to data sections
    uint64_t* offsets_raw = reinterpret_cast<uint64_t*>(data + 3 * 8);
    uint32_t* neighbors_raw = reinterpret_cast<uint32_t*>(data + 3 * 8 + (n + 1) * 8);
    
    cout << "Converting to flow graph format..." << endl;
    
    // Set up vertex structure
    for (intT i = 0; i < n; ++i) {
      intT start = offsets_raw[i];
      intT end = offsets_raw[i + 1];
      intT degree = end - start;
      
      v[i].degree = degree;
      v[i].Neighbors = adj + start;
      v[i].nghWeights = weights + start;
      
      // Copy neighbors and generate random weights
      for (intT j = 0; j < degree; j++) {
        adj[start + j] = neighbors_raw[start + j];
        // Generate deterministic weight based on edge ID
        weights[start + j] = generate_random_weight(start + j);
      }
    }
    
    // Unmap the file
    munmap(data, sb.st_size);
    close(fd);
    
    cout << "Binary graph loaded successfully, added random weights" << endl;
    return FlowGraph<intT>(wghGraph<intT>(v, n, m, adj, weights), 0, (n > 1) ? 1 : 0);
  }

  template<typename intT>
  FlowGraph<intT> readFlowGraphFromAdj(const char* filename) {
    cout << "Reading .adj file: " << filename << endl;
    
    _seq<char> file_content = readStringFromFile(const_cast<char*>(filename));
    words W = stringToWords(file_content.A, file_content.n);
    
    if (W.m < 3) {
      errorOut("Bad .adj file: too few tokens");
    }
    
    string header(W.Strings[0]);
    bool weighted_input = false;
    
    if (header == "WeightedAdjacencyGraph") {
      weighted_input = true;
    } else if (header == "AdjacencyGraph") {
      weighted_input = false;
    } else {
      errorOut("Bad .adj file: invalid header");
    }
    
    intT n = atol(W.Strings[1]);
    intT m = atol(W.Strings[2]);
    
    cout << "Graph: n=" << n << ", m=" << m << ", weighted=" << weighted_input << endl;
    
    // Expected file format:
    // AdjacencyGraph (or WeightedAdjacencyGraph)
    // n
    // m
    // offset[0] offset[1] ... offset[n-1]
    // neighbor[0] neighbor[1] ... neighbor[m-1]
    // [weight[0] weight[1] ... weight[m-1]]  // only for WeightedAdjacencyGraph
    
    intT expected_tokens = n + m + 3;
    if (weighted_input) {
      expected_tokens += m;
    }
    
    if (W.m != expected_tokens) {
      cout << "Expected " << expected_tokens << " tokens, got " << W.m << endl;
      errorOut("Bad .adj file: incorrect number of tokens");
    }
    
    // Read offsets
    intT* offsets = newA(intT, n + 1);
    for (intT i = 0; i < n; i++) {
      offsets[i] = atol(W.Strings[i + 3]);
    }
    offsets[n] = m;
    
    // Read neighbors
    intT* adj = newA(intT, m);
    for (intT i = 0; i < m; i++) {
      adj[i] = atol(W.Strings[i + n + 3]);
    }
    cout << "finished reading" << endl;
    
    // Read or generate weights
    intT* weights = newA(intT, m);
    if (weighted_input) {
      // Read weights from file
      for (intT i = 0; i < m; i++) {
        weights[i] = atol(W.Strings[i + n + m + 3]);
      }
      cout << "Read weights from file" << endl;
    } else {
      // Generate deterministic weights
      for (intT i = 0; i < m; i++) {
        // Generate deterministic weight based on edge ID
        weights[i] = generate_random_weight(i);
      }
      cout << "Added random weights in range [1, 10000]" << endl;
    }
    
    // Create vertex array
    wghVertex<intT>* v = newA(wghVertex<intT>, n);
    for (intT i = 0; i < n; i++) {
      v[i].Neighbors = adj + offsets[i];
      v[i].nghWeights = weights + offsets[i];
      v[i].degree = offsets[i + 1] - offsets[i];
    }
    
    free(offsets);
    
    // Default source and sink
    intT source = 0;
    intT sink = (n > 1) ? 1 : 0;
    
    cout << "Successfully loaded .adj file with " << n << " nodes, " << m << " edges" << endl;
    cout << "Using default source=" << source << ", sink=" << sink << endl;
    
    return FlowGraph<intT>(wghGraph<intT>(v, n, m, adj, weights), source, sink);
  }
};

#endif // _BENCH_GRAPH_IO

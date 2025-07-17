# Makefile for FLOWFLOW graph reader test

ifdef GCC
CXX = g++
else
CXX = clang++
endif

CXXFLAGS = -std=c++20 -Wall -Wextra -Werror

INCLUDE_PATH = -I../../external/parlaylib/include/ -I../ -I../../external/ 

ifdef CILKPLUS
CXX = clang++
CXXFLAGS += -DPARLAY_CILKPLUS -DCILK -fcilkplus
else ifdef OPENCILK
CXXFLAGS += -DPARLAY_OPENCILK -DCILK -fopencilk
else ifdef SERIAL
CXXFLAGS += -DPARLAY_SEQUENTIAL
else
CXXFLAGS += -pthread
endif

ifdef DEBUG
CXXFLAGS += -DDEBUG -Og -g
else ifdef PERF
CXX = g++
CXXFLAGS += -Og -mcx16 -march=native -g
else ifdef MEMCHECK
CXXFLAGS += -Og -mcx16 -DPARLAY_SEQUENTIAL -g
else
CXXFLAGS += -O3 -mcx16 -march=native
endif

ifdef STDALLOC
CXXFLAGS += -DPARLAY_USE_STD_ALLOC
endif

TARGET = maxflow
SOURCE = maxflow.C

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -o $(TARGET) $(SOURCE)

clean:
	rm -f $(TARGET) test_graph.flowflow

.PHONY: clean test

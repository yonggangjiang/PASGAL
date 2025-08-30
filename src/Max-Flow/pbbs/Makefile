ifdef FORCE_CILK
  CILKP := 1
else
  OPENMP := 1
endif

include parallelDefs

# required files
GLOBAL_REQUIRE = utils.h parallel.h sequence.h graph.h blockRadixSort.h transpose.h stlParallelSort.h
LOCAL_REQUIRE = 
OBJS = maxFlow.o
PCFLAGS := $(PCFLAGS) -std=c++11 -Wall -g -static-libgcc
PLFLAGS := $(PLFLAGS) -lpthread -static-libgcc -static-libstdc++

include common/timeRequiredFiles
include MakeBench

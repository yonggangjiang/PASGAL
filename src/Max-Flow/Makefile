ifdef GCC
CC = g++
else
CC = clang++
endif

CPPFLAGS = -std=c++20 -Wall -Wextra -Werror

INCLUDE_PATH = -I../../external/parlaylib/include/ -I..

ifdef CILKPLUS
CC = clang++
CPPFLAGS += -DPARLAY_CILKPLUS -DCILK -fcilkplus
else ifdef OPENCILK
CPPFLAGS += -DPARLAY_OPENCILK -DCILK -fopencilk
else ifdef SERIAL
CPPFLAGS += -DPARLAY_SEQUENTIAL
else
CPPFLAGS += -pthread
endif

ifdef DEBUG
CPPFLAGS += -DDEBUG -Og -g
else ifdef PERF
CC = g++
CPPFLAGS += -Og -mcx16 -march=native -g
else ifdef MEMCHECK
CPPFLAGS += -Og -mcx16 -DPARLAY_SEQUENTIAL -g
else
CPPFLAGS += -O3 -mcx16 -march=native
endif

ifdef STDALLOC
CPPFLAGS += -DPARLAY_USE_STD_ALLOC
endif

ifdef FLOAT
CPPFLAGS += -DFLOAT
endif

all: dinic push-relabel

dinic:	dinic.cpp dinic.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) dinic.cpp -o dinic

push-relabel:	push-relabel.cpp push-relabel.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) push-relabel.cpp -o push-relabel

clean:
	rm dinic push-relabel

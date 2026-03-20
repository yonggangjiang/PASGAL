####################################################################################
#                                COMPILER FLAGS                                    #
####################################################################################
#                                                                                  #
# -DSTATS        : Display operation counts.                                       #
# -DDISPLAY_CUT  : Output includes the nodes in the source set of the minimum cut. #
# -DDISPLAY_FLOW : Displays flow values for each arc.                              #
#                                                                                  #
# Build modes:                                                                     #
#   make              - Normal build (optimized)                                   #
#   make DEBUG=1      - Debug build (with sanitizers and debug flags)              #
#                                                                                  #
####################################################################################

# Base flags
CFLAGS_BASE=-O4

# Debug flags
CFLAGS_DEBUG=-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fsanitize=undefined -fsanitize=address -g -fno-sanitize-recover

# Set CFLAGS based on DEBUG variable
ifdef DEBUG
    CFLAGS=$(CFLAGS_BASE) $(CFLAGS_DEBUG)
else
    CFLAGS=$(CFLAGS_BASE)
endif

BINDIR=bin
CC=g++
INCLUDE_PATH = -I../../../external/parlaylib/include

all: ${BINDIR}/ pseudo_fifo

pseudo_fifo: ${BINDIR}/pseudo_fifo

${BINDIR}/pseudo_fifo: src/3.23/pseudo.cpp makefile
	${CC} ${CFLAGS} ${INCLUDE_PATH} -DFIFO_BUCKET src/3.23/pseudo.cpp -o bin/pseudo_fifo

serial: ${BINDIR}/pseudo_fifo_serial

${BINDIR}/pseudo_fifo_serial: src/3.23/pseudo.cpp makefile
	${CC} ${CFLAGS} ${INCLUDE_PATH} -DFIFO_BUCKET -DPARLAY_SEQUENTIAL src/3.23/pseudo.cpp -o bin/pseudo_fifo_serial

clean:
	rm -f ${BINDIR}/*

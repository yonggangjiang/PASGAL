#!/bin/bash

# Usage function
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -g, --graph GRAPH        Specify input graph file (required)"
    echo "  -a, --algorithm ALGO     Algorithm: basic, nondet, dinic, or pbbs (default: run all)"
    echo "  -t, --threads THREADS    Thread counts (comma-separated, e.g., 1,2,4,8)"
    echo "  -p, --path PATH          Graph path directory (default: /data/graphs/)"
    echo "  -s, --symmetrized        Use symmetrized graph flag"
    echo "  -o, --output FILE        Output TSV file prefix (default: max-flow-results)"
    echo "  --timeout SECONDS        Timeout for each run in seconds (default: no timeout)"
    echo "  -h, --help               Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 -g twitter_sym.bin -a nondet -t 1,2,4,8"
    echo "  $0 -g asia_sym.bin -a basic -t 4,16,64 -s"
    echo "  $0 -g soc-LiveJournal1_sym.bin -a pbbs -t 1,4,8,16 --timeout 300"
    echo "  $0 -g soc-LiveJournal1_sym.bin -t 1,4,8,16  # Run all algorithms"
}

# Default values
GRAPH=""
ALGORITHM=""  # Empty means run all algorithms
THREAD_COUNTS="1,4,16,64"
GRAPH_PATH="/data/graphs/"
SYMMETRIZED=""
OUTPUT_PREFIX="max-flow-results"
TIMEOUT="600"  # Empty means no timeout

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -g|--graph)
            GRAPH="$2"
            shift 2
            ;;
        -a|--algorithm)
            ALGORITHM="$2"
            shift 2
            ;;
        -t|--threads)
            THREAD_COUNTS="$2"
            shift 2
            ;;
        -p|--path)
            GRAPH_PATH="$2"
            shift 2
            ;;
        -s|--symmetrized)
            SYMMETRIZED="-s"
            shift
            ;;
        -o|--output)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        --timeout)
            TIMEOUT="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Check if graph is specified
if [[ -z "$GRAPH" ]]; then
    echo "Error: Graph file must be specified with -g or --graph"
    usage
    exit 1
fi

# Convert thread counts to array
IFS=',' read -ra THREADS <<< "$THREAD_COUNTS"

# Determine which algorithms to run
if [[ -z "$ALGORITHM" ]]; then
    ALGORITHMS=("basic" "nondet" "dinic" "pbbs")
    echo "No algorithm specified. Running all algorithms: ${ALGORITHMS[*]}"
else
    # Validate algorithm
    if [[ "$ALGORITHM" != "basic" && "$ALGORITHM" != "nondet" && "$ALGORITHM" != "dinic" && "$ALGORITHM" != "pbbs" ]]; then
        echo "Error: Algorithm must be 'basic', 'nondet', 'dinic', or 'pbbs'"
        exit 1
    fi
    ALGORITHMS=("$ALGORITHM")
fi

# Full path to graph
FULL_GRAPH_PATH="${GRAPH_PATH}${GRAPH}"

echo "=== Max Flow Thread Scaling Test ==="
echo "Graph: $FULL_GRAPH_PATH"
echo "Algorithms: ${ALGORITHMS[*]}"
echo "Thread counts: ${THREADS[*]}"
echo "Symmetrized: ${SYMMETRIZED:-"No"}"
echo "Timeout: ${TIMEOUT:-"No timeout"}"
echo "========================================"

# Function to get clean output filename
get_output_filename() {
    local algo=$1
    case $algo in
        "basic") echo "Result_push-relabel-basic.tsv" ;;
        "nondet") echo "Result_push-relabel-nondet.tsv" ;;
        "dinic") echo "Result_dinic.tsv" ;;
        "pbbs") echo "Result_pbbs.tsv" ;;
    esac
}

# Create combined results file if running all algorithms
if [[ ${#ALGORITHMS[@]} -gt 1 ]]; then
    COMBINED_OUTPUT="Results.tsv"
    echo -e "Algorithm\tGraph\tThreads\tMaxFlow\tTime(s)" > "$COMBINED_OUTPUT"
    echo "Combined results will be saved to: $COMBINED_OUTPUT"
fi

# Function to get taskset command based on thread count
get_taskset_cmd() {
    local threads=$1
    case $threads in
        1) echo "taskset -c 0-3:4" ;;
        2) echo "taskset -c 0-7:4" ;;
        4) echo "taskset -c 0-15:4" ;;
        8) echo "taskset -c 0-31:4" ;;
        16) echo "taskset -c 0-63:4" ;;
        24) echo "taskset -c 0-95:4" ;;
        48) echo "taskset -c 0-95:2" ;;
        96) echo "taskset -c 0-95" ;;
        *) echo "" ;;  # No taskset for other thread counts
    esac
}

# Main execution loop for each algorithm
for current_algo in "${ALGORITHMS[@]}"; do
    echo ""
    echo "=========================================="
    echo "Running algorithm: $current_algo"
    echo "=========================================="
    
    # Create individual output file
    OUTPUT_FILE=$(get_output_filename "$current_algo")
    echo -e "Graph\tThreads\tMaxFlow\tTime(s)" > "$OUTPUT_FILE"
    echo "Results for $current_algo will be saved to: $OUTPUT_FILE"

    # Run tests for each thread count
    for thread_count in "${THREADS[@]}"; do
        echo ""
        echo "Testing $current_algo with $thread_count threads..."
        
        # Drop caches if running as root
        if [[ $EUID -eq 0 ]]; then
            echo "Dropping caches..."
            sync && echo 3 > /proc/sys/vm/drop_caches
        fi
        
        # Build appropriate version
        if [[ "$current_algo" == "dinic" ]]; then
            # Dinic has its own binary
            if [[ $thread_count -eq 1 ]]; then
                echo "Building dinic (serial version)..."
                make dinic SERIAL=1 -B > /dev/null 2>&1
            else
                echo "Building dinic (parallel version)..."
                make dinic -B > /dev/null 2>&1
            fi
        elif [[ "$current_algo" == "pbbs" ]]; then
            # PBBS maxFlow
            echo "Building pbbs/maxFlow..."
            make -C pbbs maxFlow -B > /dev/null 2>&1
        else
            # Push-relabel algorithms
            if [[ $thread_count -eq 1 ]]; then
                echo "Building push-relabel (serial version)..."
                make push-relabel SERIAL=1 -B > /dev/null 2>&1
            else
                echo "Building push-relabel (parallel version)..."
                make push-relabel -B > /dev/null 2>&1
            fi
        fi
        
        # Set thread count
        if [[ "$current_algo" == "pbbs" ]]; then
            export OMP_NUM_THREADS=$thread_count
        else
            export PARLAY_NUM_THREADS=$thread_count
        fi
        
        # Get taskset command
        taskset_cmd=$(get_taskset_cmd $thread_count)
        
        # Prepare command based on algorithm
        if [[ "$current_algo" == "dinic" ]]; then
            # Dinic command
            if [[ -z "$taskset_cmd" ]]; then
                cmd="numactl -i all ./dinic -i $FULL_GRAPH_PATH $SYMMETRIZED"
            else
                cmd="$taskset_cmd numactl -i all ./dinic -i $FULL_GRAPH_PATH $SYMMETRIZED"
            fi
        elif [[ "$current_algo" == "pbbs" ]]; then
            # PBBS maxFlow command (no -i flag)
            if [[ -z "$taskset_cmd" ]]; then
                cmd="numactl -i all ./pbbs/maxFlow $FULL_GRAPH_PATH $SYMMETRIZED"
            else
                cmd="$taskset_cmd numactl -i all ./pbbs/maxFlow $FULL_GRAPH_PATH $SYMMETRIZED"
            fi
        else
            # Push-relabel commands  
            if [[ -z "$taskset_cmd" ]]; then
                cmd="numactl -i all ./push-relabel -i $FULL_GRAPH_PATH $SYMMETRIZED -a $current_algo"
            else
                cmd="$taskset_cmd numactl -i all ./push-relabel -i $FULL_GRAPH_PATH $SYMMETRIZED -a $current_algo"
            fi
        fi
        
        echo "Running: $cmd"
        if [[ -n "$TIMEOUT" ]]; then
            echo "Timeout: $TIMEOUT seconds"
        else
            echo "Timeout: No timeout"
        fi
        
        # Run the command with or without timeout and capture output
        if [[ -n "$TIMEOUT" ]]; then
            output=$(timeout "${TIMEOUT}s" $cmd 2>&1)
            exit_code=$?
        else
            output=$($cmd 2>&1)
            exit_code=$?
        fi
        # echo "$output"  # Uncomment this line if you want to see full output
        
        # Check if command timed out
        if [[ $exit_code -eq 124 && -n "$TIMEOUT" ]]; then
            echo "Warning: Command timed out after $TIMEOUT seconds"
            output="TIMEOUT: Command exceeded $TIMEOUT seconds"
        fi
        
        # Extract average time and max flow from output
        if [[ $exit_code -eq 124 && -n "$TIMEOUT" ]]; then
            # Command timed out
            avg_time="TIMEOUT"
            max_flow="TIMEOUT"
        elif [[ "$current_algo" == "pbbs" ]]; then
            # PBBS format: "PBBS-time:" for time and "flow=" for flow value
            avg_time=$(echo "$output" | grep "PBBS-time:" | tail -1 | awk '{print $2}')
            max_flow=$(echo "$output" | grep "flow=" | tail -1 | sed 's/.*flow=//' | awk '{print $1}')
        else
            # Standard format for other algorithms
            avg_time=$(echo "$output" | grep "Average time:" | tail -1 | awk '{print $3}')
            max_flow=$(echo "$output" | grep "Max flow:" | tail -1 | awk '{print $3}')
        fi
        
        # Write results to individual output file
        if [[ "$avg_time" == "TIMEOUT" ]]; then
            echo -e "${GRAPH}\t${thread_count}\tTIMEOUT\tTIMEOUT" >> "$OUTPUT_FILE"
            echo "Results: TIMEOUT (exceeded ${TIMEOUT} seconds)"
            
            # Also write to combined file if running multiple algorithms
            if [[ ${#ALGORITHMS[@]} -gt 1 ]]; then
                echo -e "${current_algo}\t${GRAPH}\t${thread_count}\tTIMEOUT\tTIMEOUT" >> "$COMBINED_OUTPUT"
            fi
        elif [[ -n "$avg_time" && -n "$max_flow" ]]; then
            echo -e "${GRAPH}\t${thread_count}\t${max_flow}\t${avg_time}" >> "$OUTPUT_FILE"
            echo "Results: Max Flow = $max_flow, Average Time = $avg_time seconds"
            
            # Also write to combined file if running multiple algorithms
            if [[ ${#ALGORITHMS[@]} -gt 1 ]]; then
                echo -e "${current_algo}\t${GRAPH}\t${thread_count}\t${max_flow}\t${avg_time}" >> "$COMBINED_OUTPUT"
            fi
        else
            echo "Warning: Could not extract results from output"
            echo -e "${GRAPH}\t${thread_count}\tERROR\tERROR" >> "$OUTPUT_FILE"
            
            # Also write error to combined file if running multiple algorithms
            if [[ ${#ALGORITHMS[@]} -gt 1 ]]; then
                echo -e "${current_algo}\t${GRAPH}\t${thread_count}\tERROR\tERROR" >> "$COMBINED_OUTPUT"
            fi
        fi
        
        echo "Completed $current_algo with $thread_count threads"
    done
    
    echo "Completed all tests for algorithm: $current_algo"
    echo "Individual results saved to: $OUTPUT_FILE"
done

echo ""
echo "=========================================="
echo "All tests completed!"
echo "=========================================="

if [[ ${#ALGORITHMS[@]} -gt 1 ]]; then
    echo "Combined results saved to: $COMBINED_OUTPUT"
    echo "Individual results saved to:"
    for algo in "${ALGORITHMS[@]}"; do
        echo "  - $(get_output_filename "$algo")"
    done
else
    echo "Results saved to: $(get_output_filename "${ALGORITHMS[0]}")"
fi

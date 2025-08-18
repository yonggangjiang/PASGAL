#!/bin/bash

# PASGAL Comprehensive Scalability Testing Script
# Tests maxflow algorithms with different thread counts
# Computes graph properties and measures performance
#
# Usage: ./test_scalability.sh [input_directory]
# If no directory specified, uses "inputs" as default

# Check for input directory parameter
if [ $# -gt 0 ]; then
    INPUTS_DIR="$1"
else
    INPUTS_DIR="inputs"
fi

# Configuration
RESULTS_FILE="scalability_results_$(basename "$INPUTS_DIR").out"
PBBS_MAXFLOW="src/pbbs/maxFlow"
PR_SEQ_MAXFLOW="src/PR_seq/maxflow"
COMPUTE_DIAMETER="./compute_diameter"
THREADS=(1 2 4 8 16 32 64)
TIMEOUT=60  # 1 minute timeout per test
NUM_RUNS=3   # Number of runs for averaging

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to extract timing from pbbs output
extract_pbbs_time() {
    local output="$1"
    # Look for "PBBS-time:" entry (main algorithm time)
    echo "$output" | grep "PBBS-time:" | grep -o "[0-9.]\+" | head -1
}

# Function to extract flow value from pbbs output
extract_pbbs_flow() {
    local output="$1"
    echo "$output" | grep -o "flow=[0-9]*" | cut -d'=' -f2 | head -1
}

# Function to extract timing from PR_seq output
extract_pr_seq_time() {
    local output="$1"
    echo "$output" | grep -o "Elapsed time (s): [0-9]*\.[0-9e\-]*" | cut -d' ' -f4 | head -1
}

# Function to extract flow value from PR_seq output
extract_pr_seq_flow() {
    local output="$1"
    echo "$output" | grep -o "Max flow value: [0-9]*" | grep -o "[0-9]*"
}

# Function to run multiple times and get average
run_pbbs_multiple() {
    local graph="$1"
    local threads="$2"
    local times=()
    local flow=""
    
    for ((i=1; i<=NUM_RUNS; i++)); do
        local output
        output=$(timeout ${TIMEOUT}s env OMP_NUM_THREADS=$threads $PBBS_MAXFLOW "$graph" 2>&1)
        if [ $? -eq 0 ]; then
            local time=$(extract_pbbs_time "$output")
            local curr_flow=$(extract_pbbs_flow "$output")
            if [ -n "$time" ] && [ -n "$curr_flow" ]; then
                times+=($time)
                flow="$curr_flow"
            else
                echo "ERROR: Could not extract results from run $i"
                return 1
            fi
        else
            echo "ERROR: PBBS run $i failed or timed out"
            return 1
        fi
    done
    
    # Calculate average time
    local sum=0
    for time in "${times[@]}"; do
        sum=$(echo "$sum + $time" | bc -l)
    done
    local avg=$(echo "scale=6; $sum / ${#times[@]}" | bc -l)
    
    echo "$flow:$avg"
    return 0
}

# Function to run PR_seq multiple times and get average
run_pr_seq_multiple() {
    local graph="$1"
    local times=()
    local flow=""
    
    for ((i=1; i<=NUM_RUNS; i++)); do
        local output
        output=$(timeout ${TIMEOUT}s $PR_SEQ_MAXFLOW "$graph" 2>&1)
        if [ $? -eq 0 ]; then
            local time=$(extract_pr_seq_time "$output")
            local curr_flow=$(extract_pr_seq_flow "$output")
            if [ -n "$time" ] && [ -n "$curr_flow" ]; then
                times+=($time)
                flow="$curr_flow"
            else
                echo "ERROR: Could not extract results from run $i"
                return 1
            fi
        else
            echo "ERROR: PR_seq run $i failed or timed out"
            return 1
        fi
    done
    
    # Calculate average time
    local sum=0
    for time in "${times[@]}"; do
        sum=$(echo "$sum + $time" | bc -l)
    done
    local avg=$(echo "scale=6; $sum / ${#times[@]}" | bc -l)
    
    echo "$flow:$avg"
    return 0
}

# Function to calculate speedup
calculate_speedup() {
    local baseline="$1"
    local current="$2"
    if [ "$baseline" != "0" ] && [ "$current" != "0" ]; then
        echo "scale=2; $baseline / $current" | bc -l
    else
        echo "N/A"
    fi
}

# Function to calculate efficiency
calculate_efficiency() {
    local speedup="$1"
    local threads="$2"
    if [ "$speedup" != "N/A" ] && [ "$threads" != "0" ]; then
        echo "scale=2; $speedup / $threads * 100" | bc -l
    else
        echo "N/A"
    fi
}

# Check if required tools exist
for tool in "$PBBS_MAXFLOW" "$PR_SEQ_MAXFLOW" "$COMPUTE_DIAMETER"; do
    if [ ! -f "$tool" ] && [ ! -x "$tool" ]; then
        echo -e "${RED}ERROR: $tool not found or not executable${NC}"
        echo "Please build the required tools first."
        exit 1
    fi
done

# Check if bc is available for calculations
if ! command -v bc &> /dev/null; then
    echo -e "${RED}ERROR: 'bc' calculator not found. Please install bc.${NC}"
    exit 1
fi

echo -e "${BLUE}=== PASGAL Comprehensive Scalability Testing ===${NC}"
echo "Starting at: $(date)"
echo "Number of runs per test: $NUM_RUNS"
echo "Thread counts: ${THREADS[*]}"
echo "Timeout per run: ${TIMEOUT}s"

# Find graph files
echo -e "\n${YELLOW}Looking for graph files in $INPUTS_DIR...${NC}"
graph_files=()

# Look for files with .adj and .bin extensions
for ext in adj bin; do
    while IFS= read -r -d '' file; do
        graph_files+=("$file")
    done < <(find "$INPUTS_DIR" -name "*.$ext" -print0 2>/dev/null)
done

# Also look for files without extensions that could be flow graph format
while IFS= read -r -d '' file; do
    # Skip directories and files with extensions
    if [ -f "$file" ] && [[ "$(basename "$file")" != *.* ]]; then
        graph_files+=("$file")
    fi
done < <(find "$INPUTS_DIR" -maxdepth 1 -type f -print0 2>/dev/null)

if [ ${#graph_files[@]} -eq 0 ]; then
    echo -e "${RED}No graph files found in $INPUTS_DIR${NC}"
    exit 1
fi

echo "Found ${#graph_files[@]} graph files:"
for file in "${graph_files[@]}"; do
    echo "  - $(basename "$file")"
done

# Initialize results file
cat > "$RESULTS_FILE" << EOF
# PASGAL Comprehensive Scalability Testing Results
# Generated on: $(date)
# Format: Graph | Nodes | Edges | Density | Src_Sink_Dist | Undirected_Diam | PBBS_Results | PR_Seq_Results | Speedup_Analysis
# PBBS_Results format: threads:flow:avg_time;...
# Speedup_Analysis format: threads:speedup:efficiency;...

EOF

# Process each graph
for graph in "${graph_files[@]}"; do
    graph_name=$(basename "$graph")
    echo -e "\n${GREEN}=== Testing $graph_name ===${NC}"
    
    # Compute graph properties
    echo "Computing graph properties..."
    props_output=$(timeout 60s $COMPUTE_DIAMETER "$graph" 2>&1)
    if [ $? -eq 0 ]; then
        echo "Graph properties: $props_output"
        # Parse properties
        N=$(echo "$props_output" | grep -o "Nodes: [0-9]*" | grep -o "[0-9]*")
        M=$(echo "$props_output" | grep -o "Edges: [0-9]*" | grep -o "[0-9]*")
        density=$(echo "$props_output" | grep -o "Density: [0-9.]*" | grep -o "[0-9.]*")
        src_sink=$(echo "$props_output" | grep -o "Source-to-sink distance: [-0-9]*" | cut -d':' -f2 | tr -d ' ')
        diameter=$(echo "$props_output" | grep -o "Undirected diameter: [0-9]*" | cut -d':' -f2 | tr -d ' ')
    else
        echo -e "${RED}WARNING: Could not compute graph properties (timeout or error)${NC}"
        N="unknown"; M="unknown"; density="unknown"; src_sink="unknown"; diameter="unknown"
    fi
    
    # Test PBBS with different thread counts
    echo "Testing PBBS scalability..."
    pbbs_results=()
    pbbs_baseline_time=""
    
    for threads in "${THREADS[@]}"; do
        echo "  Testing with $threads threads..."
        result=$(run_pbbs_multiple "$graph" "$threads")
        if [ $? -eq 0 ]; then
            flow=$(echo "$result" | cut -d':' -f1)
            time=$(echo "$result" | cut -d':' -f2)
            echo -e "    ${GREEN}SUCCESS: flow=$flow, avg_time=${time}s${NC}"
            pbbs_results+=("$threads:$flow:$time")
            
            # Store baseline (1 thread) time for speedup calculation
            if [ "$threads" -eq 1 ]; then
                pbbs_baseline_time="$time"
            fi
        else
            echo -e "    ${RED}FAILED${NC}"
            pbbs_results+=("$threads:ERROR:ERROR")
        fi
    done
    
    # Test PR_seq
    echo "Testing PR_seq..."
    pr_seq_result=$(run_pr_seq_multiple "$graph")
    if [ $? -eq 0 ]; then
        pr_flow=$(echo "$pr_seq_result" | cut -d':' -f1)
        pr_time=$(echo "$pr_seq_result" | cut -d':' -f2)
        echo -e "  ${GREEN}SUCCESS: flow=$pr_flow, avg_time=${pr_time}s${NC}"
    else
        echo -e "  ${RED}FAILED${NC}"
        pr_flow="ERROR"
        pr_time="ERROR"
    fi
    
    # Calculate speedup and efficiency
    speedup_results=()
    for result in "${pbbs_results[@]}"; do
        threads=$(echo "$result" | cut -d':' -f1)
        time=$(echo "$result" | cut -d':' -f3)
        
        if [ "$time" != "ERROR" ] && [ -n "$pbbs_baseline_time" ] && [ "$pbbs_baseline_time" != "ERROR" ]; then
            speedup=$(calculate_speedup "$pbbs_baseline_time" "$time")
            efficiency=$(calculate_efficiency "$speedup" "$threads")
            speedup_results+=("$threads:$speedup:$efficiency")
        else
            speedup_results+=("$threads:N/A:N/A")
        fi
    done
    
    # Format results for file
    pbbs_formatted=$(IFS=";"; echo "${pbbs_results[*]}")
    speedup_formatted=$(IFS=";"; echo "${speedup_results[*]}")
    
    # Save to results file
    echo "$graph_name | $N | $M | $density | $src_sink | $diameter | $pbbs_formatted | $pr_flow:$pr_time | $speedup_formatted" >> "$RESULTS_FILE"
    
    echo -e "${GREEN}Completed $graph_name${NC}"
done

echo -e "\n${BLUE}=== Scalability Testing Completed ===${NC}"
echo "Results saved to: $RESULTS_FILE"

# Display summary
echo -e "\n${YELLOW}=== Results Summary ===${NC}"
cat "$RESULTS_FILE"

# Generate performance analysis
echo -e "\n${YELLOW}=== Performance Analysis ===${NC}"
echo "Comprehensive Scalability Analysis"
echo "=================================="
printf "%-25s | %8s | %10s | %8s | %8s | %10s | %12s | %15s | %s\n" "Graph" "Nodes" "Edges" "Density" "S-T_Dist" "Diameter" "PR_seq(s)" "Best_PBBS(s)" "Scalability"
echo "-------------------------|----------|------------|----------|----------|------------|-------------|----------------|------------"

# Also write performance analysis to file
echo "" >> "$RESULTS_FILE"
echo "=== Performance Analysis ===" >> "$RESULTS_FILE"
echo "Comprehensive Scalability Analysis" >> "$RESULTS_FILE"
echo "==================================" >> "$RESULTS_FILE"
printf "%-25s | %8s | %10s | %8s | %8s | %10s | %12s | %15s | %s\n" "Graph" "Nodes" "Edges" "Density" "S-T_Dist" "Diameter" "PR_seq(s)" "Best_PBBS(s)" "Scalability" >> "$RESULTS_FILE"
echo "-------------------------|----------|------------|----------|----------|------------|-------------|----------------|------------" >> "$RESULTS_FILE"

# Extract and analyze performance for each graph
for graph in "${graph_files[@]}"; do
    graph_name=$(basename "$graph")
    line=$(grep "^$graph_name" "$RESULTS_FILE")
    if [ -n "$line" ]; then
        # Parse graph properties
        N=$(echo "$line" | cut -d'|' -f2 | tr -d ' ')
        M=$(echo "$line" | cut -d'|' -f3 | tr -d ' ')
        density=$(echo "$line" | cut -d'|' -f4 | tr -d ' ')
        src_sink_dist=$(echo "$line" | cut -d'|' -f5 | tr -d ' ')
        diameter=$(echo "$line" | cut -d'|' -f6 | tr -d ' ')
        
        # Parse PR_seq results
        pr_seq_data=$(echo "$line" | cut -d'|' -f8 | tr -d ' ')
        pr_seq_time=$(echo "$pr_seq_data" | cut -d':' -f2)
        
        # Parse PBBS results to find best time
        pbbs_data=$(echo "$line" | cut -d'|' -f7 | tr -d ' ')
        IFS=';' read -ra PBBS_RESULTS <<< "$pbbs_data"
        
        best_pbbs_time="N/A"
        best_threads="N/A"
        speedup_1_to_best="N/A"
        
        # Find single-thread time and best time
        single_thread_time=""
        for result in "${PBBS_RESULTS[@]}"; do
            threads=$(echo "$result" | cut -d':' -f1)
            time=$(echo "$result" | cut -d':' -f3)
            
            if [ "$threads" = "1" ]; then
                single_thread_time="$time"
            fi
            
            if [ "$time" != "ERROR" ] && [ "$time" != "TIMEOUT" ] && [ -n "$time" ]; then
                if [ "$best_pbbs_time" = "N/A" ] || [ $(echo "$time < $best_pbbs_time" | bc -l 2>/dev/null || echo 0) -eq 1 ]; then
                    best_pbbs_time="$time"
                    best_threads="$threads"
                fi
            fi
        done
        
        # Calculate speedup from 1 thread to best
        if [ -n "$single_thread_time" ] && [ "$single_thread_time" != "0.00" ] && [ "$best_pbbs_time" != "N/A" ] && [ "$best_pbbs_time" != "0.00" ]; then
            speedup_1_to_best=$(echo "scale=2; $single_thread_time / $best_pbbs_time" | bc -l 2>/dev/null || echo "N/A")
        fi
        
        # Format scalability info
        if [ "$speedup_1_to_best" != "N/A" ] && [ "$best_threads" != "N/A" ]; then
            scalability="${speedup_1_to_best}x@${best_threads}T"
        else
            scalability="N/A"
        fi
        
        # Truncate graph name if too long
        short_name=$(echo "$graph_name" | cut -c1-24)
        
        printf "%-25s | %8s | %10s | %8s | %8s | %10s | %12s | %15s | %s\n" \
            "$short_name" "$N" "$M" "$density" "$src_sink_dist" "$diameter" "$pr_seq_time" "$best_pbbs_time" "$scalability"
        
        # Also write to file
        printf "%-25s | %8s | %10s | %8s | %8s | %10s | %12s | %15s | %s\n" \
            "$short_name" "$N" "$M" "$density" "$src_sink_dist" "$diameter" "$pr_seq_time" "$best_pbbs_time" "$scalability" >> "$RESULTS_FILE"
    fi
done

echo ""
echo "Detailed Thread Scalability:"
echo "----------------------------"

# Also write to file
echo "" >> "$RESULTS_FILE"
echo "Detailed Thread Scalability:" >> "$RESULTS_FILE"
echo "----------------------------" >> "$RESULTS_FILE"

# Detailed per-graph analysis
for graph in "${graph_files[@]}"; do
    graph_name=$(basename "$graph")
    line=$(grep "^$graph_name" "$RESULTS_FILE")
    if [ -n "$line" ]; then
        echo ""
        echo "Graph: $graph_name"
        
        # Also write to file
        echo "" >> "$RESULTS_FILE"
        echo "Graph: $graph_name" >> "$RESULTS_FILE"
        
        # Parse graph info
        N=$(echo "$line" | cut -d'|' -f2 | tr -d ' ')
        M=$(echo "$line" | cut -d'|' -f3 | tr -d ' ')
        density=$(echo "$line" | cut -d'|' -f4 | tr -d ' ')
        src_sink_dist=$(echo "$line" | cut -d'|' -f5 | tr -d ' ')
        diameter=$(echo "$line" | cut -d'|' -f6 | tr -d ' ')
        pr_seq_data=$(echo "$line" | cut -d'|' -f8 | tr -d ' ')
        pr_seq_flow=$(echo "$pr_seq_data" | cut -d':' -f1)
        pr_seq_time=$(echo "$pr_seq_data" | cut -d':' -f2)
        
        echo "  Properties: N=$N, M=$M, Density=$density, S-T_Dist=$src_sink_dist, Diameter=$diameter"
        echo "  PR_seq: Flow=$pr_seq_flow, Time=${pr_seq_time}s"
        
        # Also write to file
        echo "  Properties: N=$N, M=$M, Density=$density, S-T_Dist=$src_sink_dist, Diameter=$diameter" >> "$RESULTS_FILE"
        echo "  PR_seq: Flow=$pr_seq_flow, Time=${pr_seq_time}s" >> "$RESULTS_FILE"
        
        # Parse thread scalability
        pbbs_data=$(echo "$line" | cut -d'|' -f7 | tr -d ' ')
        echo "  PBBS Thread Scalability:"
        echo "    Threads | Flow | Time(s) | Speedup | Efficiency(%)"
        echo "    --------|------|---------|---------|-------------"
        
        # Also write to file
        echo "  PBBS Thread Scalability:" >> "$RESULTS_FILE"
        echo "    Threads | Flow | Time(s) | Speedup | Efficiency(%)" >> "$RESULTS_FILE"
        echo "    --------|------|---------|---------|-------------" >> "$RESULTS_FILE"
        
        IFS=';' read -ra RESULTS <<< "$pbbs_data"
        IFS=';' read -ra SPEEDUPS <<< "$(echo "$line" | cut -d'|' -f9 | tr -d ' ')"
        
        for i in "${!RESULTS[@]}"; do
            threads=$(echo "${RESULTS[$i]}" | cut -d':' -f1)
            flow=$(echo "${RESULTS[$i]}" | cut -d':' -f2)
            time=$(echo "${RESULTS[$i]}" | cut -d':' -f3)
            
            if [ ${#SPEEDUPS[@]} -gt $i ]; then
                speedup=$(echo "${SPEEDUPS[$i]}" | cut -d':' -f2)
                efficiency=$(echo "${SPEEDUPS[$i]}" | cut -d':' -f3)
            else
                speedup="N/A"
                efficiency="N/A"
            fi
            
            printf "    %7s | %4s | %7s | %7s | %11s\n" "$threads" "$flow" "$time" "$speedup" "$efficiency"
            
            # Also write to file
            printf "    %7s | %4s | %7s | %7s | %11s\n" "$threads" "$flow" "$time" "$speedup" "$efficiency" >> "$RESULTS_FILE"
        done
    fi
done

echo ""
echo -e "${GREEN}Testing completed successfully!${NC}"
echo -e "${BLUE}Input directory: $INPUTS_DIR${NC}"
echo -e "${BLUE}Results saved to: $RESULTS_FILE${NC}"

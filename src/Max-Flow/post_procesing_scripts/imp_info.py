#!/usr/bin/env python3
"""
Script to filter log files and extract only tests with at least a minimum node count.
Also computes relative speedup between serial and parallel versions for each phase.
"""

import argparse
import re
import os

MIN_NODES = 10**6

phase_name = [
    "", 
	"Get highest label roots", 
	"Process roots to find merges", 
	"Perform merges and push excess", 
	"Add new strong roots to buckets", 
	"Update highest strong label"
]

def extract_phase_times(block_lines):
    """
    Extract phase times for hpf-p (parallel) and hpf-p-serial from a test block.
    Returns dict with 'parallel' and 'serial' phase times.
    """
    parallel_phases = {}
    serial_phases = {}
    current_algo = None
    
    for line in block_lines:
        if 'Algorithm hpf-p got:' in line:
            current_algo = 'parallel'
        elif 'Algorithm hpf-p-serial got:' in line:
            current_algo = 'serial'
        elif 'Algorithm' in line and 'got:' in line:
            current_algo = None
        
        # Extract phase times.
        # Supports both:
        #   Total time for phase 1 = 0.123 seconds
        #   Total time for phase 1(Some phase name) = 0.123 seconds
        match = re.search(
            r'Total time for phase\s+(\d+)(?:\([^)]*\))?\s*=\s*([\d.]+)\s+seconds',
            line,
        )
        if match and current_algo:
            phase_num = int(match.group(1))
            phase_time = float(match.group(2))
            if current_algo == 'parallel':
                parallel_phases[phase_num] = phase_time
            else:
                serial_phases[phase_num] = phase_time
    
    return parallel_phases, serial_phases

def compute_speedups(parallel_phases, serial_phases):
    """
    Compute speedup (serial/parallel) for each phase.
    """
    speedups = {}
    all_phases = set(parallel_phases.keys()) | set(serial_phases.keys())
    for phase in sorted(all_phases):
        p_time = parallel_phases.get(phase, 0)
        s_time = serial_phases.get(phase, 0)
        if p_time > 0:
            speedups[phase] = s_time / p_time
        else:
            speedups[phase] = float('inf') if s_time > 0 else 0
    return speedups

def parse_log_file(log_path):
    """
    Parse the log file and return test blocks where Num nodes >= MIN_NODES,
    along with phase speedup information.
    """
    with open(log_path, 'r') as f:
        lines = f.readlines()
    
    filtered_results = []
    current_block = []
    num_nodes = None
    test_name = None
    
    for line in lines:
        # Check if this is a new test
        if 'Running test' in line:
            # Save previous block if node count is above threshold
            if current_block and num_nodes is not None and num_nodes >= MIN_NODES:
                parallel_phases, serial_phases = extract_phase_times(current_block)
                speedups = compute_speedups(parallel_phases, serial_phases)
                filtered_results.append({
                    'name': test_name,
                    'block': ''.join(current_block),
                    'num_nodes': num_nodes,
                    'parallel_phases': parallel_phases,
                    'serial_phases': serial_phases,
                    'speedups': speedups
                })
            # Start new block
            current_block = [line]
            num_nodes = None
            match = re.search(r'Running test (.+)$', line)
            test_name = match.group(1) if match else 'Unknown'
        else:
            current_block.append(line)
            # Check for node count
            if 'Num nodes =' in line:
                match = re.search(r'Num nodes\s*=\s*(\d+)', line)
                if match:
                    num_nodes = int(match.group(1))
    
    # Don't forget the last block
    if current_block and num_nodes is not None and num_nodes >= MIN_NODES:
        parallel_phases, serial_phases = extract_phase_times(current_block)
        speedups = compute_speedups(parallel_phases, serial_phases)
        filtered_results.append({
            'name': test_name,
            'block': ''.join(current_block),
            'num_nodes': num_nodes,
            'parallel_phases': parallel_phases,
            'serial_phases': serial_phases,
            'speedups': speedups
        })
    
    return filtered_results

def main():
    parser = argparse.ArgumentParser(
        description='Filter log files to extract tests with at least 1,000,000 nodes'
    )
    parser.add_argument('log_file', help='Path to the log file to process')
    parser.add_argument(
        '-o', '--output',
        help='Output file name (default: filtered_<input_name>)',
        default=None
    )
    args = parser.parse_args()
    
    # Determine output path - always output to processed_logs folder
    output_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(args.log_file))),
        'processed_logs'
    )
    if args.output:
        output_path = os.path.join(output_dir, args.output)
    else:
        input_basename = os.path.basename(args.log_file)
        output_path = os.path.join(output_dir, f'filtered_{input_basename}')
    
    # Parse and filter
    filtered_results = parse_log_file(args.log_file)
    
    # Write output
    with open(output_path, 'w') as f:
        f.write(f'# Filtered results from: {args.log_file}\n')
        f.write(f'# Tests where Num nodes >= {MIN_NODES}\n')
        f.write(f'# Total matching tests: {len(filtered_results)}\n\n')
        
        for result in filtered_results:
            f.write(result['block'])
            
            # Add speedup summary
            if result['speedups']:
                f.write('\n--- Phase Speedups (Serial / Parallel) ---\n')
                for phase in sorted(result['speedups'].keys()):
                    speedup = result['speedups'][phase]
                    p_time = result['parallel_phases'].get(phase, 0)
                    s_time = result['serial_phases'].get(phase, 0)
                    if speedup == float('inf'):
                        f.write(f'Phase {phase}: N/A (parallel time = 0)\n')
                    else:
                        f.write(f'Phase {phase} ({phase_name[phase]}): {speedup:.3f}x (serial: {s_time:.4f}s, parallel: {p_time:.4f}s)\n')
                f.write('-------------------------------------------\n')
            f.write('\n')
    
    print(f'Filtered {len(filtered_results)} tests to: {output_path}')

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Script to generate figures from log data.
"""

import argparse
import re
import os
import matplotlib.pyplot as plt

def parse_test_data(log_path, test_name):
    """
    Parse the log file and extract data for a specific test.
    Returns dict with algorithm -> bucket -> highest_label_count
    """
    with open(log_path, 'r') as f:
        content = f.read()
    
    # Find the test block
    pattern = rf'Running test {re.escape(test_name)}(.*?)(?=Running test |$)'
    match = re.search(pattern, content, re.DOTALL)
    
    if not match:
        raise ValueError(f"Test '{test_name}' not found in log file")
    
    block = match.group(1)
    lines = block.split('\n')
    
    # Extract highest label counts per algorithm
    data = {}
    current_algo = None
    
    for line in lines:
        # Detect algorithm
        algo_match = re.search(r'Algorithm (\S+) got:', line)
        if algo_match:
            current_algo = algo_match.group(1)
            data[current_algo] = {'highest_label': {}, 'root_count': {}}
            continue
        
        # Extract highest label count
        hl_match = re.search(r'Bucket (\d+): Average highest label count = ([\d.]+)', line)
        if hl_match and current_algo:
            bucket = int(hl_match.group(1))
            count = float(hl_match.group(2))
            data[current_algo]['highest_label'][bucket] = count
        
        # Extract root count
        rc_match = re.search(r'Bucket (\d+): Average root count = ([\d.]+)', line)
        if rc_match and current_algo:
            bucket = int(rc_match.group(1))
            count = float(rc_match.group(2))
            data[current_algo]['root_count'][bucket] = count
    
    return data

def plot_highest_label_counts(data, test_name, output_path):
    """
    Generate a bar chart of highest label counts per bucket for each algorithm.
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Get algorithms that have highest_label data
    algos_with_data = [algo for algo, d in data.items() if d['highest_label']]
    
    if not algos_with_data:
        raise ValueError("No highest label count data found")
    
    # Get all bucket numbers
    all_buckets = set()
    for algo in algos_with_data:
        all_buckets.update(data[algo]['highest_label'].keys())
    buckets = sorted(all_buckets)
    
    # Set up bar positions
    n_algos = len(algos_with_data)
    bar_width = 0.8 / n_algos
    x = range(len(buckets))
    
    # Plot bars for each algorithm
    colors = plt.cm.tab10.colors
    for i, algo in enumerate(algos_with_data):
        values = [data[algo]['highest_label'].get(b, 0) for b in buckets]
        offset = (i - n_algos / 2 + 0.5) * bar_width
        bars = ax.bar([xi + offset for xi in x], values, bar_width, 
                      label=algo, color=colors[i % len(colors)])
    
    ax.set_xlabel('Bucket')
    ax.set_ylabel('Average Highest Label Count')
    ax.set_title(f'Highest Label Count per Bucket - {test_name}')
    ax.set_xticks(x)
    ax.set_xticklabels(buckets)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    
    print(f'Figure saved to: {output_path}')

def main():
    parser = argparse.ArgumentParser(
        description='Generate figure of highest label counts per bucket for a test'
    )
    parser.add_argument('log_file', help='Path to the log file')
    parser.add_argument('test_name', help='Name of the test to plot')
    parser.add_argument(
        '-o', '--output',
        help='Output file name (default: <test_name>_highest_label.png)',
        default=None
    )
    args = parser.parse_args()
    
    # Determine output path
    output_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(args.log_file))),
        'processed_logs'
    )
    if args.output:
        output_path = os.path.join(output_dir, args.output)
    else:
        safe_name = args.test_name.replace('/', '_').replace('.', '_')
        output_path = os.path.join(output_dir, f'{safe_name}_highest_label.png')
    
    # Parse and plot
    data = parse_test_data(args.log_file, args.test_name)
    plot_highest_label_counts(data, args.test_name, output_path)

if __name__ == '__main__':
    main()

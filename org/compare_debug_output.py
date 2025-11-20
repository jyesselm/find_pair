#!/usr/bin/env python3
"""
Compare debug output from original x3dna with our modernized code.

Usage:
    python compare_debug_output.py <original_debug.txt> <our_debug.txt>
"""

import sys
import re

def parse_debug_output(filename):
    """Parse debug output and extract key values."""
    data = {
        'calculate_more_bppars': [],
        'check_pair': [],
        'parameters': []
    }
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # Extract calculate_more_bppars calls
    for match in re.finditer(r'\[DEBUG calculate_more_bppars\] i=(\d+) j=(\d+) dir_z=([\d\.-]+)', content):
        data['calculate_more_bppars'].append({
            'i': int(match.group(1)),
            'j': int(match.group(2)),
            'dir_z': float(match.group(3))
        })
    
    # Extract orien values
    for match in re.finditer(r'\[DEBUG\] orien\[(\d+)\]\[7-9\]=([\d\.-]+) ([\d\.-]+) ([\d\.-]+)', content):
        idx = int(match.group(1))
        z_axis = [float(match.group(2)), float(match.group(3)), float(match.group(4))]
        if 'orien' not in data:
            data['orien'] = {}
        data['orien'][idx] = z_axis
    
    # Extract parameters
    for match in re.finditer(r'\[DEBUG\] Result: Shift=([\d\.-]+) Slide=([\d\.-]+) Rise=([\d\.-]+) Tilt=([\d\.-]+) Roll=([\d\.-]+) Twist=([\d\.-]+)', content):
        data['parameters'].append({
            'shift': float(match.group(1)),
            'slide': float(match.group(2)),
            'rise': float(match.group(3)),
            'tilt': float(match.group(4)),
            'roll': float(match.group(5)),
            'twist': float(match.group(6))
        })
    
    return data

def compare_outputs(orig_file, our_file):
    """Compare two debug output files."""
    orig_data = parse_debug_output(orig_file)
    our_data = parse_debug_output(our_file)
    
    print("=" * 80)
    print("COMPARISON: Original vs Our Code")
    print("=" * 80)
    
    # Compare first pair parameters
    if orig_data['parameters'] and our_data['parameters']:
        orig_params = orig_data['parameters'][0]
        our_params = our_data['parameters'][0]
        
        print("\nFirst Base Pair Parameters:")
        print("-" * 80)
        params = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
        for param in params:
            orig_val = orig_params.get(param, 0)
            our_val = our_params.get(param, 0)
            diff = abs(orig_val - our_val)
            match = "✓" if diff < 0.001 else "✗"
            print(f"  {param:6s}: Original={orig_val:10.6f}, Our={our_val:10.6f}, Diff={diff:10.6f} {match}")
    
    # Compare dir_z values
    if orig_data['calculate_more_bppars'] and our_data['calculate_more_bppars']:
        orig_call = orig_data['calculate_more_bppars'][0]
        our_call = our_data['calculate_more_bppars'][0]
        
        print("\ndir_z values:")
        print("-" * 80)
        orig_dir_z = orig_call['dir_z']
        our_dir_z = our_call['dir_z']
        diff = abs(orig_dir_z - our_dir_z)
        match = "✓" if diff < 0.001 else "✗"
        print(f"  Original: {orig_dir_z:.6f}")
        print(f"  Our:      {our_dir_z:.6f}")
        print(f"  Diff:     {diff:.6f} {match}")
    
    # Compare orien z-axes
    if 'orien' in orig_data and 'orien' in our_data:
        print("\nOrientation z-axes:")
        print("-" * 80)
        for idx in sorted(set(list(orig_data['orien'].keys()) + list(our_data['orien'].keys()))):
            if idx in orig_data['orien'] and idx in our_data['orien']:
                orig_z = orig_data['orien'][idx]
                our_z = our_data['orien'][idx]
                diff = [abs(orig_z[i] - our_z[i]) for i in range(3)]
                match = "✓" if all(d < 0.001 for d in diff) else "✗"
                print(f"  orien[{idx}][7-9]:")
                print(f"    Original: [{orig_z[0]:8.6f}, {orig_z[1]:8.6f}, {orig_z[2]:8.6f}]")
                print(f"    Our:      [{our_z[0]:8.6f}, {our_z[1]:8.6f}, {our_z[2]:8.6f}]")
                print(f"    Diff:     [{diff[0]:8.6f}, {diff[1]:8.6f}, {diff[2]:8.6f}] {match}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python compare_debug_output.py <original_debug.txt> <our_debug.txt>")
        sys.exit(1)
    
    compare_outputs(sys.argv[1], sys.argv[2])


#!/usr/bin/env python3
"""
Parse legacy debug output from calculate_more_bppars to extract step parameters.

Legacy code outputs debug information like:
[DEBUG] Calling bpstep_par(r2, org[j], r1, org[i], ...)
[DEBUG] Result: Shift=%.6f Slide=%.6f Rise=%.6f Tilt=%.6f Roll=%.6f Twist=%.6f

This script parses that output and compares with modern calculations.
"""

import re
import sys
from typing import Dict, Optional, Tuple

def parse_legacy_debug_line(line: str) -> Optional[Dict]:
    """Parse a legacy debug line to extract step parameters."""
    # Pattern: [DEBUG] Result: Shift=X.XXXXXX Slide=X.XXXXXX Rise=X.XXXXXX Tilt=X.XXXXXX Roll=X.XXXXXX Twist=X.XXXXXX
    pattern = r'\[DEBUG\] Result: Shift=([-\d.]+)\s+Slide=([-\d.]+)\s+Rise=([-\d.]+)\s+Tilt=([-\d.]+)\s+Roll=([-\d.]+)\s+Twist=([-\d.]+)'
    match = re.search(pattern, line)
    if match:
        return {
            'shift': float(match.group(1)),
            'slide': float(match.group(2)),
            'rise': float(match.group(3)),
            'tilt': float(match.group(4)),
            'roll': float(match.group(5)),
            'twist': float(match.group(6))
        }
    return None

def parse_legacy_pair_info(line: str) -> Optional[Tuple[int, int]]:
    """Parse pair indices from legacy debug output."""
    # Pattern: [DEBUG calculate_more_bppars] i=XXX j=YYY dir_z=Z.ZZZZZZ
    pattern = r'\[DEBUG calculate_more_bppars\] i=(\d+)\s+j=(\d+)'
    match = re.search(pattern, line)
    if match:
        return (int(match.group(1)), int(match.group(2)))
    return None

def parse_legacy_debug_file(file_path: str) -> Dict[Tuple[int, int], Dict]:
    """Parse legacy debug output file and extract step parameters for each pair."""
    results = {}
    current_pair = None
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                # Check for pair info
                pair_info = parse_legacy_pair_info(line)
                if pair_info:
                    current_pair = pair_info
                
                # Check for step parameters
                params = parse_legacy_debug_line(line)
                if params and current_pair:
                    results[current_pair] = params
                    current_pair = None  # Reset after finding parameters
    except FileNotFoundError:
        print(f"Error: File not found: {file_path}")
    except Exception as e:
        print(f"Error parsing file: {e}")
    
    return results

def compare_with_modern(legacy_params: Dict, modern_params: Dict, pair: Tuple[int, int]):
    """Compare legacy and modern step parameters."""
    print(f"\n{'='*70}")
    print(f"Comparing pair {pair}")
    print(f"{'='*70}\n")
    
    print(f"{'Parameter':<15} {'Legacy':<15} {'Modern':<15} {'Difference':<15} {'Match':<10}")
    print("-" * 70)
    
    # Parameters used for bp_type_id
    key_params = ['slide', 'rise', 'twist']
    
    for param in ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']:
        leg_val = legacy_params.get(param, 0.0)
        mod_val = modern_params.get(param, 0.0)
        diff = abs(leg_val - mod_val)
        match = "✅" if diff < 0.0001 else "❌"
        
        # Highlight key parameters
        marker = " ⭐" if param in key_params else ""
        print(f"{param.capitalize():<15}{leg_val:>15.6f}{mod_val:>15.6f}{diff:>15.6f}{match:<10}{marker}")
    
    # Check bp_type_id thresholds
    print("\nbp_type_id Threshold Checks:")
    for param_name, param_key, threshold in [('Stretch', 'rise', 2.0), ('Opening', 'twist', 60.0), ('Shear', 'slide', 1.8)]:
        leg_abs = abs(legacy_params.get(param_key, 0.0))
        mod_abs = abs(modern_params.get(param_key, 0.0))
        leg_ok = leg_abs <= threshold
        mod_ok = mod_abs <= threshold
        
        print(f"  {param_name} (fabs({param_key})):")
        print(f"    Legacy: {leg_abs:.6f} {'✅' if leg_ok else '❌'} (threshold: {threshold})")
        print(f"    Modern: {mod_abs:.6f} {'✅' if mod_ok else '❌'} (threshold: {threshold})")
        if leg_ok != mod_ok:
            print(f"    ⚠️  THRESHOLD MISMATCH - This could explain bp_type_id difference!")

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 parse_legacy_debug_output.py <legacy_debug_file> [modern_json_file]")
        print("\nExample:")
        print("  # Parse legacy stderr output")
        print("  python3 parse_legacy_debug_output.py legacy_stderr.log")
        print("\n  # Compare with modern JSON")
        print("  python3 parse_legacy_debug_output.py legacy_stderr.log data/json/pair_validation/6CAQ.json")
        sys.exit(1)
    
    legacy_file = sys.argv[1]
    legacy_results = parse_legacy_debug_file(legacy_file)
    
    if not legacy_results:
        print(f"No step parameters found in {legacy_file}")
        print("\nExpected format:")
        print("  [DEBUG calculate_more_bppars] i=1024 j=1188 dir_z=-0.973039")
        print("  [DEBUG] Result: Shift=-5.397193 Slide=-1.719321 Rise=-0.035735 Tilt=-9.237754 Roll=9.616516 Twist=-50.870333")
        sys.exit(1)
    
    print(f"Found step parameters for {len(legacy_results)} pairs:")
    for pair in legacy_results:
        print(f"  {pair}")
    
    # If modern JSON provided, compare
    if len(sys.argv) >= 3:
        modern_file = sys.argv[2]
        # TODO: Parse modern JSON and compare
        print(f"\nModern comparison not yet implemented for {modern_file}")
        print("Use tools/compare_bp_type_id_calculation for modern values")
    else:
        print("\nTo compare with modern, provide modern JSON file as second argument")
        print("Or use: build/compare_bp_type_id_calculation data/pdb/6CAQ.pdb <idx1> <idx2>")

if __name__ == "__main__":
    main()


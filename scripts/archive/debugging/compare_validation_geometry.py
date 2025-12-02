#!/usr/bin/env python3
"""
Compare validation geometry results side-by-side for a specific pair.

Shows all calculated values (dorg, d_v, plane_angle, dNN, quality_score)
and validation checks for both legacy and modern.

Usage:
    python3 scripts/compare_validation_geometry.py <PDB_ID> <res_i> <res_j>
    python3 scripts/compare_validation_geometry.py 9CF3 25 27
"""

import sys
import json
import argparse
from pathlib import Path
from typing import Optional, Dict

def find_validation_record(json_file: Path, res_i: int, res_j: int) -> Optional[Dict]:
    """Find validation record for a specific pair."""
    if not json_file.exists():
        return None
    
    with open(json_file) as f:
        data = json.load(f)
    
    if not isinstance(data, list):
        return None
    
    for record in data:
        if record.get('type') == 'pair_validation':
            base_i = record.get('base_i')
            base_j = record.get('base_j')
            if (base_i == res_i and base_j == res_j) or (base_i == res_j and base_j == res_i):
                return record
    
    return None

def format_value(val, fmt='.6f'):
    """Format a value, handling None."""
    if val is None:
        return "N/A"
    return f"{val:{fmt}}"

def compare_validation(pdb_id: str, res_i: int, res_j: int):
    """Compare validation results for a pair."""
    legacy_file = Path(f"data/json_legacy/pair_validation/{pdb_id}.json")
    modern_file = Path(f"data/json/pair_validation/{pdb_id}.json")
    
    leg_record = find_validation_record(legacy_file, res_i, res_j)
    mod_record = find_validation_record(modern_file, res_i, res_j)
    
    if not leg_record and not mod_record:
        print(f"ERROR: No validation records found for pair ({res_i}, {res_j})")
        return
    
    print(f"Validation Comparison for Pair ({res_i}, {res_j}) in PDB: {pdb_id}")
    print("=" * 100)
    print()
    
    leg_calc = leg_record.get('calculated_values', {}) if leg_record else {}
    mod_calc = mod_record.get('calculated_values', {}) if mod_record else {}
    
    leg_dir = leg_record.get('direction_vectors', {}) if leg_record else {}
    mod_dir = mod_record.get('direction_vectors', {}) if mod_record else {}
    
    # Validation status
    print("Validation Status:")
    print(f"  Legacy: is_valid={leg_record.get('is_valid') if leg_record else 'N/A'}, "
          f"bp_type_id={leg_record.get('bp_type_id') if leg_record else 'N/A'}")
    print(f"  Modern: is_valid={mod_record.get('is_valid') if mod_record else 'N/A'}, "
          f"bp_type_id={mod_record.get('bp_type_id') if mod_record else 'N/A'}")
    print()
    
    # Geometry values
    print("Geometry Values:")
    print(f"{'Parameter':<15} {'Legacy':<15} {'Modern':<15} {'Difference':<15} {'Status'}")
    print("-" * 100)
    
    params = [
        ('dorg', 'Distance between origins (Å)'),
        ('d_v', 'Vertical distance (Å)'),
        ('plane_angle', 'Plane angle (degrees)'),
        ('dNN', 'N-N distance (Å)'),
        ('quality_score', 'Quality score'),
    ]
    
    for param, desc in params:
        leg_val = leg_calc.get(param)
        mod_val = mod_calc.get(param)
        
        if leg_val is not None and mod_val is not None:
            diff = abs(mod_val - leg_val)
            status = "✅" if diff < 0.01 else "❌"
        else:
            diff = None
            status = "?"
        
        print(f"{desc:<15} {format_value(leg_val, '>15.6f'):<15} "
              f"{format_value(mod_val, '>15.6f'):<15} {format_value(diff, '>15.6f'):<15} {status}")
    
    print()
    
    # Direction vectors
    print("Direction Vectors:")
    if leg_dir and mod_dir:
        for axis in ['dir_x', 'dir_y', 'dir_z']:
            leg_val = leg_dir.get(axis)
            mod_val = mod_dir.get(axis)
            if leg_val is not None and mod_val is not None:
                diff = abs(mod_val - leg_val)
                status = "✅" if diff < 0.01 else "❌"
            else:
                diff = None
                status = "?"
            print(f"  {axis:7s}: Legacy={format_value(leg_val, '>10.6f'):<12} "
                  f"Modern={format_value(mod_val, '>10.6f'):<12} "
                  f"Diff={format_value(diff, '>10.6f'):<12} {status}")
    print()
    
    # Validation thresholds
    print("Validation Thresholds (from code defaults):")
    print("  dorg:     0.0 - 15.0 Å")
    print("  d_v:      0.0 - 2.5 Å")
    print("  plane_angle: 0.0 - 65.0 degrees")
    print("  dNN:      4.5 - 1e18 Å")
    print()
    
    # Check which validation tests fail
    print("Validation Checks:")
    if leg_record:
        leg_dorg = leg_calc.get('dorg')
        leg_dv = leg_calc.get('d_v')
        leg_plane = leg_calc.get('plane_angle')
        leg_dnn = leg_calc.get('dNN')
        
        leg_passes = []
        leg_fails = []
        
        if leg_dorg is not None:
            if 0 <= leg_dorg <= 15.0:
                leg_passes.append('dorg')
            else:
                leg_fails.append(f'dorg={leg_dorg:.2f}')
        
        if leg_dv is not None:
            if 0 <= leg_dv <= 2.5:
                leg_passes.append('d_v')
            else:
                leg_fails.append(f'd_v={leg_dv:.2f}')
        
        if leg_plane is not None:
            if 0 <= leg_plane <= 65.0:
                leg_passes.append('plane_angle')
            else:
                leg_fails.append(f'plane_angle={leg_plane:.2f}')
        
        if leg_dnn is not None:
            if 4.5 <= leg_dnn <= 1e18:
                leg_passes.append('dNN')
            else:
                leg_fails.append(f'dNN={leg_dnn:.2f}')
        
        print(f"  Legacy: Passes=[{', '.join(leg_passes)}], Fails=[{', '.join(leg_fails)}]")
    
    if mod_record:
        mod_dorg = mod_calc.get('dorg')
        mod_dv = mod_calc.get('d_v')
        mod_plane = mod_calc.get('plane_angle')
        mod_dnn = mod_calc.get('dNN')
        
        mod_passes = []
        mod_fails = []
        
        if mod_dorg is not None:
            if 0 <= mod_dorg <= 15.0:
                mod_passes.append('dorg')
            else:
                mod_fails.append(f'dorg={mod_dorg:.2f}')
        
        if mod_dv is not None:
            if 0 <= mod_dv <= 2.5:
                mod_passes.append('d_v')
            else:
                mod_fails.append(f'd_v={mod_dv:.2f}')
        
        if mod_plane is not None:
            if 0 <= mod_plane <= 65.0:
                mod_passes.append('plane_angle')
            else:
                mod_fails.append(f'plane_angle={mod_plane:.2f}')
        
        if mod_dnn is not None:
            if 4.5 <= mod_dnn <= 1e18:
                mod_passes.append('dNN')
            else:
                mod_fails.append(f'dNN={mod_dnn:.2f}')
        
        print(f"  Modern: Passes=[{', '.join(mod_passes)}], Fails=[{', '.join(mod_fails)}]")

def main():
    parser = argparse.ArgumentParser(
        description='Compare validation geometry between legacy and modern')
    parser.add_argument('pdb_id', help='PDB ID (e.g., 9CF3)')
    parser.add_argument('res_i', type=int, help='First residue index')
    parser.add_argument('res_j', type=int, help='Second residue index')
    
    args = parser.parse_args()
    
    compare_validation(args.pdb_id, args.res_i, args.res_j)

if __name__ == '__main__':
    main()


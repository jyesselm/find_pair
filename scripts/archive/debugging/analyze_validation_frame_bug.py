#!/usr/bin/env python3
"""
Analyze validation frame bug - compares frame_calc frames vs frames used in validation.

This tool helps identify when validation uses wrong frames by comparing:
1. Frame origins from frame_calc JSON (calculated frames)
2. Frame origins implied by validation dorg values

Usage:
    python3 scripts/analyze_validation_frame_bug.py <PDB_ID> <res_i> <res_j>
    python3 scripts/analyze_validation_frame_bug.py 9CF3 25 27
"""

import sys
import json
import math
import argparse
from pathlib import Path
from typing import Dict, Tuple, Optional

def load_frame_origins(json_file: Path) -> Dict[int, Tuple[float, float, float]]:
    """Load frame origins from frame_calc JSON."""
    origins = {}
    
    if not json_file.exists():
        return origins
    
    with open(json_file) as f:
        data = json.load(f)
    
    if not isinstance(data, list):
        return origins
    
    for record in data:
        if isinstance(record, dict):
            idx = record.get('residue_idx') or record.get('legacy_residue_idx')
            translation = record.get('translation', [])
            
            if idx is not None and translation and len(translation) >= 3:
                origins[idx] = tuple(translation[:3])
    
    return origins

def find_validation_dorg(legacy_json: Path, modern_json: Path, 
                        res_i: int, res_j: int) -> Tuple[Optional[float], Optional[float]]:
    """Find dorg values from validation JSON."""
    leg_dorg = None
    mod_dorg = None
    
    # Legacy
    if legacy_json.exists():
        with open(legacy_json) as f:
            leg_data = json.load(f)
        
        if isinstance(leg_data, list):
            for record in leg_data:
                if record.get('type') == 'pair_validation':
                    base_i = record.get('base_i')
                    base_j = record.get('base_j')
                    if (base_i == res_i and base_j == res_j) or (base_i == res_j and base_j == res_i):
                        calc_vals = record.get('calculated_values', {})
                        leg_dorg = calc_vals.get('dorg')
                        break
    
    # Modern
    if modern_json.exists():
        with open(modern_json) as f:
            mod_data = json.load(f)
        
        if isinstance(mod_data, list):
            for record in mod_data:
                if record.get('type') == 'pair_validation':
                    base_i = record.get('base_i')
                    base_j = record.get('base_j')
                    if (base_i == res_i and base_j == res_j) or (base_i == res_j and base_j == res_i):
                        calc_vals = record.get('calculated_values', {})
                        mod_dorg = calc_vals.get('dorg')
                        break
    
    return leg_dorg, mod_dorg

def calculate_dorg_from_origins(origins: Dict[int, Tuple[float, float, float]],
                                res_i: int, res_j: int) -> Optional[float]:
    """Calculate dorg from frame origins."""
    if res_i not in origins or res_j not in origins:
        return None
    
    origin_i = origins[res_i]
    origin_j = origins[res_j]
    
    return math.sqrt(sum((a - b)**2 for a, b in zip(origin_i, origin_j)))

def analyze_pair(pdb_id: str, res_i: int, res_j: int):
    """Analyze a specific pair to detect frame retrieval bug."""
    print(f"Analyzing pair ({res_i}, {res_j}) for PDB: {pdb_id}")
    print("=" * 80)
    print()
    
    # Load frame origins
    legacy_frames_file = Path(f"data/json_legacy/frame_calc/{pdb_id}.json")
    modern_frames_file = Path(f"data/json/frame_calc/{pdb_id}.json")
    
    legacy_origins = load_frame_origins(legacy_frames_file)
    modern_origins = load_frame_origins(modern_frames_file)
    
    # Load validation dorg
    legacy_val_file = Path(f"data/json_legacy/pair_validation/{pdb_id}.json")
    modern_val_file = Path(f"data/json/pair_validation/{pdb_id}.json")
    
    leg_val_dorg, mod_val_dorg = find_validation_dorg(legacy_val_file, modern_val_file, res_i, res_j)
    
    # Calculate expected dorg from frame origins
    leg_expected_dorg = calculate_dorg_from_origins(legacy_origins, res_i, res_j)
    mod_expected_dorg = calculate_dorg_from_origins(modern_origins, res_i, res_j)
    
    print("Frame Origins:")
    if res_i in legacy_origins and res_i in modern_origins:
        leg_orig_i = legacy_origins[res_i]
        mod_orig_i = modern_origins[res_i]
        print(f"  Residue {res_i}:")
        print(f"    Legacy: ({leg_orig_i[0]:10.6f}, {leg_orig_i[1]:10.6f}, {leg_orig_i[2]:10.6f})")
        print(f"    Modern: ({mod_orig_i[0]:10.6f}, {mod_orig_i[1]:10.6f}, {mod_orig_i[2]:10.6f})")
        orig_dist_i = math.sqrt(sum((a - b)**2 for a, b in zip(leg_orig_i, mod_orig_i)))
        print(f"    Distance: {orig_dist_i:.6f} Å")
    
    if res_j in legacy_origins and res_j in modern_origins:
        leg_orig_j = legacy_origins[res_j]
        mod_orig_j = modern_origins[res_j]
        print(f"  Residue {res_j}:")
        print(f"    Legacy: ({leg_orig_j[0]:10.6f}, {leg_orig_j[1]:10.6f}, {leg_orig_j[2]:10.6f})")
        print(f"    Modern: ({mod_orig_j[0]:10.6f}, {mod_orig_j[1]:10.6f}, {mod_orig_j[2]:10.6f})")
        orig_dist_j = math.sqrt(sum((a - b)**2 for a, b in zip(leg_orig_j, mod_orig_j)))
        print(f"    Distance: {orig_dist_j:.6f} Å")
    
    print()
    print("dorg Comparison:")
    print(f"  Legacy validation dorg: {leg_val_dorg:.6f}" if leg_val_dorg else "  Legacy validation dorg: Not found")
    print(f"  Modern validation dorg: {mod_val_dorg:.6f}" if mod_val_dorg else "  Modern validation dorg: Not found")
    print()
    print("Expected dorg (from frame origins):")
    print(f"  Legacy expected dorg: {leg_expected_dorg:.6f}" if leg_expected_dorg else "  Legacy expected dorg: Not found")
    print(f"  Modern expected dorg: {mod_expected_dorg:.6f}" if mod_expected_dorg else "  Modern expected dorg: Not found")
    print()
    
    # Check for frame retrieval bug
    if mod_val_dorg and mod_expected_dorg:
        diff = abs(mod_val_dorg - mod_expected_dorg)
        if diff > 0.01:
            print("🚨 FRAME RETRIEVAL BUG DETECTED!")
            print(f"   Validation uses dorg={mod_val_dorg:.6f}")
            print(f"   But frame origins yield dorg={mod_expected_dorg:.6f}")
            print(f"   Difference: {diff:.6f} Å")
            print()
            print("   This suggests validation is using wrong frames!")
        else:
            print("✅ Validation dorg matches expected (frames are correct)")
    
    if leg_val_dorg and leg_expected_dorg:
        diff = abs(leg_val_dorg - leg_expected_dorg)
        if diff > 0.01:
            print(f"⚠️  Legacy validation dorg differs from expected by {diff:.6f} Å")
        else:
            print("✅ Legacy validation dorg matches expected")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze validation frame bug')
    parser.add_argument('pdb_id', help='PDB ID (e.g., 9CF3)')
    parser.add_argument('res_i', type=int, help='First residue index')
    parser.add_argument('res_j', type=int, help='Second residue index')
    
    args = parser.parse_args()
    
    analyze_pair(args.pdb_id, args.res_i, args.res_j)

if __name__ == '__main__':
    main()


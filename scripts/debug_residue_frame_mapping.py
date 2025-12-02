#!/usr/bin/env python3
"""
Debug residue-to-frame mapping to identify frame retrieval bugs.

This tool checks if residues retrieved by legacy index have the correct frames
set on them, helping identify when frames are lost or wrong residues are used.

Usage:
    python3 scripts/debug_residue_frame_mapping.py <PDB_ID> <res_i> [<res_j>]
    python3 scripts/debug_residue_frame_mapping.py 9CF3 25 27
"""

import sys
import json
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

def find_validation_dorg(json_file: Path, res_i: int, res_j: int) -> Optional[float]:
    """Find dorg from validation JSON."""
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
                calc_vals = record.get('calculated_values', {})
                return calc_vals.get('dorg')
    
    return None

def calculate_dorg_from_origins(origins: Dict[int, Tuple[float, float, float]],
                                res_i: int, res_j: int) -> Optional[float]:
    """Calculate dorg from frame origins."""
    if res_i not in origins or res_j not in origins:
        return None
    
    import math
    origin_i = origins[res_i]
    origin_j = origins[res_j]
    
    return math.sqrt(sum((a - b)**2 for a, b in zip(origin_i, origin_j)))

def analyze_mapping(pdb_id: str, res_i: int, res_j: Optional[int] = None):
    """Analyze residue-frame mapping for debugging."""
    print(f"Residue-Frame Mapping Analysis for PDB: {pdb_id}")
    print("=" * 80)
    print()
    
    # Load frame origins from frame_calc (the "correct" frames)
    modern_frames_file = Path(f"data/json/frame_calc/{pdb_id}.json")
    modern_origins = load_frame_origins(modern_frames_file)
    
    # Load validation dorg (what validation actually calculates)
    modern_val_file = Path(f"data/json/pair_validation/{pdb_id}.json")
    
    print("Frame Origins (from frame_calc JSON):")
    if res_i in modern_origins:
        orig_i = modern_origins[res_i]
        print(f"  Residue {res_i}: ({orig_i[0]:.6f}, {orig_i[1]:.6f}, {orig_i[2]:.6f})")
    else:
        print(f"  Residue {res_i}: NOT FOUND in frame_calc")
    
    if res_j and res_j in modern_origins:
        orig_j = modern_origins[res_j]
        print(f"  Residue {res_j}: ({orig_j[0]:.6f}, {orig_j[1]:.6f}, {orig_j[2]:.6f})")
    elif res_j:
        print(f"  Residue {res_j}: NOT FOUND in frame_calc")
    
    print()
    
    # Calculate expected dorg
    if res_j:
        expected_dorg = calculate_dorg_from_origins(modern_origins, res_i, res_j)
        validation_dorg = find_validation_dorg(modern_val_file, res_i, res_j)
        
        print("dorg Comparison:")
        print(f"  Expected (from frame origins): {expected_dorg:.6f}" if expected_dorg else "  Expected: Not calculable")
        print(f"  Validation (from pair_validation): {validation_dorg:.6f}" if validation_dorg else "  Validation: Not found")
        
        if expected_dorg and validation_dorg:
            diff = abs(validation_dorg - expected_dorg)
            print(f"  Difference: {diff:.6f} Å")
            if diff > 0.01:
                print()
                print("🚨 BUG DETECTED: Validation uses wrong frames!")
                print()
                print("Possible causes:")
                print("  1. Wrong residue objects retrieved by legacy index")
                print("  2. Frames not set correctly on residue objects")
                print("  3. Frames overwritten between calculation and validation")
                print("  4. Residue index mapping issue (wrong residues)")
                print()
                print("Next steps:")
                print("  - Check if residue-by-legacy_idx mapping is correct")
                print("  - Verify frames are set after calculate_all_frames")
                print("  - Check if Structure object is passed correctly")
                print("  - Verify no frame recalculation happens between phases")
            else:
                print("  ✅ Frames match correctly!")
    
    print()
    print("How base_pair_finder gets residues:")
    print("  1. Iterate through structure.chains() and chain.residues()")
    print("  2. Get legacy_residue_idx from residue.atoms()[0].legacy_residue_idx()")
    print("  3. Store in residue_by_legacy_idx[legacy_idx] = &residue")
    print("  4. Retrieve via residue_by_legacy_idx.find(legacy_idx)->second")
    print("  5. Validate using validator_.validate(*res1, *res2)")
    print("  6. Validate gets frames via res1.reference_frame() and res2.reference_frame()")
    print()
    print("Potential issues:")
    print("  - If residue pointers point to wrong residue objects")
    print("  - If frames aren't set on the correct residue objects")
    print("  - If Structure object is modified/copied between phases")

def main():
    parser = argparse.ArgumentParser(
        description='Debug residue-to-frame mapping')
    parser.add_argument('pdb_id', help='PDB ID (e.g., 9CF3)')
    parser.add_argument('res_i', type=int, help='First residue index')
    parser.add_argument('res_j', type=int, nargs='?', help='Second residue index (optional)')
    
    args = parser.parse_args()
    
    analyze_mapping(args.pdb_id, args.res_i, args.res_j)

if __name__ == '__main__':
    main()


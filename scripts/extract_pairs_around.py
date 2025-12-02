#!/usr/bin/env python3
"""
Extract minimal PDB around specific residue pairs.
"""

import sys
import json
from pathlib import Path
from extract_minimal_pairs import extract_minimal_pdb

if __name__ == '__main__':
    pdb_file = "data/pdb/1TTT.pdb"
    legacy_json = "data/json_legacy/find_bestpair_selection/1TTT.json"
    output_dir = "data/pdb/minimal"
    
    # Find pairs around residue 162 or 177
    with open(legacy_json) as f:
        data = json.load(f)
    
    pairs = data['pairs']
    
    # Find pairs that include 162 or 177
    target_pairs = []
    for i, (r1, r2) in enumerate(pairs):
        if r1 == 162 or r2 == 162 or r1 == 177 or r2 == 177:
            target_pairs.append((i, (r1, r2)))
            if len(target_pairs) >= 5:
                break
    
    print(f"Found {len(target_pairs)} pairs involving residues 162 or 177:")
    for idx, (r1, r2) in target_pairs:
        print(f"  Index {idx}: ({r1}, {r2})")
    
    # Try to find a consecutive pair that includes one of these
    print("\nLooking for consecutive pairs...")
    for i in range(len(pairs) - 1):
        r1a, r2a = pairs[i]
        r1b, r2b = pairs[i+1]
        # Check if this pair involves 162 or 177
        if 162 in (r1a, r2a, r1b, r2b) or 177 in (r1a, r2a, r1b, r2b):
            print(f"  Found consecutive pair at index {i}: ({r1a}, {r2a}), ({r1b}, {r2b})")
            # Extract this pair
            try:
                extract_minimal_pdb(
                    pdb_file=pdb_file,
                    legacy_json_file=legacy_json,
                    output_dir=output_dir,
                    pair_indices=[i, i+1],
                    fragment_name=f"1TTT_around_162_177"
                )
                print(f"  ✅ Extracted to {output_dir}/1TTT_around_162_177.pdb")
            except Exception as e:
                print(f"  ❌ Error: {e}")
            break


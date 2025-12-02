#!/usr/bin/env python3
"""Extract minimal case around a specific residue."""

import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from extract_minimal_pairs import *

if __name__ == '__main__':
    target_residue = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    
    pdb_file = Path("data/pdb/1TTT.pdb")
    legacy_json = Path("data/json_legacy/find_bestpair_selection/1TTT.json")
    output_dir = Path("data/pdb/minimal")
    json_dir = legacy_json.parent.parent
    
    # Load pairs
    pairs = load_legacy_find_bestpair_json(legacy_json)
    
    # Find pair closest to target residue
    print(f"Looking for pairs involving residue {target_residue}...")
    closest_idx = None
    min_dist = float('inf')
    
    for i, (r1, r2) in enumerate(pairs):
        dist = min(abs(r1 - target_residue), abs(r2 - target_residue))
        if dist < min_dist:
            min_dist = dist
            closest_idx = i
    
    if closest_idx is None:
        print("No pairs found")
        sys.exit(1)
    
    print(f"Found closest pair at index {closest_idx}: {pairs[closest_idx]} (distance: {min_dist})")
    
    # Extract this pair and one neighbor
    if closest_idx > 0:
        target_indices = [closest_idx - 1, closest_idx]
    elif closest_idx < len(pairs) - 1:
        target_indices = [closest_idx, closest_idx + 1]
    else:
        target_indices = [closest_idx]
    
    target_pairs = [pairs[i] for i in target_indices]
    print(f"Extracting pairs at indices {target_indices}: {target_pairs}")
    
    # Load residue mapping
    pdb_id = pdb_file.stem.upper()
    residue_map = load_legacy_residue_mapping(json_dir, pdb_id)
    
    # Get residue range
    result = get_residue_range_for_pairs(target_pairs, residue_map)
    if result is None:
        print("ERROR: Could not map residues")
        sys.exit(1)
    
    residue_set, legacy_indices = result
    
    # Extract fragment
    output_file = output_dir / f"1TTT_res_{target_residue}.pdb"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if extract_pdb_fragment(pdb_file, residue_set, output_file, target_pairs):
        print(f"✅ Successfully extracted to {output_file}")


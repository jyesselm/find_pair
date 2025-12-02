#!/usr/bin/env python3
"""Extract specific pair indices from legacy JSON."""

import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from extract_minimal_pairs import *

if __name__ == '__main__':
    pdb_file = Path("data/pdb/1TTT.pdb")
    legacy_json = Path("data/json_legacy/find_bestpair_selection/1TTT.json")
    output_dir = Path("data/pdb/minimal")
    json_dir = legacy_json.parent.parent
    
    # Load pairs
    pairs = load_legacy_find_bestpair_json(legacy_json)
    
    # Extract pairs at indices 68-69: (160, 166), (162, 177)
    target_indices = [68, 69]
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
    output_file = output_dir / "1TTT_162_177.pdb"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if extract_pdb_fragment(pdb_file, residue_set, output_file, target_pairs):
        print(f"✅ Successfully extracted to {output_file}")


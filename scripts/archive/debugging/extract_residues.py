#!/usr/bin/env python3
"""Extract PDB fragment containing specific residues."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from extract_minimal_pairs import *

if __name__ == '__main__':
    # Residues to include: 16, 59, and nearby paired residues
    target_residues = [13, 15, 16, 18, 19, 20, 22, 25, 48, 55, 56, 57, 59]
    
    pdb_file = Path("data/pdb/1TTT.pdb")
    legacy_json = Path("data/json_legacy/find_bestpair_selection/1TTT.json")
    output_dir = Path("data/pdb/minimal")
    json_dir = legacy_json.parent.parent
    
    # Load residue mapping
    pdb_id = pdb_file.stem.upper()
    residue_map = load_legacy_residue_mapping(json_dir, pdb_id)
    
    # Get residue set for these residues
    residue_set = set()
    legacy_indices = set()
    for res_idx in target_residues:
        if res_idx in residue_map:
            residue_set.add(residue_map[res_idx])
            legacy_indices.add(res_idx)
    
    print(f"Extracting {len(residue_set)} residues: {sorted(target_residues)}")
    
    # Extract fragment
    output_file = output_dir / "1TTT_16_59.pdb"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create dummy pair_info for the remark
    pair_info = [(16, 59)]
    
    if extract_pdb_fragment(pdb_file, residue_set, output_file, pair_info):
        print(f"✅ Successfully extracted to {output_file}")


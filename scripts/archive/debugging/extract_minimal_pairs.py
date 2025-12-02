#!/usr/bin/env python3
"""
Extract minimal PDB fragments containing exactly N consecutive base pairs.

This script reads legacy JSON output to identify consecutive base pairs,
then extracts the corresponding residues from the PDB file to create
minimal test cases for easier debugging.

Usage:
    python3 scripts/extract_minimal_pairs.py \
        --pdb data/pdb/1AQ4.pdb \
        --legacy-json data/json_legacy/find_bestpair_selection/1AQ4.json \
        --output-dir data/pdb/minimal \
        --num-pairs 2
"""

import sys
import json
import argparse
from pathlib import Path
from typing import List, Tuple, Dict, Set, Optional

def load_legacy_find_bestpair_json(json_file: Path) -> List[Tuple[int, int]]:
    """Load legacy find_bestpair_selection JSON and return list of pairs."""
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Find the find_bestpair_selection record
    for record in data:
        if isinstance(record, dict) and record.get("type") == "find_bestpair_selection":
            pairs = record.get("pairs", [])
            return [(pair[0], pair[1]) for pair in pairs if len(pair) >= 2]
    
    return []

def load_legacy_residue_mapping(json_dir: Path, pdb_id: str) -> Dict[int, Tuple[str, str, int, str]]:
    """
    Load legacy residue index to PDB properties mapping from base_frame_calc JSON.
    
    Returns: Dict mapping legacy_residue_idx -> (residue_name, chain_id, residue_seq, insertion)
    """
    # Try to find base_frame_calc JSON
    json_file = json_dir / "base_frame_calc" / f"{pdb_id}.json"
    if not json_file.exists():
        # Try main JSON file
        json_file = json_dir / f"{pdb_id}_globals.json"
    
    if not json_file.exists():
        print(f"WARNING: Could not find legacy JSON file for {pdb_id}")
        return {}
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    residue_map = {}
    
    # Look for base_frame_calc records
    for record in data:
        if isinstance(record, dict):
            # Check if this is a base_frame_calc record
            is_base_frame = False
            if record.get("type") == "base_frame_calc":
                is_base_frame = True
            elif "residue_idx" in record and "residue_name" in record:
                is_base_frame = True
            
            if is_base_frame:
                residue_idx = record.get("residue_idx", 0)
                residue_name = record.get("residue_name", "")
                if not residue_name and "base_type" in record:
                    # Convert base_type to residue_name
                    base_type = record.get("base_type", "")
                    base_to_resname = {
                        "A": "  A", "C": "  C", "G": "  G", 
                        "U": "  U", "T": "  T"
                    }
                    residue_name = base_to_resname.get(base_type, "")
                
                chain_str = record.get("chain_id", "")
                chain_id = chain_str[0] if chain_str else ' '
                
                residue_seq = record.get("residue_seq", 0)
                
                ins_str = record.get("insertion", "")
                insertion = ins_str[0] if ins_str else ' '
                
                if residue_idx > 0 and residue_name:
                    residue_map[residue_idx] = (residue_name.strip(), chain_id, residue_seq, insertion)
    
    return residue_map

def find_consecutive_pairs(pairs: List[Tuple[int, int]], num_pairs: int = 2) -> List[List[Tuple[int, int]]]:
    """
    Find consecutive base pairs in the pair list.
    
    Two pairs are consecutive if they appear next to each other in the pair list.
    Returns list of groups, each containing num_pairs consecutive pairs.
    """
    if len(pairs) < num_pairs:
        return []
    
    groups = []
    for i in range(len(pairs) - num_pairs + 1):
        group = pairs[i:i+num_pairs]
        groups.append(group)
    
    return groups

def get_residue_range_for_pairs(pairs: List[Tuple[int, int]], 
                                residue_map: Dict[int, Tuple[str, str, int, str]]) -> Optional[Tuple[Set[Tuple[str, str, int, str]], Set[int]]]:
    """
    Get all residues (by PDB properties) needed for a set of base pairs.
    
    Returns: (set of (residue_name, chain_id, residue_seq, insertion), set of legacy_indices)
    """
    residues = set()
    legacy_indices = set()
    
    for res1_idx, res2_idx in pairs:
        legacy_indices.add(res1_idx)
        legacy_indices.add(res2_idx)
        
        if res1_idx in residue_map:
            residues.add(residue_map[res1_idx])
        if res2_idx in residue_map:
            residues.add(residue_map[res2_idx])
    
    if not residues:
        return None
    
    return (residues, legacy_indices)

def extract_pdb_fragment(pdb_file: Path, 
                        residue_set: Set[Tuple[str, str, int, str]],
                        output_file: Path,
                        pair_info: List[Tuple[int, int]]) -> bool:
    """
    Extract PDB fragment containing only the specified residues.
    
    Args:
        pdb_file: Input PDB file
        residue_set: Set of (residue_name, chain_id, residue_seq, insertion)
        output_file: Output PDB file path
        pair_info: List of base pairs (for remarks)
    """
    atoms = []
    
    # Build lookup set for faster checking
    # Normalize residue names (strip spaces)
    residue_lookup = set()
    for res_name, chain_id, res_seq, insertion in residue_set:
        normalized_name = res_name.strip()
        residue_lookup.add((normalized_name, chain_id, res_seq, insertion))
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Parse PDB line
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21] if len(line) > 21 else ' '
                res_seq_str = line[22:26].strip()
                insertion = line[26] if len(line) > 26 else ' '
                
                try:
                    res_seq = int(res_seq_str)
                except ValueError:
                    continue
                
                # Check if this residue is in our set
                if (res_name, chain_id, res_seq, insertion) in residue_lookup:
                    atoms.append(line.rstrip())
    
    if not atoms:
        print(f"ERROR: No atoms found for specified residues")
        return False
    
    # Write output PDB
    with open(output_file, 'w') as f:
        f.write("REMARK   Minimal test case extracted from " + str(pdb_file.name) + "\n")
        f.write(f"REMARK   Contains {len(pair_info)} base pair(s): {pair_info}\n")
        f.write("REMARK   Residues extracted by legacy indices\n")
        for atom in atoms:
            f.write(atom + "\n")
        f.write("END\n")
    
    print(f"Created minimal PDB: {output_file}")
    print(f"  Pairs: {pair_info}")
    print(f"  Atoms: {len(atoms)}")
    return True

def main():
    parser = argparse.ArgumentParser(
        description='Extract minimal PDB fragments with N consecutive base pairs')
    parser.add_argument('--pdb', type=Path, required=True,
                       help='Input PDB file')
    parser.add_argument('--legacy-json', type=Path, required=True,
                       help='Legacy find_bestpair_selection JSON file')
    parser.add_argument('--legacy-json-dir', type=Path,
                       help='Directory containing legacy JSON files (for residue mapping)')
    parser.add_argument('--output-dir', type=Path, required=True,
                       help='Output directory for minimal PDB files')
    parser.add_argument('--num-pairs', type=int, default=2,
                       help='Number of consecutive pairs to extract (default: 2)')
    parser.add_argument('--max-fragments', type=int,
                       help='Maximum number of fragments to create (default: all)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.pdb.exists():
        print(f"ERROR: PDB file not found: {args.pdb}")
        sys.exit(1)
    
    if not args.legacy_json.exists():
        print(f"ERROR: Legacy JSON file not found: {args.legacy_json}")
        sys.exit(1)
    
    # Determine legacy JSON directory
    json_dir = args.legacy_json_dir
    if json_dir is None:
        # Try to infer from legacy_json path
        json_dir = args.legacy_json.parent.parent
        if not json_dir.exists():
            json_dir = args.legacy_json.parent
    
    # Get PDB ID
    pdb_id = args.pdb.stem.upper()
    
    # Load base pairs from legacy JSON
    print(f"Loading base pairs from: {args.legacy_json}")
    pairs = load_legacy_find_bestpair_json(args.legacy_json)
    if not pairs:
        print(f"ERROR: No base pairs found in {args.legacy_json}")
        sys.exit(1)
    
    print(f"Found {len(pairs)} base pairs")
    
    # Load residue mapping
    print(f"Loading residue mapping from: {json_dir}")
    residue_map = load_legacy_residue_mapping(json_dir, pdb_id)
    if not residue_map:
        print(f"WARNING: Could not load residue mapping. Extraction may fail.")
    else:
        print(f"Loaded mapping for {len(residue_map)} residues")
    
    # Find consecutive pair groups
    print(f"\nFinding groups of {args.num_pairs} consecutive pairs...")
    pair_groups = find_consecutive_pairs(pairs, args.num_pairs)
    if not pair_groups:
        print(f"ERROR: No groups of {args.num_pairs} consecutive pairs found")
        sys.exit(1)
    
    print(f"Found {len(pair_groups)} group(s) of {args.num_pairs} consecutive pairs")
    
    # Limit number of fragments if specified
    if args.max_fragments:
        pair_groups = pair_groups[:args.max_fragments]
        print(f"Limited to {len(pair_groups)} fragment(s)")
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract fragments
    print(f"\nExtracting fragments to: {args.output_dir}")
    success_count = 0
    
    for i, group in enumerate(pair_groups, 1):
        # Get residue range for this group
        result = get_residue_range_for_pairs(group, residue_map)
        if result is None:
            print(f"WARNING: Could not map residues for group {i}: {group}")
            continue
        
        residue_set, legacy_indices = result
        
        # Create output filename
        output_file = args.output_dir / f"{pdb_id}_minimal_pairs_{i}_{i+args.num_pairs-1}.pdb"
        
        # Extract fragment
        if extract_pdb_fragment(args.pdb, residue_set, output_file, group):
            success_count += 1
            print(f"  Group {i}: pairs {group} -> {output_file.name}")
    
    print(f"\nSuccessfully created {success_count} minimal PDB fragment(s)")
    
    if success_count > 0:
        print(f"\nNext steps:")
        print(f"1. Run legacy find_pair on each minimal PDB")
        print(f"2. Run modern find_pair on each minimal PDB")
        print(f"3. Compare results using: python3 scripts/compare_json.py compare ...")

if __name__ == '__main__':
    main()


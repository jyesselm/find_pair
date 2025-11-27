#!/usr/bin/env python3
"""
Compare step parameters for specific pairs between legacy and modern implementations.
This script helps debug why legacy assigns bp_type_id=-1 while modern assigns bp_type_id=2.
"""

import json
import sys
from pathlib import Path
from typing import Optional, Dict, Tuple

def load_json_records(pdb_id: str, record_type: str, legacy: bool = False) -> list:
    """Load JSON records of a specific type."""
    base_dir = Path("data/json_legacy" if legacy else "data/json")
    
    # Try split file format first
    split_file = base_dir / record_type / f"{pdb_id}.json"
    if split_file.exists():
        with open(split_file) as f:
            data = json.load(f)
            if isinstance(data, list):
                return data
            elif isinstance(data, dict) and "calculations" in data:
                return data["calculations"]
    
    # Try main file
    main_file = base_dir / f"{pdb_id}.json"
    if main_file.exists():
        with open(main_file) as f:
            data = json.load(f)
            if isinstance(data, dict) and "calculations" in data:
                calcs = data["calculations"]
                if isinstance(calcs, list):
                    return [c for c in calcs if isinstance(c, dict) and c.get("type") == record_type]
                elif isinstance(calcs, dict):
                    return calcs.get(record_type, [])
    
    return []


def find_pair_validation(pdb_id: str, idx1: int, idx2: int, legacy: bool = False) -> Optional[Dict]:
    """Find pair validation record for a specific pair."""
    records = load_json_records(pdb_id, "pair_validation", legacy)
    
    for record in records:
        # Handle both 0-based and 1-based indexing
        r1 = record.get("residue_idx1")
        r2 = record.get("residue_idx2")
        
        # Try 0-based match
        if r1 == idx1 - 1 and r2 == idx2 - 1:
            return record
        # Try 1-based match
        if r1 == idx1 and r2 == idx2:
            return record
    
    return None


def find_base_pair_record(pdb_id: str, idx1: int, idx2: int, legacy: bool = False) -> Optional[Dict]:
    """Find base_pair record for a specific pair."""
    records = load_json_records(pdb_id, "base_pair", legacy)
    
    for record in records:
        # Handle both 0-based and 1-based indexing
        base_i = record.get("base_i")
        base_j = record.get("base_j")
        
        # Try 1-based match (legacy uses 1-based)
        if base_i == idx1 and base_j == idx2:
            return record
        # Try 0-based match
        if base_i == idx1 - 1 and base_j == idx2 - 1:
            return record
    
    return None


def compare_pair_step_params(pdb_id: str, idx1: int, idx2: int):
    """Compare step parameters and bp_type_id for a specific pair."""
    print(f"\n{'='*70}")
    print(f"Comparing pair ({idx1}, {idx2}) in {pdb_id}")
    print(f"{'='*70}\n")
    
    # Get validation records
    legacy_val = find_pair_validation(pdb_id, idx1, idx2, legacy=True)
    modern_val = find_pair_validation(pdb_id, idx1, idx2, legacy=False)
    
    # Get base_pair records (may contain frame information)
    legacy_bp = find_base_pair_record(pdb_id, idx1, idx2, legacy=True)
    modern_bp = find_base_pair_record(pdb_id, idx1, idx2, legacy=False)
    
    print("LEGACY:")
    if legacy_val:
        print(f"  bp_type_id: {legacy_val.get('bp_type_id')}")
        print(f"  bp_type: {legacy_val.get('bp_type')}")
        print(f"  dir_x: {legacy_val.get('dir_x')}")
        print(f"  dir_y: {legacy_val.get('dir_y')}")
        print(f"  dir_z: {legacy_val.get('dir_z')}")
        print(f"  quality_score: {legacy_val.get('quality_score')}")
    else:
        print("  Validation record not found")
    
    if legacy_bp:
        print(f"  Base pair record found")
        if "orien1" in legacy_bp:
            print(f"  Has frame1 orientation")
        if "orien2" in legacy_bp:
            print(f"  Has frame2 orientation")
    
    print("\nMODERN:")
    if modern_val:
        print(f"  bp_type_id: {modern_val.get('bp_type_id')}")
        print(f"  bp_type: {modern_val.get('bp_type')}")
        print(f"  dir_x: {modern_val.get('dir_x')}")
        print(f"  dir_y: {modern_val.get('dir_y')}")
        print(f"  dir_z: {modern_val.get('dir_z')}")
        print(f"  quality_score: {modern_val.get('quality_score')}")
    else:
        print("  Validation record not found")
    
    if modern_bp:
        print(f"  Base pair record found")
    
    # Check if step parameters are stored anywhere
    print("\nSTEP PARAMETERS:")
    print("  Note: Step parameters are typically only stored for consecutive pairs in helices")
    print("  For single pairs, they are calculated on-the-fly during bp_type_id calculation")
    
    # Check if we can extract from rtn_val array (legacy stores pars in rtn_val[26+k])
    if legacy_val and "rtn_val" in legacy_val:
        rtn_val = legacy_val["rtn_val"]
        if len(rtn_val) > 31:
            print(f"\n  Legacy rtn_val array (if available):")
            print(f"    pars[1] (Slide/Shear): {rtn_val[26] if len(rtn_val) > 26 else 'N/A'}")
            print(f"    pars[2] (Rise/Stretch): {rtn_val[27] if len(rtn_val) > 27 else 'N/A'}")
            print(f"    pars[6] (Twist/Opening): {rtn_val[31] if len(rtn_val) > 31 else 'N/A'}")


def main():
    if len(sys.argv) < 4:
        print("Usage: python3 compare_step_params_for_pairs.py <pdb_id> <idx1> <idx2> [idx3] [idx4]")
        print("Example: python3 compare_step_params_for_pairs.py 6CAQ 1024 1188 980 997")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    pairs = []
    
    # Parse pairs
    if len(sys.argv) >= 4:
        idx1 = int(sys.argv[2])
        idx2 = int(sys.argv[3])
        pairs.append((idx1, idx2))
    
    if len(sys.argv) >= 6:
        idx3 = int(sys.argv[4])
        idx4 = int(sys.argv[5])
        pairs.append((idx3, idx4))
    
    for idx1, idx2 in pairs:
        compare_pair_step_params(pdb_id, idx1, idx2)


if __name__ == "__main__":
    main()


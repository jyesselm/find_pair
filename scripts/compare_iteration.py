#!/usr/bin/env python3
"""
Compare iteration states between legacy and modern code.

Usage:
    python3 scripts/compare_iteration.py <PDB_ID> [<iteration_num>] [--verbose]
"""

import sys
import json
import argparse
from pathlib import Path

def load_json_file(json_file: Path):
    """Load JSON file."""
    if not json_file.exists():
        return None
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    if isinstance(data, list):
        return data
    return [data]

def compare_iterations(legacy_data, modern_data, iteration_num=None, verbose=False):
    """Compare iteration states."""
    if legacy_data is None and modern_data is None:
        return [], "Both missing"
    
    if legacy_data is None:
        return [], "Legacy missing"
    
    if modern_data is None:
        return [], "Modern missing"
    
    # Filter by iteration number if specified
    if iteration_num is not None:
        legacy_data = [r for r in legacy_data if r.get("iteration_num") == iteration_num]
        modern_data = [r for r in modern_data if r.get("iteration_num") == iteration_num]
    
    # Build maps by iteration number
    # Note: Legacy uses "iteration_states" (plural) as type
    legacy_map = {r["iteration_num"]: r for r in legacy_data if r.get("type") in ["iteration_state", "iteration_states"]}
    modern_map = {r["iteration_num"]: r for r in modern_data if r.get("type") in ["iteration_state", "iteration_states"]}
    
    all_iterations = set(legacy_map.keys()) | set(modern_map.keys())
    
    differences = []
    
    for iter_num in sorted(all_iterations):
        leg_iter = legacy_map.get(iter_num)
        mod_iter = modern_map.get(iter_num)
        
        if leg_iter is None:
            differences.append(f"Iteration {iter_num}: Missing in legacy")
            continue
        
        if mod_iter is None:
            differences.append(f"Iteration {iter_num}: Missing in modern")
            continue
        
        # Compare fields
        leg_matched = leg_iter.get("num_matched", 0)
        mod_matched = mod_iter.get("num_matched", 0)
        
        if leg_matched != mod_matched:
            differences.append(f"Iteration {iter_num}: num_matched differs (legacy={leg_matched}, modern={mod_matched})")
        
        # Compare pairs
        leg_pairs = set(tuple(p) for p in leg_iter.get("pairs_found_in_iteration", []))
        mod_pairs = set(tuple(p) for p in mod_iter.get("pairs_found_in_iteration", []))
        
        if leg_pairs != mod_pairs:
            leg_only = leg_pairs - mod_pairs
            mod_only = mod_pairs - leg_pairs
            if leg_only:
                differences.append(f"Iteration {iter_num}: Pairs only in legacy: {sorted(leg_only)}")
            if mod_only:
                differences.append(f"Iteration {iter_num}: Pairs only in modern: {sorted(mod_only)}")
        
        # Compare matched residues
        leg_matched_res = set(leg_iter.get("matched_residues", []))
        mod_matched_res = set(mod_iter.get("matched_residues", []))
        
        if leg_matched_res != mod_matched_res:
            leg_only = leg_matched_res - mod_matched_res
            mod_only = mod_matched_res - leg_matched_res
            if leg_only:
                differences.append(f"Iteration {iter_num}: Matched residues only in legacy: {sorted(leg_only)}")
            if mod_only:
                differences.append(f"Iteration {iter_num}: Matched residues only in modern: {sorted(mod_only)}")
    
    return differences

def main():
    parser = argparse.ArgumentParser(description='Compare iteration states')
    parser.add_argument('pdb_id', help='PDB ID')
    parser.add_argument('iteration_num', type=int, nargs='?',
                       help='Specific iteration number (optional)')
    parser.add_argument('--legacy-dir', default='data/json_legacy',
                       help='Legacy JSON directory')
    parser.add_argument('--modern-dir', default='data/json',
                       help='Modern JSON directory')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Load JSON files
    # Legacy: data/json_legacy/iteration_states/<PDB_ID>.json
    # Modern: data/json/<PDB_ID>.json/iteration_states/<PDB_ID>.json
    legacy_file = Path(args.legacy_dir) / "iteration_states" / f"{args.pdb_id}.json"
    # Modern dir might be the PDB-specific directory, so check if it already includes the PDB name
    if Path(args.modern_dir).name == f"{args.pdb_id}.json":
        modern_file = Path(args.modern_dir) / "iteration_states" / f"{args.pdb_id}.json"
    else:
        modern_file = Path(args.modern_dir) / "iteration_states" / f"{args.pdb_id}.json"
    
    legacy_data = load_json_file(legacy_file)
    modern_data = load_json_file(modern_file)
    
    if legacy_data is None:
        print(f"ERROR: Legacy JSON not found: {legacy_file}")
        sys.exit(1)
    
    if modern_data is None:
        print(f"ERROR: Modern JSON not found: {modern_file}")
        print(f"       Note: Modern code needs to be updated to output iteration_states JSON")
        sys.exit(1)
    
    # Compare
    title = f"Comparing iteration states for {args.pdb_id}"
    if args.iteration_num is not None:
        title += f" (iteration {args.iteration_num})"
    print(title)
    print("=" * 70)
    
    differences = compare_iterations(legacy_data, modern_data, args.iteration_num, args.verbose)
    
    if not differences:
        print("✅ All iterations match")
    else:
        print(f"❌ {len(differences)} difference(s) found:")
        for diff in differences:
            print(f"  {diff}")
    
    if args.verbose:
        print("\n--- Legacy Iterations ---")
        for record in sorted([r for r in legacy_data if r.get("type") == "iteration_state"], 
                            key=lambda x: x.get("iteration_num", 0)):
            print(f"  Iteration {record.get('iteration_num')}: "
                  f"matched={record.get('num_matched')}, "
                  f"pairs={record.get('pairs_found_in_iteration', [])}")
        
        print("\n--- Modern Iterations ---")
        for record in sorted([r for r in modern_data if r.get("type") == "iteration_state"], 
                            key=lambda x: x.get("iteration_num", 0)):
            print(f"  Iteration {record.get('iteration_num')}: "
                  f"matched={record.get('num_matched')}, "
                  f"pairs={record.get('pairs_found_in_iteration', [])}")

if __name__ == '__main__':
    main()


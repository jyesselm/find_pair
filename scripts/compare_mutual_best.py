#!/usr/bin/env python3
"""
Compare mutual best decisions between legacy and modern code.

Usage:
    python3 scripts/compare_mutual_best.py <PDB_ID> [--verbose]
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

def compare_mutual_best(legacy_data, modern_data, verbose=False):
    """Compare mutual best decisions."""
    if legacy_data is None and modern_data is None:
        return [], "Both missing"
    
    if legacy_data is None:
        return [], "Legacy missing"
    
    if modern_data is None:
        return [], "Modern missing"
    
    # Filter by type
    legacy_decisions = [r for r in legacy_data if r.get("type") == "mutual_best_decision"]
    modern_decisions = [r for r in modern_data if r.get("type") == "mutual_best_decision"]
    
    # Build maps by (res1, res2) pair (normalized to smaller first)
    def normalize_pair(r1, r2):
        return (min(r1, r2), max(r1, r2))
    
    legacy_map = {}
    for dec in legacy_decisions:
        key = normalize_pair(dec["res1"], dec["res2"])
        legacy_map[key] = dec
    
    modern_map = {}
    for dec in modern_decisions:
        key = normalize_pair(dec["res1"], dec["res2"])
        modern_map[key] = dec
    
    all_pairs = set(legacy_map.keys()) | set(modern_map.keys())
    
    differences = []
    matches = 0
    
    for pair_key in sorted(all_pairs):
        leg_dec = legacy_map.get(pair_key)
        mod_dec = modern_map.get(pair_key)
        
        if leg_dec is None:
            differences.append(f"Pair {pair_key}: Missing in legacy")
            continue
        
        if mod_dec is None:
            differences.append(f"Pair {pair_key}: Missing in modern")
            continue
        
        # Compare fields
        leg_mutual = leg_dec.get("is_mutual", 0)
        mod_mutual = mod_dec.get("is_mutual", 0)
        leg_selected = leg_dec.get("was_selected", 0)
        mod_selected = mod_dec.get("was_selected", 0)
        leg_best_for_res1 = leg_dec.get("best_partner_for_res1", 0)
        mod_best_for_res1 = mod_dec.get("best_partner_for_res1", 0)
        leg_best_for_res2 = leg_dec.get("best_partner_for_res2", 0)
        mod_best_for_res2 = mod_dec.get("best_partner_for_res2", 0)
        
        if leg_mutual != mod_mutual:
            differences.append(f"Pair {pair_key}: is_mutual differs (legacy={leg_mutual}, modern={mod_mutual})")
        
        if leg_selected != mod_selected:
            differences.append(f"Pair {pair_key}: was_selected differs (legacy={leg_selected}, modern={mod_selected})")
        
        if leg_best_for_res1 != mod_best_for_res1:
            differences.append(f"Pair {pair_key}: best_partner_for_res1 differs (legacy={leg_best_for_res1}, modern={mod_best_for_res1})")
        
        if leg_best_for_res2 != mod_best_for_res2:
            differences.append(f"Pair {pair_key}: best_partner_for_res2 differs (legacy={leg_best_for_res2}, modern={mod_best_for_res2})")
        
        if (leg_mutual == mod_mutual and leg_selected == mod_selected and 
            leg_best_for_res1 == mod_best_for_res1 and leg_best_for_res2 == mod_best_for_res2):
            matches += 1
    
    return differences, matches

def main():
    parser = argparse.ArgumentParser(description='Compare mutual best decisions')
    parser.add_argument('pdb_id', help='PDB ID')
    parser.add_argument('--legacy-dir', default='data/json_legacy',
                       help='Legacy JSON directory')
    parser.add_argument('--modern-dir', default='data/json',
                       help='Modern JSON directory')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Load JSON files
    legacy_file = Path(args.legacy_dir) / "mutual_best_decisions" / f"{args.pdb_id}.json"
    # Modern dir might be the PDB-specific directory
    if Path(args.modern_dir).name == f"{args.pdb_id}.json":
        modern_file = Path(args.modern_dir) / "mutual_best_decisions" / f"{args.pdb_id}.json"
    else:
        modern_file = Path(args.modern_dir) / "mutual_best_decisions" / f"{args.pdb_id}.json"
    
    legacy_data = load_json_file(legacy_file)
    modern_data = load_json_file(modern_file)
    
    if legacy_data is None:
        print(f"ERROR: Legacy JSON not found: {legacy_file}")
        sys.exit(1)
    
    if modern_data is None:
        print(f"ERROR: Modern JSON not found: {modern_file}")
        print(f"       Note: Modern code needs to output mutual_best_decisions JSON")
        sys.exit(1)
    
    # Compare
    print(f"Comparing mutual best decisions for {args.pdb_id}")
    print("=" * 70)
    
    differences, matches = compare_mutual_best(legacy_data, modern_data, args.verbose)
    
    print(f"Total decisions: Legacy={len([r for r in legacy_data if r.get('type') == 'mutual_best_decision'])}, "
          f"Modern={len([r for r in modern_data if r.get('type') == 'mutual_best_decision'])})")
    print(f"Matches: {matches}")
    
    if not differences:
        print("✅ All mutual best decisions match")
    else:
        print(f"❌ {len(differences)} difference(s) found:")
        for diff in differences:
            print(f"  {diff}")
    
    if args.verbose:
        print("\n--- Legacy Decisions ---")
        for dec in sorted([r for r in legacy_data if r.get("type") == "mutual_best_decision"], 
                         key=lambda x: (x.get("res1", 0), x.get("res2", 0))):
            print(f"  ({dec['res1']}, {dec['res2']}): mutual={dec.get('is_mutual', 0)}, "
                  f"selected={dec.get('was_selected', 0)}, "
                  f"best_for_res1={dec.get('best_partner_for_res1', 0)}, "
                  f"best_for_res2={dec.get('best_partner_for_res2', 0)}")
        
        print("\n--- Modern Decisions ---")
        for dec in sorted([r for r in modern_data if r.get("type") == "mutual_best_decision"], 
                         key=lambda x: (x.get("res1", 0), x.get("res2", 0))):
            print(f"  ({dec['res1']}, {dec['res2']}): mutual={dec.get('is_mutual', 0)}, "
                  f"selected={dec.get('was_selected', 0)}, "
                  f"best_for_res1={dec.get('best_partner_for_res1', 0)}, "
                  f"best_for_res2={dec.get('best_partner_for_res2', 0)}")

if __name__ == '__main__':
    main()


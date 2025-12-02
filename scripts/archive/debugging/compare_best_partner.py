#!/usr/bin/env python3
"""
Compare best partner finding between legacy and modern code.

Usage:
    python3 scripts/compare_best_partner.py <PDB_ID> <res_i> [--verbose]
"""

import sys
import json
import argparse
from pathlib import Path

def load_json_file(json_file: Path):
    """Load JSON file, handling both array and single object formats."""
    if not json_file.exists():
        return None
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # If it's an array, return first matching record or all records
    if isinstance(data, list):
        return data
    
    # If it's a single object, wrap in list
    return [data]

def find_best_partner_for_residue(legacy_data, modern_data, res_i):
    """Find best partner candidate records for a specific residue."""
    legacy_record = None
    modern_record = None
    
    for record in legacy_data:
        if record.get("type") == "best_partner_candidates" and record.get("res_i") == res_i:
            legacy_record = record
            break
    
    for record in modern_data:
        if record.get("type") == "best_partner_candidates" and record.get("res_i") == res_i:
            modern_record = record
            break
    
    return legacy_record, modern_record

def compare_candidates(legacy_candidates, modern_candidates, verbose=False):
    """Compare candidate lists and report differences."""
    if legacy_candidates is None and modern_candidates is None:
        return True, "Both missing"
    
    if legacy_candidates is None:
        return False, "Legacy missing"
    
    if modern_candidates is None:
        return False, "Modern missing"
    
    # Build maps by res_j for easier comparison
    legacy_map = {c["res_j"]: c for c in legacy_candidates}
    modern_map = {c["res_j"]: c for c in modern_candidates}
    
    all_res_j = set(legacy_map.keys()) | set(modern_map.keys())
    
    matches = 0
    differences = []
    
    for res_j in sorted(all_res_j):
        leg_cand = legacy_map.get(res_j)
        mod_cand = modern_map.get(res_j)
        
        if leg_cand is None:
            differences.append(f"  Candidate {res_j}: Missing in legacy")
            continue
        
        if mod_cand is None:
            differences.append(f"  Candidate {res_j}: Missing in modern")
            continue
        
        # Compare fields
        leg_eligible = leg_cand.get("is_eligible", 0)
        mod_eligible = mod_cand.get("is_eligible", 0)
        leg_score = leg_cand.get("score", 0.0)
        mod_score = mod_cand.get("score", 0.0)
        leg_bp_type = leg_cand.get("bp_type_id", 0)
        mod_bp_type = mod_cand.get("bp_type_id", 0)
        
        if leg_eligible != mod_eligible:
            differences.append(f"  Candidate {res_j}: Eligibility differs (legacy={leg_eligible}, modern={mod_eligible})")
        
        if leg_eligible and mod_eligible:  # Only compare scores if both eligible
            score_diff = abs(leg_score - mod_score)
            if score_diff > 1e-6:
                differences.append(f"  Candidate {res_j}: Score differs (legacy={leg_score:.6f}, modern={mod_score:.6f}, diff={score_diff:.6e})")
        
        if leg_bp_type != mod_bp_type:
            differences.append(f"  Candidate {res_j}: bp_type_id differs (legacy={leg_bp_type}, modern={mod_bp_type})")
        
        if leg_eligible == mod_eligible and abs(leg_score - mod_score) < 1e-6 and leg_bp_type == mod_bp_type:
            matches += 1
    
    return len(differences) == 0, differences

def main():
    parser = argparse.ArgumentParser(description='Compare best partner finding')
    parser.add_argument('pdb_id', help='PDB ID')
    parser.add_argument('res_i', type=int, help='Residue index to check')
    parser.add_argument('--legacy-dir', default='data/json_legacy',
                       help='Legacy JSON directory')
    parser.add_argument('--modern-dir', default='data/json',
                       help='Modern JSON directory')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Load JSON files
    legacy_file = Path(args.legacy_dir) / "best_partner_candidates" / f"{args.pdb_id}.json"
    modern_file = Path(args.modern_dir) / "best_partner_candidates" / f"{args.pdb_id}.json"
    
    legacy_data = load_json_file(legacy_file)
    modern_data = load_json_file(modern_file)
    
    if legacy_data is None:
        print(f"ERROR: Legacy JSON not found: {legacy_file}")
        sys.exit(1)
    
    if modern_data is None:
        print(f"ERROR: Modern JSON not found: {modern_file}")
        print(f"       Note: Modern code needs to be updated to output best_partner_candidates JSON")
        sys.exit(1)
    
    # Find records for this residue
    legacy_record, modern_record = find_best_partner_for_residue(legacy_data, modern_data, args.res_i)
    
    if legacy_record is None:
        print(f"ERROR: No legacy record found for residue {args.res_i}")
        sys.exit(1)
    
    if modern_record is None:
        print(f"ERROR: No modern record found for residue {args.res_i}")
        sys.exit(1)
    
    # Compare best partner and score
    print(f"Comparing best partner for residue {args.res_i} in {args.pdb_id}")
    print("=" * 70)
    
    leg_best = legacy_record.get("best_partner")
    mod_best = modern_record.get("best_partner")
    leg_score = legacy_record.get("best_score", 0.0)
    mod_score = modern_record.get("best_score", 0.0)
    
    print(f"Best Partner: Legacy={leg_best}, Modern={mod_best}", end="")
    if leg_best == mod_best:
        print(" ✅ MATCH")
    else:
        print(" ❌ DIFFER")
    
    print(f"Best Score: Legacy={leg_score:.6f}, Modern={mod_score:.6f}", end="")
    if abs(leg_score - mod_score) < 1e-6:
        print(" ✅ MATCH")
    else:
        print(f" ❌ DIFFER (diff={abs(leg_score - mod_score):.6e})")
    
    # Compare candidates
    print(f"\nCandidates: Legacy={len(legacy_record.get('candidates', []))}, Modern={len(modern_record.get('candidates', []))}")
    
    is_match, differences = compare_candidates(
        legacy_record.get("candidates"),
        modern_record.get("candidates"),
        args.verbose
    )
    
    if is_match:
        print("✅ All candidates match")
    else:
        print(f"❌ {len(differences)} difference(s) found:")
        for diff in differences:
            print(diff)
    
    if args.verbose:
        print("\n--- Legacy Candidates ---")
        for cand in legacy_record.get("candidates", []):
            print(f"  res_j={cand['res_j']}, eligible={cand.get('is_eligible', 0)}, "
                  f"score={cand.get('score', 0.0):.6f}, bp_type={cand.get('bp_type_id', 0)}, "
                  f"is_best={cand.get('is_best', 0)}")
        
        print("\n--- Modern Candidates ---")
        for cand in modern_record.get("candidates", []):
            print(f"  res_j={cand['res_j']}, eligible={cand.get('is_eligible', 0)}, "
                  f"score={cand.get('score', 0.0):.6f}, bp_type={cand.get('bp_type_id', 0)}, "
                  f"is_best={cand.get('is_best', 0)}")

if __name__ == '__main__':
    main()


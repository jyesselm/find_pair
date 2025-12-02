#!/usr/bin/env python3
"""
Simple, direct check: Do all residue frames match between legacy and modern?

For each PDB:
1. Load base_pair JSON from both
2. For each pair, compare org_i/org_j and orien_i/orien_j
3. Report max differences

Usage:
    python3 scripts/verify_all_frames_match.py [--all] [--show-details]
"""

import sys
import json
import csv
import math
from pathlib import Path
import argparse


def compare_one_pdb(pdb_id: str, project_root: Path) -> dict:
    """Compare frames for one PDB - returns max differences."""
    result = {
        'pdb_id': pdb_id,
        'status': 'UNKNOWN',
        'pairs_compared': 0,
        'max_org_diff': 0.0,
        'max_orien_diff': 0.0
    }
    
    modern_file = project_root / "data" / "json" / "base_pair" / f"{pdb_id}.json"
    legacy_file = project_root / "data" / "json_legacy" / "base_pair" / f"{pdb_id}.json"
    
    if not modern_file.exists() or not legacy_file.exists():
        result['status'] = 'SKIP'
        return result
    
    try:
        with open(modern_file) as f:
            modern = json.load(f)
        with open(legacy_file) as f:
            legacy = json.load(f)
    except Exception as e:
        result['status'] = 'ERROR'
        result['error'] = str(e)
        return result
    
    # Match pairs by (base_i, base_j)
    legacy_by_pair = {(p['base_i'], p['base_j']): p for p in legacy}
    
    for mod_pair in modern:
        idx_i = mod_pair.get('base_i')
        idx_j = mod_pair.get('base_j')
        key = (idx_i, idx_j)
        
        if key not in legacy_by_pair:
            continue
        
        leg_pair = legacy_by_pair[key]
        result['pairs_compared'] += 1
        
        # Compare org_i
        if 'org_i' in mod_pair and 'org_i' in leg_pair:
            mod_org = mod_pair['org_i']
            leg_org = leg_pair['org_i']
            if len(mod_org) == 3 and len(leg_org) == 3:
                diff = math.sqrt(sum((mod_org[i] - leg_org[i])**2 for i in range(3)))
                result['max_org_diff'] = max(result['max_org_diff'], diff)
        
        # Compare org_j
        if 'org_j' in mod_pair and 'org_j' in leg_pair:
            mod_org = mod_pair['org_j']
            leg_org = leg_pair['org_j']
            if len(mod_org) == 3 and len(leg_org) == 3:
                diff = math.sqrt(sum((mod_org[i] - leg_org[i])**2 for i in range(3)))
                result['max_org_diff'] = max(result['max_org_diff'], diff)
        
        # Compare orien_i (nested 3x3)
        if 'orien_i' in mod_pair and 'orien_i' in leg_pair:
            mod_orien = mod_pair['orien_i']
            leg_orien = leg_pair['orien_i']
            
            # Both are nested [[],[],[]]
            if len(mod_orien) == 3 and len(leg_orien) == 3:
                for i in range(3):
                    if len(mod_orien[i]) == 3 and len(leg_orien[i]) == 3:
                        for j in range(3):
                            diff = abs(mod_orien[i][j] - leg_orien[i][j])
                            result['max_orien_diff'] = max(result['max_orien_diff'], diff)
        
        # Compare orien_j
        if 'orien_j' in mod_pair and 'orien_j' in leg_pair:
            mod_orien = mod_pair['orien_j']
            leg_orien = leg_pair['orien_j']
            
            if len(mod_orien) == 3 and len(leg_orien) == 3:
                for i in range(3):
                    if len(mod_orien[i]) == 3 and len(leg_orien[i]) == 3:
                        for j in range(3):
                            diff = abs(mod_orien[i][j] - leg_orien[i][j])
                            result['max_orien_diff'] = max(result['max_orien_diff'], diff)
    
    # Determine status
    if result['max_org_diff'] < 0.001 and result['max_orien_diff'] < 0.0001:
        result['status'] = 'PERFECT'
    elif result['max_org_diff'] < 0.01 and result['max_orien_diff'] < 0.001:
        result['status'] = 'EXCELLENT'
    elif result['max_org_diff'] < 0.1 and result['max_orien_diff'] < 0.01:
        result['status'] = 'ACCEPTABLE'
    else:
        result['status'] = 'FAIL'
    
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--all', action='store_true', help='Check all validated PDBs')
    parser.add_argument('--show-details', action='store_true', help='Show detailed output')
    parser.add_argument('--stop-on-fail', action='store_true', help='Stop on first failure')
    args = parser.parse_args()
    
    project_root = Path(__file__).parent.parent
    
    if not args.all:
        # Just test 1EHZ
        result = compare_one_pdb('1EHZ', project_root)
        print(f"\n1EHZ Frame Comparison:")
        print(f"  Status: {result['status']}")
        print(f"  Pairs compared: {result['pairs_compared']}")
        print(f"  Max origin diff: {result['max_org_diff']:.10f} Å")
        print(f"  Max orientation diff: {result['max_orien_diff']:.10f}")
        
        if result['status'] in ['PERFECT', 'EXCELLENT']:
            print(f"\n✅ Frames match perfectly!")
        return 0
    
    # Check all validated PDBs
    status_csv = project_root / "data" / "index_validation_status.csv"
    
    pass_pdbs = []
    with open(status_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['match_status'] == 'PASS':
                pass_pdbs.append(row['pdb_id'])
    
    print(f"Checking frames for {len(pass_pdbs)} validated PDBs...\n")
    
    perfect = 0
    excellent = 0
    acceptable = 0
    fail = 0
    skip = 0
    
    for i, pdb_id in enumerate(pass_pdbs):
        result = compare_one_pdb(pdb_id, project_root)
        
        if result['status'] == 'PERFECT':
            perfect += 1
        elif result['status'] == 'EXCELLENT':
            excellent += 1
            if args.show_details:
                print(f"  ⚠️  {pdb_id}: {result['status']} (org={result['max_org_diff']:.6f}, orien={result['max_orien_diff']:.6f})")
        elif result['status'] == 'ACCEPTABLE':
            acceptable += 1
            print(f"  ⚠️  {pdb_id}: Small diffs (org={result['max_org_diff']:.4f}, orien={result['max_orien_diff']:.4f})")
        elif result['status'] == 'FAIL':
            fail += 1
            print(f"  ❌ {pdb_id}: LARGE diffs (org={result['max_org_diff']:.4f}, orien={result['max_orien_diff']:.4f})")
            if args.stop_on_fail:
                print(f"\n⚠️  Stopping on first failure")
                return 1
        else:
            skip += 1
        
        if (i + 1) % 500 == 0:
            print(f"  Progress: {i+1}/{len(pass_pdbs)} - {perfect} perfect, {excellent} excellent")
    
    print(f"\n{'='*60}")
    print("FRAME VERIFICATION COMPLETE")
    print(f"{'='*60}")
    print(f"Total: {len(pass_pdbs)}")
    print(f"✅ Perfect (< 0.001 Å, < 0.0001): {perfect}")
    print(f"✅ Excellent (< 0.01 Å, < 0.001): {excellent}")
    print(f"⚠️  Acceptable (< 0.1 Å, < 0.01): {acceptable}")
    print(f"❌ Failed: {fail}")
    print(f"⏭️  Skipped: {skip}")
    
    if fail == 0:
        print(f"\n✅ All frames match within tolerance!")
        return 0
    else:
        print(f"\n❌ {fail} PDBs have large frame differences")
        return 1


if __name__ == '__main__':
    sys.exit(main())


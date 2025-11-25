#!/usr/bin/env python3
"""
Compare all PDBs between legacy and modern JSON outputs
"""

import json
import os
import glob
import sys
from pathlib import Path

def compare_base_pairs(pdb_id):
    """Compare base_pair records for a single PDB"""
    legacy_file = f'data/json_legacy/{pdb_id}_base_pair.json'
    modern_file = f'data/json/{pdb_id}_base_pair.json'
    
    if not os.path.exists(legacy_file):
        return None, f"No legacy file"
    if not os.path.exists(modern_file):
        return None, f"No modern file"
    
    try:
        legacy = json.load(open(legacy_file))
        modern = json.load(open(modern_file))
    except Exception as e:
        return None, f"JSON error: {e}"
    
    # Normalize pairs (min, max) to handle order differences
    legacy_pairs = set((min(r.get('base_i'), r.get('base_j')), max(r.get('base_i'), r.get('base_j'))) 
                       for r in legacy if r.get('base_i') and r.get('base_j'))
    modern_pairs = set((min(r.get('base_i'), r.get('base_j')), max(r.get('base_i'), r.get('base_j'))) 
                       for r in modern if r.get('base_i') and r.get('base_j'))
    
    common = legacy_pairs & modern_pairs
    missing = legacy_pairs - modern_pairs
    extra = modern_pairs - legacy_pairs
    
    # Build maps for detailed comparison
    legacy_map = {}
    for r in legacy:
        base_i, base_j = r.get('base_i'), r.get('base_j')
        if base_i and base_j:
            key = (min(base_i, base_j), max(base_i, base_j))
            legacy_map[key] = r
    
    modern_map = {}
    for r in modern:
        base_i, base_j = r.get('base_i'), r.get('base_j')
        if base_i and base_j:
            key = (min(base_i, base_j), max(base_i, base_j))
            modern_map[key] = r
    
    # Check dir_xyz, hbonds, orien_i, org_i for common pairs
    dir_xyz_mismatches = 0
    hbonds_mismatches = 0
    orien_mismatches = 0
    org_mismatches = 0
    
    for pair in list(common)[:10]:  # Check first 10 for detailed comparison
        leg = legacy_map[pair]
        mod = modern_map[pair]
        
        # Check dir_xyz
        leg_dir = leg.get('dir_xyz', [])
        mod_dir = mod.get('dir_xyz', [])
        if leg_dir and mod_dir and len(leg_dir) == 3 and len(mod_dir) == 3:
            diff = [abs(leg_dir[i] - mod_dir[i]) for i in range(3)]
            if max(diff) >= 0.0001:
                dir_xyz_mismatches += 1
        
        # Check hbonds
        leg_hb = leg.get('hbonds')
        mod_hb = mod.get('hbonds')
        if leg_hb != mod_hb:
            hbonds_mismatches += 1
        
        # Check orien_i first row
        leg_orien_i = leg.get('orien_i', [])
        mod_orien_i = mod.get('orien_i', [])
        if leg_orien_i and mod_orien_i and len(leg_orien_i) == 3 and len(mod_orien_i) == 3:
            if len(leg_orien_i[0]) == 3 and len(mod_orien_i[0]) == 3:
                diff = [abs(leg_orien_i[0][i] - mod_orien_i[0][i]) for i in range(3)]
                if max(diff) >= 0.0001:
                    orien_mismatches += 1
        
        # Check org_i
        leg_org_i = leg.get('org_i', [])
        mod_org_i = mod.get('org_i', [])
        if leg_org_i and mod_org_i and len(leg_org_i) == 3 and len(mod_org_i) == 3:
            diff = [abs(leg_org_i[i] - mod_org_i[i]) for i in range(3)]
            if max(diff) >= 0.0001:
                org_mismatches += 1
    
    match_rate = len(common) / len(legacy_pairs) * 100 if len(legacy_pairs) > 0 else 0
    
    return {
        'pdb': pdb_id,
        'legacy_count': len(legacy_pairs),
        'modern_count': len(modern_pairs),
        'common': len(common),
        'missing': len(missing),
        'extra': len(extra),
        'match_rate': match_rate,
        'dir_xyz_mismatches': dir_xyz_mismatches,
        'hbonds_mismatches': hbonds_mismatches,
        'orien_mismatches': orien_mismatches,
        'org_mismatches': org_mismatches,
    }, None

def main():
    # Get all PDBs with legacy base_pair files
    legacy_dir = 'data/json_legacy'
    legacy_files = glob.glob(f'{legacy_dir}/*_base_pair.json')
    pdbs = sorted([os.path.basename(f).replace('_base_pair.json', '') for f in legacy_files])
    
    print(f"Comparing {len(pdbs)} PDBs...")
    print()
    
    results = []
    errors = []
    
    for pdb in pdbs:
        result, error = compare_base_pairs(pdb)
        if result:
            results.append(result)
        else:
            errors.append((pdb, error))
    
    # Print summary
    print(f"{'PDB':<8} {'Legacy':<8} {'Modern':<8} {'Common':<8} {'Missing':<8} {'Extra':<8} {'Match %':<8} {'Issues':<10}")
    print('-' * 80)
    
    total_legacy = 0
    total_modern = 0
    total_common = 0
    total_missing = 0
    total_extra = 0
    
    for r in results:
        issues = []
        if r['dir_xyz_mismatches'] > 0:
            issues.append(f"dir_xyz({r['dir_xyz_mismatches']})")
        if r['hbonds_mismatches'] > 0:
            issues.append(f"hbonds({r['hbonds_mismatches']})")
        if r['orien_mismatches'] > 0:
            issues.append(f"orien({r['orien_mismatches']})")
        if r['org_mismatches'] > 0:
            issues.append(f"org({r['org_mismatches']})")
        
        issues_str = ", ".join(issues) if issues else "OK"
        
        print(f"{r['pdb']:<8} {r['legacy_count']:<8} {r['modern_count']:<8} {r['common']:<8} "
              f"{r['missing']:<8} {r['extra']:<8} {r['match_rate']:<8.1f} {issues_str:<10}")
        
        total_legacy += r['legacy_count']
        total_modern += r['modern_count']
        total_common += r['common']
        total_missing += r['missing']
        total_extra += r['extra']
    
    print('-' * 80)
    overall_match = total_common / total_legacy * 100 if total_legacy > 0 else 0
    print(f"{'TOTAL':<8} {total_legacy:<8} {total_modern:<8} {total_common:<8} "
          f"{total_missing:<8} {total_extra:<8} {overall_match:<8.1f}%")
    
    if errors:
        print()
        print(f"Errors ({len(errors)}):")
        for pdb, error in errors[:10]:
            print(f"  {pdb}: {error}")
    
    # Count PDBs with issues
    pdbs_with_issues = sum(1 for r in results if r['dir_xyz_mismatches'] > 0 or 
                          r['hbonds_mismatches'] > 0 or r['orien_mismatches'] > 0 or r['org_mismatches'] > 0)
    pdbs_perfect_match = len(results) - pdbs_with_issues
    
    print()
    print(f"Summary:")
    print(f"  Total PDBs: {len(results)}")
    print(f"  Perfect match (100% pairs + no data issues): {pdbs_perfect_match}")
    print(f"  PDBs with issues: {pdbs_with_issues}")
    print(f"  Overall match rate: {overall_match:.1f}%")

if __name__ == '__main__':
    main()


#!/usr/bin/env python3
"""
Analyze differences between legacy and modern JSON for all PDBs in test set
"""

import json
import subprocess
import sys
from pathlib import Path
from collections import defaultdict

def load_json_records(json_file, record_type):
    """Load all records of a specific type from JSON file."""
    records = []
    if not json_file.exists():
        return records
    
    with open(json_file) as f:
        data = json.load(f)
    
    if isinstance(data, list):
        for record in data:
            if record.get('type') == record_type:
                records.append(record)
    elif 'calculations' in data:
        for record in data['calculations']:
            if record.get('type') == record_type:
                records.append(record)
    
    return records

def normalize_pair_key(base_i, base_j):
    """Normalize pair key (order-independent)."""
    return tuple(sorted([base_i, base_j]))

def analyze_pdb(pdb_id, project_root):
    """Analyze a single PDB file."""
    print(f"\n{'='*80}")
    print(f"Analyzing {pdb_id}")
    print(f"{'='*80}")
    
    # Generate modern JSON
    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    json_file = project_root / "data" / "json" / f"{pdb_id}.json"
    
    if not pdb_file.exists():
        print(f"  ⚠️  PDB file not found: {pdb_file}")
        return None
    
    print(f"  Generating modern JSON...")
    result = subprocess.run(
        [str(project_root / "build" / "generate_modern_json"),
         str(pdb_file), str(json_file)],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(f"  ❌ Error generating JSON: {result.stderr}")
        return None
    
    # Load records
    legacy_base_pair = project_root / "data" / "json_legacy" / f"{pdb_id}_base_pair.json"
    modern_base_pair = project_root / "data" / "json" / f"{pdb_id}_base_pair.json"
    if not legacy_base_pair.exists():
        legacy_base_pair = project_root / "data" / "json_legacy" / f"{pdb_id}.json"
    if not modern_base_pair.exists():
        modern_base_pair = project_root / "data" / "json" / f"{pdb_id}.json"
    
    legacy_validation = project_root / "data" / "json_legacy" / f"{pdb_id}_pair_validation.json"
    modern_validation = project_root / "data" / "json" / f"{pdb_id}_pair_validation.json"
    if not legacy_validation.exists():
        legacy_validation = project_root / "data" / "json_legacy" / f"{pdb_id}.json"
    if not modern_validation.exists():
        modern_validation = project_root / "data" / "json" / f"{pdb_id}.json"
    
    if not legacy_base_pair.exists() or not modern_base_pair.exists():
        print(f"  ⚠️  JSON files not found")
        return None
    
    # Compare base pairs
    legacy_bp = load_json_records(legacy_base_pair, 'base_pair')
    modern_bp = load_json_records(modern_base_pair, 'base_pair')
    
    legacy_map = {}
    modern_map = {}
    
    for rec in legacy_bp:
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            key = normalize_pair_key(base_i, base_j)
            if key not in legacy_map:
                legacy_map[key] = []
            legacy_map[key].append(rec)
    
    for rec in modern_bp:
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            key = normalize_pair_key(base_i, base_j)
            if key not in modern_map:
                modern_map[key] = []
            modern_map[key].append(rec)
    
    legacy_keys = set(legacy_map.keys())
    modern_keys = set(modern_map.keys())
    common_keys = legacy_keys & modern_keys
    missing_keys = legacy_keys - modern_keys
    extra_keys = modern_keys - legacy_keys
    
    print(f"\n  Base Pairs:")
    print(f"    Legacy: {len(legacy_bp)} records ({len(legacy_keys)} unique)")
    print(f"    Modern: {len(modern_bp)} records ({len(modern_keys)} unique)")
    print(f"    Common: {len(common_keys)}")
    print(f"    Missing in modern: {len(missing_keys)}")
    print(f"    Extra in modern: {len(extra_keys)}")
    
    # Compare pair validations
    legacy_val = load_json_records(legacy_validation, 'pair_validation')
    modern_val = load_json_records(modern_validation, 'pair_validation')
    
    legacy_val_map = {}
    modern_val_map = {}
    
    for rec in legacy_val:
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            key = normalize_pair_key(base_i, base_j)
            legacy_val_map[key] = rec
    
    for rec in modern_val:
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            key = normalize_pair_key(base_i, base_j)
            modern_val_map[key] = rec
    
    legacy_val_keys = set(legacy_val_map.keys())
    modern_val_keys = set(modern_val_map.keys())
    common_val_keys = legacy_val_keys & modern_val_keys
    
    # Count validation mismatches
    mismatches = []
    for key in common_val_keys:
        leg = legacy_val_map[key]
        mod = modern_val_map[key]
        if leg.get('is_valid') != mod.get('is_valid'):
            mismatches.append((key, leg.get('is_valid'), mod.get('is_valid')))
    
    legacy_valid = sum(1 for r in legacy_val if r.get('is_valid') == 1)
    modern_valid = sum(1 for r in modern_val if r.get('is_valid') == 1)
    
    print(f"\n  Pair Validations:")
    print(f"    Legacy: {len(legacy_val)} records ({legacy_valid} valid)")
    print(f"    Modern: {len(modern_val)} records ({modern_valid} valid)")
    print(f"    Common: {len(common_val_keys)}")
    print(f"    Validation mismatches: {len(mismatches)}")
    
    # Analyze missing pairs
    if missing_keys:
        print(f"\n  Missing in modern ({len(missing_keys)} pairs):")
        for key in sorted(missing_keys)[:5]:
            rec = legacy_map[key][0]
            print(f"    {key}: {rec.get('bp_type', 'N/A')}")
        if len(missing_keys) > 5:
            print(f"    ... and {len(missing_keys) - 5} more")
    
    # Analyze extra pairs
    if extra_keys:
        print(f"\n  Extra in modern ({len(extra_keys)} pairs):")
        for key in sorted(extra_keys)[:5]:
            rec = modern_map[key][0]
            # Check if this pair is valid in modern
            is_valid = False
            if key in modern_val_map:
                is_valid = modern_val_map[key].get('is_valid') == 1
            print(f"    {key}: {rec.get('bp_type', 'N/A')} (valid={is_valid})")
        if len(extra_keys) > 5:
            print(f"    ... and {len(extra_keys) - 5} more")
    
    # Show validation mismatches
    if mismatches:
        print(f"\n  Validation Mismatches ({len(mismatches)} pairs):")
        for key, leg_val, mod_val in mismatches[:5]:
            print(f"    {key}: legacy={leg_val}, modern={mod_val}")
        if len(mismatches) > 5:
            print(f"    ... and {len(mismatches) - 5} more")
    
    return {
        'pdb_id': pdb_id,
        'base_pairs': {
            'legacy_total': len(legacy_bp),
            'legacy_unique': len(legacy_keys),
            'modern_total': len(modern_bp),
            'modern_unique': len(modern_keys),
            'common': len(common_keys),
            'missing': len(missing_keys),
            'extra': len(extra_keys),
            'missing_pairs': list(missing_keys),
            'extra_pairs': list(extra_keys)
        },
        'validations': {
            'legacy_total': len(legacy_val),
            'legacy_valid': legacy_valid,
            'modern_total': len(modern_val),
            'modern_valid': modern_valid,
            'common': len(common_val_keys),
            'mismatches': len(mismatches),
            'mismatch_pairs': [k for k, _, _ in mismatches]
        }
    }

def main():
    project_root = Path(__file__).parent.parent
    
    # Load test set
    test_set_file = project_root / "data" / "test_sets" / "test_set_10.json"
    if not test_set_file.exists():
        print(f"Test set file not found: {test_set_file}")
        sys.exit(1)
    
    with open(test_set_file) as f:
        test_set = json.load(f)
    
    pdb_ids = test_set.get('pdb_ids', [])
    
    print(f"Analyzing {len(pdb_ids)} PDB files from test set")
    print(f"PDB IDs: {', '.join(pdb_ids)}")
    
    results = []
    for pdb_id in pdb_ids:
        result = analyze_pdb(pdb_id, project_root)
        if result:
            results.append(result)
    
    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    
    print(f"\nBase Pairs:")
    total_legacy = sum(r['base_pairs']['legacy_unique'] for r in results)
    total_modern = sum(r['base_pairs']['modern_unique'] for r in results)
    total_common = sum(r['base_pairs']['common'] for r in results)
    total_missing = sum(r['base_pairs']['missing'] for r in results)
    total_extra = sum(r['base_pairs']['extra'] for r in results)
    
    print(f"  Total legacy unique: {total_legacy}")
    print(f"  Total modern unique: {total_modern}")
    print(f"  Total common: {total_common}")
    print(f"  Total missing in modern: {total_missing}")
    print(f"  Total extra in modern: {total_extra}")
    
    print(f"\nPair Validations:")
    total_val_legacy = sum(r['validations']['legacy_total'] for r in results)
    total_val_modern = sum(r['validations']['modern_total'] for r in results)
    total_val_mismatches = sum(r['validations']['mismatches'] for r in results)
    
    print(f"  Total legacy validations: {total_val_legacy}")
    print(f"  Total modern validations: {total_val_modern}")
    print(f"  Total validation mismatches: {total_val_mismatches}")
    
    # Find patterns in missing/extra pairs
    print(f"\nPatterns:")
    
    # PDBs with most missing pairs
    missing_by_pdb = [(r['pdb_id'], r['base_pairs']['missing']) for r in results]
    missing_by_pdb.sort(key=lambda x: x[1], reverse=True)
    if missing_by_pdb and missing_by_pdb[0][1] > 0:
        print(f"  PDBs with most missing pairs:")
        for pdb_id, count in missing_by_pdb[:5]:
            if count > 0:
                print(f"    {pdb_id}: {count} missing")
    
    # PDBs with most extra pairs
    extra_by_pdb = [(r['pdb_id'], r['base_pairs']['extra']) for r in results]
    extra_by_pdb.sort(key=lambda x: x[1], reverse=True)
    if extra_by_pdb and extra_by_pdb[0][1] > 0:
        print(f"  PDBs with most extra pairs:")
        for pdb_id, count in extra_by_pdb[:5]:
            if count > 0:
                print(f"    {pdb_id}: {count} extra")
    
    # PDBs with validation mismatches
    mismatch_by_pdb = [(r['pdb_id'], r['validations']['mismatches']) for r in results]
    mismatch_by_pdb.sort(key=lambda x: x[1], reverse=True)
    if mismatch_by_pdb and mismatch_by_pdb[0][1] > 0:
        print(f"  PDBs with validation mismatches:")
        for pdb_id, count in mismatch_by_pdb:
            if count > 0:
                print(f"    {pdb_id}: {count} mismatches")

if __name__ == '__main__':
    main()


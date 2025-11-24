#!/usr/bin/env python3
"""
Quick comparison script to analyze differences between legacy and modern JSON
"""

import json
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
    
    # Check if it's a split file format
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

def compare_base_pairs(legacy_file, modern_file):
    """Compare base_pair records."""
    legacy_records = load_json_records(legacy_file, 'base_pair')
    modern_records = load_json_records(modern_file, 'base_pair')
    
    print(f"\n=== Base Pair Comparison ===")
    print(f"Legacy: {len(legacy_records)} records")
    print(f"Modern: {len(modern_records)} records")
    
    # Build maps
    legacy_map = {}
    modern_map = {}
    
    for rec in legacy_records:
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            key = normalize_pair_key(base_i, base_j)
            legacy_map[key] = rec
    
    for rec in modern_records:
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        basepair_idx = rec.get('basepair_idx')
        if base_i is not None and base_j is not None:
            key = normalize_pair_key(base_i, base_j)
            modern_map[key] = rec
    
    legacy_keys = set(legacy_map.keys())
    modern_keys = set(modern_map.keys())
    common_keys = legacy_keys & modern_keys
    
    print(f"Common pairs: {len(common_keys)}")
    print(f"Missing in modern: {len(legacy_keys - modern_keys)}")
    print(f"Extra in modern: {len(modern_keys - legacy_keys)}")
    
    # Show missing pairs
    missing = legacy_keys - modern_keys
    if missing:
        print(f"\nMissing in modern ({len(missing)} pairs):")
        for key in sorted(missing)[:10]:  # Show first 10
            rec = legacy_map[key]
            print(f"  ({rec['base_i']}, {rec['base_j']}) - {rec.get('bp_type', 'N/A')}")
        if len(missing) > 10:
            print(f"  ... and {len(missing) - 10} more")
    
    # Show extra pairs
    extra = modern_keys - legacy_keys
    if extra:
        print(f"\nExtra in modern ({len(extra)} pairs):")
        for key in sorted(extra)[:10]:  # Show first 10
            rec = modern_map[key]
            print(f"  ({rec['base_i']}, {rec['base_j']}) - {rec.get('bp_type', 'N/A')} - idx={rec.get('basepair_idx', 'N/A')}")
        if len(extra) > 10:
            print(f"  ... and {len(extra) - 10} more")
    
    # Check for index assignment
    modern_with_idx = sum(1 for rec in modern_records if 'basepair_idx' in rec)
    print(f"\nModern records with basepair_idx: {modern_with_idx}/{len(modern_records)}")
    
    # Check hbond indices
    hbond_count = 0
    hbond_with_idx = 0
    for rec in modern_records:
        hbonds = rec.get('hbonds', [])
        hbond_count += len(hbonds)
        hbond_with_idx += sum(1 for hb in hbonds if 'hbond_idx' in hb)
    
    print(f"Modern hbonds with hbond_idx: {hbond_with_idx}/{hbond_count}")
    
    return len(common_keys), len(missing), len(extra)

def compare_pair_validations(legacy_file, modern_file):
    """Compare pair_validation records."""
    legacy_records = load_json_records(legacy_file, 'pair_validation')
    modern_records = load_json_records(modern_file, 'pair_validation')
    
    print(f"\n=== Pair Validation Comparison ===")
    print(f"Legacy: {len(legacy_records)} records")
    print(f"Modern: {len(modern_records)} records")
    
    # Count validations
    legacy_valid = sum(1 for r in legacy_records if r.get('is_valid') == 1)
    modern_valid = sum(1 for r in modern_records if r.get('is_valid') == 1)
    
    print(f"Legacy valid: {legacy_valid}/{len(legacy_records)}")
    print(f"Modern valid: {modern_valid}/{len(modern_records)}")
    
    # Build maps
    legacy_map = {}
    modern_map = {}
    
    for rec in legacy_records:
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            key = normalize_pair_key(base_i, base_j)
            legacy_map[key] = rec
    
    for rec in modern_records:
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            key = normalize_pair_key(base_i, base_j)
            modern_map[key] = rec
    
    legacy_keys = set(legacy_map.keys())
    modern_keys = set(modern_map.keys())
    common_keys = legacy_keys & modern_keys
    
    print(f"Common pairs: {len(common_keys)}")
    
    # Check validation mismatches
    mismatches = []
    for key in common_keys:
        leg = legacy_map[key]
        mod = modern_map[key]
        if leg.get('is_valid') != mod.get('is_valid'):
            mismatches.append((key, leg.get('is_valid'), mod.get('is_valid')))
    
    print(f"Validation mismatches: {len(mismatches)}")
    if mismatches:
        print(f"\nFirst 10 validation mismatches:")
        for key, leg_val, mod_val in mismatches[:10]:
            print(f"  {key}: legacy={leg_val}, modern={mod_val}")
    
    return len(common_keys), len(mismatches)

def compare_find_bestpair(legacy_file, modern_file):
    """Compare find_bestpair_selection records."""
    legacy_records = load_json_records(legacy_file, 'find_bestpair_selection')
    modern_records = load_json_records(modern_file, 'find_bestpair_selection')
    
    print(f"\n=== Find Bestpair Selection Comparison ===")
    
    if not legacy_records and not modern_records:
        print("No records found in either file")
        return 0, 0, 0
    
    # Extract pairs from records
    legacy_pairs = set()
    modern_pairs = set()
    
    for rec in legacy_records:
        pairs = rec.get('selected_pairs', [])
        for pair in pairs:
            if len(pair) >= 2:
                key = normalize_pair_key(pair[0], pair[1])
                legacy_pairs.add(key)
    
    for rec in modern_records:
        pairs = rec.get('selected_pairs', [])
        for pair in pairs:
            if len(pair) >= 2:
                key = normalize_pair_key(pair[0], pair[1])
                modern_pairs.add(key)
    
    print(f"Legacy selected: {len(legacy_pairs)} pairs")
    print(f"Modern selected: {len(modern_pairs)} pairs")
    
    common = legacy_pairs & modern_pairs
    missing = legacy_pairs - modern_pairs
    extra = modern_pairs - legacy_pairs
    
    print(f"Common pairs: {len(common)}")
    print(f"Missing in modern: {len(missing)}")
    print(f"Extra in modern: {len(extra)}")
    
    if missing:
        print(f"\nMissing in modern ({len(missing)} pairs):")
        for key in sorted(missing)[:10]:
            print(f"  {key}")
        if len(missing) > 10:
            print(f"  ... and {len(missing) - 10} more")
    
    if extra:
        print(f"\nExtra in modern ({len(extra)} pairs):")
        for key in sorted(extra)[:10]:
            print(f"  {key}")
        if len(extra) > 10:
            print(f"  ... and {len(extra) - 10} more")
    
    return len(common), len(missing), len(extra)

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/quick_compare.py <PDB_ID>")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    project_root = Path(__file__).parent.parent
    
    # Try split files first
    legacy_base_pair = project_root / "data" / "json_legacy" / f"{pdb_id}_base_pair.json"
    modern_base_pair = project_root / "data" / "json" / f"{pdb_id}_base_pair.json"
    
    # Fallback to main JSON files
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
    
    legacy_bestpair = project_root / "data" / "json_legacy" / f"{pdb_id}_find_bestpair_selection.json"
    modern_bestpair = project_root / "data" / "json" / f"{pdb_id}_find_bestpair_selection.json"
    if not legacy_bestpair.exists():
        legacy_bestpair = project_root / "data" / "json_legacy" / f"{pdb_id}.json"
    if not modern_bestpair.exists():
        modern_bestpair = project_root / "data" / "json" / f"{pdb_id}.json"
    
    print(f"Comparing {pdb_id}")
    print(f"Legacy: {legacy_base_pair}")
    print(f"Modern: {modern_base_pair}")
    
    if not legacy_base_pair.exists():
        print(f"ERROR: Legacy file not found: {legacy_base_pair}")
        sys.exit(1)
    if not modern_base_pair.exists():
        print(f"ERROR: Modern file not found: {modern_base_pair}")
        sys.exit(1)
    
    # Run comparisons
    compare_base_pairs(legacy_base_pair, modern_base_pair)
    compare_pair_validations(legacy_validation, modern_validation)
    compare_find_bestpair(legacy_bestpair, modern_bestpair)
    
    print("\n=== Summary ===")
    print("Comparison complete!")

if __name__ == '__main__':
    main()


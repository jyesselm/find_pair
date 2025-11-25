#!/usr/bin/env python3
"""
Compare H-bond detection between modern and legacy JSON files
"""

import json
import sys
from pathlib import Path
from collections import defaultdict

def load_json_file(filepath):
    """Load JSON file, handling both single objects and arrays, and newline-delimited JSON"""
    records = []
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Try to parse as single JSON object first
    decoder = json.JSONDecoder()
    idx = 0
    
    while idx < len(content):
        # Skip whitespace
        while idx < len(content) and content[idx].isspace():
            idx += 1
        
        if idx >= len(content):
            break
        
        # Try to parse a JSON object
        try:
            obj, new_idx = decoder.raw_decode(content, idx)
            records.append(obj)
            idx = new_idx
        except json.JSONDecodeError:
            # If parsing fails, try to find next object start
            next_brace = content.find('{', idx + 1)
            if next_brace == -1:
                break
            idx = next_brace
    
    # Handle different JSON structures
    if len(records) == 1:
        data = records[0]
        if isinstance(data, list):
            return data
        elif isinstance(data, dict):
            # Check for 'calculations' array (legacy format)
            if 'calculations' in data:
                # Merge with any additional records
                all_records = data['calculations'] + [r for r in records[1:] if isinstance(r, dict)]
                return all_records
            # Check if it's a single record
            elif 'type' in data:
                return records
        return records
    else:
        # Multiple JSON objects - return all
        return records

def extract_hbonds_from_json(json_data, record_type='hbond_list'):
    """Extract H-bond records from JSON data"""
    hbonds_by_pair = {}
    
    for record in json_data:
        # Check for hbond_list records (both legacy and modern format)
        if record.get('type') == record_type:
            base_i = record.get('base_i')
            base_j = record.get('base_j')
            if base_i is not None and base_j is not None:
                # Normalize pair order to (min, max) for consistent comparison
                pair_key = (min(base_i, base_j), max(base_i, base_j))
                hbonds = record.get('hbonds', [])
                hbonds_by_pair[pair_key] = hbonds
    
    return hbonds_by_pair

def normalize_hbond(hbond):
    """Normalize H-bond for comparison (handle atom name variations)"""
    donor = hbond.get('donor_atom', '').strip()
    acceptor = hbond.get('acceptor_atom', '').strip()
    distance = round(hbond.get('distance', 0.0), 3)
    hb_type = hbond.get('type', ' ')
    
    return {
        'donor': donor,
        'acceptor': acceptor,
        'distance': distance,
        'type': hb_type
    }

def compare_hbonds(legacy_hbonds, modern_hbonds):
    """Compare H-bonds between legacy and modern"""
    legacy_norm = [normalize_hbond(hb) for hb in legacy_hbonds]
    modern_norm = [normalize_hbond(hb) for hb in modern_hbonds]
    
    # Create sets for comparison (using tuple of key fields)
    legacy_set = {tuple(sorted([hb['donor'], hb['acceptor']])) + (hb['distance'], hb['type']) 
                  for hb in legacy_norm}
    modern_set = {tuple(sorted([hb['donor'], hb['acceptor']])) + (hb['distance'], hb['type']) 
                  for hb in modern_norm}
    
    missing_in_modern = legacy_set - modern_set
    extra_in_modern = modern_set - legacy_set
    common = legacy_set & modern_set
    
    return {
        'legacy_count': len(legacy_hbonds),
        'modern_count': len(modern_hbonds),
        'common': len(common),
        'missing_in_modern': len(missing_in_modern),
        'extra_in_modern': len(extra_in_modern),
        'missing_details': missing_in_modern,
        'extra_details': extra_in_modern,
        'legacy_hbonds': legacy_norm,
        'modern_hbonds': modern_norm
    }

def main():
    if len(sys.argv) < 3:
        print("Usage: compare_hbond_detection.py <pdb_id> <legacy_json_dir> <modern_json_dir>")
        print("Example: compare_hbond_detection.py 3G8T data/json_legacy data/json")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    legacy_dir = Path(sys.argv[2])
    modern_dir = Path(sys.argv[3])
    
    # Find legacy H-bond file (try split file first, then main file)
    legacy_file = legacy_dir / f"{pdb_id}_hbond_list.json"
    if not legacy_file.exists():
        # Try main JSON file
        legacy_main = legacy_dir / f"{pdb_id}.json"
        if legacy_main.exists():
            legacy_data = load_json_file(legacy_main)
        else:
            print(f"ERROR: Legacy H-bond file not found: {legacy_file} or {legacy_main}")
            sys.exit(1)
    else:
        legacy_data = load_json_file(legacy_file)
    
    # Find modern H-bond file (try split file first, then main file)
    modern_file = modern_dir / f"{pdb_id}_hbond_list.json"
    if not modern_file.exists():
        modern_main = modern_dir / f"{pdb_id}.json"
        if modern_main.exists():
            modern_data = load_json_file(modern_main)
        else:
            print(f"ERROR: Modern H-bond file not found: {modern_file} or {modern_main}")
            sys.exit(1)
    else:
        modern_data = load_json_file(modern_file)
    
    # Extract H-bonds from hbond_list records
    legacy_hbonds = extract_hbonds_from_json(legacy_data, 'hbond_list')
    modern_hbonds = extract_hbonds_from_json(modern_data, 'hbond_list')
    
    print(f"\n=== H-bond Comparison for {pdb_id} ===\n")
    print(f"Legacy pairs with H-bonds: {len(legacy_hbonds)}")
    print(f"Modern pairs with H-bonds: {len(modern_hbonds)}\n")
    
    # Compare each pair
    all_pairs = set(legacy_hbonds.keys()) | set(modern_hbonds.keys())
    
    total_legacy = 0
    total_modern = 0
    total_common = 0
    total_missing = 0
    total_extra = 0
    
    differences = []
    
    for pair in sorted(all_pairs):
        legacy_pair_hbonds = legacy_hbonds.get(pair, [])
        modern_pair_hbonds = modern_hbonds.get(pair, [])
        
        if not legacy_pair_hbonds and not modern_pair_hbonds:
            continue
        
        comparison = compare_hbonds(legacy_pair_hbonds, modern_pair_hbonds)
        
        total_legacy += comparison['legacy_count']
        total_modern += comparison['modern_count']
        total_common += comparison['common']
        total_missing += comparison['missing_in_modern']
        total_extra += comparison['extra_in_modern']
        
        if comparison['missing_in_modern'] > 0 or comparison['extra_in_modern'] > 0:
            differences.append((pair, comparison))
            print(f"Pair ({pair[0]}, {pair[1]}):")
            print(f"  Legacy: {comparison['legacy_count']} H-bonds")
            print(f"  Modern: {comparison['modern_count']} H-bonds")
            print(f"  Common: {comparison['common']}")
            
            if comparison['missing_in_modern'] > 0:
                print(f"  Missing in modern: {comparison['missing_in_modern']}")
                for missing in comparison['missing_details']:
                    print(f"    - {missing}")
            
            if comparison['extra_in_modern'] > 0:
                print(f"  Extra in modern: {comparison['extra_in_modern']}")
                for extra in comparison['extra_details']:
                    print(f"    + {extra}")
            print()
    
    print(f"\n=== Summary ===")
    print(f"Total legacy H-bonds: {total_legacy}")
    print(f"Total modern H-bonds: {total_modern}")
    print(f"Common H-bonds: {total_common}")
    print(f"Missing in modern: {total_missing}")
    print(f"Extra in modern: {total_extra}")
    print(f"Match rate: {total_common / total_legacy * 100:.1f}%" if total_legacy > 0 else "N/A")
    print(f"Pairs with differences: {len(differences)}")

if __name__ == '__main__':
    main()


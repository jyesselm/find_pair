#!/usr/bin/env python3
"""
Debug script to investigate validation mismatches between legacy and modern implementations.
This script will help identify exactly where the differences occur.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

def load_json_records(json_file: Path) -> List[Dict]:
    """Load JSON records from file."""
    try:
        with open(json_file) as f:
            data = json.load(f)
        if isinstance(data, list):
            return data
        else:
            return data.get('calculations', [])
    except Exception as e:
        print(f"Error loading {json_file}: {e}")
        return []

def find_validation_record(records: List[Dict], base_i: int, base_j: int) -> Optional[Dict]:
    """Find pair_validation record for given pair."""
    for r in records:
        if (r.get('type') == 'pair_validation' and 
            r.get('base_i') == base_i and r.get('base_j') == base_j):
            return r
    return None

def find_hbond_record(records: List[Dict], base_i: int, base_j: int) -> Optional[Dict]:
    """Find hbond_list record for given pair."""
    for r in records:
        if (r.get('type') == 'hbond_list' and 
            r.get('base_i') == base_i and r.get('base_j') == base_j):
            return r
    return None

def analyze_pair(pdb_id: str, base_i: int, base_j: int, project_root: Path):
    """Analyze a specific pair to find validation differences."""
    print(f"\n{'='*80}")
    print(f"Analyzing {pdb_id} pair ({base_i}, {base_j})")
    print(f"{'='*80}\n")
    
    legacy_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
    modern_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
    
    if not legacy_file.exists():
        print(f"ERROR: Legacy file not found: {legacy_file}")
        return
    if not modern_file.exists():
        print(f"ERROR: Modern file not found: {modern_file}")
        return
    
    legacy_records = load_json_records(legacy_file)
    modern_records = load_json_records(modern_file)
    
    # Get validation records
    legacy_val = find_validation_record(legacy_records, base_i, base_j)
    modern_val = find_validation_record(modern_records, base_i, base_j)
    
    print("=== VALIDATION RECORDS ===")
    if legacy_val:
        print(f"Legacy:")
        print(f"  is_valid: {legacy_val.get('is_valid')}")
        print(f"  num_base_hb: {legacy_val.get('num_base_hb', 'N/A')}")
        print(f"  num_o2_hb: {legacy_val.get('num_o2_hb', 'N/A')}")
        print(f"  min_base_hb: {legacy_val.get('min_base_hb', 'N/A')}")
        print(f"  bp_type_id: {legacy_val.get('bp_type_id', 'N/A')}")
    else:
        print("Legacy: No validation record found")
    
    if modern_val:
        print(f"\nModern:")
        print(f"  is_valid: {modern_val.get('is_valid')}")
        print(f"  num_base_hb: {modern_val.get('num_base_hb', 'N/A')}")
        print(f"  num_o2_hb: {modern_val.get('num_o2_hb', 'N/A')}")
        print(f"  min_base_hb: {modern_val.get('min_base_hb', 'N/A')}")
        print(f"  bp_type_id: {modern_val.get('bp_type_id', 'N/A')}")
    else:
        print("Modern: No validation record found")
    
    # Get H-bond records
    legacy_hb = find_hbond_record(legacy_records, base_i, base_j)
    modern_hb = find_hbond_record(modern_records, base_i, base_j)
    
    print(f"\n=== H-BOND RECORDS ===")
    if legacy_hb:
        print(f"Legacy: {legacy_hb.get('num_hbonds', 0)} H-bonds")
        for i, hb in enumerate(legacy_hb.get('hbonds', [])[:10]):
            donor = hb.get('donor_atom', '')
            acceptor = hb.get('acceptor_atom', '')
            dist = hb.get('distance', 0)
            print(f"  {i+1}. {donor} - {acceptor} ({dist:.3f} Å)")
    else:
        print("Legacy: No hbond_list record (pair rejected)")
    
    if modern_hb:
        print(f"\nModern: {modern_hb.get('num_hbonds', 0)} H-bonds")
        for i, hb in enumerate(modern_hb.get('hbonds', [])[:10]):
            donor = hb.get('donor_atom', '')
            acceptor = hb.get('acceptor_atom', '')
            dist = hb.get('distance', 0)
            print(f"  {i+1}. {donor} - {acceptor} ({dist:.3f} Å)")
            
            # Analyze each H-bond
            print(f"      Analysis:")
            PO = [' O1P', ' O2P', " O3'", " O4'", " O5'", ' N7 ']
            donor_in_po = donor in PO
            acceptor_in_po = acceptor in PO
            both_in_po = donor_in_po and acceptor_in_po
            print(f"        Donor in PO: {donor_in_po}, Acceptor in PO: {acceptor_in_po}")
            print(f"        Both in PO: {both_in_po}")
            
            # Check is_baseatom
            def is_baseatom(atomname):
                if atomname == ' C5M':
                    return True
                if (len(atomname) >= 4 and atomname[0] == ' ' and 
                    atomname[1] not in 'HP' and 
                    atomname[2].isdigit() and 
                    atomname[3] == ' '):
                    return True
                return False
            
            donor_is_base = is_baseatom(donor)
            acceptor_is_base = is_baseatom(acceptor)
            print(f"        Donor is_baseatom: {donor_is_base}")
            print(f"        Acceptor is_baseatom: {acceptor_is_base}")
            
            # Simulate pattern creation
            def create_pattern(atom_name):
                pattern = list(atom_name)
                for i in range(min(4, len(pattern))):
                    if not pattern[i].isalpha():
                        pattern[i] = '.'
                return ''.join(pattern[:4])
            
            donor_pattern = create_pattern(donor)
            acceptor_pattern = create_pattern(acceptor)
            print(f"        Donor pattern: '{donor_pattern}'")
            print(f"        Acceptor pattern: '{acceptor_pattern}'")
    else:
        print("Modern: No hbond_list record")
    
    # Summary
    print(f"\n=== SUMMARY ===")
    if legacy_val and modern_val:
        legacy_valid = legacy_val.get('is_valid', 0)
        modern_valid = modern_val.get('is_valid', 0)
        if legacy_valid != modern_valid:
            print(f"❌ MISMATCH: legacy={legacy_valid}, modern={modern_valid}")
            if legacy_valid == 0 and modern_valid == 1:
                print(f"   Modern accepts but legacy rejects")
                print(f"   Legacy num_base_hb: {legacy_val.get('num_base_hb', 'N/A')}")
                print(f"   Modern num_base_hb: {modern_val.get('num_base_hb', 'N/A')}")
                print(f"   Legacy num_o2_hb: {legacy_val.get('num_o2_hb', 'N/A')}")
                print(f"   Modern num_o2_hb: {modern_val.get('num_o2_hb', 'N/A')}")
        else:
            print(f"✓ Match: both is_valid={legacy_valid}")

def main():
    if len(sys.argv) < 4:
        print("Usage: python3 debug_validation_mismatch.py <pdb_id> <base_i> <base_j>")
        print("Example: python3 debug_validation_mismatch.py 1VBY 20 21")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    base_i = int(sys.argv[2])
    base_j = int(sys.argv[3])
    
    project_root = Path(__file__).parent.parent
    analyze_pair(pdb_id, base_i, base_j, project_root)

if __name__ == '__main__':
    main()


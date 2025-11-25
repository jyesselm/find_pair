#!/usr/bin/env python3
"""Analyze missing pairs in base pair finding comparison"""

import json
import sys
from pathlib import Path

def load_pairs(json_file, key='pairs'):
    """Load pairs from JSON file"""
    with open(json_file) as f:
        data = json.load(f)
    if isinstance(data, list) and len(data) > 0:
        if key in data[0]:
            return set(tuple(sorted(p)) for p in data[0][key])
    return set()

def load_validation(json_file):
    """Load validation records"""
    with open(json_file) as f:
        return json.load(f)

def find_validation(pair, validations):
    """Find validation record for a pair"""
    for v in validations:
        v_pair = tuple(sorted([v.get('base_i'), v.get('base_j')]))
        if v_pair == tuple(sorted(pair)):
            return v
    return None

def analyze_pdb(pdb_id):
    """Analyze missing pairs for a PDB"""
    json_dir = Path('data/json')
    legacy_dir = Path('data/json_legacy')
    
    modern_pairs = load_pairs(json_dir / f'{pdb_id}_find_bestpair_selection.json')
    legacy_pairs = load_pairs(legacy_dir / f'{pdb_id}_find_bestpair_selection.json')
    
    missing = sorted(legacy_pairs - modern_pairs)
    extra = sorted(modern_pairs - legacy_pairs)
    
    print(f"\n{pdb_id}:")
    print(f"  Legacy pairs: {len(legacy_pairs)}")
    print(f"  Modern pairs: {len(modern_pairs)}")
    print(f"  Common: {len(modern_pairs & legacy_pairs)}")
    print(f"  Missing: {len(missing)}")
    print(f"  Extra: {len(extra)}")
    print(f"  Match rate: {100*len(modern_pairs & legacy_pairs)/len(legacy_pairs):.1f}%")
    
    if missing:
        print(f"\n  Missing pairs: {missing}")
        validations = load_validation(json_dir / f'{pdb_id}_pair_validation.json')
        for pair in missing:
            val = find_validation(pair, validations)
            if val:
                print(f"    {pair}: is_valid={val.get('is_valid')}, quality_score={val.get('calculated_values', {}).get('quality_score', 'N/A')}")
            else:
                print(f"    {pair}: NO VALIDATION RECORD")
    
    if extra:
        print(f"\n  Extra pairs: {extra}")

if __name__ == '__main__':
    pdbs = ['3G8T', '6CAQ', '3AVY']
    for pdb in pdbs:
        try:
            analyze_pdb(pdb)
        except FileNotFoundError as e:
            print(f"Error: {e}")
        except Exception as e:
            print(f"Error analyzing {pdb}: {e}")


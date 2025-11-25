#!/usr/bin/env python3
"""Investigate why specific pairs passed validation but weren't selected"""

import json
from pathlib import Path

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

def find_selected_pairs(json_file):
    """Load selected pairs"""
    with open(json_file) as f:
        data = json.load(f)
    if isinstance(data, list) and len(data) > 0:
        if 'pairs' in data[0]:
            return set(tuple(sorted(p)) for p in data[0]['pairs'])
    return set()

def analyze_pair(pdb_id, pair):
    """Analyze why a pair wasn't selected"""
    json_dir = Path('data/json')
    legacy_dir = Path('data/json_legacy')
    
    print(f"\n{pdb_id} Pair {pair}:")
    
    # Check if pair is in modern selection
    modern_pairs = find_selected_pairs(json_dir / f'{pdb_id}_find_bestpair_selection.json')
    legacy_pairs = find_selected_pairs(legacy_dir / f'{pdb_id}_find_bestpair_selection.json')
    
    is_in_modern = tuple(sorted(pair)) in modern_pairs
    is_in_legacy = tuple(sorted(pair)) in legacy_pairs
    
    print(f"  In modern selection: {is_in_modern}")
    print(f"  In legacy selection: {is_in_legacy}")
    
    # Get validation record
    validations = load_validation(json_dir / f'{pdb_id}_pair_validation.json')
    val = find_validation(pair, validations)
    
    if val:
        print(f"  Validation: is_valid={val.get('is_valid')}")
        calc_vals = val.get('calculated_values', {})
        print(f"  Quality score: {calc_vals.get('quality_score', 'N/A')}")
        print(f"  dorg: {calc_vals.get('dorg', 'N/A')}")
        print(f"  d_v: {calc_vals.get('d_v', 'N/A')}")
        print(f"  plane_angle: {calc_vals.get('plane_angle', 'N/A')}")
        print(f"  dNN: {calc_vals.get('dNN', 'N/A')}")
        print(f"  bp_type_id: {val.get('bp_type_id', 'N/A')}")
        
        # Check what pairs were selected for these residues
        i, j = pair
        i_pairs = [p for p in modern_pairs if i in p]
        j_pairs = [p for p in modern_pairs if j in p]
        
        print(f"\n  Residue {i} selected pairs: {i_pairs}")
        print(f"  Residue {j} selected pairs: {j_pairs}")
        
        # Check legacy selections
        legacy_i_pairs = [p for p in legacy_pairs if i in p]
        legacy_j_pairs = [p for p in legacy_pairs if j in p]
        
        print(f"  Legacy residue {i} selected pairs: {legacy_i_pairs}")
        print(f"  Legacy residue {j} selected pairs: {legacy_j_pairs}")
        
        # Get quality scores for competing pairs
        print(f"\n  Quality scores for competing pairs:")
        for comp_pair in i_pairs + j_pairs:
            if comp_pair != tuple(sorted(pair)):
                comp_val = find_validation(comp_pair, validations)
                if comp_val:
                    comp_qs = comp_val.get('calculated_values', {}).get('quality_score', 'N/A')
                    print(f"    {comp_pair}: {comp_qs}")
    else:
        print(f"  NO VALIDATION RECORD")

# Analyze the 4 missing pairs
missing_pairs = [
    ('3G8T', (92, 160)),
    ('3G8T', (946, 947)),
    ('6CAQ', (75, 78)),
    ('6CAQ', (968, 1024)),
]

for pdb_id, pair in missing_pairs:
    analyze_pair(pdb_id, pair)


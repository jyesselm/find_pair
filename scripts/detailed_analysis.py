#!/usr/bin/env python3
"""
Detailed analysis of specific PDB differences
"""

import json
import sys
from pathlib import Path

def analyze_3knc():
    """Detailed analysis of 3KNC - why no pairs found?"""
    project_root = Path(__file__).parent.parent
    
    print("="*80)
    print("DETAILED ANALYSIS: 3KNC (0 pairs found in modern)")
    print("="*80)
    
    # Check legacy base pairs
    legacy_bp_file = project_root / "data" / "json_legacy" / "3KNC_base_pair.json"
    if not legacy_bp_file.exists():
        legacy_bp_file = project_root / "data" / "json_legacy" / "3KNC.json"
    
    with open(legacy_bp_file) as f:
        data = json.load(f)
    
    records = data if isinstance(data, list) else data.get('calculations', [])
    legacy_bp = [r for r in records if r.get('type') == 'base_pair']
    
    print(f"\nLegacy found {len(legacy_bp)} base pairs:")
    for rec in legacy_bp[:5]:
        print(f"  ({rec['base_i']}, {rec['base_j']}): {rec.get('bp_type', 'N/A')}")
    
    # Check legacy validations
    legacy_val_file = project_root / "data" / "json_legacy" / "3KNC_pair_validation.json"
    if not legacy_val_file.exists():
        legacy_val_file = project_root / "data" / "json_legacy" / "3KNC.json"
    
    with open(legacy_val_file) as f:
        data = json.load(f)
    
    records = data if isinstance(data, list) else data.get('calculations', [])
    legacy_val = [r for r in records if r.get('type') == 'pair_validation']
    valid_ones = [r for r in legacy_val if r.get('is_valid') == 1]
    
    print(f"\nLegacy validations: {len(legacy_val)} total, {len(valid_ones)} valid")
    
    if valid_ones:
        print(f"\nFirst valid pair details:")
        v = valid_ones[0]
        print(f"  Pair: ({v.get('base_i')}, {v.get('base_j')})")
        calc = v.get('calculated_values', {})
        print(f"  dorg: {calc.get('dorg')}")
        print(f"  d_v: {calc.get('d_v')}")
        print(f"  plane_angle: {calc.get('plane_angle')}")
        print(f"  dNN: {calc.get('dNN')}")
        checks = v.get('validation_checks', {})
        print(f"  Checks:")
        print(f"    distance_check: {checks.get('distance_check')}")
        print(f"    d_v_check: {checks.get('d_v_check')}")
        print(f"    plane_angle_check: {checks.get('plane_angle_check')}")
        print(f"    dNN_check: {checks.get('dNN_check')}")
        thresholds = v.get('thresholds', {})
        print(f"  Thresholds:")
        print(f"    min_dorg: {thresholds.get('min_dorg')}, max_dorg: {thresholds.get('max_dorg')}")
        print(f"    min_dv: {thresholds.get('min_dv')}, max_dv: {thresholds.get('max_dv')}")
        print(f"    min_dNN: {thresholds.get('min_dNN')}, max_dNN: {thresholds.get('max_dNN')}")
        print(f"    min_plane_angle: {thresholds.get('min_plane_angle')}, max_plane_angle: {thresholds.get('max_plane_angle')}")
    
    # Check modern - should have 0
    modern_bp_file = project_root / "data" / "json" / "3KNC_base_pair.json"
    if not modern_bp_file.exists():
        modern_bp_file = project_root / "data" / "json" / "3KNC.json"
    
    if modern_bp_file.exists():
        with open(modern_bp_file) as f:
            data = json.load(f)
        
        records = data if isinstance(data, list) else data.get('calculations', [])
        modern_bp = [r for r in records if r.get('type') == 'base_pair']
        print(f"\nModern found {len(modern_bp)} base pairs (expected 0)")
        
        modern_val_file = project_root / "data" / "json" / "3KNC_pair_validation.json"
        if not modern_val_file.exists():
            modern_val_file = project_root / "data" / "json" / "3KNC.json"
        
        if modern_val_file.exists():
            with open(modern_val_file) as f:
                data = json.load(f)
            
            records = data if isinstance(data, list) else data.get('calculations', [])
            modern_val = [r for r in records if r.get('type') == 'pair_validation']
            print(f"Modern validations: {len(modern_val)} total")
            
            if modern_val:
                print(f"\nFirst modern validation (if any):")
                v = modern_val[0]
                print(f"  Pair: ({v.get('base_i')}, {v.get('base_j')})")
                print(f"  is_valid: {v.get('is_valid')}")
                calc = v.get('calculated_values', {})
                print(f"  dorg: {calc.get('dorg')}")
                print(f"  d_v: {calc.get('d_v')}")
                print(f"  plane_angle: {calc.get('plane_angle')}")
                print(f"  dNN: {calc.get('dNN')}")

def analyze_validation_mismatches():
    """Analyze validation mismatches across all PDBs"""
    project_root = Path(__file__).parent.parent
    
    test_set_file = project_root / "data" / "test_sets" / "test_set_10.json"
    with open(test_set_file) as f:
        test_set = json.load(f)
    
    pdb_ids = test_set.get('pdb_ids', [])
    
    print("\n" + "="*80)
    print("VALIDATION MISMATCHES ANALYSIS")
    print("="*80)
    
    for pdb_id in pdb_ids:
        legacy_val_file = project_root / "data" / "json_legacy" / f"{pdb_id}_pair_validation.json"
        modern_val_file = project_root / "data" / "json" / f"{pdb_id}_pair_validation.json"
        
        if not legacy_val_file.exists():
            legacy_val_file = project_root / "data" / "json_legacy" / f"{pdb_id}.json"
        if not modern_val_file.exists():
            modern_val_file = project_root / "data" / "json" / f"{pdb_id}.json"
        
        if not legacy_val_file.exists() or not modern_val_file.exists():
            continue
        
        # Load and compare
        with open(legacy_val_file) as f:
            data = json.load(f)
        records = data if isinstance(data, list) else data.get('calculations', [])
        legacy_val = {tuple(sorted([r['base_i'], r['base_j']])): r 
                     for r in records if r.get('type') == 'pair_validation'}
        
        with open(modern_val_file) as f:
            data = json.load(f)
        records = data if isinstance(data, list) else data.get('calculations', [])
        modern_val = {tuple(sorted([r['base_i'], r['base_j']])): r 
                     for r in records if r.get('type') == 'pair_validation'}
        
        mismatches = []
        for key in set(legacy_val.keys()) & set(modern_val.keys()):
            leg = legacy_val[key]
            mod = modern_val[key]
            if leg.get('is_valid') != mod.get('is_valid'):
                mismatches.append((key, leg, mod))
        
        if mismatches:
            print(f"\n{pdb_id}: {len(mismatches)} mismatches")
            for key, leg, mod in mismatches[:3]:
                print(f"  Pair {key}:")
                print(f"    Legacy: valid={leg.get('is_valid')}, bp_type_id={leg.get('bp_type_id')}")
                print(f"    Modern: valid={mod.get('is_valid')}, bp_type_id={mod.get('bp_type_id')}")
                leg_calc = leg.get('calculated_values', {})
                mod_calc = mod.get('calculated_values', {})
                print(f"    Legacy dorg={leg_calc.get('dorg')}, dNN={leg_calc.get('dNN')}")
                print(f"    Modern dorg={mod_calc.get('dorg')}, dNN={mod_calc.get('dNN')}")

def main():
    if len(sys.argv) > 1:
        if sys.argv[1] == '3knc':
            analyze_3knc()
        elif sys.argv[1] == 'mismatches':
            analyze_validation_mismatches()
        else:
            print(f"Unknown analysis: {sys.argv[1]}")
    else:
        analyze_3knc()
        analyze_validation_mismatches()

if __name__ == '__main__':
    main()


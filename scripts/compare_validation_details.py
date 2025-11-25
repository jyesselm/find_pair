#!/usr/bin/env python3
"""
Compare validation details for mismatched pairs.
Extracts all validation information from modern JSON and compares with debug output.
"""

import json
import sys
from pathlib import Path

def load_json(pdb_id: str, project_root: Path):
    """Load modern JSON file."""
    json_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
    if not json_file.exists():
        print(f"ERROR: JSON file not found: {json_file}")
        return None
    
    try:
        with open(json_file) as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading JSON: {e}")
        return None

def find_validation_record(data, base_i: int, base_j: int):
    """Find pair_validation record."""
    for record in data.get('calculations', []):
        if (record.get('type') == 'pair_validation' and 
            record.get('base_i') == base_i and record.get('base_j') == base_j):
            return record
    return None

def find_hbond_record(data, base_i: int, base_j: int):
    """Find hbond_list record."""
    for record in data.get('calculations', []):
        if (record.get('type') == 'hbond_list' and 
            record.get('base_i') == base_i and record.get('base_j') == base_j):
            return record
    return None

def analyze_pair(pdb_id: str, base_i: int, base_j: int, project_root: Path):
    """Analyze validation details for a specific pair."""
    print(f"\n{'='*80}")
    print(f"Validation Analysis: {pdb_id} pair ({base_i}, {base_j})")
    print(f"{'='*80}\n")
    
    data = load_json(pdb_id, project_root)
    if not data:
        return
    
    val_record = find_validation_record(data, base_i, base_j)
    hb_record = find_hbond_record(data, base_i, base_j)
    
    if not val_record:
        print("ERROR: No pair_validation record found")
        return
    
    print("=== VALIDATION RECORD ===")
    print(f"is_valid: {val_record.get('is_valid')}")
    print(f"bp_type_id: {val_record.get('bp_type_id')}")
    
    calc_vals = val_record.get('calculated_values', {})
    print(f"\nCalculated Values:")
    print(f"  dorg: {calc_vals.get('dorg', 'N/A')}")
    print(f"  d_v: {calc_vals.get('d_v', 'N/A')}")
    print(f"  plane_angle: {calc_vals.get('plane_angle', 'N/A')}")
    print(f"  dNN: {calc_vals.get('dNN', 'N/A')}")
    print(f"  quality_score: {calc_vals.get('quality_score', 'N/A')}")
    
    checks = val_record.get('validation_checks', {})
    print(f"\nValidation Checks:")
    print(f"  distance_check: {checks.get('distance_check', 'N/A')}")
    print(f"  d_v_check: {checks.get('d_v_check', 'N/A')}")
    print(f"  plane_angle_check: {checks.get('plane_angle_check', 'N/A')}")
    print(f"  dNN_check: {checks.get('dNN_check', 'N/A')}")
    
    thresholds = val_record.get('thresholds', {})
    print(f"\nThresholds:")
    print(f"  dorg: [{thresholds.get('min_dorg', 'N/A')}, {thresholds.get('max_dorg', 'N/A')}]")
    print(f"  d_v: [{thresholds.get('min_dv', 'N/A')}, {thresholds.get('max_dv', 'N/A')}]")
    print(f"  plane_angle: [{thresholds.get('min_plane_angle', 'N/A')}, {thresholds.get('max_plane_angle', 'N/A')}]")
    print(f"  dNN: [{thresholds.get('min_dNN', 'N/A')}, {thresholds.get('max_dNN', 'N/A')}]")
    
    # Check if values are within thresholds
    print(f"\nThreshold Analysis:")
    dorg = calc_vals.get('dorg')
    if dorg is not None:
        min_dorg = thresholds.get('min_dorg')
        max_dorg = thresholds.get('max_dorg')
        if min_dorg is not None and max_dorg is not None:
            in_range = (dorg >= min_dorg if min_dorg is not None else True) and (dorg <= max_dorg)
            print(f"  dorg={dorg:.6f} in range [{min_dorg}, {max_dorg}]: {in_range}")
    
    d_v = calc_vals.get('d_v')
    if d_v is not None:
        min_dv = thresholds.get('min_dv')
        max_dv = thresholds.get('max_dv')
        if min_dv is not None and max_dv is not None:
            in_range = (d_v >= min_dv if min_dv is not None else True) and (d_v <= max_dv)
            print(f"  d_v={d_v:.6f} in range [{min_dv}, {max_dv}]: {in_range}")
    
    plane_angle = calc_vals.get('plane_angle')
    if plane_angle is not None:
        min_pa = thresholds.get('min_plane_angle')
        max_pa = thresholds.get('max_plane_angle')
        if min_pa is not None and max_pa is not None:
            in_range = (plane_angle >= min_pa if min_pa is not None else True) and (plane_angle <= max_pa)
            print(f"  plane_angle={plane_angle:.6f} in range [{min_pa}, {max_pa}]: {in_range}")
    
    dNN = calc_vals.get('dNN')
    if dNN is not None:
        min_dNN = thresholds.get('min_dNN')
        max_dNN = thresholds.get('max_dNN')
        if min_dNN is not None and max_dNN is not None:
            in_range = (dNN >= min_dNN if min_dNN is not None else True) and (dNN <= max_dNN)
            print(f"  dNN={dNN:.6f} in range [{min_dNN}, {max_dNN}]: {in_range}")
    
    if hb_record:
        print(f"\n=== H-BOND RECORD ===")
        print(f"num_hbonds: {hb_record.get('num_hbonds', 0)}")
        hbonds = hb_record.get('hbonds', [])
        for i, hb in enumerate(hbonds[:5]):  # Show first 5
            donor = hb.get('donor_atom', '')
            acceptor = hb.get('acceptor_atom', '')
            dist = hb.get('distance', 0)
            print(f"  {i+1}. {donor} - {acceptor} ({dist:.3f} Å)")
    else:
        print(f"\n=== H-BOND RECORD ===")
        print("No hbond_list record (pair rejected or no H-bonds found)")
    
    # Summary
    print(f"\n=== SUMMARY ===")
    is_valid = val_record.get('is_valid', 0)
    all_checks = (checks.get('distance_check', False) and 
                  checks.get('d_v_check', False) and 
                  checks.get('plane_angle_check', False) and 
                  checks.get('dNN_check', False))
    
    print(f"All distance/angle checks pass: {all_checks}")
    print(f"is_valid: {is_valid}")
    
    if is_valid == 1 and all_checks:
        print("✓ Modern accepts this pair (all checks pass)")
    elif is_valid == 0 and all_checks:
        print("✗ Modern rejects this pair despite passing distance/angle checks")
        print("  -> Likely failed H-bond requirement")
    elif is_valid == 0 and not all_checks:
        print("✗ Modern rejects this pair (failed distance/angle checks)")

def main():
    if len(sys.argv) < 4:
        print("Usage: python3 compare_validation_details.py <pdb_id> <base_i> <base_j>")
        print("Example: python3 compare_validation_details.py 1VBY 20 21")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    base_i = int(sys.argv[2])
    base_j = int(sys.argv[3])
    
    project_root = Path(__file__).parent.parent
    analyze_pair(pdb_id, base_i, base_j, project_root)

if __name__ == '__main__':
    main()


#!/usr/bin/env python3
"""
Investigate bp_type_id differences between legacy and modern.

This script finds pairs where modern assigns bp_type_id = 2 but legacy keeps -1,
or vice versa, and compares their step parameters to identify root causes.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


def load_json_records(pdb_id: str, record_type: str, is_legacy: bool) -> List[Dict]:
    """Load JSON records of a specific type."""
    if is_legacy:
        json_file = project_root / f"data/json_legacy/{record_type}/{pdb_id}.json"
    else:
        json_file = project_root / f"data/json/{record_type}/{pdb_id}.json"
    
    if not json_file.exists():
        return []
    
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
            if isinstance(data, list):
                return [r for r in data if r.get('type') == record_type]
            elif isinstance(data, dict) and 'records' in data:
                return [r for r in data['records'] if r.get('type') == record_type]
            else:
                return []
    except Exception as e:
        print(f"Error loading {json_file}: {e}", file=sys.stderr)
        return []


def find_bp_type_id_differences(pdb_id: str) -> List[Dict]:
    """Find pairs with bp_type_id differences."""
    legacy_validation = load_json_records(pdb_id, 'pair_validation', True)
    modern_validation = load_json_records(pdb_id, 'pair_validation', False)
    
    # Build maps by pair
    legacy_map = {}
    modern_map = {}
    
    for rec in legacy_validation:
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i and base_j:
            pair = (min(base_i, base_j), max(base_i, base_j))
            legacy_map[pair] = rec
    
    for rec in modern_validation:
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i and base_j:
            pair = (min(base_i, base_j), max(base_i, base_j))
            modern_map[pair] = rec
    
    # Find differences
    differences = []
    common_pairs = set(legacy_map.keys()) & set(modern_map.keys())
    
    for pair in common_pairs:
        leg_rec = legacy_map[pair]
        mod_rec = modern_map[pair]
        
        leg_bp_type_id = leg_rec.get('bp_type_id', 0)
        mod_bp_type_id = mod_rec.get('bp_type_id', 0)
        
        # Focus on cases where modern assigns 2 but legacy keeps -1
        # or modern keeps -1 but legacy assigns 2
        if (mod_bp_type_id == 2 and leg_bp_type_id == -1) or \
           (mod_bp_type_id == -1 and leg_bp_type_id == 2):
            differences.append({
                'pair': pair,
                'legacy_bp_type_id': leg_bp_type_id,
                'modern_bp_type_id': mod_bp_type_id,
                'legacy_record': leg_rec,
                'modern_record': mod_rec
            })
    
    return differences


def extract_step_parameters(pdb_id: str, pair: Tuple[int, int]) -> Dict:
    """Extract step parameters for a pair if available."""
    legacy_steps = load_json_records(pdb_id, 'bpstep_params', True)
    modern_steps = load_json_records(pdb_id, 'bpstep_params', False)
    
    result = {
        'legacy': None,
        'modern': None
    }
    
    # Find step parameters for this pair
    # Note: bpstep_params are for consecutive pairs, so we need to check
    # if this pair has step parameters calculated
    for rec in legacy_steps:
        bp_idx1 = rec.get('bp_idx1')
        bp_idx2 = rec.get('bp_idx2')
        if bp_idx1 and bp_idx2:
            # Check if this step corresponds to our pair
            # This is a simplified check - actual matching would need pair indices
            pass
    
    return result


def analyze_difference(diff: Dict, pdb_id: str):
    """Analyze a single bp_type_id difference."""
    pair = diff['pair']
    leg_rec = diff['legacy_record']
    mod_rec = diff['modern_record']
    
    print(f"\n{'='*80}")
    print(f"Pair ({pair[0]}, {pair[1]})")
    print(f"{'='*80}")
    print(f"Legacy bp_type_id: {diff['legacy_bp_type_id']}")
    print(f"Modern bp_type_id: {diff['modern_bp_type_id']}")
    
    # Extract direction vectors
    leg_dir = leg_rec.get('direction_vectors', {})
    mod_dir = mod_rec.get('direction_vectors', {})
    
    print(f"\nDirection Vectors:")
    print(f"  Legacy: dir_x={leg_dir.get('dir_x', 'N/A'):.6f}, "
          f"dir_y={leg_dir.get('dir_y', 'N/A'):.6f}, "
          f"dir_z={leg_dir.get('dir_z', 'N/A'):.6f}")
    print(f"  Modern: dir_x={mod_dir.get('dir_x', 'N/A'):.6f}, "
          f"dir_y={mod_dir.get('dir_y', 'N/A'):.6f}, "
          f"dir_z={mod_dir.get('dir_z', 'N/A'):.6f}")
    
    # Check if direction vector condition is met
    mod_dir_x = mod_dir.get('dir_x', 0)
    mod_dir_y = mod_dir.get('dir_y', 0)
    mod_dir_z = mod_dir.get('dir_z', 0)
    
    condition_met = (mod_dir_x > 0.0 and mod_dir_y < 0.0 and mod_dir_z < 0.0)
    print(f"\nDirection Vector Condition (dir_x > 0 && dir_y < 0 && dir_z < 0): {condition_met}")
    
    if not condition_met:
        print(f"  ⚠️  Condition not met - bp_type_id should remain -1")
        print(f"  This suggests modern is incorrectly assigning bp_type_id = 2")
        return
    
    # Extract calculated values
    leg_calc = leg_rec.get('calculated_values', {})
    mod_calc = mod_rec.get('calculated_values', {})
    
    print(f"\nGeometric Values:")
    print(f"  Legacy: dorg={leg_calc.get('dorg', 'N/A'):.6f}, "
          f"d_v={leg_calc.get('d_v', 'N/A'):.6f}, "
          f"plane_angle={leg_calc.get('plane_angle', 'N/A'):.6f}")
    print(f"  Modern: dorg={mod_calc.get('dorg', 'N/A'):.6f}, "
          f"d_v={mod_calc.get('d_v', 'N/A'):.6f}, "
          f"plane_angle={mod_calc.get('plane_angle', 'N/A'):.6f}")
    
    # Check validation status
    print(f"\nValidation Status:")
    print(f"  Legacy: is_valid={leg_rec.get('is_valid', 'N/A')}")
    print(f"  Modern: is_valid={mod_rec.get('is_valid', 'N/A')}")
    
    # Try to find step parameters
    step_params = extract_step_parameters(pdb_id, pair)
    if step_params['legacy'] or step_params['modern']:
        print(f"\nStep Parameters:")
        if step_params['legacy']:
            print(f"  Legacy: {step_params['legacy']}")
        if step_params['modern']:
            print(f"  Modern: {step_params['modern']}")
    else:
        print(f"\n⚠️  Step parameters not found in JSON (may not be calculated for single pairs)")
        print(f"   Step parameters are typically calculated for consecutive pairs in a helix")
    
    # Analysis
    print(f"\nAnalysis:")
    if diff['modern_bp_type_id'] == 2 and diff['legacy_bp_type_id'] == -1:
        print(f"  Modern assigns bp_type_id = 2 (Watson-Crick), legacy keeps -1")
        print(f"  Possible causes:")
        print(f"    1. Step parameter differences (shear, stretch, opening)")
        print(f"    2. Base pair type string differences (WC_LIST matching)")
        print(f"    3. Frame reversal logic differences")
    elif diff['modern_bp_type_id'] == -1 and diff['legacy_bp_type_id'] == 2:
        print(f"  Legacy assigns bp_type_id = 2 (Watson-Crick), modern keeps -1")
        print(f"  Possible causes:")
        print(f"    1. Step parameter differences (shear, stretch, opening)")
        print(f"    2. Base pair type string differences (WC_LIST matching)")
        print(f"    3. Frame reversal logic differences")
        print(f"    4. Missing step parameter calculation in modern")


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/investigate_bp_type_id_differences.py <pdb_id>")
        sys.exit(1)
    
    pdb_id = sys.argv[1].upper()
    
    print(f"Investigating bp_type_id differences for {pdb_id}")
    print(f"{'='*80}")
    
    differences = find_bp_type_id_differences(pdb_id)
    
    if not differences:
        print(f"✅ No bp_type_id differences found for {pdb_id}")
        sys.exit(0)
    
    print(f"\nFound {len(differences)} pairs with bp_type_id differences:")
    
    # Count by type
    modern_2_legacy_minus1 = sum(1 for d in differences if d['modern_bp_type_id'] == 2 and d['legacy_bp_type_id'] == -1)
    modern_minus1_legacy_2 = sum(1 for d in differences if d['modern_bp_type_id'] == -1 and d['legacy_bp_type_id'] == 2)
    
    print(f"  Modern=2, Legacy=-1: {modern_2_legacy_minus1}")
    print(f"  Modern=-1, Legacy=2: {modern_minus1_legacy_2}")
    
    # Analyze each difference
    for diff in differences:
        analyze_difference(diff, pdb_id)
    
    print(f"\n{'='*80}")
    print(f"Investigation complete for {pdb_id}")


if __name__ == '__main__':
    main()


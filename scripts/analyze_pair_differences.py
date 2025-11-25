#!/usr/bin/env python3
"""
Analyze pair validation and distance checks differences between legacy and modern JSON.

Provides detailed per-PDB breakdown of differences.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List
from collections import defaultdict

from x3dna_json_compare import JsonComparator, ComparisonResult


def load_json_file(json_file: Path):
    """Load JSON file, handling both objects and arrays."""
    if not json_file.exists():
        return None
    
    try:
        with open(json_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        return None


def count_records_in_file(json_file: Path, record_type: str) -> int:
    """Count records of a specific type in a JSON file."""
    data = load_json_file(json_file)
    if data is None:
        return 0
    
    count = 0
    if isinstance(data, list):
        for record in data:
            if record.get('type') == record_type:
                count += 1
    elif isinstance(data, dict):
        if 'calculation_records' in data:
            for record in data['calculation_records']:
                if record.get('type') == record_type:
                    count += 1
    
    return count


def analyze_pdb_pair_data(pdb_id: str, project_root: Path) -> Dict:
    """Analyze pair validation and distance checks for a single PDB."""
    result = {
        'pdb_id': pdb_id,
        'legacy': {},
        'modern': {},
        'comparison': {}
    }
    
    legacy_dir = project_root / 'data' / 'json_legacy'
    modern_dir = project_root / 'data' / 'json'
    
    # Count legacy records
    legacy_pair_val_file = legacy_dir / f"{pdb_id}_pair_validation.json"
    legacy_dist_file = legacy_dir / f"{pdb_id}_distance_checks.json"
    
    result['legacy']['pair_validation_count'] = count_records_in_file(legacy_pair_val_file, 'pair_validation')
    result['legacy']['distance_checks_count'] = count_records_in_file(legacy_dist_file, 'distance_checks')
    
    # Count modern records
    modern_pair_val_file = modern_dir / f"{pdb_id}_pair_validation.json"
    modern_dist_file = modern_dir / f"{pdb_id}_distance_checks.json"
    
    result['modern']['pair_validation_count'] = count_records_in_file(modern_pair_val_file, 'pair_validation')
    result['modern']['distance_checks_count'] = count_records_in_file(modern_dist_file, 'distance_checks')
    
    # Run comparison
    legacy_file = legacy_dir / f"{pdb_id}.json"
    modern_file = modern_dir / f"{pdb_id}.json"
    pdb_file = project_root / 'data' / 'pdb' / f"{pdb_id}.pdb"
    
    comparator = JsonComparator(compare_pairs=True, compare_atoms=False, compare_frames=False, compare_steps=False)
    comparison_result = comparator.compare_files(legacy_file, modern_file, pdb_file, pdb_id)
    
    if comparison_result.pair_validation_comparison:
        pvc = comparison_result.pair_validation_comparison
        result['comparison']['pair_validation'] = {
            'total_legacy': pvc.total_legacy,
            'total_modern': pvc.total_modern,
            'common': pvc.common_count,
            'missing_in_modern': len(pvc.missing_in_modern),
            'extra_in_modern': len(pvc.extra_in_modern),
            'mismatched': len(pvc.mismatched_validations)
        }
    
    if comparison_result.distance_checks_comparison:
        dcc = comparison_result.distance_checks_comparison
        result['comparison']['distance_checks'] = {
            'total_legacy': dcc.total_legacy,
            'total_modern': dcc.total_modern,
            'common': dcc.common_count,
            'missing_in_modern': len(dcc.missing_in_modern),
            'extra_in_modern': len(dcc.extra_in_modern),
            'mismatched': len(dcc.mismatched_checks)
        }
    
    # Analyze common pairs for mismatches
    if comparison_result.pair_validation_comparison and comparison_result.pair_validation_comparison.mismatched_validations:
        mismatches = comparison_result.pair_validation_comparison.mismatched_validations
        result['comparison']['pair_validation']['mismatch_details'] = []
        
        # Sample first few mismatches
        for mismatch in mismatches[:5]:
            mismatch_info = {
                'base_i': mismatch['base_i'],
                'base_j': mismatch['base_j'],
                'mismatched_fields': list(mismatch['mismatches'].keys())
            }
            result['comparison']['pair_validation']['mismatch_details'].append(mismatch_info)
    
    return result


def main():
    project_root = Path(__file__).parent.parent
    
    # Load test set
    test_set_file = project_root / 'data' / 'test_sets' / 'test_set_10.json'
    if not test_set_file.exists():
        print(f"Error: Test set file not found: {test_set_file}")
        sys.exit(1)
    
    with open(test_set_file) as f:
        test_set = json.load(f)
    
    pdb_ids = test_set['pdb_ids']
    
    print("=" * 100)
    print("PAIR VALIDATION & DISTANCE CHECKS ANALYSIS - 10 PDB TEST SET")
    print("=" * 100)
    print()
    
    all_results = []
    total_stats = {
        'legacy_pair_val': 0,
        'legacy_dist_checks': 0,
        'modern_pair_val': 0,
        'modern_dist_checks': 0,
        'common_pair_val': 0,
        'common_dist_checks': 0,
        'missing_pair_val': 0,
        'extra_pair_val': 0,
        'mismatched_pair_val': 0,
        'missing_dist_checks': 0,
        'extra_dist_checks': 0,
        'mismatched_dist_checks': 0,
    }
    
    for pdb_id in pdb_ids:
        print(f"Analyzing {pdb_id}...")
        result = analyze_pdb_pair_data(pdb_id, project_root)
        all_results.append(result)
        
        # Accumulate stats
        if 'pair_validation' in result['comparison']:
            pv = result['comparison']['pair_validation']
            total_stats['legacy_pair_val'] += pv['total_legacy']
            total_stats['modern_pair_val'] += pv['total_modern']
            total_stats['common_pair_val'] += pv['common']
            total_stats['missing_pair_val'] += pv['missing_in_modern']
            total_stats['extra_pair_val'] += pv['extra_in_modern']
            total_stats['mismatched_pair_val'] += pv['mismatched']
        
        if 'distance_checks' in result['comparison']:
            dc = result['comparison']['distance_checks']
            total_stats['legacy_dist_checks'] += dc['total_legacy']
            total_stats['modern_dist_checks'] += dc['total_modern']
            total_stats['common_dist_checks'] += dc['common']
            total_stats['missing_dist_checks'] += dc['missing_in_modern']
            total_stats['extra_dist_checks'] += dc['extra_in_modern']
            total_stats['mismatched_dist_checks'] += dc['mismatched']
    
    print()
    print("=" * 100)
    print("SUMMARY STATISTICS")
    print("=" * 100)
    print()
    
    print("PAIR VALIDATION:")
    print(f"  Legacy total:     {total_stats['legacy_pair_val']:>10,}")
    print(f"  Modern total:     {total_stats['modern_pair_val']:>10,}")
    print(f"  Common pairs:     {total_stats['common_pair_val']:>10,}")
    print(f"  Missing in modern: {total_stats['missing_pair_val']:>8,}")
    print(f"  Extra in modern:   {total_stats['extra_pair_val']:>8,}")
    print(f"  Mismatched:        {total_stats['mismatched_pair_val']:>8,}")
    print()
    
    print("DISTANCE CHECKS:")
    print(f"  Legacy total:     {total_stats['legacy_dist_checks']:>10,}")
    print(f"  Modern total:     {total_stats['modern_dist_checks']:>10,}")
    print(f"  Common pairs:     {total_stats['common_dist_checks']:>10,}")
    print(f"  Missing in modern: {total_stats['missing_dist_checks']:>8,}")
    print(f"  Extra in modern:   {total_stats['extra_dist_checks']:>8,}")
    print(f"  Mismatched:        {total_stats['mismatched_dist_checks']:>8,}")
    print()
    
    print("=" * 100)
    print("PER-PDB BREAKDOWN")
    print("=" * 100)
    print()
    
    # Sort by total differences
    all_results.sort(key=lambda x: (
        x['comparison'].get('pair_validation', {}).get('mismatched', 0) +
        x['comparison'].get('distance_checks', {}).get('mismatched', 0)
    ), reverse=True)
    
    for result in all_results:
        pdb_id = result['pdb_id']
        print(f"{pdb_id}:")
        
        if 'pair_validation' in result['comparison']:
            pv = result['comparison']['pair_validation']
            print(f"  Pair Validation:")
            print(f"    Legacy: {pv['total_legacy']:>6,}  Modern: {pv['total_modern']:>6,}  Common: {pv['common']:>6,}")
            print(f"    Missing: {pv['missing_in_modern']:>5,}  Extra: {pv['extra_in_modern']:>7,}  Mismatched: {pv['mismatched']:>4,}")
            
            if pv['mismatched'] > 0 and 'mismatch_details' in pv:
                print(f"    Sample mismatches:")
                for mm in pv['mismatch_details'][:3]:
                    print(f"      Base pair ({mm['base_i']}, {mm['base_j']}): {', '.join(mm['mismatched_fields'][:3])}")
        
        if 'distance_checks' in result['comparison']:
            dc = result['comparison']['distance_checks']
            print(f"  Distance Checks:")
            print(f"    Legacy: {dc['total_legacy']:>6,}  Modern: {dc['total_modern']:>6,}  Common: {dc['common']:>6,}")
            print(f"    Missing: {dc['missing_in_modern']:>5,}  Extra: {dc['extra_in_modern']:>7,}  Mismatched: {dc['mismatched']:>4,}")
        
        print()
    
    print("=" * 100)
    print("ANALYSIS")
    print("=" * 100)
    print()
    print("Key Observations:")
    print()
    
    # Calculate ratios
    if total_stats['legacy_pair_val'] > 0:
        ratio = total_stats['modern_pair_val'] / total_stats['legacy_pair_val']
        print(f"1. Modern code checks {ratio:.1f}x more pairs than legacy")
        print(f"   - This is because modern code validates ALL residue pairs during finding")
        print(f"   - Legacy code only validates pairs that pass initial distance/angle checks")
        print()
    
    if total_stats['common_pair_val'] > 0:
        match_rate = (total_stats['common_pair_val'] / total_stats['legacy_pair_val']) * 100
        print(f"2. Common pairs: {match_rate:.1f}% of legacy pairs are also in modern")
        print()
    
    if total_stats['mismatched_pair_val'] > 0:
        mismatch_rate = (total_stats['mismatched_pair_val'] / total_stats['common_pair_val']) * 100 if total_stats['common_pair_val'] > 0 else 0
        print(f"3. Mismatch rate: {mismatch_rate:.2f}% of common pairs have differences")
        print(f"   - These differences may be due to:")
        print(f"     * Floating point precision differences")
        print(f"     * Different calculation order")
        print(f"     * Slight algorithm differences")
        print()
    
    print("4. Missing in modern:")
    print(f"   - {total_stats['missing_pair_val']} pair validations from legacy not found in modern")
    print(f"   - These are pairs that legacy checked but modern didn't (or vice versa)")
    print()
    
    print("5. Extra in modern:")
    print(f"   - {total_stats['extra_pair_val']} pair validations in modern not in legacy")
    print(f"   - Modern code checks more pairs during the finding process")
    print()


if __name__ == '__main__':
    main()


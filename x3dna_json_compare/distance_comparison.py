"""
Comparison functions for distance_checks JSON records.

Compares geometric distance and angle measurements between base pairs.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional


@dataclass
class DistanceComparison:
    """Result of comparing distance_checks records."""
    total_legacy: int = 0
    total_modern: int = 0
    matched: int = 0
    mismatched_values: List[Dict] = field(default_factory=list)
    missing_in_modern: List[Tuple[int, int]] = field(default_factory=list)
    extra_in_modern: List[Tuple[int, int]] = field(default_factory=list)


def normalize_pair_key(base_i: int, base_j: int) -> Tuple[int, int]:
    """Normalize pair key to (min, max) for consistent comparison."""
    return (min(base_i, base_j), max(base_i, base_j))


def compare_distance_checks(
    legacy_records: List[Dict],
    modern_records: List[Dict],
    tolerance: float = 1e-6
) -> DistanceComparison:
    """
    Compare distance_checks records between legacy and modern.
    
    Matching Strategy:
    - Match by (base_i, base_j) - order normalized to (min, max)
    - Compare: dorg, dNN, plane_angle, d_v, overlap_area
    - All values must match within tolerance
    
    Args:
        legacy_records: List of legacy distance_checks records
        modern_records: List of modern distance_checks records
        tolerance: Numerical tolerance for floating point comparisons (default: 1e-6)
    
    Returns:
        DistanceComparison object with detailed results
    """
    result = DistanceComparison()
    
    result.total_legacy = len(legacy_records)
    result.total_modern = len(modern_records)
    
    # Build dictionaries keyed by normalized (base_i, base_j)
    legacy_dict = {}
    for record in legacy_records:
        base_i = record.get('base_i')
        base_j = record.get('base_j')
        if base_i is not None and base_j is not None:
            key = normalize_pair_key(base_i, base_j)
            legacy_dict[key] = record
    
    modern_dict = {}
    for record in modern_records:
        base_i = record.get('base_i')
        base_j = record.get('base_j')
        if base_i is not None and base_j is not None:
            key = normalize_pair_key(base_i, base_j)
            modern_dict[key] = record
    
    # Find common, missing, and extra pairs
    legacy_keys = set(legacy_dict.keys())
    modern_keys = set(modern_dict.keys())
    
    common_keys = legacy_keys & modern_keys
    result.missing_in_modern = list(legacy_keys - modern_keys)
    result.extra_in_modern = list(modern_keys - legacy_keys)
    
    # Compare common pairs
    for key in common_keys:
        legacy_rec = legacy_dict[key]
        modern_rec = modern_dict[key]
        
        # Get values dictionaries
        legacy_values = legacy_rec.get('values', {})
        modern_values = modern_rec.get('values', {})
        
        # Fields to compare
        fields_to_compare = ['dorg', 'dNN', 'plane_angle', 'd_v', 'overlap_area']
        
        pair_matched = True
        mismatches = {}
        
        for field in fields_to_compare:
            legacy_val = legacy_values.get(field)
            modern_val = modern_values.get(field)
            
            if legacy_val is None and modern_val is None:
                continue  # Both missing, OK
            
            if legacy_val is None or modern_val is None:
                mismatches[field] = {
                    'legacy': legacy_val,
                    'modern': modern_val,
                    'reason': 'One value missing'
                }
                pair_matched = False
                continue
            
            # Compare with tolerance
            diff = abs(legacy_val - modern_val)
            if diff > tolerance:
                mismatches[field] = {
                    'legacy': legacy_val,
                    'modern': modern_val,
                    'diff': diff,
                    'tolerance': tolerance
                }
                pair_matched = False
        
        if pair_matched:
            result.matched += 1
        else:
            result.mismatched_values.append({
                'pair': key,
                'base_i': key[0],
                'base_j': key[1],
                'mismatches': mismatches
            })
    
    return result


def print_distance_comparison_summary(comparison: DistanceComparison, verbose: bool = False):
    """
    Print a summary of the distance_checks comparison.
    
    Args:
        comparison: DistanceComparison result
        verbose: If True, print detailed mismatches
    """
    print(f"\nDistance Checks Comparison Summary:")
    print(f"  Total legacy records: {comparison.total_legacy}")
    print(f"  Total modern records: {comparison.total_modern}")
    print(f"  Matched pairs: {comparison.matched}")
    print(f"  Missing in modern: {len(comparison.missing_in_modern)}")
    print(f"  Extra in modern: {len(comparison.extra_in_modern)}")
    print(f"  Mismatched values: {len(comparison.mismatched_values)}")
    
    if comparison.total_legacy > 0:
        match_rate = 100.0 * comparison.matched / comparison.total_legacy
        print(f"  Match rate: {match_rate:.2f}%")
    
    if verbose and comparison.mismatched_values:
        print(f"\nDetailed mismatches:")
        for mismatch in comparison.mismatched_values[:10]:  # Show first 10
            pair = mismatch['pair']
            print(f"\n  Pair ({pair[0]}, {pair[1]}):")
            for field, details in mismatch['mismatches'].items():
                print(f"    {field}:")
                print(f"      Legacy: {details.get('legacy')}")
                print(f"      Modern: {details.get('modern')}")
                if 'diff' in details:
                    print(f"      Diff: {details['diff']}")
    
    if verbose and comparison.missing_in_modern:
        print(f"\nMissing in modern (first 10):")
        for pair in comparison.missing_in_modern[:10]:
            print(f"  ({pair[0]}, {pair[1]})")
    
    if verbose and comparison.extra_in_modern:
        print(f"\nExtra in modern (first 10):")
        for pair in comparison.extra_in_modern[:10]:
            print(f"  ({pair[0]}, {pair[1]})")


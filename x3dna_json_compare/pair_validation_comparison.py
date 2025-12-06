"""
Comparison functions for pair_validation JSON records.

Compares validation results for each tested base pair.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional


@dataclass
class PairValidationComparison:
    """Result of comparing pair_validation records."""
    total_legacy: int = 0
    total_modern: int = 0
    matched: int = 0
    validation_mismatches: List[Dict] = field(default_factory=list)
    value_mismatches: List[Dict] = field(default_factory=list)
    missing_in_modern: List[Tuple[int, int]] = field(default_factory=list)
    extra_in_modern: List[Tuple[int, int]] = field(default_factory=list)


def normalize_pair_key(base_i: int, base_j: int) -> Tuple[int, int]:
    """Normalize pair key to (min, max) for consistent comparison."""
    return (min(base_i, base_j), max(base_i, base_j))


def compare_pair_validation(
    legacy_records: List[Dict],
    modern_records: List[Dict],
    tolerance: float = 1e-6
) -> PairValidationComparison:
    """
    Compare pair_validation records between legacy and modern.
    
    IMPORTANT: Legacy uses 1-based indices, modern uses 0-based indices.
    This function converts modern indices to 1-based for comparison.
    
    Fields compared:
    1. is_valid: 0 or 1 (exact)
    2. bp_type_id: -1, 0, 1, 2 (exact)
    3. Direction vectors: dir_x, dir_y, dir_z (±tolerance)
    4. Calculated values: dorg, d_v, plane_angle, dNN, quality_score (±tolerance)
    5. Validation checks: distance_check, d_v_check, plane_angle_check, dNN_check (exact)
    
    Args:
        legacy_records: List of legacy pair_validation records (1-based indices)
        modern_records: List of modern pair_validation records (0-based indices)
        tolerance: Numerical tolerance for floating point comparisons (default: 1e-6)
    
    Returns:
        PairValidationComparison object with detailed results
    """
    result = PairValidationComparison()
    
    result.total_legacy = len(legacy_records)
    result.total_modern = len(modern_records)
    
    # Build dictionaries keyed by normalized (base_i, base_j)
    # Legacy uses 1-based indices
    legacy_dict = {}
    for record in legacy_records:
        base_i = record.get('base_i')
        base_j = record.get('base_j')
        if base_i is not None and base_j is not None:
            key = normalize_pair_key(base_i, base_j)
            legacy_dict[key] = record
    
    # Modern now uses 1-based indices (same as legacy)
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
        
        pair_matched = True
        validation_diffs = {}
        value_diffs = {}
        
        # Compare is_valid (exact)
        if legacy_rec.get('is_valid') != modern_rec.get('is_valid'):
            validation_diffs['is_valid'] = {
                'legacy': legacy_rec.get('is_valid'),
                'modern': modern_rec.get('is_valid')
            }
            pair_matched = False
        
        # Compare bp_type_id (exact)
        if legacy_rec.get('bp_type_id') != modern_rec.get('bp_type_id'):
            validation_diffs['bp_type_id'] = {
                'legacy': legacy_rec.get('bp_type_id'),
                'modern': modern_rec.get('bp_type_id')
            }
            pair_matched = False
        
        # Compare direction vectors (with tolerance)
        dir_vecs = legacy_rec.get('direction_vectors', {})
        mod_dir_vecs = modern_rec.get('direction_vectors', {})
        
        for dir_key in ['dir_x', 'dir_y', 'dir_z']:
            leg_val = dir_vecs.get(dir_key, 0.0)
            mod_val = mod_dir_vecs.get(dir_key, 0.0)
            if abs(leg_val - mod_val) > tolerance:
                value_diffs[f'direction_vectors.{dir_key}'] = {
                    'legacy': leg_val,
                    'modern': mod_val,
                    'diff': abs(leg_val - mod_val)
                }
                pair_matched = False
        
        # Compare calculated values (with tolerance)
        calc_vals = legacy_rec.get('calculated_values', {})
        mod_calc_vals = modern_rec.get('calculated_values', {})
        
        for calc_key in ['dorg', 'd_v', 'plane_angle', 'dNN', 'quality_score']:
            leg_val = calc_vals.get(calc_key, 0.0)
            mod_val = mod_calc_vals.get(calc_key, 0.0)
            if abs(leg_val - mod_val) > tolerance:
                value_diffs[f'calculated_values.{calc_key}'] = {
                    'legacy': leg_val,
                    'modern': mod_val,
                    'diff': abs(leg_val - mod_val)
                }
                pair_matched = False
        
        # Compare validation checks (exact)
        val_checks = legacy_rec.get('validation_checks', {})
        mod_val_checks = modern_rec.get('validation_checks', {})
        
        for check_key in ['distance_check', 'd_v_check', 'plane_angle_check', 'dNN_check']:
            if val_checks.get(check_key) != mod_val_checks.get(check_key):
                validation_diffs[f'validation_checks.{check_key}'] = {
                    'legacy': val_checks.get(check_key),
                    'modern': mod_val_checks.get(check_key)
                }
                pair_matched = False
        
        if validation_diffs:
            result.validation_mismatches.append({
                'pair': key,
                'base_i': key[0],
                'base_j': key[1],
                'mismatches': validation_diffs
            })
        
        if value_diffs:
            result.value_mismatches.append({
                'pair': key,
                'base_i': key[0],
                'base_j': key[1],
                'mismatches': value_diffs
            })
        
        if pair_matched:
            result.matched += 1
    
    return result


def print_pair_validation_summary(comparison: PairValidationComparison, verbose: bool = False):
    """
    Print a summary of the pair_validation comparison.
    
    Args:
        comparison: PairValidationComparison result
        verbose: If True, print detailed mismatches
    """
    print(f"\nPair Validation Comparison Summary:")
    print(f"  Total legacy records: {comparison.total_legacy}")
    print(f"  Total modern records: {comparison.total_modern}")
    print(f"  Matched pairs: {comparison.matched}")
    print(f"  Missing in modern: {len(comparison.missing_in_modern)}")
    print(f"  Extra in modern: {len(comparison.extra_in_modern)}")
    print(f"  Validation mismatches: {len(comparison.validation_mismatches)}")
    print(f"  Value mismatches: {len(comparison.value_mismatches)}")
    
    if comparison.total_legacy > 0:
        match_rate = 100.0 * comparison.matched / comparison.total_legacy
        print(f"  Match rate: {match_rate:.2f}%")
    
    if verbose and comparison.validation_mismatches:
        print(f"\nValidation mismatches (first 10):")
        for item in comparison.validation_mismatches[:10]:
            pair = item['pair']
            print(f"  Pair ({pair[0]}, {pair[1]}):")
            for field, details in item['mismatches'].items():
                print(f"    {field}: legacy={details['legacy']}, modern={details['modern']}")
    
    if verbose and comparison.value_mismatches:
        print(f"\nValue mismatches (first 10):")
        for item in comparison.value_mismatches[:10]:
            pair = item['pair']
            print(f"  Pair ({pair[0]}, {pair[1]}):")
            for field, details in item['mismatches'].items():
                print(f"    {field}: {details['legacy']:.6f} vs {details['modern']:.6f} (diff={details['diff']:.2e})")


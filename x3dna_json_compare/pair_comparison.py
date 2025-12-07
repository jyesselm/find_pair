"""
Base pair validation and distance checks comparison utilities.

Provides functions to compare pair_validation and distance_checks records
between legacy and modern JSON outputs.
"""

from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field


def normalize_pair_key(base_i: int, base_j: int) -> Tuple[int, int]:
    """Normalize pair key to (min, max) for consistent comparison."""
    return (min(base_i, base_j), max(base_i, base_j))


@dataclass
class PairValidationComparison:
    """Result of pair_validation comparison."""
    missing_in_modern: List[Dict] = field(default_factory=list)
    extra_in_modern: List[Dict] = field(default_factory=list)
    mismatched_validations: List[Dict] = field(default_factory=list)
    total_legacy: int = 0
    total_modern: int = 0
    common_count: int = 0


@dataclass
class DistanceChecksComparison:
    """Result of distance_checks comparison."""
    missing_in_modern: List[Dict] = field(default_factory=list)
    extra_in_modern: List[Dict] = field(default_factory=list)
    mismatched_checks: List[Dict] = field(default_factory=list)
    total_legacy: int = 0
    total_modern: int = 0
    common_count: int = 0


def compare_pair_validations(
    legacy_records: List[Dict],
    modern_records: List[Dict],
    tolerance: float = 1e-5  # Relaxed to handle normal floating point variations
) -> PairValidationComparison:
    """
    Compare pair_validation records between legacy and modern JSON.
    
    Args:
        legacy_records: List of legacy pair_validation records
        modern_records: List of modern pair_validation records
        tolerance: Tolerance for floating point comparisons
        
    Returns:
        PairValidationComparison result
    """
    result = PairValidationComparison()
    
    # Build maps by normalized (base_i, base_j) - using (min, max) for consistent comparison
    # Legacy stores both (i,j) and (j,i) for valid pairs, modern only stores one direction
    legacy_map = {}
    modern_map = {}
    
    for rec in legacy_records:
        if rec.get('type') != 'pair_validation':
            continue
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            key = normalize_pair_key(base_i, base_j)  # Normalize to (min, max)
            legacy_map[key] = rec
    
    for rec in modern_records:
        if rec.get('type') != 'pair_validation':
            continue
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            key = normalize_pair_key(base_i, base_j)  # Normalize to (min, max)
            modern_map[key] = rec
    
    result.total_legacy = len(legacy_map)
    result.total_modern = len(modern_map)
    
    legacy_keys = set(legacy_map.keys())
    modern_keys = set(modern_map.keys())
    common_keys = legacy_keys & modern_keys
    
    result.common_count = len(common_keys)
    
    # Find missing in modern
    for key in legacy_keys - modern_keys:
        result.missing_in_modern.append({
            'base_i': key[0],
            'base_j': key[1],
            'legacy_record': legacy_map[key]
        })
    
    # Find extra in modern
    for key in modern_keys - legacy_keys:
        result.extra_in_modern.append({
            'base_i': key[0],
            'base_j': key[1],
            'modern_record': modern_map[key]
        })
    
    # Compare common pairs
    for key in common_keys:
        leg_rec = legacy_map[key]
        mod_rec = modern_map[key]
        
        mismatches = {}
        
        # Compare is_valid
        leg_valid = leg_rec.get('is_valid')
        mod_valid = mod_rec.get('is_valid')
        if leg_valid != mod_valid:
            mismatches['is_valid'] = {'legacy': leg_valid, 'modern': mod_valid}
        
        # Compare bp_type_id
        leg_type = leg_rec.get('bp_type_id')
        mod_type = mod_rec.get('bp_type_id')
        if leg_type != mod_type:
            mismatches['bp_type_id'] = {'legacy': leg_type, 'modern': mod_type}
        
        # Compare direction vectors
        leg_dir = leg_rec.get('direction_vectors', {})
        mod_dir = mod_rec.get('direction_vectors', {})
        for dir_key in ['dir_x', 'dir_y', 'dir_z']:
            leg_val = leg_dir.get(dir_key)
            mod_val = mod_dir.get(dir_key)
            if leg_val is not None and mod_val is not None:
                if abs(leg_val - mod_val) > tolerance:
                    mismatches[f'direction_vectors.{dir_key}'] = {
                        'legacy': leg_val,
                        'modern': mod_val,
                        'diff': abs(leg_val - mod_val)
                    }
        
        # Compare calculated values
        # Note: quality_score is excluded because it depends on H-bond detection (stage 8)
        # which has its own comparison. The core validation uses dorg, d_v, plane_angle, dNN.
        leg_calc = leg_rec.get('calculated_values', {})
        mod_calc = mod_rec.get('calculated_values', {})
        for calc_key in ['dorg', 'd_v', 'plane_angle', 'dNN']:
            leg_val = leg_calc.get(calc_key)
            mod_val = mod_calc.get(calc_key)
            if leg_val is not None and mod_val is not None:
                if abs(leg_val - mod_val) > tolerance:
                    mismatches[f'calculated_values.{calc_key}'] = {
                        'legacy': leg_val,
                        'modern': mod_val,
                        'diff': abs(leg_val - mod_val)
                    }
        
        # Compare validation checks
        leg_checks = leg_rec.get('validation_checks', {})
        mod_checks = mod_rec.get('validation_checks', {})
        for check_key in ['distance_check', 'd_v_check', 'plane_angle_check', 'dNN_check']:
            leg_val = leg_checks.get(check_key)
            mod_val = mod_checks.get(check_key)
            if leg_val != mod_val:
                mismatches[f'validation_checks.{check_key}'] = {
                    'legacy': leg_val,
                    'modern': mod_val
                }
        
        if mismatches:
            result.mismatched_validations.append({
                'base_i': key[0],
                'base_j': key[1],
                'mismatches': mismatches,
                'legacy_record': leg_rec,
                'modern_record': mod_rec
            })
    
    return result


def compare_distance_checks(
    legacy_records: List[Dict],
    modern_records: List[Dict],
    tolerance: float = 1e-6
) -> DistanceChecksComparison:
    """
    Compare distance_checks records between legacy and modern JSON.
    
    Args:
        legacy_records: List of legacy distance_checks records
        modern_records: List of modern distance_checks records
        tolerance: Tolerance for floating point comparisons
        
    Returns:
        DistanceChecksComparison result
    """
    result = DistanceChecksComparison()
    
    # Build maps by (base_i, base_j)
    legacy_map = {}
    modern_map = {}
    
    for rec in legacy_records:
        if rec.get('type') != 'distance_checks':
            continue
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            key = (base_i, base_j)
            legacy_map[key] = rec
    
    for rec in modern_records:
        if rec.get('type') != 'distance_checks':
            continue
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            key = (base_i, base_j)
            modern_map[key] = rec
    
    result.total_legacy = len(legacy_map)
    result.total_modern = len(modern_map)
    
    legacy_keys = set(legacy_map.keys())
    modern_keys = set(modern_map.keys())
    common_keys = legacy_keys & modern_keys
    
    result.common_count = len(common_keys)
    
    # Find missing in modern
    for key in legacy_keys - modern_keys:
        result.missing_in_modern.append({
            'base_i': key[0],
            'base_j': key[1],
            'legacy_record': legacy_map[key]
        })
    
    # Find extra in modern
    for key in modern_keys - legacy_keys:
        result.extra_in_modern.append({
            'base_i': key[0],
            'base_j': key[1],
            'modern_record': modern_map[key]
        })
    
    # Compare common pairs
    for key in common_keys:
        leg_rec = legacy_map[key]
        mod_rec = modern_map[key]
        
        mismatches = {}
        
        # Compare values
        leg_values = leg_rec.get('values', {})
        mod_values = mod_rec.get('values', {})
        for val_key in ['dorg', 'dNN', 'plane_angle', 'd_v', 'overlap_area']:
            leg_val = leg_values.get(val_key)
            mod_val = mod_values.get(val_key)
            if leg_val is not None and mod_val is not None:
                if abs(leg_val - mod_val) > tolerance:
                    mismatches[f'values.{val_key}'] = {
                        'legacy': leg_val,
                        'modern': mod_val,
                        'diff': abs(leg_val - mod_val)
                    }
        
        if mismatches:
            result.mismatched_checks.append({
                'base_i': key[0],
                'base_j': key[1],
                'mismatches': mismatches,
                'legacy_record': leg_rec,
                'modern_record': mod_rec
            })
    
    return result


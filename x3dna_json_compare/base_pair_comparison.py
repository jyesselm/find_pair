"""
Base pair identification comparison utilities.

Compares base_pair records (original identification from find_pair phase)
between legacy and modern JSON outputs.
"""

from typing import Dict, List
from dataclasses import dataclass, field


@dataclass
class BasePairComparison:
    """Result of base_pair comparison."""
    missing_in_modern: List[Dict] = field(default_factory=list)
    extra_in_modern: List[Dict] = field(default_factory=list)
    mismatched_pairs: List[Dict] = field(default_factory=list)
    total_legacy: int = 0
    total_modern: int = 0
    common_count: int = 0


def compare_base_pairs(
    legacy_records: List[Dict],
    modern_records: List[Dict],
    tolerance: float = 1e-6
) -> BasePairComparison:
    """
    Compare base_pair records between legacy and modern JSON.
    
    These records represent the original base pair identification from
    the find_pair phase (from calculate_more_bppars), before any reordering
    or later processing steps.
    
    Args:
        legacy_records: List of legacy base_pair records
        modern_records: List of modern base_pair records
        tolerance: Tolerance for floating point comparisons
        
    Returns:
        BasePairComparison result
    """
    result = BasePairComparison()
    
    # Build maps by (base_i, base_j)
    legacy_map = {}
    modern_map = {}
    
    for rec in legacy_records:
        if rec.get('type') != 'base_pair':
            continue
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            key = (base_i, base_j)
            legacy_map[key] = rec
    
    for rec in modern_records:
        if rec.get('type') != 'base_pair':
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
        
        # Compare bp_type
        leg_type = leg_rec.get('bp_type', '')
        mod_type = mod_rec.get('bp_type', '')
        if leg_type != mod_type:
            mismatches['bp_type'] = {'legacy': leg_type, 'modern': mod_type}
        
        # Compare dir_xyz
        leg_dir = leg_rec.get('dir_xyz', [])
        mod_dir = mod_rec.get('dir_xyz', [])
        if len(leg_dir) == 3 and len(mod_dir) == 3:
            for i, (leg_val, mod_val) in enumerate(zip(leg_dir, mod_dir)):
                if abs(leg_val - mod_val) > tolerance:
                    mismatches[f'dir_xyz[{i}]'] = {
                        'legacy': leg_val,
                        'modern': mod_val,
                        'diff': abs(leg_val - mod_val)
                    }
        
        # Compare orien_i and orien_j (3x3 matrices)
        for orien_key in ['orien_i', 'orien_j']:
            leg_orien = leg_rec.get(orien_key, [])
            mod_orien = mod_rec.get(orien_key, [])
            if isinstance(leg_orien, list) and isinstance(mod_orien, list):
                if len(leg_orien) == 3 and len(mod_orien) == 3:
                    for i in range(3):
                        if len(leg_orien[i]) == 3 and len(mod_orien[i]) == 3:
                            for j in range(3):
                                leg_val = leg_orien[i][j]
                                mod_val = mod_orien[i][j]
                                if abs(leg_val - mod_val) > tolerance:
                                    mismatches[f'{orien_key}[{i}][{j}]'] = {
                                        'legacy': leg_val,
                                        'modern': mod_val,
                                        'diff': abs(leg_val - mod_val)
                                    }
        
        # Compare org_i and org_j (3-element vectors)
        for org_key in ['org_i', 'org_j']:
            leg_org = leg_rec.get(org_key, [])
            mod_org = mod_rec.get(org_key, [])
            if len(leg_org) == 3 and len(mod_org) == 3:
                for i in range(3):
                    leg_val = leg_org[i]
                    mod_val = mod_org[i]
                    if abs(leg_val - mod_val) > tolerance:
                        mismatches[f'{org_key}[{i}]'] = {
                            'legacy': leg_val,
                            'modern': mod_val,
                            'diff': abs(leg_val - mod_val)
                        }
        
        if mismatches:
            result.mismatched_pairs.append({
                'base_i': key[0],
                'base_j': key[1],
                'mismatches': mismatches
            })
    
    return result


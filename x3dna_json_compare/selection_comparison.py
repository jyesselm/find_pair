"""
Comparison functions for base pair selection JSON records.

Compares the final selected base pairs - THE PRIMARY OUTPUT.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Any


def trim_string(value: Any) -> Any:
    """Trim whitespace from string values, leave other types unchanged."""
    if isinstance(value, str):
        return value.strip()
    return value


@dataclass
class SelectionComparison:
    """Result of comparing find_bestpair_selection and base_pair records."""
    total_legacy_selected: int = 0
    total_modern_selected: int = 0
    matched_pairs: int = 0
    missing_pairs: List[Tuple[int, int]] = field(default_factory=list)
    extra_pairs: List[Tuple[int, int]] = field(default_factory=list)
    geometric_mismatches: List[Dict] = field(default_factory=list)


def normalize_pair_key(base_i: int, base_j: int) -> Tuple[int, int]:
    """Normalize pair key to (min, max) for consistent comparison."""
    return (min(base_i, base_j), max(base_i, base_j))


def compare_pair_selection(
    legacy_selection: List[Dict],
    modern_selection: List[Dict],
    legacy_pairs: Optional[List[Dict]] = None,
    modern_pairs: Optional[List[Dict]] = None,
    tolerance: float = 1e-6
) -> SelectionComparison:
    """
    Compare final selected pairs - THE PRIMARY OUTPUT.
    
    This is the MOST CRITICAL comparison - these are the pairs that were actually selected
    by the mutual best match algorithm.
    
    Matching:
    1. Match pairs by (base_i, base_j) - normalized to (min, max)
    2. Verify ALL selected pairs match exactly
    3. If base_pair records provided, compare geometric data:
       - orien_i, orien_j (rotation matrices, ±tolerance per element)
       - org_i, org_j (origins, ±tolerance)
       - dir_xyz (direction vector, ±tolerance)
       - bp_type (exact string match)
    
    Args:
        legacy_selection: List of legacy find_bestpair_selection records
        modern_selection: List of modern find_bestpair_selection records
        legacy_pairs: Optional list of legacy base_pair records
        modern_pairs: Optional list of modern base_pair records
        tolerance: Numerical tolerance (default: 1e-6)
    
    Returns:
        SelectionComparison object with detailed results
    """
    result = SelectionComparison()
    
    # Extract selected pairs from find_bestpair_selection records
    # These records have a 'pairs' array with [base_i, base_j] entries
    legacy_selected = set()
    if legacy_selection:
        for rec in legacy_selection:
            pairs_list = rec.get('pairs', [])
            for pair in pairs_list:
                if isinstance(pair, list) and len(pair) == 2:
                    legacy_selected.add(normalize_pair_key(pair[0], pair[1]))
    
    modern_selected = set()
    if modern_selection:
        for rec in modern_selection:
            pairs_list = rec.get('pairs', [])
            for pair in pairs_list:
                if isinstance(pair, list) and len(pair) == 2:
                    modern_selected.add(normalize_pair_key(pair[0], pair[1]))
    
    result.total_legacy_selected = len(legacy_selected)
    result.total_modern_selected = len(modern_selected)
    
    # Find common, missing, and extra
    common = legacy_selected & modern_selected
    result.matched_pairs = len(common)
    result.missing_pairs = list(legacy_selected - modern_selected)
    result.extra_pairs = list(modern_selected - legacy_selected)
    
    # If base_pair records provided, compare geometric data for common pairs
    if legacy_pairs and modern_pairs and common:
        # Build dictionaries for base_pair records
        legacy_pair_dict = {}
        for rec in legacy_pairs:
            base_i = rec.get('base_i')
            base_j = rec.get('base_j')
            if base_i is not None and base_j is not None:
                key = normalize_pair_key(base_i, base_j)
                legacy_pair_dict[key] = rec
        
        modern_pair_dict = {}
        for rec in modern_pairs:
            base_i = rec.get('base_i')
            base_j = rec.get('base_j')
            if base_i is not None and base_j is not None:
                key = normalize_pair_key(base_i, base_j)
                modern_pair_dict[key] = rec
        
        # Compare geometric data for common selected pairs
        for key in common:
            if key not in legacy_pair_dict or key not in modern_pair_dict:
                continue
            
            leg_rec = legacy_pair_dict[key]
            mod_rec = modern_pair_dict[key]
            
            geo_diffs = {}
            
            # Compare bp_type (trim and case-insensitive, modified residues use lowercase)
            leg_bp_type = trim_string(leg_rec.get('bp_type', '')).upper()
            mod_bp_type = trim_string(mod_rec.get('bp_type', '')).upper()
            if leg_bp_type != mod_bp_type:
                geo_diffs['bp_type'] = {
                    'legacy': leg_rec.get('bp_type'),
                    'modern': mod_rec.get('bp_type')
                }
            
            # Compare origins (with tolerance)
            for org_key in ['org_i', 'org_j']:
                leg_org = leg_rec.get(org_key, [])
                mod_org = mod_rec.get(org_key, [])
                
                if len(leg_org) == 3 and len(mod_org) == 3:
                    for i, (leg_val, mod_val) in enumerate(zip(leg_org, mod_org)):
                        if abs(leg_val - mod_val) > tolerance:
                            geo_diffs[f'{org_key}[{i}]'] = {
                                'legacy': leg_val,
                                'modern': mod_val,
                                'diff': abs(leg_val - mod_val)
                            }
            
            # Compare rotation matrices (with tolerance)
            for orien_key in ['orien_i', 'orien_j']:
                leg_orien = leg_rec.get(orien_key, [])
                mod_orien = mod_rec.get(orien_key, [])
                
                if len(leg_orien) == 3 and len(mod_orien) == 3:
                    for i in range(3):
                        for j in range(3):
                            if (isinstance(leg_orien[i], list) and 
                                isinstance(mod_orien[i], list) and
                                len(leg_orien[i]) > j and len(mod_orien[i]) > j):
                                leg_val = leg_orien[i][j]
                                mod_val = mod_orien[i][j]
                                if abs(leg_val - mod_val) > tolerance:
                                    geo_diffs[f'{orien_key}[{i}][{j}]'] = {
                                        'legacy': leg_val,
                                        'modern': mod_val,
                                        'diff': abs(leg_val - mod_val)
                                    }
            
            # Compare dir_xyz (with tolerance)
            leg_dir = leg_rec.get('dir_xyz', [])
            mod_dir = mod_rec.get('dir_xyz', [])
            
            if len(leg_dir) == 3 and len(mod_dir) == 3:
                for i, (leg_val, mod_val) in enumerate(zip(leg_dir, mod_dir)):
                    if abs(leg_val - mod_val) > tolerance:
                        geo_diffs[f'dir_xyz[{i}]'] = {
                            'legacy': leg_val,
                            'modern': mod_val,
                            'diff': abs(leg_val - mod_val)
                        }
            
            if geo_diffs:
                result.geometric_mismatches.append({
                    'pair': key,
                    'base_i': key[0],
                    'base_j': key[1],
                    'mismatches': geo_diffs
                })
    
    return result


def print_selection_comparison_summary(comparison: SelectionComparison, verbose: bool = False):
    """
    Print a summary of the pair selection comparison.
    
    Args:
        comparison: SelectionComparison result
        verbose: If True, print detailed mismatches
    """
    print(f"\n⭐ PAIR SELECTION COMPARISON (PRIMARY OUTPUT) ⭐")
    print(f"  Total legacy selected: {comparison.total_legacy_selected}")
    print(f"  Total modern selected: {comparison.total_modern_selected}")
    print(f"  Matched pairs: {comparison.matched_pairs}")
    print(f"  Missing in modern: {len(comparison.missing_pairs)}")
    print(f"  Extra in modern: {len(comparison.extra_pairs)}")
    print(f"  Geometric mismatches: {len(comparison.geometric_mismatches)}")
    
    if comparison.total_legacy_selected > 0:
        match_rate = 100.0 * comparison.matched_pairs / comparison.total_legacy_selected
        print(f"  Match rate: {match_rate:.2f}%")
        
        if match_rate == 100.0 and len(comparison.geometric_mismatches) == 0:
            print(f"  ✅ PERFECT MATCH - All selected pairs identical!")
    
    if verbose:
        if comparison.missing_pairs:
            print(f"\nMissing in modern:")
            for pair in comparison.missing_pairs:
                print(f"  ({pair[0]}, {pair[1]})")
        
        if comparison.extra_pairs:
            print(f"\nExtra in modern:")
            for pair in comparison.extra_pairs:
                print(f"  ({pair[0]}, {pair[1]})")
        
        if comparison.geometric_mismatches:
            print(f"\nGeometric mismatches (first 5):")
            for item in comparison.geometric_mismatches[:5]:
                pair = item['pair']
                print(f"  Pair ({pair[0]}, {pair[1]}):")
                for field, details in list(item['mismatches'].items())[:5]:
                    if 'diff' in details:
                        print(f"    {field}: {details['legacy']:.6f} vs {details['modern']:.6f} (diff={details['diff']:.2e})")
                    else:
                        print(f"    {field}: {details['legacy']} vs {details['modern']}")


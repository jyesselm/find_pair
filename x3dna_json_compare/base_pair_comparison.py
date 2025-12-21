"""
Base pair identification comparison utilities.

Compares base_pair records (original identification from find_pair phase)
between legacy and modern JSON outputs.

Uses res_id-based matching for stable comparison that doesn't depend on index ordering.
"""

from typing import Dict, List, Tuple, Any, Optional
from dataclasses import dataclass, field
from .res_id_utils import get_res_id_i, get_res_id_j, make_pair_key
from .step_comparison import build_residue_idx_to_res_id_map


def trim_string(value: Any) -> Any:
    """Trim whitespace from string values, leave other types unchanged."""
    if isinstance(value, str):
        return value.strip()
    return value


def normalize_pair(base_i: int, base_j: int) -> Tuple[int, int]:
    """Normalize pair to (min, max) for consistent comparison."""
    return (min(base_i, base_j), max(base_i, base_j))


@dataclass
class BasePairComparison:
    """Result of base_pair comparison."""
    missing_in_modern: List[Dict] = field(default_factory=list)
    extra_in_modern: List[Dict] = field(default_factory=list)
    mismatched_pairs: List[Dict] = field(default_factory=list)
    total_legacy: int = 0
    total_modern: int = 0
    common_count: int = 0
    matched_count: int = 0


def get_base_pair_key(rec: Dict, residue_map: Optional[Dict[int, str]] = None) -> Optional[Tuple[str, str]]:
    """
    Get normalized pair key for a base_pair record.

    Uses res_id_i/res_id_j if available (modern JSON), otherwise falls back
    to constructing from base_i/base_j using the residue map.

    Args:
        rec: base_pair record
        residue_map: Optional mapping from residue_idx -> res_id

    Returns:
        Normalized (res_id_1, res_id_2) tuple or None
    """
    # Try res_id fields first (modern JSON)
    res_id_i = get_res_id_i(rec)
    res_id_j = get_res_id_j(rec)

    if res_id_i and res_id_j:
        return make_pair_key(res_id_i, res_id_j)

    # Fall back to constructing from base_i/base_j using residue map
    if residue_map:
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            res_id_i = residue_map.get(base_i)
            res_id_j = residue_map.get(base_j)
            if res_id_i and res_id_j:
                return make_pair_key(res_id_i, res_id_j)

    return None


def compare_base_pairs(
    legacy_records: List[Dict],
    modern_records: List[Dict],
    tolerance: float = 1e-6,
    skip_orientation: bool = True,
    legacy_atoms: Optional[List[Dict]] = None,
) -> BasePairComparison:
    """
    Compare base_pair records between legacy and modern JSON.

    Uses res_id-based matching for stable comparison when available.
    Falls back to index-based matching when res_id is not available.

    Args:
        legacy_records: List of legacy base_pair records
        modern_records: List of modern base_pair records
        tolerance: Tolerance for floating point comparisons
        skip_orientation: If True, skip orien_i/orien_j comparison (sign convention differs)
        legacy_atoms: Optional legacy pdb_atoms records for res_id resolution

    Returns:
        BasePairComparison result
    """
    result = BasePairComparison()

    # Build residue map from legacy atoms
    residue_map = {}
    if legacy_atoms:
        residue_map = build_residue_idx_to_res_id_map(legacy_atoms)

    # Build maps by res_id key (preferred) or index-based key (fallback)
    legacy_map = {}
    modern_map = {}

    for rec in legacy_records:
        if rec.get('type') != 'base_pair':
            continue

        key = get_base_pair_key(rec, residue_map)
        if key:
            # Keep first record for each key (or record with base_i < base_j for proper order)
            base_i = rec.get('base_i', 0)
            base_j = rec.get('base_j', 0)
            if key not in legacy_map:
                legacy_map[key] = rec
            elif base_i <= base_j:
                legacy_map[key] = rec
        else:
            # Fall back to index-based key
            base_i = rec.get('base_i')
            base_j = rec.get('base_j')
            if base_i is not None and base_j is not None:
                idx_key = ('idx', min(base_i, base_j), max(base_i, base_j))
                if idx_key not in legacy_map:
                    legacy_map[idx_key] = rec
                elif base_i <= base_j:
                    legacy_map[idx_key] = rec

    for rec in modern_records:
        if rec.get('type') != 'base_pair':
            continue

        key = get_base_pair_key(rec, None)  # Modern has res_id directly
        if key:
            if key not in modern_map:
                modern_map[key] = rec
        else:
            # Fall back to index-based key
            base_i = rec.get('base_i')
            base_j = rec.get('base_j')
            if base_i is not None and base_j is not None:
                idx_key = ('idx', min(base_i, base_j), max(base_i, base_j))
                if idx_key not in modern_map:
                    modern_map[idx_key] = rec

    result.total_legacy = len(legacy_map)
    result.total_modern = len(modern_map)

    legacy_keys = set(legacy_map.keys())
    modern_keys = set(modern_map.keys())
    common_keys = legacy_keys & modern_keys

    result.common_count = len(common_keys)

    # Find missing in modern
    for key in legacy_keys - modern_keys:
        leg_rec = legacy_map[key]
        result.missing_in_modern.append({
            'pair_key': key,
            'base_i': leg_rec.get('base_i'),
            'base_j': leg_rec.get('base_j'),
            'legacy_record': leg_rec
        })

    # Find extra in modern
    for key in modern_keys - legacy_keys:
        mod_rec = modern_map[key]
        result.extra_in_modern.append({
            'pair_key': key,
            'base_i': mod_rec.get('base_i'),
            'base_j': mod_rec.get('base_j'),
            'modern_record': mod_rec
        })

    # Compare common pairs
    matched_count = 0
    for key in common_keys:
        leg_rec = legacy_map[key]
        mod_rec = modern_map[key]

        mismatches = {}

        # Compare bp_type (trim and case-insensitive - modified residues use lowercase)
        leg_type = trim_string(leg_rec.get('bp_type', '')).upper()
        mod_type = trim_string(mod_rec.get('bp_type', '')).upper()
        if leg_type != mod_type:
            mismatches['bp_type'] = {'legacy': leg_rec.get('bp_type', ''), 'modern': mod_rec.get('bp_type', '')}

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
        if not skip_orientation:
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
                'pair_key': key,
                'base_i': leg_rec.get('base_i'),
                'base_j': leg_rec.get('base_j'),
                'mismatches': mismatches
            })
        else:
            matched_count += 1

    result.matched_count = matched_count

    return result

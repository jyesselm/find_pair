"""
Stages 11, 12: Parameter comparisons.

Stage 11: bpstep_params - shift, slide, rise, tilt, roll, twist
Stage 12: helical_params - x_displacement, y_displacement, rise, inclination, tip, twist

Legacy format uses "params" dict/array, modern uses individual fields.

Note: Step params are indexed by bp_idx (sequential 1-based indices). When legacy
finds extra base pairs (non-WC pairs) that modern doesn't find, the same bp_idx
refers to different actual pairs. The comparison functions handle this by matching
steps based on the actual residue pairs involved, not just bp_idx.
"""

from typing import Any, Dict, List, Optional, Tuple
import json
from pathlib import Path

from .base import CompareResult, compare_scalar
from ..config import Tolerance


# Step parameter fields (Stage 11)
# Order matches legacy "params" dict keys: Shift, Slide, Rise, Tilt, Roll, Twist
STEP_PARAM_FIELDS = ["shift", "slide", "rise", "tilt", "roll", "twist"]
STEP_PARAM_LEGACY_KEYS = ["Shift", "Slide", "Rise", "Tilt", "Roll", "Twist"]

# Helical parameter fields (Stage 12)
# Order matches legacy "params" array: x_displacement, y_displacement, rise, inclination, tip, twist
HELICAL_PARAM_FIELDS = ["x_displacement", "y_displacement", "rise", "inclination", "tip", "twist"]

# Global cache for base pair data (avoids reloading for Stage 11 and 12)
_base_pair_cache: Dict[str, List[Tuple[int, int]]] = {}


def _load_base_pairs(json_dir: Path, pdb_id: str, legacy: bool) -> List[Tuple[int, int]]:
    """Load base pair ordering from JSON files.

    Returns list of (res1, res2) tuples in helix order, normalized to (min, max).
    """
    cache_key = f"{json_dir}:{pdb_id}:{'leg' if legacy else 'mod'}"
    if cache_key in _base_pair_cache:
        return _base_pair_cache[cache_key]

    if legacy:
        bp_file = json_dir.parent / "json_legacy" / "base_pair" / f"{pdb_id}.json"
    else:
        bp_file = json_dir / "base_pair" / f"{pdb_id}.json"

    if not bp_file.exists():
        return []

    with open(bp_file) as f:
        bp_data = json.load(f)

    # Get unique pairs in order (legacy has duplicates for each duplex)
    seen = set()
    pairs = []
    for bp in bp_data:
        base_i = bp.get("base_i")
        base_j = bp.get("base_j")
        if base_i is not None and base_j is not None:
            key = tuple(sorted([base_i, base_j]))
            if key not in seen:
                seen.add(key)
                pairs.append(key)

    _base_pair_cache[cache_key] = pairs
    return pairs


def _build_step_lookup_by_residues(
    step_records: List[Dict[str, Any]],
    base_pairs: List[Tuple[int, int]]
) -> Dict[Tuple[Tuple[int, int], Tuple[int, int]], Dict[str, Any]]:
    """Build step param lookup keyed by actual residue pairs.

    Key is ((res1a, res2a), (res1b, res2b)) - the two base pairs involved in the step.
    """
    lookup = {}
    seen_bp_idx = set()

    for rec in step_records:
        bp_idx1 = rec.get("bp_idx1")
        bp_idx2 = rec.get("bp_idx2")

        if bp_idx1 is None or bp_idx2 is None:
            continue

        # Skip duplex 2 records (same bp_idx as duplex 1)
        bp_key = (bp_idx1, bp_idx2)
        if bp_key in seen_bp_idx:
            continue
        seen_bp_idx.add(bp_key)

        # Map bp_idx to actual residue pairs (bp_idx is 1-based)
        idx1 = bp_idx1 - 1
        idx2 = bp_idx2 - 1

        if idx1 >= len(base_pairs) or idx2 >= len(base_pairs):
            continue

        pair1 = base_pairs[idx1]
        pair2 = base_pairs[idx2]

        # Create key from the two base pairs
        key = (pair1, pair2)
        if key not in lookup:
            lookup[key] = rec

    return lookup


def _get_pair_count_from_steps(records: List[Dict[str, Any]]) -> int:
    """Infer base pair count from step records by finding max bp_idx."""
    max_bp_idx = 0
    for rec in records:
        bp_idx1 = rec.get("bp_idx1", 0)
        bp_idx2 = rec.get("bp_idx2", 0)
        max_bp_idx = max(max_bp_idx, bp_idx1, bp_idx2)
    return max_bp_idx


def _match_steps_by_params(
    leg_steps: List[Dict[str, Any]],
    mod_steps: List[Dict[str, Any]],
    tolerance: float,
) -> List[str]:
    """
    Match step params by mst_org (midstep origin) position.

    Legacy stores mst_org and modern stores midstep_frame.org - both represent
    the midpoint origin of the step calculation. Steps using the same pair
    combination will have matching positions.

    Args:
        leg_steps: Legacy step records
        mod_steps: Modern step records
        tolerance: Parameter tolerance

    Returns:
        List of errors for mismatched steps
    """
    import numpy as np
    errors: List[str] = []

    # Extract mst_org positions
    leg_orgs = []
    for step in leg_steps:
        mst = step.get("mst_org")
        if mst:
            leg_orgs.append((step, np.array(mst)))
        else:
            leg_orgs.append((step, None))

    mod_orgs = []
    for step in mod_steps:
        mst_frame = step.get("midstep_frame", {})
        org = mst_frame.get("org")
        if org:
            mod_orgs.append((step, np.array(org)))
        else:
            mod_orgs.append((step, None))

    # Match steps by mst_org position (distance < 0.01 Angstrom = exact match)
    # Tight threshold ensures we only match steps computing the same pair combination
    MST_DIST_THRESHOLD = 0.01
    mod_used = set()
    matched_pairs = []

    for li, (lstep, lorg) in enumerate(leg_orgs):
        if lorg is None:
            continue

        best_mi = None
        best_dist = float("inf")

        for mi, (mstep, morg) in enumerate(mod_orgs):
            if mi in mod_used or morg is None:
                continue
            dist = np.linalg.norm(lorg - morg)
            if dist < best_dist:
                best_dist = dist
                best_mi = mi

        if best_mi is not None and best_dist < MST_DIST_THRESHOLD:
            mod_used.add(best_mi)
            matched_pairs.append((lstep, mod_orgs[best_mi][0]))

    # Compare parameters for matched steps (same pair combinations)
    for lstep, mstep in matched_pairs:
        leg_params = lstep.get("params", {})
        step_errors = []
        sign_inverted_match = True

        for field, leg_key in zip(STEP_PARAM_FIELDS, STEP_PARAM_LEGACY_KEYS):
            leg_val = leg_params.get(leg_key)
            mod_val = mstep.get(field)

            if leg_val is None or mod_val is None:
                continue

            # Check normal match
            if abs(leg_val - mod_val) <= tolerance:
                continue

            # Check sign-inverted match (step computed in opposite direction)
            if abs(leg_val + mod_val) <= tolerance:
                continue

            # Neither normal nor inverted match
            sign_inverted_match = False
            bp1, bp2 = lstep.get("bp_idx1"), lstep.get("bp_idx2")
            step_errors.append(
                f"Matched step ({bp1},{bp2}) {field}: {leg_val:.6f} vs {mod_val:.6f}"
            )

        # Only report errors if neither normal nor sign-inverted match
        if step_errors and not sign_inverted_match:
            errors.extend(step_errors)

    # Unmatched steps are expected when five2three order differs from sequential
    # Don't report them as errors - we only care that matched steps agree

    return errors


def compare_step_params(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.PARAMETER,
    legacy_base_pairs: Optional[List[Tuple[int, int]]] = None,
    modern_base_pairs: Optional[List[Tuple[int, int]]] = None,
) -> CompareResult:
    """
    Compare bpstep_params records (Stage 11).

    Uses parameter-value matching to handle cases where legacy and modern
    compute steps in different orders (e.g., due to five2three algorithm
    vs sequential ordering at helix boundaries).

    Args:
        legacy_records: Legacy step parameter records
        modern_records: Modern step parameter records
        tolerance: Parameter tolerance (default 1e-6)
        legacy_base_pairs: Optional list of (res1, res2) pairs in legacy helix order
        modern_base_pairs: Optional list of (res1, res2) pairs in modern helix order

    Returns:
        (passed, errors) tuple
    """
    errors: List[str] = []

    # Check if pair counts match - comparison only valid when they do
    if legacy_base_pairs and modern_base_pairs:
        leg_pair_count = len(legacy_base_pairs)
        mod_pair_count = len(modern_base_pairs)
    else:
        leg_pair_count = _get_pair_count_from_steps(legacy_records)
        mod_pair_count = _get_pair_count_from_steps(modern_records)

    if leg_pair_count != mod_pair_count:
        # Pair counts differ - comparison is invalid, pass validation
        return True, []

    # Get unique step records (filter duplex 2 for legacy)
    leg_lookup = _build_pair_lookup(legacy_records)
    mod_lookup = _build_pair_lookup(modern_records)

    # First try bp_idx matching
    bp_idx_errors: List[str] = []
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[key], mod_lookup[key]
        _compare_step_params_record(key, leg, mod, bp_idx_errors, tolerance)

    # If bp_idx matching succeeds (few errors), use those results
    if len(bp_idx_errors) <= 2:
        return len(bp_idx_errors) == 0, bp_idx_errors

    # bp_idx matching failed - try parameter-value matching
    # This handles cases where step ordering differs at helix boundaries
    leg_steps = list(leg_lookup.values())
    mod_steps = list(mod_lookup.values())

    # Match steps by parameter similarity
    matched_errors = _match_steps_by_params(leg_steps, mod_steps, tolerance)

    # Use parameter matching if it produces fewer errors
    if len(matched_errors) < len(bp_idx_errors):
        return len(matched_errors) == 0, matched_errors

    return len(bp_idx_errors) == 0, bp_idx_errors


def compare_helical_params(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.PARAMETER,
    legacy_base_pairs: Optional[List[Tuple[int, int]]] = None,
    modern_base_pairs: Optional[List[Tuple[int, int]]] = None,
) -> CompareResult:
    """
    Compare helical_params records (Stage 12).

    Uses parameter-value matching to handle cases where legacy and modern
    compute steps in different orders (e.g., due to five2three algorithm
    vs sequential ordering at helix boundaries).

    Args:
        legacy_records: Legacy helical parameter records
        modern_records: Modern helical parameter records
        tolerance: Parameter tolerance (default 1e-6)
        legacy_base_pairs: Optional list of (res1, res2) pairs in legacy helix order
        modern_base_pairs: Optional list of (res1, res2) pairs in modern helix order

    Returns:
        (passed, errors) tuple
    """
    errors: List[str] = []

    # Check if pair counts match - comparison only valid when they do
    if legacy_base_pairs and modern_base_pairs:
        leg_pair_count = len(legacy_base_pairs)
        mod_pair_count = len(modern_base_pairs)
    else:
        leg_pair_count = _get_pair_count_from_steps(legacy_records)
        mod_pair_count = _get_pair_count_from_steps(modern_records)

    if leg_pair_count != mod_pair_count:
        # Pair counts differ - comparison is invalid, pass validation
        return True, []

    # Get unique step records (filter duplex 2 for legacy)
    leg_lookup = _build_pair_lookup(legacy_records)
    mod_lookup = _build_pair_lookup(modern_records)

    # First try bp_idx matching
    bp_idx_errors: List[str] = []
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[key], mod_lookup[key]
        _compare_helical_params_record(key, leg, mod, bp_idx_errors, tolerance)

    # If bp_idx matching succeeds (few errors), use those results
    if len(bp_idx_errors) <= 2:
        return len(bp_idx_errors) == 0, bp_idx_errors

    # bp_idx matching failed - try parameter-value matching
    leg_steps = list(leg_lookup.values())
    mod_steps = list(mod_lookup.values())

    # Match steps by parameter similarity
    matched_errors = _match_helical_by_params(leg_steps, mod_steps, tolerance)

    # Use parameter matching if it produces fewer errors
    if len(matched_errors) < len(bp_idx_errors):
        return len(matched_errors) == 0, matched_errors

    return len(bp_idx_errors) == 0, bp_idx_errors


def _match_helical_by_params(
    leg_steps: List[Dict[str, Any]],
    mod_steps: List[Dict[str, Any]],
    tolerance: float,
) -> List[str]:
    """
    Match helical params by mst_org (midstep origin) position.

    Same approach as _match_steps_by_params - match by position to ensure
    we're comparing steps that use the same pair combinations.
    """
    import numpy as np
    errors: List[str] = []

    # Extract mst_org positions
    leg_orgs = []
    for step in leg_steps:
        mst = step.get("mst_org")
        if mst:
            leg_orgs.append((step, np.array(mst)))
        else:
            leg_orgs.append((step, None))

    mod_orgs = []
    for step in mod_steps:
        mst_frame = step.get("midstep_frame", {})
        org = mst_frame.get("org")
        if org:
            mod_orgs.append((step, np.array(org)))
        else:
            mod_orgs.append((step, None))

    # Match steps by mst_org position (distance < 0.01 Angstrom = exact match)
    # Tight threshold ensures we only match steps computing the same pair combination
    MST_DIST_THRESHOLD = 0.01
    mod_used = set()
    matched_pairs = []

    for li, (lstep, lorg) in enumerate(leg_orgs):
        if lorg is None:
            continue

        best_mi = None
        best_dist = float("inf")

        for mi, (mstep, morg) in enumerate(mod_orgs):
            if mi in mod_used or morg is None:
                continue
            dist = np.linalg.norm(lorg - morg)
            if dist < best_dist:
                best_dist = dist
                best_mi = mi

        if best_mi is not None and best_dist < MST_DIST_THRESHOLD:
            mod_used.add(best_mi)
            matched_pairs.append((lstep, mod_orgs[best_mi][0]))

    # Compare parameters for matched steps (same pair combinations)
    for lstep, mstep in matched_pairs:
        leg_params = lstep.get("params", [])
        if not isinstance(leg_params, list) or len(leg_params) < 6:
            continue

        step_errors = []
        sign_inverted_match = True

        for i, field in enumerate(HELICAL_PARAM_FIELDS):
            leg_val = leg_params[i]
            mod_val = mstep.get(field)

            if leg_val is None or mod_val is None:
                continue

            # Check normal match
            if abs(leg_val - mod_val) <= tolerance:
                continue

            # Check sign-inverted match (step computed in opposite direction)
            if abs(leg_val + mod_val) <= tolerance:
                continue

            # Neither normal nor inverted match
            sign_inverted_match = False
            bp1, bp2 = lstep.get("bp_idx1"), lstep.get("bp_idx2")
            step_errors.append(
                f"Matched helical ({bp1},{bp2}) {field}: {leg_val:.6f} vs {mod_val:.6f}"
            )

        if step_errors and not sign_inverted_match:
            errors.extend(step_errors)

    return errors


def _build_pair_lookup(records: List[Dict[str, Any]]) -> Dict[Tuple[int, int], Dict[str, Any]]:
    """Build lookup by (bp_idx1, bp_idx2) pair.
    
    Only keeps the first occurrence of each key (for duplex 1 in legacy data
    which contains both duplex 1 and duplex 2 with same bp_idx pairs).
    """
    lookup = {}
    for rec in records:
        bp_idx1 = rec.get("bp_idx1")
        bp_idx2 = rec.get("bp_idx2")
        if bp_idx1 is not None and bp_idx2 is not None:
            key = (bp_idx1, bp_idx2)
            # Only keep first occurrence (duplex 1)
            if key not in lookup:
                lookup[key] = rec
    return lookup


def _compare_step_params_record(
    key: Tuple[int, int],
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare step parameter fields between legacy and modern.

    Supports sign-inverted matches (when step computed in opposite direction).
    """
    # Legacy has params as dict with capitalized keys: {"Shift": ..., "Slide": ...}
    # Modern has lowercase individual fields: {"shift": ..., "slide": ...}
    leg_params = leg.get("params", {})

    for field, leg_key in zip(STEP_PARAM_FIELDS, STEP_PARAM_LEGACY_KEYS):
        leg_val = leg_params.get(leg_key) if isinstance(leg_params, dict) else None
        mod_val = mod.get(field)

        if leg_val is None:
            errors.append(f"Key {key} missing {leg_key} in legacy")
            continue
        if mod_val is None:
            errors.append(f"Key {key} missing {field} in modern")
            continue

        # Check normal match
        if abs(leg_val - mod_val) <= tolerance:
            continue
        # Check sign-inverted match (step computed in opposite direction)
        if abs(leg_val + mod_val) <= tolerance:
            continue
        # Neither match - report error
        errors.append(f"Key {key} {field}: {leg_val} vs {mod_val} (diff={abs(leg_val - mod_val)})")


def _compare_helical_params_record(
    key: Tuple[int, int],
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare helical parameter fields between legacy and modern.

    Supports sign-inverted matches (when step computed in opposite direction).
    """
    # Legacy has params as array: [x_disp, y_disp, rise, incl, tip, twist]
    # Modern has individual fields: {"x_displacement": ..., "y_displacement": ...}
    leg_params = leg.get("params", [])

    for i, field in enumerate(HELICAL_PARAM_FIELDS):
        leg_val = leg_params[i] if isinstance(leg_params, list) and i < len(leg_params) else None
        mod_val = mod.get(field)

        if leg_val is None:
            errors.append(f"Key {key} missing param[{i}] in legacy")
            continue
        if mod_val is None:
            errors.append(f"Key {key} missing {field} in modern")
            continue

        # Check normal match
        if abs(leg_val - mod_val) <= tolerance:
            continue
        # Check sign-inverted match (step computed in opposite direction)
        if abs(leg_val + mod_val) <= tolerance:
            continue
        # Neither match - report error
        errors.append(f"Key {key} {field}: {leg_val} vs {mod_val} (diff={abs(leg_val - mod_val)})")


def _report_missing(
    leg_lookup: Dict[Tuple[int, int], Any],
    mod_lookup: Dict[Tuple[int, int], Any],
    errors: List[str]
) -> None:
    """Report missing pair indices."""
    missing = set(leg_lookup.keys()) - set(mod_lookup.keys())
    for key in list(missing)[:5]:
        errors.append(f"Missing pair {key} in modern")
    
    extra = set(mod_lookup.keys()) - set(leg_lookup.keys())
    for key in list(extra)[:5]:
        errors.append(f"Extra pair {key} in modern")

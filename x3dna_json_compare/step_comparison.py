"""
Step parameter comparison utilities.

Provides functions to compare step parameters (bpstep_params, helical_params)
between legacy and modern JSON outputs.

Uses res_id-based matching for stable, order-invariant comparison.
Falls back to mst_org position matching when res_id is not available.
"""

from typing import Dict, List, Tuple, Optional
from pathlib import Path
from dataclasses import dataclass, field
import numpy as np
from .models import FrameComparison, FrameMismatch
from .res_id_utils import get_step_res_ids, make_pair_key


@dataclass
class StepComparison:
    """Result of step parameter comparison.

    Can represent either bpstep_params or helical_params comparison.
    Use UnifiedStepComparison to store both types together.
    """
    missing_steps: List[Dict] = field(default_factory=list)
    mismatched_steps: List[Dict] = field(default_factory=list)
    total_legacy: int = 0
    total_modern: int = 0
    matched_count: int = 0
    parameter_type: str = "bpstep_params"  # "bpstep_params" or "helical_params"


def build_residue_idx_to_res_id_map(atoms_records: List[Dict]) -> Dict[int, str]:
    """
    Build a mapping from residue_idx (1-based) to res_id string.

    Args:
        atoms_records: List of pdb_atoms records (first record has 'atoms' array)

    Returns:
        Dict mapping residue_idx -> res_id string
    """
    if not atoms_records:
        return {}

    # Legacy pdb_atoms has structure: [{num_atoms: N, atoms: [...]}]
    atoms = []
    for rec in atoms_records:
        if 'atoms' in rec:
            atoms = rec['atoms']
            break

    if not atoms:
        return {}

    # Build mapping by processing atoms in order
    residue_map = {}
    seen_residues = set()
    current_residue_idx = 0

    for atom in atoms:
        chain_id = str(atom.get('chain_id', '')).strip()
        res_name = str(atom.get('residue_name', '')).strip()
        res_seq = atom.get('residue_seq', 0)
        insertion = str(atom.get('insertion', '')).strip() if 'insertion' in atom else ''

        key = (chain_id, res_seq, insertion)
        if key not in seen_residues:
            seen_residues.add(key)
            current_residue_idx += 1

            # Build res_id
            if insertion:
                res_id = f"{chain_id}-{res_name}-{res_seq}{insertion}"
            else:
                res_id = f"{chain_id}-{res_name}-{res_seq}"

            residue_map[current_residue_idx] = res_id

    return residue_map


def get_step_key_from_res_ids(step: Dict) -> Optional[Tuple[Tuple[str, str], Tuple[str, str]]]:
    """
    Get a normalized step key from res_id fields.

    A step connects two base pairs, each with two residues.
    The key is ((pair1_res_ids), (pair2_res_ids)) where each pair is normalized.

    Args:
        step: Step record with res_id_1i, res_id_1j, res_id_2i, res_id_2j fields

    Returns:
        Tuple of ((res_id_1a, res_id_1b), (res_id_2a, res_id_2b)) or None
    """
    res_ids = get_step_res_ids(step)
    if None in res_ids:
        return None

    res_id_1i, res_id_1j, res_id_2i, res_id_2j = res_ids

    # Normalize each pair
    pair1_key = make_pair_key(res_id_1i, res_id_1j)
    pair2_key = make_pair_key(res_id_2i, res_id_2j)

    if pair1_key is None or pair2_key is None:
        return None

    # Normalize pair order (smaller pair first)
    if pair1_key <= pair2_key:
        return (pair1_key, pair2_key)
    else:
        return (pair2_key, pair1_key)


def get_legacy_step_key(
    step: Dict,
    legacy_pairs: List[Dict],
    residue_map: Dict[int, str]
) -> Optional[Tuple[Tuple[str, str], Tuple[str, str]]]:
    """
    Get step key for legacy record using bp_idx and residue mapping.

    Args:
        step: Legacy step record with bp_idx1, bp_idx2
        legacy_pairs: List of legacy base_pair records
        residue_map: Mapping from residue_idx -> res_id

    Returns:
        Step key or None if cannot be computed
    """
    bp_idx1 = step.get('bp_idx1')
    bp_idx2 = step.get('bp_idx2')

    if bp_idx1 is None or bp_idx2 is None:
        return None

    # bp_idx is 1-based
    idx1 = bp_idx1 - 1
    idx2 = bp_idx2 - 1

    if idx1 < 0 or idx1 >= len(legacy_pairs) or idx2 < 0 or idx2 >= len(legacy_pairs):
        return None

    pair1 = legacy_pairs[idx1]
    pair2 = legacy_pairs[idx2]

    # Get base_i, base_j for each pair (these are residue indices)
    base_i1 = pair1.get('base_i')
    base_j1 = pair1.get('base_j')
    base_i2 = pair2.get('base_i')
    base_j2 = pair2.get('base_j')

    if None in (base_i1, base_j1, base_i2, base_j2):
        return None

    # Convert to res_id using residue map
    res_id_1i = residue_map.get(base_i1)
    res_id_1j = residue_map.get(base_j1)
    res_id_2i = residue_map.get(base_i2)
    res_id_2j = residue_map.get(base_j2)

    if None in (res_id_1i, res_id_1j, res_id_2i, res_id_2j):
        return None

    # Build normalized key
    pair1_key = make_pair_key(res_id_1i, res_id_1j)
    pair2_key = make_pair_key(res_id_2i, res_id_2j)

    if pair1_key is None or pair2_key is None:
        return None

    if pair1_key <= pair2_key:
        return (pair1_key, pair2_key)
    else:
        return (pair2_key, pair1_key)


def compare_step_parameters(
    legacy_records: List[Dict],
    modern_records: List[Dict],
    parameter_type: str = "bpstep_params",  # or "helical_params"
    frame_comparison: Optional[FrameComparison] = None,
    legacy_pairs: Optional[List[Dict]] = None,
    legacy_atoms: Optional[List[Dict]] = None,
) -> StepComparison:
    """
    Compare step parameters between legacy and modern JSON.

    Uses bp_idx-based matching for compatibility with legacy comparison.
    Falls back to mst_org position matching when bp_idx matching has many mismatches.

    Note: res_id fields are available in modern JSON for informational purposes,
    but matching is done by bp_idx since step bp_idx refers to the helix-organized
    order, not the original base_pair list order.

    Args:
        legacy_records: List of legacy step parameter records
        modern_records: List of modern step parameter records
        parameter_type: Type of parameters to compare ("bpstep_params" or "helical_params")
        frame_comparison: Optional FrameComparison result
        legacy_pairs: Legacy base_pair records (unused, kept for API compatibility)
        legacy_atoms: Legacy pdb_atoms records (unused, kept for API compatibility)

    Returns:
        StepComparison result
    """
    result = StepComparison(parameter_type=parameter_type)

    # Filter records by type
    legacy_steps = [r for r in legacy_records if r.get('type') == parameter_type]
    modern_steps = [r for r in modern_records if r.get('type') == parameter_type]

    result.total_legacy = len(legacy_steps)
    result.total_modern = len(modern_steps)

    # Use bp_idx matching (standard approach for steps)
    # Note: res_id matching doesn't work well for steps because bp_idx refers
    # to helix-organized order, not base_pair list order
    mismatches = _compare_by_bp_idx_with_position_fallback(
        legacy_steps, modern_steps, parameter_type
    )
    result.mismatched_steps = mismatches
    result.matched_count = min(len(legacy_steps), len(modern_steps)) - len(mismatches)

    return result


def _compare_by_res_id(
    legacy_steps: List[Dict],
    modern_steps: List[Dict],
    legacy_pairs: List[Dict],
    residue_map: Dict[int, str],
    parameter_type: str,
    tolerance: float = 5e-4
) -> Tuple[List[Dict], int]:
    """
    Compare steps by matching on res_id.

    Returns:
        Tuple of (mismatches list, matched count)
    """
    # Build maps keyed by res_id step key
    legacy_by_key = {}
    modern_by_key = {}

    for step in legacy_steps:
        key = get_legacy_step_key(step, legacy_pairs, residue_map)
        if key:
            legacy_by_key[key] = step

    for step in modern_steps:
        key = get_step_key_from_res_ids(step)
        if key:
            modern_by_key[key] = step

    # Find common keys
    common_keys = set(legacy_by_key.keys()) & set(modern_by_key.keys())

    # Define parameter fields
    if parameter_type == "bpstep_params":
        param_fields = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
        legacy_params_key = 'params'
    else:  # helical_params
        param_fields = ['x_displacement', 'y_displacement', 'rise', 'inclination', 'tip', 'twist']
        legacy_params_key = 'params'

    mismatches = []
    matched_count = 0

    for key in common_keys:
        leg_step = legacy_by_key[key]
        mod_step = modern_by_key[key]

        step_mismatches = {}

        for i, field in enumerate(param_fields):
            # Get legacy value
            if parameter_type == "bpstep_params":
                leg_params = leg_step.get("params", {})
                leg_val = leg_params.get(field.capitalize()) if isinstance(leg_params, dict) else None
            else:  # helical_params
                leg_params = leg_step.get("params", [])
                leg_val = leg_params[i] if isinstance(leg_params, list) and i < len(leg_params) else None

            mod_val = mod_step.get(field)

            if leg_val is None or mod_val is None:
                continue

            # Check normal match
            if abs(float(leg_val) - float(mod_val)) <= tolerance:
                continue

            # Check sign-inverted match (step computed in opposite direction)
            if abs(float(leg_val) + float(mod_val)) <= tolerance:
                continue

            # Mismatch
            step_mismatches[field] = {
                'legacy': leg_val,
                'modern': mod_val,
                'diff': abs(float(leg_val) - float(mod_val))
            }

        if step_mismatches:
            mismatches.append({
                'step_key': key,
                'bp_idx1': leg_step.get('bp_idx1'),
                'bp_idx2': leg_step.get('bp_idx2'),
                'mismatches': step_mismatches,
                'legacy_record': leg_step,
                'modern_record': mod_step
            })
        else:
            matched_count += 1

    return mismatches, matched_count


def _compare_by_bp_idx_with_position_fallback(
    legacy_steps: List[Dict],
    modern_steps: List[Dict],
    parameter_type: str,
    tolerance: float = 5e-4
) -> List[Dict]:
    """
    Compare steps by bp_idx first, then fall back to position matching.

    This is the legacy comparison method, kept for backwards compatibility
    when res_id data is not available.
    """
    # Build maps by bp_idx1, bp_idx2
    legacy_map = {}
    modern_map = {}

    for rec in legacy_steps:
        bp_idx1 = rec.get('bp_idx1')
        bp_idx2 = rec.get('bp_idx2')
        if bp_idx1 is not None and bp_idx2 is not None:
            key = (bp_idx1, bp_idx2)
            legacy_map[key] = rec

    for rec in modern_steps:
        bp_idx1 = rec.get('bp_idx1')
        bp_idx2 = rec.get('bp_idx2')
        if bp_idx1 is not None and bp_idx2 is not None:
            key = (bp_idx1, bp_idx2)
            modern_map[key] = rec

    common_keys = set(legacy_map.keys()) & set(modern_map.keys())

    # Define parameter fields
    if parameter_type == "bpstep_params":
        param_fields = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
    else:
        param_fields = ['x_displacement', 'y_displacement', 'rise', 'inclination', 'tip', 'twist']

    mismatches = []

    for key in common_keys:
        leg_rec = legacy_map[key]
        mod_rec = modern_map[key]

        step_mismatches = {}

        for i, field in enumerate(param_fields):
            if parameter_type == "bpstep_params":
                leg_val = None
                if 'params' in leg_rec and isinstance(leg_rec['params'], dict):
                    leg_val = leg_rec['params'].get(field.capitalize())
                mod_val = mod_rec.get(field)
            else:
                leg_val = None
                if 'params' in leg_rec and isinstance(leg_rec['params'], list):
                    leg_val = leg_rec['params'][i] if i < len(leg_rec['params']) else None
                mod_val = mod_rec.get(field)

            if leg_val is not None and mod_val is not None:
                if abs(float(leg_val) - float(mod_val)) > tolerance:
                    if abs(float(leg_val) + float(mod_val)) > tolerance:
                        step_mismatches[field] = {
                            'legacy': leg_val,
                            'modern': mod_val,
                            'diff': abs(float(leg_val) - float(mod_val))
                        }

        if step_mismatches:
            mismatches.append({
                'bp_idx1': key[0],
                'bp_idx2': key[1],
                'mismatches': step_mismatches,
                'legacy_record': leg_rec,
                'modern_record': mod_rec
            })

    # If bp_idx matching has errors, try mst_org position matching
    if len(mismatches) > 0:
        position_mismatches, _ = _match_steps_by_position(
            list(legacy_map.values()),
            list(modern_map.values()),
            parameter_type,
            tolerance
        )
        if len(position_mismatches) < len(mismatches):
            mismatches = position_mismatches

    return mismatches


def _match_steps_by_position(
    legacy_steps: List[Dict],
    modern_steps: List[Dict],
    parameter_type: str,
    tolerance: float
) -> Tuple[List[Dict], int]:
    """
    Match steps by mst_org (midstep origin) position.

    This is a fallback matching method when res_id and bp_idx matching fail.
    """
    mismatches = []

    # Extract mst_org positions
    leg_orgs = []
    for step in legacy_steps:
        mst = step.get("mst_org")
        if mst:
            leg_orgs.append((step, np.array(mst)))
        else:
            leg_orgs.append((step, None))

    mod_orgs = []
    for step in modern_steps:
        mst_frame = step.get("midstep_frame", {})
        org = mst_frame.get("org")
        if org:
            mod_orgs.append((step, np.array(org)))
        else:
            mod_orgs.append((step, None))

    # Match steps by mst_org position
    MST_DIST_THRESHOLD = 1.0
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

    # Define parameter fields
    if parameter_type == "bpstep_params":
        param_fields = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
        legacy_keys = ['Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist']
    else:
        param_fields = ['x_displacement', 'y_displacement', 'rise', 'inclination', 'tip', 'twist']
        legacy_keys = None

    # Compare parameters for matched steps
    for lstep, mstep in matched_pairs:
        step_mismatches = {}

        for i, field in enumerate(param_fields):
            if parameter_type == "bpstep_params":
                leg_params = lstep.get("params", {})
                leg_val = leg_params.get(legacy_keys[i]) if isinstance(leg_params, dict) else None
            else:
                leg_params = lstep.get("params", [])
                leg_val = leg_params[i] if isinstance(leg_params, list) and i < len(leg_params) else None

            mod_val = mstep.get(field)

            if leg_val is None or mod_val is None:
                continue

            # Check normal match
            if abs(leg_val - mod_val) <= tolerance:
                continue

            # Check sign-inverted match
            if abs(leg_val + mod_val) <= tolerance:
                continue

            step_mismatches[field] = {
                'legacy': leg_val,
                'modern': mod_val,
                'diff': abs(leg_val - mod_val)
            }

        if step_mismatches:
            mismatches.append({
                'bp_idx1': lstep.get('bp_idx1'),
                'bp_idx2': lstep.get('bp_idx2'),
                'mismatches': step_mismatches,
                'legacy_record': lstep,
                'modern_record': mstep
            })

    return mismatches, len(matched_pairs)


def compare_steps_by_residue_pairs(
    legacy_steps: List[Dict],
    modern_steps: List[Dict],
    legacy_pairs: List[Dict],
    modern_pairs: List[Dict],
    parameter_type: str = "bpstep_params",
    tolerance: float = 1e-4
) -> Tuple[List[Dict], int, int]:
    """
    Compare step parameters by matching steps that use the same residue pairs.

    This is a legacy function kept for backwards compatibility.
    Prefer compare_step_parameters with res_id matching.
    """
    def get_step_residue_key(step: Dict, pairs: List[Dict]) -> Optional[Tuple]:
        bp_idx1 = step.get('bp_idx1', 0)
        bp_idx2 = step.get('bp_idx2', 0)

        idx1 = bp_idx1 - 1
        idx2 = bp_idx2 - 1

        if idx1 < 0 or idx1 >= len(pairs) or idx2 < 0 or idx2 >= len(pairs):
            return None

        pair1 = pairs[idx1]
        pair2 = pairs[idx2]

        res1 = (pair1.get('base_i'), pair1.get('base_j'))
        res2 = (pair2.get('base_i'), pair2.get('base_j'))

        if None in res1 or None in res2:
            return None

        return (res1, res2)

    legacy_by_residues = {}
    for step in legacy_steps:
        key = get_step_residue_key(step, legacy_pairs)
        if key:
            legacy_by_residues[key] = step

    modern_by_residues = {}
    for step in modern_steps:
        key = get_step_residue_key(step, modern_pairs)
        if key:
            modern_by_residues[key] = step

    common_keys = set(legacy_by_residues.keys()) & set(modern_by_residues.keys())

    if parameter_type == "bpstep_params":
        param_fields = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
        legacy_param_keys = ['Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist']
    else:
        param_fields = ['x_displacement', 'y_displacement', 'rise', 'inclination', 'tip', 'twist']
        legacy_param_keys = None

    mismatches = []
    matched_count = 0

    for key in common_keys:
        lstep = legacy_by_residues[key]
        mstep = modern_by_residues[key]
        step_mismatches = {}

        for i, field in enumerate(param_fields):
            if parameter_type == "bpstep_params":
                leg_params = lstep.get("params", {})
                leg_val = leg_params.get(legacy_param_keys[i]) if isinstance(leg_params, dict) else None
            else:
                leg_params = lstep.get("params", [])
                leg_val = leg_params[i] if isinstance(leg_params, list) and i < len(leg_params) else None

            mod_val = mstep.get(field)

            if leg_val is None or mod_val is None:
                continue

            if abs(float(leg_val) - float(mod_val)) > tolerance:
                step_mismatches[field] = {
                    'legacy': leg_val,
                    'modern': mod_val,
                    'diff': abs(float(leg_val) - float(mod_val))
                }

        if step_mismatches:
            mismatches.append({
                'residue_key': key,
                'legacy_bp_idx': (lstep.get('bp_idx1'), lstep.get('bp_idx2')),
                'modern_bp_idx': (mstep.get('bp_idx1'), mstep.get('bp_idx2')),
                'mismatches': step_mismatches
            })
        else:
            matched_count += 1

    return mismatches, matched_count, len(common_keys)

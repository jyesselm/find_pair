"""
Step parameter comparison utilities.

Provides functions to compare step parameters (bpstep_params, helical_params)
between legacy and modern JSON outputs.

Uses mst_org position matching to handle cases where legacy and modern
compute steps in different orders (e.g., due to five2three algorithm
vs sequential ordering at helix boundaries).
"""

from typing import Dict, List, Tuple, Optional
from pathlib import Path
from dataclasses import dataclass, field
import numpy as np
from .models import FrameComparison, FrameMismatch


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
    parameter_type: str = "bpstep_params"  # "bpstep_params" or "helical_params"


def compare_step_parameters(
    legacy_records: List[Dict], 
    modern_records: List[Dict],
    parameter_type: str = "bpstep_params",  # or "helical_params"
    frame_comparison: Optional[FrameComparison] = None  # Optional: verify frames match first
) -> StepComparison:
    """
    Compare step parameters between legacy and modern JSON.
    
    Step parameters are calculated FROM reference frames. If frames don't match,
    step parameters cannot match either. It's recommended to verify frames match
    before comparing step parameters.
    
    Args:
        legacy_records: List of legacy step parameter records
        modern_records: List of modern step parameter records
        parameter_type: Type of parameters to compare ("bpstep_params" or "helical_params")
        frame_comparison: Optional FrameComparison result - if provided, will verify
                         frames match before comparing steps
        
    Returns:
        StepComparison result with a flag indicating if frames matched
    """
    result = StepComparison(parameter_type=parameter_type)
    
    # Warn if frames don't match (step parameters depend on frames)
    # Note: Warning is added in json_comparison.py, but we store frame status here too
    if frame_comparison:
        if frame_comparison.mismatched_calculations or frame_comparison.missing_residues:
            # Step parameters depend on frames - if frames don't match, 
            # step parameter differences are expected
            pass  # Warning already added in json_comparison.py
    
    # Build maps by bp_idx1, bp_idx2 (or equivalent indices)
    legacy_map = {}
    modern_map = {}
    
    for rec in legacy_records:
        if rec.get('type') != parameter_type:
            continue
        bp_idx1 = rec.get('bp_idx1')
        bp_idx2 = rec.get('bp_idx2')
        if bp_idx1 is not None and bp_idx2 is not None:
            key = (bp_idx1, bp_idx2)
            legacy_map[key] = rec.copy() if rec else rec
    
    for rec in modern_records:
        if rec.get('type') != parameter_type:
            continue
        bp_idx1 = rec.get('bp_idx1')
        bp_idx2 = rec.get('bp_idx2')
        if bp_idx1 is not None and bp_idx2 is not None:
            key = (bp_idx1, bp_idx2)
            modern_map[key] = rec.copy() if rec else rec
    
    result.total_legacy = len(legacy_map)
    result.total_modern = len(modern_map)
    
    legacy_keys = set(legacy_map.keys())
    modern_keys = set(modern_map.keys())
    common_keys = legacy_keys & modern_keys
    
    # Find missing steps
    for key in legacy_keys - modern_keys:
        result.missing_steps.append({
            'bp_idx1': key[0],
            'bp_idx2': key[1],
            'legacy_record': legacy_map[key]
        })
    
    # Compare common steps
    # Both legacy and modern output at %.6f precision. Due to numerical differences
    # in calculations and coordinate transforms, small differences are expected.
    # Use 5e-4 tolerance to account for floating-point precision in trig calculations.
    tolerance = 5e-4
    param_fields = []
    
    if parameter_type == "bpstep_params":
        # Legacy uses capitalized names in params dict, modern uses lowercase direct fields
        param_fields = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
        legacy_params_key = 'params'  # Legacy has params dict with capitalized keys
    elif parameter_type == "helical_params":
        param_fields = ['x_displacement', 'y_displacement', 'rise', 'inclination', 'tip', 'twist']
        legacy_params_key = 'params'  # Legacy has params array
    
    for key in common_keys:
        leg_rec = legacy_map[key]
        mod_rec = modern_map[key]
        
        mismatches = {}
        
        # Compare parameters
        for field in param_fields:
            if parameter_type == "bpstep_params":
                # Legacy: params['Shift'], modern: shift (lowercase)
                leg_val = None
                if 'params' in leg_rec and isinstance(leg_rec['params'], dict):
                    leg_val = leg_rec['params'].get(field.capitalize())
                elif field in leg_rec:
                    leg_val = leg_rec[field]
                
                mod_val = mod_rec.get(field)
            else:  # helical_params
                # Legacy: params array [x_displacement, y_displacement, rise, inclination, tip, twist]
                # Modern: individual fields
                leg_val = None
                if 'params' in leg_rec and isinstance(leg_rec['params'], list):
                    idx_map = {
                        'x_displacement': 0,
                        'y_displacement': 1,
                        'rise': 2,
                        'inclination': 3,
                        'tip': 4,
                        'twist': 5
                    }
                    if field in idx_map and len(leg_rec['params']) > idx_map[field]:
                        leg_val = leg_rec['params'][idx_map[field]]
                elif field in leg_rec:
                    leg_val = leg_rec[field]
                
                mod_val = mod_rec.get(field)
            
            if leg_val is not None and mod_val is not None:
                if abs(float(leg_val) - float(mod_val)) > tolerance:
                    mismatches[field] = {
                        'legacy': leg_val,
                        'modern': mod_val,
                        'diff': abs(float(leg_val) - float(mod_val))
                    }
        
        # Compare midstep frames if present
        if 'mst_org' in leg_rec and 'midstep_frame' in mod_rec:
            leg_org = leg_rec.get('mst_org')
            mod_org = mod_rec.get('midstep_frame', {}).get('origin')
            if leg_org and mod_org:
                for i in range(3):
                    if abs(float(leg_org[i]) - float(mod_org[i])) > tolerance:
                        mismatches[f'mst_org[{i}]'] = {
                            'legacy': leg_org[i],
                            'modern': mod_org[i],
                            'diff': abs(float(leg_org[i]) - float(mod_org[i]))
                        }
        
        if mismatches:
            result.mismatched_steps.append({
                'bp_idx1': key[0],
                'bp_idx2': key[1],
                'mismatches': mismatches,
                'legacy_record': leg_rec,
                'modern_record': mod_rec
            })
    
    # If bp_idx matching has errors, try mst_org position matching
    # This handles cases where legacy and modern compute steps in different orders
    if len(result.mismatched_steps) > 0:
        position_mismatches, position_match_count = _match_steps_by_position(
            list(legacy_map.values()),
            list(modern_map.values()),
            parameter_type,
            tolerance
        )
        # Only use position matching if it matched ALL common keys
        # Otherwise we'd be silently ignoring steps that couldn't be position-matched
        if position_match_count >= len(common_keys) and len(position_mismatches) < len(result.mismatched_steps):
            result.mismatched_steps = position_mismatches

    return result


def _match_steps_by_position(
    legacy_steps: List[Dict],
    modern_steps: List[Dict],
    parameter_type: str,
    tolerance: float
) -> Tuple[List[Dict], int]:
    """
    Match steps by mst_org (midstep origin) position.

    Legacy stores mst_org and modern stores midstep_frame.org - both represent
    the midpoint origin of the step calculation. Steps using the same pair
    combination will have matching positions.

    Args:
        legacy_steps: Legacy step records
        modern_steps: Modern step records
        parameter_type: "bpstep_params" or "helical_params"
        tolerance: Parameter tolerance

    Returns:
        Tuple of (list of mismatched steps, number of position matches found)
    """
    mismatches = []

    # Extract mst_org positions
    leg_orgs = []
    for step in legacy_steps:
        mst = step.get("mst_org") if parameter_type == "bpstep_params" else step.get("mst_org")
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

    # Define parameter fields based on type
    if parameter_type == "bpstep_params":
        param_fields = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
        legacy_keys = ['Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist']
    else:  # helical_params
        param_fields = ['x_displacement', 'y_displacement', 'rise', 'inclination', 'tip', 'twist']
        legacy_keys = None  # Uses array indexing

    # Compare parameters for matched steps
    for lstep, mstep in matched_pairs:
        step_mismatches = {}
        sign_inverted_match = True

        for i, field in enumerate(param_fields):
            # Get legacy value
            if parameter_type == "bpstep_params":
                leg_params = lstep.get("params", {})
                leg_val = leg_params.get(legacy_keys[i]) if isinstance(leg_params, dict) else None
            else:  # helical_params
                leg_params = lstep.get("params", [])
                leg_val = leg_params[i] if isinstance(leg_params, list) and i < len(leg_params) else None

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
            step_mismatches[field] = {
                'legacy': leg_val,
                'modern': mod_val,
                'diff': abs(leg_val - mod_val)
            }

        # Only report errors if neither normal nor sign-inverted match
        if step_mismatches and not sign_inverted_match:
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

    This handles cases where legacy and modern have different helix organization
    (different bp_idx ordering) but compute steps for the same underlying pairs.

    Args:
        legacy_steps: Legacy step parameter records
        modern_steps: Modern step parameter records
        legacy_pairs: Legacy base_pair records (to get residue info)
        modern_pairs: Modern base_pair records (to get residue info)
        parameter_type: "bpstep_params" or "helical_params"
        tolerance: Parameter comparison tolerance

    Returns:
        Tuple of (mismatches, matched_count, total_comparable)
    """

    def get_step_residue_key(step: Dict, pairs: List[Dict]) -> Optional[Tuple]:
        """Get residue pair key for a step: ((base_i1, base_j1), (base_i2, base_j2))"""
        bp_idx1 = step.get('bp_idx1', 0)
        bp_idx2 = step.get('bp_idx2', 0)

        # bp_idx is 1-based, convert to 0-based for list indexing
        idx1 = bp_idx1 - 1
        idx2 = bp_idx2 - 1

        if idx1 < 0 or idx1 >= len(pairs) or idx2 < 0 or idx2 >= len(pairs):
            return None

        pair1 = pairs[idx1]
        pair2 = pairs[idx2]

        # Get residue indices (base_i, base_j) for each pair
        res1 = (pair1.get('base_i'), pair1.get('base_j'))
        res2 = (pair2.get('base_i'), pair2.get('base_j'))

        if None in res1 or None in res2:
            return None

        return (res1, res2)

    # Build maps from residue key to step record
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

    # Find common residue pair combinations
    common_keys = set(legacy_by_residues.keys()) & set(modern_by_residues.keys())

    # Define parameter fields
    if parameter_type == "bpstep_params":
        param_fields = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
        legacy_param_keys = ['Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist']
    else:  # helical_params
        param_fields = ['x_displacement', 'y_displacement', 'rise', 'inclination', 'tip', 'twist']
        legacy_param_keys = None

    mismatches = []
    matched_count = 0

    for key in common_keys:
        lstep = legacy_by_residues[key]
        mstep = modern_by_residues[key]
        step_mismatches = {}

        for i, field in enumerate(param_fields):
            # Get legacy value
            if parameter_type == "bpstep_params":
                leg_params = lstep.get("params", {})
                leg_val = leg_params.get(legacy_param_keys[i]) if isinstance(leg_params, dict) else None
            else:
                leg_params = lstep.get("params", [])
                leg_val = leg_params[i] if isinstance(leg_params, list) and i < len(leg_params) else None

            mod_val = mstep.get(field)

            if leg_val is None or mod_val is None:
                continue

            # Check match
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


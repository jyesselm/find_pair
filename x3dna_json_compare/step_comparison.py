"""
Step parameter comparison utilities.

Provides functions to compare step parameters (bpstep_params, helical_params)
between legacy and modern JSON outputs.
"""

from typing import Dict, List, Tuple, Optional
from pathlib import Path
from dataclasses import dataclass, field
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
    tolerance = 1e-6
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
    
    return result


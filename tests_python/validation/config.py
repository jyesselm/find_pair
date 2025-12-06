"""
Stage configuration and tolerance definitions.

This module defines:
- StageConfig: Configuration for each validation stage
- STAGES: Registry of all 12 stages
- STAGE_GROUPS: Convenience groups (atoms, frames, pairs, etc.)
- Tolerance: Standard tolerance values
"""

from dataclasses import dataclass, field
from typing import List, Dict


class Tolerance:
    """Standard tolerance values for comparisons."""
    
    COORDINATE = 1e-6      # xyz values
    DISTANCE = 1e-6        # dorg, dNN, d_v
    ANGLE = 1e-6           # plane_angle
    MATRIX = 1e-4          # rotation matrix elements
    RMS = 0.001            # RMS fit values


@dataclass
class StageConfig:
    """Configuration for a validation stage."""
    
    stage_num: int
    stage_id: str
    name: str
    json_type: str
    required_fields: List[str]
    dependencies: List[int] = field(default_factory=list)
    tolerance: float = Tolerance.COORDINATE


# Stage definitions aligned with docs/JSON_DATA_TYPES_AND_COMPARISONS.md
STAGES: Dict[int, StageConfig] = {
    1: StageConfig(
        stage_num=1,
        stage_id="pdb_atoms",
        name="Atom Parsing",
        json_type="pdb_atoms",
        required_fields=["xyz", "atom_name", "residue_name", "chain_id"],
    ),
    2: StageConfig(
        stage_num=2,
        stage_id="residue_indices",
        name="Residue Indices",
        json_type="residue_indices",
        required_fields=["legacy_residue_idx", "start_atom_idx", "end_atom_idx"],
        dependencies=[1],
    ),
    3: StageConfig(
        stage_num=3,
        stage_id="base_frame_calc",
        name="Base Frame Calc",
        json_type="base_frame_calc",
        required_fields=["base_type", "rms_fit", "num_matched_atoms", "matched_atoms"],
        dependencies=[1, 2],
        tolerance=Tolerance.RMS,
    ),
    4: StageConfig(
        stage_num=4,
        stage_id="ls_fitting",
        name="LS Fitting",
        json_type="ls_fitting",
        required_fields=["base_type", "rms_fit", "num_points", "rotation_matrix", "translation"],
        dependencies=[1, 2],
        tolerance=Tolerance.RMS,
    ),
    5: StageConfig(
        stage_num=5,
        stage_id="frame_calc",
        name="Frame Calc",
        json_type="frame_calc",
        required_fields=["base_type", "orien", "org"],
        dependencies=[3, 4],
        tolerance=Tolerance.MATRIX,
    ),
    6: StageConfig(
        stage_num=6,
        stage_id="pair_validation",
        name="Pair Validation",
        json_type="pair_validation",
        required_fields=["is_valid", "bp_type_id"],
        dependencies=[5],
    ),
    7: StageConfig(
        stage_num=7,
        stage_id="distance_checks",
        name="Distance Checks",
        json_type="distance_checks",
        required_fields=["dorg", "dNN", "plane_angle", "d_v"],
        dependencies=[5],
    ),
    8: StageConfig(
        stage_num=8,
        stage_id="hbond_list",
        name="H-bond List",
        json_type="hbond_list",
        required_fields=["num_hbonds", "hbonds"],
        dependencies=[6, 7],
    ),
    9: StageConfig(
        stage_num=9,
        stage_id="base_pair",
        name="Base Pair",
        json_type="base_pair",
        required_fields=["bp_type", "orien_i", "org_i"],
        dependencies=[6, 7, 8],
    ),
    10: StageConfig(
        stage_num=10,
        stage_id="find_bestpair_selection",
        name="Best Pair Selection",
        json_type="find_bestpair_selection",
        required_fields=["num_bp", "pairs"],
        dependencies=[9],
    ),
    11: StageConfig(
        stage_num=11,
        stage_id="bpstep_params",
        name="Step Parameters",
        json_type="bpstep_params",
        required_fields=["shift", "slide", "rise", "tilt", "roll", "twist"],
        dependencies=[10],
    ),
    12: StageConfig(
        stage_num=12,
        stage_id="helical_params",
        name="Helical Parameters",
        json_type="helical_params",
        required_fields=["x_displacement", "y_displacement", "rise", "inclination", "tip", "twist"],
        dependencies=[10, 11],
    ),
}

# Stage groups for convenience
STAGE_GROUPS: Dict[str, List[int]] = {
    "atoms": [1],
    "residue": [2],
    "frames": [3, 4, 5],
    "pairs": [6, 7, 9, 10],
    "hbonds": [8],
    "steps": [11, 12],
    "all": list(range(1, 13)),
}

# Map stage IDs to numbers
STAGE_ID_TO_NUM: Dict[str, int] = {
    cfg.stage_id: cfg.stage_num for cfg in STAGES.values()
}


def resolve_stages(stage_args: List[str]) -> List[int]:
    """
    Resolve stage arguments to list of stage numbers.
    
    Args:
        stage_args: List of stage identifiers (numbers, names, or groups)
        
    Returns:
        Sorted list of unique stage numbers
    """
    if not stage_args:
        return list(range(1, 13))
    
    stages = []
    for arg in stage_args:
        arg_lower = arg.lower()
        
        if arg_lower in STAGE_GROUPS:
            stages.extend(STAGE_GROUPS[arg_lower])
        elif arg.isdigit() and 1 <= int(arg) <= 12:
            stages.append(int(arg))
        elif arg_lower in STAGE_ID_TO_NUM:
            stages.append(STAGE_ID_TO_NUM[arg_lower])
    
    return sorted(set(stages))


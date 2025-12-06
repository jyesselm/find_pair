"""
Comparators for each validation stage.

Each comparator function takes legacy and modern records and returns
(passed: bool, errors: List[str]).
"""

from .atoms import compare_atoms
from .residues import compare_residues
from .frames import compare_base_frame_calc, compare_ls_fitting, compare_frame_calc
from .pairs import compare_pair_validation, compare_distance_checks, compare_base_pair
from .pairs import compare_best_pair_selection
from .hbonds import compare_hbonds
from .params import compare_step_params, compare_helical_params

# Map stage numbers to their comparator functions
COMPARATORS = {
    1: compare_atoms,
    2: compare_residues,
    3: compare_base_frame_calc,
    4: compare_ls_fitting,
    5: compare_frame_calc,
    6: compare_pair_validation,
    7: compare_distance_checks,
    8: compare_hbonds,
    9: compare_base_pair,
    10: compare_best_pair_selection,
    11: compare_step_params,
    12: compare_helical_params,
}

__all__ = [
    "COMPARATORS",
    "compare_atoms",
    "compare_residues",
    "compare_base_frame_calc",
    "compare_ls_fitting",
    "compare_frame_calc",
    "compare_pair_validation",
    "compare_distance_checks",
    "compare_hbonds",
    "compare_base_pair",
    "compare_best_pair_selection",
    "compare_step_params",
    "compare_helical_params",
]


"""
Validation package for comparing legacy vs modern JSON outputs.

Usage:
    from validation import validate_stage, StageConfig, STAGES
    from validation.comparators import compare_atoms, compare_frames
"""

from .config import STAGES, STAGE_GROUPS, StageConfig, Tolerance
from .runner import validate_pdb, validate_stage, ValidationResult
from .cli import main as run_cli

__all__ = [
    "STAGES",
    "STAGE_GROUPS", 
    "StageConfig",
    "Tolerance",
    "validate_pdb",
    "validate_stage",
    "ValidationResult",
    "run_cli",
]


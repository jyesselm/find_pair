"""Geometric validation for base pairs.

This package provides clean, focused modules for validating base pair geometry:
- result: ValidationResult dataclass with all metrics
- thresholds: Configurable ValidationThresholds
- validator: GeometricValidator class

Example:
    from validation import GeometricValidator, ValidationThresholds

    validator = GeometricValidator(ValidationThresholds.strict())
    result = validator.validate(frame1, frame2, n1_pos, n9_pos)

    if result.is_valid:
        print(f"Quality: {result.quality_score:.2f}")
    else:
        print(f"Failed: {result.rejection_reason}")
"""

from prototypes.pair_identification.validation.result import ValidationResult
from prototypes.pair_identification.validation.thresholds import ValidationThresholds
from prototypes.pair_identification.validation.validator import GeometricValidator

__all__ = [
    "ValidationResult",
    "ValidationThresholds",
    "GeometricValidator",
]

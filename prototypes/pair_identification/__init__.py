"""Pair identification prototype.

Python prototype for base pair identification algorithms, using pre-computed
reference frames from the modern C++ implementation.

Modules:
    frame_loader: Load reference frames from JSON output
    geometric_validator: Geometric validation for base pairs
    pair_cache: Cache pre-computed frames and validation results
"""

from prototypes.pair_identification.frame_loader import (
    ReferenceFrame,
    FrameLoader,
)
from prototypes.pair_identification.geometric_validator import (
    GeometricValidator,
    ValidationResult,
)
from prototypes.pair_identification.pair_cache import (
    AtomCoords,
    CachedPair,
    PairCache,
)

__all__ = [
    "ReferenceFrame",
    "FrameLoader",
    "GeometricValidator",
    "ValidationResult",
    "AtomCoords",
    "CachedPair",
    "PairCache",
]

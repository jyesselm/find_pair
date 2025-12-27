"""Pair identification prototype.

Python prototype for base pair identification algorithms, using pre-computed
reference frames from the modern C++ implementation.

High-level API:
    from pair_identification import PairFinder, FinderConfig
    finder = PairFinder(Path("data/json"))
    result = finder.find_pairs("1EHZ")

Low-level API:
    from pair_identification.validation import GeometricValidator
    from pair_identification.hbond import HBondFinder
    from pair_identification.templates import LWClassifier

Packages:
    core: Data structures, PDB parsing, alignment
    validation: Geometric validation
    hbond: H-bond detection and patterns
    templates: LW classification via template alignment
    finder: Pair finding pipeline
    analysis: Diagnostic and reporting tools
"""

# Frame loading
from prototypes.pair_identification.frame_loader import (
    ReferenceFrame,
    FrameLoader,
)

# Validation
from prototypes.pair_identification.validation import (
    GeometricValidator,
    ValidationResult,
    ValidationThresholds,
)

# Legacy pair cache (for backwards compatibility)
from prototypes.pair_identification.pair_cache import (
    AtomCoords,
    CachedPair,
    PairCache as LegacyPairCache,
)

# Finder - high-level API
from prototypes.pair_identification.finder import (
    PairFinder,
    FinderConfig,
    FinderResult,
    QualityScorer,
    MutualBestStrategy,
)

__all__ = [
    # Frame loading
    "ReferenceFrame",
    "FrameLoader",
    # Validation
    "GeometricValidator",
    "ValidationResult",
    "ValidationThresholds",
    # Legacy (backwards compatibility)
    "AtomCoords",
    "CachedPair",
    "LegacyPairCache",
    # Finder - high-level API
    "PairFinder",
    "FinderConfig",
    "FinderResult",
    "QualityScorer",
    "MutualBestStrategy",
]

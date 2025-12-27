"""Pair finding components.

This package provides the pair finding pipeline including:
- PairCache: Spatial indexing and caching of validated pairs
- QualityScorer: Composite scoring for pair ranking
- MutualBestStrategy: Greedy selection with mutual best criterion
- PairFinder: High-level facade combining all components
"""

from prototypes.pair_identification.finder.cache import (
    AtomCoords,
    CandidateInfo,
    PairCache,
)
from prototypes.pair_identification.finder.quality_scorer import (
    QualityScorer,
    EXPECTED_HBOND_COUNT,
)
from prototypes.pair_identification.finder.strategy import (
    MutualBestStrategy,
    GreedyBestStrategy,
    SelectionResult,
)
from prototypes.pair_identification.finder.finder import (
    PairFinder,
    FinderConfig,
    FinderResult,
)

__all__ = [
    # Cache
    "AtomCoords",
    "CandidateInfo",
    "PairCache",
    # Scoring
    "QualityScorer",
    "EXPECTED_HBOND_COUNT",
    # Strategy
    "MutualBestStrategy",
    "GreedyBestStrategy",
    "SelectionResult",
    # Finder facade
    "PairFinder",
    "FinderConfig",
    "FinderResult",
]

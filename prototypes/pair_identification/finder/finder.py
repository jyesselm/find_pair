"""High-level pair finder facade.

This module provides PairFinder, a high-level interface that combines:
- PairCache: Spatial indexing and geometric validation
- QualityScorer: Composite scoring for pair ranking
- MutualBestStrategy: Greedy selection with mutual criterion
- Optional H-bond and LW classification

Matches the C++ BasePairFinder architecture.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from prototypes.pair_identification.finder.cache import CandidateInfo, PairCache
from prototypes.pair_identification.finder.quality_scorer import QualityScorer
from prototypes.pair_identification.finder.strategy import (
    MutualBestStrategy,
    SelectionResult,
)
from prototypes.pair_identification.validation import (
    GeometricValidator,
    ValidationThresholds,
)


@dataclass
class FinderConfig:
    """Configuration for pair finding.

    Attributes:
        max_distance: Maximum distance between frame origins (Angstroms).
        min_score: Minimum quality score to accept.
        require_mutual: Require mutual best criterion.
        rmsd_weight: Weight for RMSD component in scoring.
        coverage_weight: Weight for H-bond coverage in scoring.
        quality_weight: Weight for H-bond quality in scoring.
        validation_thresholds: Geometric validation thresholds.
    """

    max_distance: float = 15.0
    min_score: float = 0.0
    require_mutual: bool = True
    rmsd_weight: float = 0.3
    coverage_weight: float = 0.4
    quality_weight: float = 0.3
    validation_thresholds: Optional[ValidationThresholds] = None

    @classmethod
    def default(cls) -> "FinderConfig":
        """Create default configuration."""
        return cls()

    @classmethod
    def strict(cls) -> "FinderConfig":
        """Create strict configuration with higher thresholds."""
        return cls(
            min_score=0.5,
            require_mutual=True,
            validation_thresholds=ValidationThresholds.strict(),
        )


@dataclass
class FinderResult:
    """Result of pair finding.

    Attributes:
        pdb_id: PDB identifier.
        pairs: Selected pairs.
        candidates_total: Total candidates found.
        candidates_valid: Candidates passing validation.
        selection_result: Detailed selection result.
    """

    pdb_id: str
    pairs: List[CandidateInfo]
    candidates_total: int
    candidates_valid: int
    selection_result: Optional[SelectionResult] = None


class PairFinder:
    """High-level pair finder combining all components.

    This class provides a simple interface for finding base pairs:
    1. Load pre-computed frames from JSON
    2. Find candidates within distance cutoff
    3. Validate geometry
    4. Score candidates
    5. Select pairs using mutual best strategy

    Example:
        finder = PairFinder(Path("data/json"))
        result = finder.find_pairs("1EHZ")
        for pair in result.pairs:
            print(f"{pair.res_id1}-{pair.res_id2}: {pair.quality_score:.2f}")

    Attributes:
        json_dir: Directory containing JSON output.
        config: Finder configuration.
        validator: Geometric validator.
        scorer: Quality scorer.
        strategy: Selection strategy.
    """

    def __init__(
        self,
        json_dir: Path | str,
        config: Optional[FinderConfig] = None,
    ):
        """Initialize pair finder.

        Args:
            json_dir: Directory containing JSON output.
            config: Finder configuration (default: FinderConfig.default()).
        """
        self.json_dir = Path(json_dir)
        self.config = config or FinderConfig.default()

        thresholds = self.config.validation_thresholds
        self.validator = GeometricValidator(thresholds)

        self.scorer = QualityScorer(
            rmsd_weight=self.config.rmsd_weight,
            coverage_weight=self.config.coverage_weight,
            quality_weight=self.config.quality_weight,
        )

        self.strategy = MutualBestStrategy(
            min_score=self.config.min_score,
            require_mutual=self.config.require_mutual,
        )

    def find_pairs(self, pdb_id: str) -> FinderResult:
        """Find base pairs for a PDB.

        Args:
            pdb_id: PDB identifier (e.g., "1EHZ").

        Returns:
            FinderResult with selected pairs and statistics.

        Raises:
            FileNotFoundError: If required JSON files don't exist.
            ValueError: If JSON format is invalid.
        """
        cache = PairCache(pdb_id, self.json_dir)
        cache.build(max_distance=self.config.max_distance, validator=self.validator)

        self._score_candidates(cache)

        valid_candidates = cache.get_valid_candidates()
        selection_result = self.strategy.select_with_details(valid_candidates)

        return FinderResult(
            pdb_id=pdb_id,
            pairs=selection_result.selected,
            candidates_total=len(cache.candidates),
            candidates_valid=len(valid_candidates),
            selection_result=selection_result,
        )

    def _score_candidates(self, cache: PairCache) -> None:
        """Score all candidates in cache.

        Args:
            cache: PairCache with candidates to score.
        """
        for candidate in cache.candidates:
            if not candidate.validation.is_valid:
                continue

            score = self.scorer.compute_score(
                validation=candidate.validation,
                sequence=candidate.sequence,
                hbonds=None,
                rmsd=None,
            )
            candidate.quality_score = score

    def find_pairs_batch(
        self, pdb_ids: List[str], max_workers: Optional[int] = None
    ) -> Dict[str, FinderResult]:
        """Find pairs for multiple PDBs in parallel.

        Args:
            pdb_ids: List of PDB identifiers.
            max_workers: Maximum parallel workers (default: all CPUs).

        Returns:
            Dict mapping pdb_id to FinderResult.
        """
        import os
        from concurrent.futures import ProcessPoolExecutor, as_completed

        if max_workers is None:
            max_workers = os.cpu_count() or 1

        results: Dict[str, FinderResult] = {}

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(self._find_pairs_single, pdb_id): pdb_id
                for pdb_id in pdb_ids
            }

            for future in as_completed(futures):
                pdb_id = futures[future]
                try:
                    results[pdb_id] = future.result()
                except Exception as e:
                    print(f"Error processing {pdb_id}: {e}")

        return results

    def _find_pairs_single(self, pdb_id: str) -> FinderResult:
        """Find pairs for single PDB (for parallel execution)."""
        return self.find_pairs(pdb_id)

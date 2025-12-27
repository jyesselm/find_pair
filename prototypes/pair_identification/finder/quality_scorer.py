"""Quality score calculation for base pairs.

Computes composite scores (0-1) based on:
- RMSD to expected template
- H-bond coverage (found vs expected)
- H-bond quality (distance and alignment)

Aligns with C++ BasePairValidator quality scoring.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np

from prototypes.pair_identification.validation import ValidationResult


EXPECTED_HBOND_COUNT = {
    "GC": 3,
    "CG": 3,
    "AU": 2,
    "UA": 2,
    "GU": 2,
    "UG": 2,
}

IDEAL_DISTANCE_MIN = 2.7
IDEAL_DISTANCE_MAX = 3.2


class QualityScorer:
    """Calculate quality scores for pair ranking.

    Scores combine geometric validation metrics with H-bond analysis to
    produce a 0-1 quality score (higher is better) for greedy selection.
    """

    def __init__(
        self,
        rmsd_weight: float = 0.3,
        coverage_weight: float = 0.4,
        quality_weight: float = 0.3,
    ):
        """Initialize quality scorer.

        Args:
            rmsd_weight: Weight for RMSD component (default: 0.3).
            coverage_weight: Weight for H-bond coverage (default: 0.4).
            quality_weight: Weight for H-bond quality (default: 0.3).
        """
        self.rmsd_weight = rmsd_weight
        self.coverage_weight = coverage_weight
        self.quality_weight = quality_weight

    def compute_score(
        self,
        validation: ValidationResult,
        sequence: str,
        hbonds: Optional[List[Dict]] = None,
        rmsd: Optional[float] = None,
    ) -> float:
        """Compute quality score from validation and H-bonds.

        Args:
            validation: Geometric validation result.
            sequence: Two-letter sequence (e.g., "GC").
            hbonds: Optional list of H-bond dictionaries.
            rmsd: Optional RMSD to template (uses validation.quality_score if None).

        Returns:
            Quality score (0-1, higher is better).
        """
        if not validation.is_valid:
            return 0.0

        if rmsd is None:
            rmsd = validation.quality_score / 10.0

        if hbonds is None:
            hbonds = []

        base_base_hbonds = [
            hb for hb in hbonds if hb.get("context") == "base_base"
        ]

        total_score, _ = self.compute_bp_score(
            sequence=sequence,
            rmsd=rmsd,
            found_hbonds=base_base_hbonds,
        )

        return total_score

    def compute_bp_score(
        self, sequence: str, rmsd: float, found_hbonds: List[Dict]
    ) -> Tuple[float, Dict[str, float]]:
        """Compute BP score (0-1) for pair grading.

        Args:
            sequence: Base pair sequence (e.g., "GC", "AU").
            rmsd: RMSD to cWW template in Angstroms.
            found_hbonds: List of detected H-bonds with distance, alignment.

        Returns:
            Tuple of (total_score, component_scores_dict).
            total_score: 0.0 (worst) to 1.0 (best).
            components: Dict with 'rmsd', 'coverage', 'quality' scores.
        """
        expected_count = EXPECTED_HBOND_COUNT.get(sequence, 2)
        found_count = len(found_hbonds)

        rmsd_score = self._compute_rmsd_score(rmsd)
        coverage_score = self._compute_coverage_score(found_count, expected_count)
        quality_score = self._compute_hbond_quality(found_hbonds, rmsd)

        total_score = (
            self.rmsd_weight * rmsd_score
            + self.coverage_weight * coverage_score
            + self.quality_weight * quality_score
        )

        components = {
            "rmsd": round(rmsd_score, 3),
            "coverage": round(coverage_score, 3),
            "quality": round(quality_score, 3),
        }

        return round(total_score, 3), components

    def _compute_rmsd_score(self, rmsd: float) -> float:
        """Compute RMSD component score.

        Args:
            rmsd: RMSD in Angstroms.

        Returns:
            Score 0-1 (1 is perfect).
        """
        if rmsd <= 0.3:
            return 1.0
        if rmsd >= 1.0:
            return 0.0
        return 1.0 - (rmsd - 0.3) / 0.7

    def _compute_coverage_score(self, found_count: int, expected_count: int) -> float:
        """Compute H-bond coverage score.

        Args:
            found_count: Number of H-bonds found.
            expected_count: Number of H-bonds expected.

        Returns:
            Score 0-1 (1 is full coverage).
        """
        if expected_count == 0:
            return 0.0
        return min(found_count / expected_count, 1.0)

    def _compute_hbond_quality(self, hbonds: List[Dict], rmsd: float) -> float:
        """Compute H-bond quality score from distance and alignment.

        Args:
            hbonds: List of H-bond dictionaries with distance, alignment.
            rmsd: RMSD for geometry leniency calculation.

        Returns:
            Average quality score 0-1.
        """
        if not hbonds:
            return 0.0

        geometry_leniency = self._compute_geometry_leniency(rmsd)
        quality_scores = []

        for hb in hbonds:
            dist = hb.get("distance", 3.0)
            alignment = hb.get("alignment", 1.0)

            dist_score = self._compute_distance_score(dist, geometry_leniency)
            align_score = self._compute_alignment_score(alignment)

            hb_quality = 0.7 * dist_score + 0.3 * align_score
            quality_scores.append(hb_quality)

        return sum(quality_scores) / len(quality_scores)

    def _compute_geometry_leniency(self, rmsd: float) -> float:
        """Compute geometry-based leniency for H-bond distances.

        Args:
            rmsd: RMSD in Angstroms.

        Returns:
            Leniency factor 0-1.
        """
        if rmsd <= 0.5:
            return 1.0
        if rmsd >= 0.8:
            return 0.0
        return 1.0 - (rmsd - 0.5) / 0.3

    def _compute_distance_score(self, dist: float, leniency: float) -> float:
        """Compute H-bond distance score with geometry leniency.

        Args:
            dist: H-bond distance in Angstroms.
            leniency: Geometry leniency factor 0-1.

        Returns:
            Distance score 0-1.
        """
        if IDEAL_DISTANCE_MIN <= dist <= IDEAL_DISTANCE_MAX:
            return 1.0

        if dist < IDEAL_DISTANCE_MIN:
            return max(0.5, 1.0 - (IDEAL_DISTANCE_MIN - dist) / 0.5)

        lenient_max = IDEAL_DISTANCE_MAX + (1.0 * leniency)
        if dist <= lenient_max:
            return 1.0

        return max(0.0, 1.0 - (dist - lenient_max) / 0.5)

    def _compute_alignment_score(self, alignment: float) -> float:
        """Compute H-bond alignment score.

        Args:
            alignment: Alignment metric (0-2 scale, lower is better).

        Returns:
            Alignment score 0-1.
        """
        if alignment <= 1.0:
            return 1.0
        if alignment >= 2.0:
            return 0.0
        return 1.0 - (alignment - 1.0)

    def score_to_grade(self, score: float) -> str:
        """Convert numeric score to letter grade.

        Args:
            score: 0.0 to 1.0 score.

        Returns:
            Letter grade: A, B, C, D, or F.
        """
        if score >= 0.9:
            return "A"
        if score >= 0.8:
            return "B"
        if score >= 0.7:
            return "C"
        if score >= 0.6:
            return "D"
        return "F"

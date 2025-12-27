"""Combined LW classifier using template RMSD + H-bond scoring."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from prototypes.pair_identification.core.residue import Residue
from prototypes.pair_identification.hbond.patterns import get_expected_hbonds
from prototypes.pair_identification.hbond_scorer import HBondScorer, HBond
from prototypes.pair_identification.templates.aligner import TemplateAligner


@dataclass
class LWClassificationResult:
    """Result of classifying a pair using RMSD + H-bond scoring.

    Attributes:
        res_id1: First residue ID.
        res_id2: Second residue ID.
        sequence: Target sequence.
        best_lw: Best matching LW class.
        best_rmsd: RMSD of best match (Angstroms).
        best_hbond_score: H-bond pattern match score [0, 1].
        best_combined_score: Combined RMSD + H-bond score.
        second_lw: Second best LW class.
        second_rmsd: RMSD of second best.
        second_hbond_score: H-bond score of second best.
        second_combined_score: Combined score of second best.
        all_results: All scored results sorted by combined score.
    """

    res_id1: str
    res_id2: str
    sequence: str
    best_lw: str
    best_rmsd: float
    best_hbond_score: float
    best_combined_score: float
    second_lw: Optional[str] = None
    second_rmsd: Optional[float] = None
    second_hbond_score: Optional[float] = None
    second_combined_score: Optional[float] = None
    all_results: List[dict] = field(default_factory=list)

    @property
    def confidence(self) -> float:
        """Compute confidence from gap between best and second best.

        Returns:
            Confidence value in [0, 1], where 1 is highest confidence.
        """
        if self.second_combined_score is None:
            return 1.0
        gap = self.best_combined_score - self.second_combined_score
        return min(1.0, gap / 0.3)

    @property
    def rmsd_gap(self) -> Optional[float]:
        """Gap between best and second best RMSD."""
        if self.second_rmsd is None:
            return None
        return self.second_rmsd - self.best_rmsd

    def __repr__(self):
        """Compact string representation."""
        conf_str = f"{self.confidence:.0%}"
        return (
            f"LWClassificationResult({self.res_id1}-{self.res_id2} "
            f"{self.sequence}: best={self.best_lw} "
            f"RMSD={self.best_rmsd:.3f}Å "
            f"hbond={self.best_hbond_score:.0%} "
            f"score={self.best_combined_score:.3f} conf={conf_str})"
        )


class LWClassifier:
    """Classify base pairs using template RMSD + H-bond scoring.

    Combines geometric alignment (RMSD to templates) with hydrogen bond
    pattern matching to determine the most likely LW classification.
    """

    DEFAULT_LW_CLASSES = [
        "cWW", "tWW", "cWH", "tWH", "cWS", "tWS",
        "cHH", "tHH", "cHS", "tHS", "cSS", "tSS"
    ]

    def __init__(
        self,
        idealized_dir: Path = Path("basepair-idealized"),
        exemplar_dir: Path = Path("basepair-exemplars"),
        rmsd_weight: float = 0.4,
        hbond_weight: float = 0.6,
        rmsd_cutoff: float = 2.0,
    ):
        """Initialize LW classifier.

        Args:
            idealized_dir: Directory with idealized templates.
            exemplar_dir: Directory with exemplar templates.
            rmsd_weight: Weight for RMSD score component [0, 1].
            hbond_weight: Weight for H-bond score component [0, 1].
            rmsd_cutoff: RMSD above this gets score of 0 (Angstroms).
        """
        self.aligner = TemplateAligner(idealized_dir, exemplar_dir)
        self.hbond_scorer = HBondScorer()
        self.rmsd_weight = rmsd_weight
        self.hbond_weight = hbond_weight
        self.rmsd_cutoff = rmsd_cutoff

    def classify(
        self,
        res1: Residue,
        res2: Residue,
        hbonds: List[HBond],
        lw_classes: Optional[List[str]] = None,
    ) -> LWClassificationResult:
        """Classify a pair by trying all LW class templates.

        Args:
            res1: First residue.
            res2: Second residue.
            hbonds: Observed H-bonds between the residues.
            lw_classes: LW classes to try (default: all 12).

        Returns:
            LWClassificationResult with best match and confidence.
        """
        if lw_classes is None:
            lw_classes = self.DEFAULT_LW_CLASSES

        sequence = res1.base_type + res2.base_type

        alignment_result = self.aligner.classify_pair(
            res1, res2, lw_classes
        )

        scored_results = []

        for align in alignment_result.all_results:
            if align.rmsd == float('inf'):
                continue

            scored = self._score_alignment(align, hbonds)
            scored_results.append(scored)

        scored_results.sort(key=lambda x: -x["combined_score"])

        if not scored_results:
            return self._create_unknown_result(res1, res2, sequence)

        best = scored_results[0]
        second = scored_results[1] if len(scored_results) > 1 else None

        return LWClassificationResult(
            res_id1=res1.res_id,
            res_id2=res2.res_id,
            sequence=sequence,
            best_lw=best["lw_class"],
            best_rmsd=best["rmsd"],
            best_hbond_score=best["hbond_score"],
            best_combined_score=best["combined_score"],
            second_lw=second["lw_class"] if second else None,
            second_rmsd=second["rmsd"] if second else None,
            second_hbond_score=second["hbond_score"] if second else None,
            second_combined_score=second["combined_score"] if second else None,
            all_results=scored_results,
        )

    def _score_alignment(
        self, align, hbonds: List[HBond]
    ) -> Dict[str, float]:
        """Score a single alignment result with RMSD + H-bonds."""
        rmsd_score = self._rmsd_to_score(align.rmsd)

        hbond_seq = align.sequence.replace("(rev)", "")
        hbond_result = self.hbond_scorer.score_pattern(
            hbonds, align.lw_class, hbond_seq
        )

        hbond_count_bonus = hbond_result.matched_hbonds * 0.05
        combined = (
            self.hbond_weight * hbond_result.score
            + self.rmsd_weight * rmsd_score
            + hbond_count_bonus
        )

        return {
            "lw_class": align.lw_class,
            "sequence": align.sequence,
            "rmsd": align.rmsd,
            "rmsd_score": rmsd_score,
            "hbond_score": hbond_result.score,
            "hbond_matched": hbond_result.matched_hbonds,
            "hbond_expected": hbond_result.expected_hbonds,
            "combined_score": combined,
            "num_atoms": align.num_atoms_aligned,
        }

    def _rmsd_to_score(self, rmsd: float) -> float:
        """Convert RMSD to a 0-1 score (lower RMSD = higher score).

        Args:
            rmsd: RMSD value in Angstroms.

        Returns:
            Score in [0, 1] where 0 is poor match and 1 is perfect match.
        """
        if rmsd >= self.rmsd_cutoff:
            return 0.0
        return 1.0 - (rmsd / self.rmsd_cutoff)

    def _create_unknown_result(
        self, res1: Residue, res2: Residue, sequence: str
    ) -> LWClassificationResult:
        """Create result for pair with no valid templates."""
        return LWClassificationResult(
            res_id1=res1.res_id,
            res_id2=res2.res_id,
            sequence=sequence,
            best_lw="unknown",
            best_rmsd=float('inf'),
            best_hbond_score=0.0,
            best_combined_score=0.0,
        )

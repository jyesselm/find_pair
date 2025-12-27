"""Classify base pairs by matching to LW templates.

This module takes validated base pairs from the pair cache and classifies them
by matching against templates with statistical geometry and idealized coordinates.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json

import numpy as np

from prototypes.pair_identification.frame_loader import ReferenceFrame
from prototypes.pair_identification.validation import ValidationResult
from prototypes.pair_identification.pair_cache import CachedPair, PairCache, AtomCoords
from prototypes.pair_identification.template_generator import (
    TemplateGenerator,
    IdealizedTemplate,
    TemplateStats,
)


@dataclass
class MatchResult:
    """Result of matching a pair against a template.

    Attributes:
        sequence: Template sequence (e.g., "GC").
        lw_class: Leontis-Westhof class (e.g., "cWW").
        score: Overall match score (lower = better).
        geometry_score: Score from geometry comparison.
        coord_rmsd: RMSD to idealized coordinates (if available).
        n1n9_zscore: Z-score for N1N9 distance.
        angle_zscore: Z-score for interbase angle.
        planarity_zscore: Z-score for planarity.
        hbond_match: Whether H-bond pattern matches.
        confidence: Confidence level (0.0-1.0).
        has_idealized_coords: Whether template has coordinate data.
    """

    sequence: str
    lw_class: str
    score: float
    geometry_score: float
    coord_rmsd: Optional[float] = None
    n1n9_zscore: float = 0.0
    angle_zscore: float = 0.0
    planarity_zscore: float = 0.0
    hbond_match: bool = False
    confidence: float = 0.0
    has_idealized_coords: bool = False


@dataclass
class ClassificationResult:
    """Classification result for a base pair.

    Attributes:
        res1_id: Residue 1 ID.
        res2_id: Residue 2 ID.
        sequence: Pair sequence (e.g., "GC").
        best_match: Best matching template result.
        all_matches: All template matches, sorted by score.
        is_classified: Whether a confident match was found.
    """

    res1_id: str
    res2_id: str
    sequence: str
    best_match: Optional[MatchResult] = None
    all_matches: List[MatchResult] = field(default_factory=list)
    is_classified: bool = False


class PairIdentifier:
    """Classify base pairs using template matching.

    This class loads templates and uses them to classify validated base pairs
    by their Leontis-Westhof type. Classification considers:
    - Statistical geometry (N1N9 distance, interbase angle, planarity)
    - RMSD to idealized coordinates (when available)
    - Hydrogen bond patterns

    Example:
        identifier = PairIdentifier()
        identifier.load_templates(Path("data/templates.json"))

        cache = PairCache.load(Path("1EHZ_cache.json"))
        for pair in cache.get_valid_pairs():
            result = identifier.classify(pair)
            print(f"{pair.res1_id}-{pair.res2_id}: {result.best_match.lw_class}")
    """

    def __init__(self):
        """Initialize pair identifier."""
        self.generator: Optional[TemplateGenerator] = None
        self.templates: Dict[Tuple[str, str], IdealizedTemplate] = {}

        # Weights for scoring
        self.geometry_weight = 0.6
        self.coord_weight = 0.3
        self.hbond_weight = 0.1

        # Z-score weights
        self.n1n9_weight = 1.0
        self.angle_weight = 0.5
        self.planarity_weight = 0.3

        # Confidence thresholds
        self.high_confidence_score = 2.0
        self.min_confidence_score = 5.0

    def load_templates(self, path: Path | str) -> None:
        """Load templates from JSON file.

        Args:
            path: Path to templates.json file.
        """
        self.generator = TemplateGenerator.load_registry(path)
        self.templates = self.generator.templates

    def classify(
        self,
        pair: CachedPair,
        n1n9_pos1: Optional[np.ndarray] = None,
        n1n9_pos2: Optional[np.ndarray] = None,
        interbase_angle: Optional[float] = None,
        planarity: Optional[float] = None,
    ) -> ClassificationResult:
        """Classify a base pair.

        Args:
            pair: CachedPair with frames and validation.
            n1n9_pos1: N1/N9 position of residue 1.
            n1n9_pos2: N1/N9 position of residue 2.
            interbase_angle: Interbase angle in degrees (optional).
            planarity: Planarity measure (optional).

        Returns:
            ClassificationResult with matches ranked by score.
        """
        sequence = pair.res1_name + pair.res2_name

        # Compute N1N9 distance from validation or positions
        if n1n9_pos1 is not None and n1n9_pos2 is not None:
            n1n9_dist = np.linalg.norm(n1n9_pos1 - n1n9_pos2)
        else:
            n1n9_dist = pair.validation.dNN

        # Use interbase angle if provided, otherwise estimate from plane_angle
        if interbase_angle is None:
            interbase_angle = pair.validation.plane_angle

        # Default planarity from validation (d_v is similar)
        if planarity is None:
            planarity = abs(pair.validation.d_v)

        # Get matching templates for this sequence
        matching_templates = self._get_templates_for_sequence(sequence)

        # Score each template
        matches = []
        for template in matching_templates:
            match = self._score_template(
                template=template,
                pair=pair,
                n1n9_dist=n1n9_dist,
                interbase_angle=interbase_angle,
                planarity=planarity,
            )
            matches.append(match)

        # Sort by score (lower is better)
        matches.sort(key=lambda m: m.score)

        # Determine best match and confidence
        result = ClassificationResult(
            res1_id=pair.res1_id,
            res2_id=pair.res2_id,
            sequence=sequence,
            all_matches=matches,
        )

        if matches:
            result.best_match = matches[0]
            result.is_classified = matches[0].confidence > 0.5

        return result

    def _get_templates_for_sequence(
        self, sequence: str
    ) -> List[IdealizedTemplate]:
        """Get all templates matching a sequence.

        Args:
            sequence: Two-letter sequence (e.g., "GC").

        Returns:
            List of matching templates.
        """
        templates = []
        for (seq, _), template in self.templates.items():
            if seq == sequence:
                templates.append(template)
        return templates

    def _score_template(
        self,
        template: IdealizedTemplate,
        pair: CachedPair,
        n1n9_dist: float,
        interbase_angle: float,
        planarity: float,
    ) -> MatchResult:
        """Score how well a pair matches a template.

        Args:
            template: Template to match against.
            pair: Candidate pair.
            n1n9_dist: N1N9 distance.
            interbase_angle: Interbase angle in degrees.
            planarity: Planarity measure.

        Returns:
            MatchResult with scores.
        """
        # Initialize result
        result = MatchResult(
            sequence=template.sequence,
            lw_class=template.lw_class,
            score=float("inf"),
            geometry_score=float("inf"),
            has_idealized_coords=bool(template.res1_atoms),
        )

        # Compute geometry score from z-scores
        if template.stats:
            stats = template.stats

            # N1N9 distance z-score
            if stats.n1n9_dist_std > 0:
                result.n1n9_zscore = abs(
                    (n1n9_dist - stats.n1n9_dist_mean) / stats.n1n9_dist_std
                )
            else:
                result.n1n9_zscore = abs(n1n9_dist - stats.n1n9_dist_mean)

            # Interbase angle z-score
            if stats.interbase_angle_std > 0:
                result.angle_zscore = abs(
                    (interbase_angle - stats.interbase_angle_mean)
                    / stats.interbase_angle_std
                )
            else:
                result.angle_zscore = abs(interbase_angle - stats.interbase_angle_mean) / 10.0

            # Planarity z-score
            if stats.planarity_std > 0:
                result.planarity_zscore = abs(
                    (planarity - stats.planarity_mean) / stats.planarity_std
                )
            else:
                result.planarity_zscore = abs(planarity - stats.planarity_mean)

            # Combined geometry score
            result.geometry_score = (
                self.n1n9_weight * result.n1n9_zscore
                + self.angle_weight * result.angle_zscore
                + self.planarity_weight * result.planarity_zscore
            )

        # Compute coordinate RMSD if template has idealized coords
        if template.res1_atoms and template.res2_atoms:
            rmsd = self._compute_template_rmsd(template, pair)
            if rmsd is not None:
                result.coord_rmsd = rmsd

        # Check H-bond match (simplified for now)
        # TODO: Compare actual H-bonds with expected pattern
        if template.stats and template.stats.expected_hbonds:
            # For now, just check if number of H-bonds matches
            if template.stats.hbonds_num_mean > 0 and len(pair.hbonds) > 0:
                result.hbond_match = True

        # Compute overall score
        result.score = result.geometry_score

        if result.coord_rmsd is not None:
            # Blend geometry and coordinate scores
            result.score = (
                self.geometry_weight * result.geometry_score
                + self.coord_weight * result.coord_rmsd
            )
            if result.hbond_match:
                result.score *= (1.0 - self.hbond_weight)

        # Compute confidence
        if result.score < self.high_confidence_score:
            result.confidence = 1.0
        elif result.score < self.min_confidence_score:
            result.confidence = 1.0 - (result.score - self.high_confidence_score) / (
                self.min_confidence_score - self.high_confidence_score
            )
        else:
            result.confidence = 0.0

        return result

    def _compute_template_rmsd(
        self,
        template: IdealizedTemplate,
        pair: CachedPair,
    ) -> Optional[float]:
        """Compute RMSD between pair and idealized template.

        This aligns the pair's residue 1 frame to the template's standard
        orientation and computes RMSD of overlapping atoms.

        Args:
            template: Template with idealized coordinates.
            pair: Pair with reference frames.

        Returns:
            RMSD in Angstroms, or None if insufficient atoms.
        """
        # Get atoms common to both template and (potentially available) pair atoms
        # For now, use a simplified approach based on frame comparison
        # A full implementation would transform template coords using pair frames

        # Compute frame-based similarity
        # Template assumes residue 1 is at origin with standard orientation
        # We compare the relative position/orientation of residue 2

        if not template.res1_atoms or not template.res2_atoms:
            return None

        # Get key atoms from template
        template_n1_1 = template.res1_atoms.get("N1")
        if template_n1_1 is None:
            template_n1_1 = template.res1_atoms.get("N9")
        template_n1_2 = template.res2_atoms.get("N1")
        if template_n1_2 is None:
            template_n1_2 = template.res2_atoms.get("N9")

        if template_n1_1 is None or template_n1_2 is None:
            return None

        # Template N1N9 distance
        template_n1n9_dist = np.linalg.norm(template_n1_2 - template_n1_1)

        # Compare with pair's frame-derived distance
        pair_n1n9_dist = pair.validation.dNN

        # RMSD approximation based on distance difference
        # This is a simplification - full RMSD would require atom alignment
        dist_diff = abs(pair_n1n9_dist - template_n1n9_dist)

        # Also consider origin distance difference
        # Template origin should be near center of base pair
        template_origin = np.mean(
            [np.mean(list(template.res1_atoms.values()), axis=0),
             np.mean(list(template.res2_atoms.values()), axis=0)],
            axis=0
        )
        pair_origin_dist = pair.validation.dorg

        # Estimate RMSD from these differences
        estimated_rmsd = np.sqrt(dist_diff**2 + (pair_origin_dist / 10.0)**2)

        return estimated_rmsd

    def classify_all(
        self,
        cache: PairCache,
        atoms: Optional[Dict[str, AtomCoords]] = None,
    ) -> List[ClassificationResult]:
        """Classify all valid pairs in a cache.

        Args:
            cache: PairCache with validated pairs.
            atoms: Optional atom coordinates for N1N9 positions.

        Returns:
            List of ClassificationResult for each valid pair.
        """
        results = []
        valid_pairs = cache.get_valid_pairs()

        for pair in valid_pairs:
            # Get N1N9 positions if atoms provided
            n1n9_pos1 = None
            n1n9_pos2 = None
            if atoms:
                if pair.res1_id in atoms and atoms[pair.res1_id].n1_or_n9 is not None:
                    n1n9_pos1 = atoms[pair.res1_id].n1_or_n9
                if pair.res2_id in atoms and atoms[pair.res2_id].n1_or_n9 is not None:
                    n1n9_pos2 = atoms[pair.res2_id].n1_or_n9

            result = self.classify(pair, n1n9_pos1, n1n9_pos2)
            results.append(result)

        return results

    def print_classification_summary(
        self, results: List[ClassificationResult]
    ) -> None:
        """Print summary of classification results.

        Args:
            results: List of classification results.
        """
        total = len(results)
        classified = sum(1 for r in results if r.is_classified)

        print(f"\nClassification Summary")
        print("=" * 60)
        print(f"Total pairs: {total}")
        print(f"Classified: {classified} ({100*classified/total:.1f}%)")

        # Count by LW class
        lw_counts: Dict[str, int] = {}
        for result in results:
            if result.best_match:
                lw_class = result.best_match.lw_class
                lw_counts[lw_class] = lw_counts.get(lw_class, 0) + 1

        print(f"\nBy LW class:")
        for lw_class in sorted(lw_counts.keys()):
            count = lw_counts[lw_class]
            pct = 100 * count / total
            print(f"  {lw_class}: {count} ({pct:.1f}%)")

        # High confidence matches
        high_conf = sum(
            1 for r in results
            if r.best_match and r.best_match.confidence > 0.8
        )
        print(f"\nHigh confidence (>0.8): {high_conf} ({100*high_conf/total:.1f}%)")


def main():
    """Test pair identifier on a PDB file."""
    import argparse

    parser = argparse.ArgumentParser(description="Classify base pairs")
    parser.add_argument("pdb_id", help="PDB ID to process")
    parser.add_argument(
        "--templates",
        type=Path,
        default=Path("prototypes/pair_identification/data/templates.json"),
        help="Path to templates.json",
    )
    parser.add_argument(
        "--json-dir",
        type=Path,
        default=Path("data/json"),
        help="Directory containing JSON output",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    # Load templates
    identifier = PairIdentifier()
    identifier.load_templates(args.templates)
    print(f"Loaded {len(identifier.templates)} templates")

    # Build pair cache
    cache = PairCache(args.pdb_id, args.json_dir)
    cache.build_cache()
    print(f"Built cache with {len(cache.pairs)} pairs")

    # Classify valid pairs
    valid_pairs = cache.get_valid_pairs()
    print(f"Valid pairs: {len(valid_pairs)}")

    results = identifier.classify_all(cache, cache.atoms)

    # Print results
    if args.verbose:
        print(f"\nClassification Results for {args.pdb_id}")
        print("-" * 80)
        for result in results:
            if result.best_match:
                match = result.best_match
                print(
                    f"{result.res1_id} - {result.res2_id} ({result.sequence}): "
                    f"{match.lw_class} (score={match.score:.2f}, conf={match.confidence:.2f})"
                )
                if len(result.all_matches) > 1:
                    alt = result.all_matches[1]
                    print(
                        f"    Alt: {alt.lw_class} (score={alt.score:.2f})"
                    )

    identifier.print_classification_summary(results)


if __name__ == "__main__":
    main()

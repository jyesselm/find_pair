"""Statistical geometry for template validation."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import json
import numpy as np


@dataclass
class GeometryStats:
    """Statistical geometry for a (sequence, LW class) template.

    Attributes:
        sequence: Two-letter sequence (e.g., "GC").
        lw_class: Leontis-Westhof class (e.g., "cWW").
        count: Number of instances observed.
        n1n9_dist_mean: Mean N1/N9 distance (Angstroms).
        n1n9_dist_std: Std of N1/N9 distance.
        interbase_angle_mean: Mean interbase angle (degrees).
        interbase_angle_std: Std of interbase angle.
        planarity_mean: Mean planarity (Angstroms).
        planarity_std: Std of planarity.
        hbonds_num_mean: Mean number of H-bonds.
        expected_hbonds: Most common H-bond pattern.
    """

    sequence: str
    lw_class: str
    count: int = 0
    n1n9_dist_mean: float = 0.0
    n1n9_dist_std: float = 0.0
    interbase_angle_mean: float = 0.0
    interbase_angle_std: float = 0.0
    planarity_mean: float = 0.0
    planarity_std: float = 0.0
    hbonds_num_mean: float = 0.0
    expected_hbonds: str = ""


class TemplateStatsLoader:
    """Load and query template statistics from DSSR data.

    Provides statistical geometry constraints for validating base pairs,
    including N1-N9 distance distributions and RMSD thresholds per LW class.
    """

    def __init__(self, stats_file: Optional[Path] = None):
        """Initialize stats loader.

        Args:
            stats_file: Path to pre-computed stats JSON file. If None,
                        stats must be computed from DSSR data.
        """
        self.stats: Dict[Tuple[str, str], GeometryStats] = {}
        if stats_file and stats_file.exists():
            self.load_from_json(stats_file)

    def load_from_json(self, stats_file: Path) -> None:
        """Load pre-computed statistics from JSON file.

        Args:
            stats_file: Path to JSON file with template statistics.
        """
        with open(stats_file) as f:
            data = json.load(f)

        for template_data in data.get("templates", []):
            if "stats" not in template_data:
                continue

            stats = self._parse_template_stats(template_data)
            if stats:
                key = (stats.sequence, stats.lw_class)
                self.stats[key] = stats

    def _parse_template_stats(self, template_data: dict) -> Optional[GeometryStats]:
        """Parse stats from template data dictionary."""
        seq = template_data["sequence"]
        lw = template_data["lw_class"]
        s = template_data["stats"]

        return GeometryStats(
            sequence=seq,
            lw_class=lw,
            count=s["count"],
            n1n9_dist_mean=s["n1n9_dist_mean"],
            n1n9_dist_std=s["n1n9_dist_std"],
            interbase_angle_mean=s["interbase_angle_mean"],
            interbase_angle_std=s["interbase_angle_std"],
            planarity_mean=s["planarity_mean"],
            planarity_std=s["planarity_std"],
            hbonds_num_mean=s["hbonds_num_mean"],
            expected_hbonds=s["expected_hbonds"],
        )

    def get_stats(
        self, lw_class: str, sequence: str
    ) -> Optional[GeometryStats]:
        """Get statistics for a specific LW class and sequence.

        Args:
            lw_class: Leontis-Westhof class (e.g., "cWW").
            sequence: Two-letter sequence (e.g., "GC").

        Returns:
            GeometryStats or None if not found.
        """
        return self.stats.get((sequence, lw_class))

    def get_n1n9_range(
        self, lw_class: str, sequence: str, n_std: float = 3.0
    ) -> Optional[Tuple[float, float]]:
        """Get expected N1-N9 distance range for validation.

        Args:
            lw_class: Leontis-Westhof class.
            sequence: Two-letter sequence.
            n_std: Number of standard deviations for range.

        Returns:
            Tuple of (min_dist, max_dist) or None if stats not available.
        """
        stats = self.get_stats(lw_class, sequence)
        if not stats or stats.n1n9_dist_std == 0:
            return None

        min_dist = stats.n1n9_dist_mean - n_std * stats.n1n9_dist_std
        max_dist = stats.n1n9_dist_mean + n_std * stats.n1n9_dist_std
        return (max(0, min_dist), max_dist)

    def get_rmsd_threshold(
        self, lw_class: str, sequence: str, percentile: float = 0.95
    ) -> float:
        """Get RMSD threshold for template matching.

        Uses statistical data to estimate a reasonable RMSD cutoff.
        Falls back to class-based defaults if stats unavailable.

        Args:
            lw_class: Leontis-Westhof class.
            sequence: Two-letter sequence.
            percentile: Percentile for threshold (unused, for future use).

        Returns:
            RMSD threshold in Angstroms.
        """
        stats = self.get_stats(lw_class, sequence)
        if stats and stats.count >= 10:
            return min(1.5, 0.5 + stats.planarity_mean)

        return self._default_rmsd_threshold(lw_class)

    def _default_rmsd_threshold(self, lw_class: str) -> float:
        """Get default RMSD threshold by LW class."""
        wc_classes = {"cWW", "tWW"}
        if lw_class in wc_classes:
            return 1.0
        return 1.5

    def get_all_lw_classes(self) -> List[str]:
        """Get list of all LW classes with statistics.

        Returns:
            Sorted list of LW class names.
        """
        classes = set()
        for _, lw_class in self.stats.keys():
            classes.add(lw_class)
        return sorted(classes)

    def get_sequences_for_lw(self, lw_class: str) -> List[str]:
        """Get all sequences with statistics for a specific LW class.

        Args:
            lw_class: Leontis-Westhof class.

        Returns:
            List of sequences with stats for this LW class.
        """
        sequences = []
        for (seq, lw), _ in self.stats.items():
            if lw == lw_class:
                sequences.append(seq)
        return sorted(sequences)

    def print_summary(self) -> None:
        """Print summary of loaded statistics."""
        print(f"\nTemplate Statistics Summary")
        print("=" * 60)
        print(f"Total entries: {len(self.stats)}")
        print(f"LW classes: {', '.join(self.get_all_lw_classes())}")

        lw_counts: Dict[str, int] = {}
        for _, lw_class in self.stats.keys():
            lw_counts[lw_class] = lw_counts.get(lw_class, 0) + 1

        print(f"\nBy LW class:")
        for lw_class in sorted(lw_counts.keys()):
            count = lw_counts[lw_class]
            print(f"  {lw_class}: {count} sequences")

        print(f"\nTop 10 by instance count:")
        stats_list = sorted(
            self.stats.values(),
            key=lambda s: s.count,
            reverse=True
        )

        for stats in stats_list[:10]:
            n1n9 = f"{stats.n1n9_dist_mean:.2f}±{stats.n1n9_dist_std:.2f}"
            print(
                f"  {stats.sequence}-{stats.lw_class}: "
                f"{stats.count:,} instances, N1N9={n1n9} Å"
            )

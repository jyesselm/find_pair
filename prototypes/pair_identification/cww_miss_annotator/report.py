"""Report generation and aggregate statistics.

This module provides functionality to aggregate statistics across multiple PDB
reports and generate summary reports for cWW identification performance.
"""

from pathlib import Path
from typing import Dict, List
from dataclasses import dataclass, field
import json

from .diagnostics import PDBReport, MissAnnotation


@dataclass
class AggregateStats:
    """Aggregate statistics across multiple PDBs.

    Attributes:
        total_pdbs: Total number of PDBs analyzed.
        total_canonical_cww: Total canonical cWW pairs across all PDBs.
        total_true_positives: Total correct cWW identifications.
        total_false_negatives: Total canonical cWW pairs missed.
        total_false_positives: Total non-cWW pairs incorrectly called cWW.
        dssr_questionable: Pairs flagged as likely DSSR errors (high RMSD + no H-bonds).
        reason_counts: Distribution of reason codes across all misses.
        reason_by_sequence: Reason distribution broken down by sequence.
        wrong_atom_counts: Distribution of wrong atom errors.
        rmsd_gap_distribution: Distribution of RMSD gaps (binned).
        sequence_stats: Per-sequence statistics.
    """

    total_pdbs: int = 0
    total_canonical_cww: int = 0
    total_true_positives: int = 0
    total_false_negatives: int = 0
    total_false_positives: int = 0
    dssr_questionable: int = 0  # Likely DSSR errors, not our false negatives

    reason_counts: Dict[str, int] = field(default_factory=dict)
    reason_by_sequence: Dict[str, Dict[str, int]] = field(default_factory=dict)
    wrong_atom_counts: Dict[str, int] = field(default_factory=dict)
    rmsd_gap_distribution: Dict[str, int] = field(default_factory=dict)
    sequence_stats: Dict[str, Dict[str, int]] = field(default_factory=dict)

    @property
    def accuracy(self) -> float:
        """Calculate overall accuracy (sensitivity).

        Returns:
            Accuracy as fraction of canonical cWW correctly identified.
        """
        total = self.total_true_positives + self.total_false_negatives
        return self.total_true_positives / total if total > 0 else 0.0

    @property
    def adjusted_accuracy(self) -> float:
        """Calculate accuracy excluding DSSR questionable pairs.

        These are pairs with high RMSD (>1Å) and no H-bonds, which are
        likely DSSR errors rather than our false negatives.

        Returns:
            Adjusted accuracy excluding likely DSSR errors.
        """
        adjusted_fn = self.total_false_negatives - self.dssr_questionable
        total = self.total_true_positives + adjusted_fn
        return self.total_true_positives / total if total > 0 else 0.0

    @property
    def sensitivity(self) -> float:
        """Calculate sensitivity (same as accuracy for FN analysis).

        Returns:
            Sensitivity as fraction of canonical cWW correctly identified.
        """
        return self.accuracy

    @property
    def precision(self) -> float:
        """Calculate precision (fraction of predictions that are correct).

        Returns:
            Precision as fraction of cWW predictions that are correct.
        """
        total = self.total_true_positives + self.total_false_positives
        return self.total_true_positives / total if total > 0 else 0.0

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization.

        Returns:
            Dictionary representation suitable for JSON output.
        """
        return {
            "summary": {
                "total_pdbs": self.total_pdbs,
                "total_canonical_cww": self.total_canonical_cww,
                "true_positives": self.total_true_positives,
                "false_negatives": self.total_false_negatives,
                "false_positives": self.total_false_positives,
                "dssr_questionable": self.dssr_questionable,
                "accuracy": round(self.accuracy, 4),
                "adjusted_accuracy": round(self.adjusted_accuracy, 4),
                "sensitivity": round(self.sensitivity, 4),
                "precision": round(self.precision, 4),
            },
            "reason_distribution": dict(
                sorted(self.reason_counts.items(), key=lambda x: -x[1])
            ),
            "reason_by_sequence": {
                seq: dict(sorted(reasons.items(), key=lambda x: -x[1]))
                for seq, reasons in sorted(self.reason_by_sequence.items())
            },
            "wrong_atom_distribution": dict(
                sorted(self.wrong_atom_counts.items(), key=lambda x: -x[1])
            ),
            "rmsd_gap_distribution": dict(
                sorted(
                    self.rmsd_gap_distribution.items(),
                    key=lambda x: self._sort_rmsd_bin(x[0]),
                )
            ),
            "sequence_breakdown": {
                seq: {
                    "total": stats["total"],
                    "misses": stats["misses"],
                    "accuracy": (
                        (stats["total"] - stats["misses"]) / stats["total"]
                        if stats["total"] > 0
                        else 0.0
                    ),
                }
                for seq, stats in sorted(self.sequence_stats.items())
            },
        }

    @staticmethod
    def _sort_rmsd_bin(bin_label: str) -> float:
        """Helper to sort RMSD bin labels numerically.

        Args:
            bin_label: RMSD bin label (e.g., "0-0.1", ">1.0").

        Returns:
            Numeric value for sorting.
        """
        if bin_label.startswith(">"):
            return 999.0
        return float(bin_label.split("-")[0])


class AggregateReporter:
    """Generates aggregate statistics from multiple PDB reports.

    This class processes a collection of PDBReport objects and computes
    aggregate statistics including reason distributions, sequence breakdowns,
    and geometric diagnostics.
    """

    RMSD_BINS = [0.0, 0.1, 0.2, 0.3, 0.5, 1.0, float("inf")]
    RMSD_LABELS = ["0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.5", "0.5-1.0", ">1.0"]

    def aggregate(self, reports: List[PDBReport]) -> AggregateStats:
        """Compute aggregate statistics from list of PDB reports.

        Args:
            reports: List of PDBReport objects to aggregate.

        Returns:
            AggregateStats object with combined statistics.
        """
        stats = AggregateStats()
        stats.total_pdbs = len(reports)

        for report in reports:
            stats.total_canonical_cww += report.total_canonical_cww
            stats.total_true_positives += report.true_positives
            stats.total_false_negatives += report.false_negatives
            stats.total_false_positives += report.false_positives

            # Update sequence totals
            for ann in report.fn_annotations:
                seq = ann.sequence
                if seq not in stats.sequence_stats:
                    stats.sequence_stats[seq] = {"total": 0, "misses": 0}
                stats.sequence_stats[seq]["total"] += 1
                stats.sequence_stats[seq]["misses"] += 1

                # Count DSSR questionable pairs
                if "dssr_questionable" in ann.reasons:
                    stats.dssr_questionable += 1

            # Process each false negative annotation
            for ann in report.fn_annotations:
                self._process_annotation(ann, stats)

            # Process false positives similarly
            for ann in report.fp_annotations:
                self._process_annotation(ann, stats)

        return stats

    def _process_annotation(self, ann: MissAnnotation, stats: AggregateStats):
        """Update stats with annotation details.

        Args:
            ann: MissAnnotation to process.
            stats: AggregateStats to update in place.
        """
        # Count reasons
        for reason in ann.reasons:
            stats.reason_counts[reason] = stats.reason_counts.get(reason, 0) + 1

            # By sequence
            if ann.sequence not in stats.reason_by_sequence:
                stats.reason_by_sequence[ann.sequence] = {}
            seq_reasons = stats.reason_by_sequence[ann.sequence]
            seq_reasons[reason] = seq_reasons.get(reason, 0) + 1

        # Wrong atoms
        for donor, desc in ann.hbond_diagnostics.wrong_atoms.items():
            key = f"{donor}: {desc}"
            stats.wrong_atom_counts[key] = stats.wrong_atom_counts.get(key, 0) + 1

        # RMSD gap binning
        gap = ann.geometric_diagnostics.rmsd_gap
        bin_label = self._bin_rmsd_gap(gap)
        stats.rmsd_gap_distribution[bin_label] = (
            stats.rmsd_gap_distribution.get(bin_label, 0) + 1
        )

    def _bin_rmsd_gap(self, gap: float) -> str:
        """Bin RMSD gap into discrete ranges.

        Args:
            gap: RMSD gap value to bin.

        Returns:
            Bin label string (e.g., "0-0.1", ">1.0").
        """
        for i, threshold in enumerate(self.RMSD_BINS[1:]):
            if gap < threshold:
                return self.RMSD_LABELS[i]
        return self.RMSD_LABELS[-1]


def save_pdb_report(report: PDBReport, output_path: Path):
    """Save single PDB report to JSON.

    Args:
        report: PDBReport to save.
        output_path: Path to output JSON file.
    """
    with open(output_path, "w") as f:
        json.dump(report.to_dict(), f, indent=2)


def save_aggregate_report(stats: AggregateStats, output_path: Path):
    """Save aggregate statistics to JSON.

    Args:
        stats: AggregateStats to save.
        output_path: Path to output JSON file.
    """
    with open(output_path, "w") as f:
        json.dump(stats.to_dict(), f, indent=2)


def load_pdb_report(json_path: Path) -> PDBReport:
    """Load a single PDB report from JSON.

    Args:
        json_path: Path to JSON file.

    Returns:
        Reconstructed PDBReport object.
    """
    with open(json_path, "r") as f:
        data = json.load(f)

    from .diagnostics import (
        HBondDiagnostics,
        GeometricDiagnostics,
        ExpectedHBond,
    )

    # Reconstruct false negative annotations
    fn_annotations = []
    for ann_dict in data.get("false_negatives", []):
        hbond_diag = HBondDiagnostics(
            expected_hbonds=[
                ExpectedHBond(**hb)
                for hb in ann_dict["hbond_diagnostics"]["expected_hbonds"]
            ],
            found_hbonds=ann_dict["hbond_diagnostics"]["found_hbonds"],
            missing_hbonds=[
                ExpectedHBond(**hb)
                for hb in ann_dict["hbond_diagnostics"]["missing_hbonds"]
            ],
            extra_hbonds=ann_dict["hbond_diagnostics"]["extra_hbonds"],
            wrong_atoms=ann_dict["hbond_diagnostics"]["wrong_atoms"],
            distance_issues=[
                (item["atom"], item["distance"], item["issue"])
                for item in ann_dict["hbond_diagnostics"]["distance_issues"]
            ],
            overloaded_acceptors=ann_dict["hbond_diagnostics"]["overloaded_acceptors"],
        )

        geom_diag = GeometricDiagnostics(**ann_dict["geometric_diagnostics"])

        ann = MissAnnotation(
            res_id1=ann_dict["res_id1"],
            res_id2=ann_dict["res_id2"],
            sequence=ann_dict["sequence"],
            our_prediction=ann_dict["our_prediction"],
            dssr_class=ann_dict["dssr_class"],
            saenger=ann_dict["saenger"],
            reasons=ann_dict["reasons"],
            hbond_diagnostics=hbond_diag,
            geometric_diagnostics=geom_diag,
        )
        fn_annotations.append(ann)

    # Reconstruct false positive annotations (same structure)
    fp_annotations = []
    for ann_dict in data.get("false_positives", []):
        hbond_diag = HBondDiagnostics(
            expected_hbonds=[
                ExpectedHBond(**hb)
                for hb in ann_dict["hbond_diagnostics"]["expected_hbonds"]
            ],
            found_hbonds=ann_dict["hbond_diagnostics"]["found_hbonds"],
            missing_hbonds=[
                ExpectedHBond(**hb)
                for hb in ann_dict["hbond_diagnostics"]["missing_hbonds"]
            ],
            extra_hbonds=ann_dict["hbond_diagnostics"]["extra_hbonds"],
            wrong_atoms=ann_dict["hbond_diagnostics"]["wrong_atoms"],
            distance_issues=[
                (item["atom"], item["distance"], item["issue"])
                for item in ann_dict["hbond_diagnostics"]["distance_issues"]
            ],
            overloaded_acceptors=ann_dict["hbond_diagnostics"]["overloaded_acceptors"],
        )

        geom_diag = GeometricDiagnostics(**ann_dict["geometric_diagnostics"])

        ann = MissAnnotation(
            res_id1=ann_dict["res_id1"],
            res_id2=ann_dict["res_id2"],
            sequence=ann_dict["sequence"],
            our_prediction=ann_dict["our_prediction"],
            dssr_class=ann_dict["dssr_class"],
            saenger=ann_dict["saenger"],
            reasons=ann_dict["reasons"],
            hbond_diagnostics=hbond_diag,
            geometric_diagnostics=geom_diag,
        )
        fp_annotations.append(ann)

    summary = data["summary"]
    return PDBReport(
        pdb_id=data["pdb_id"],
        total_canonical_cww=summary["total_canonical_cww"],
        true_positives=summary["true_positives"],
        false_negatives=summary["false_negatives"],
        false_positives=summary["false_positives"],
        fn_annotations=fn_annotations,
        fp_annotations=fp_annotations,
    )


def load_pdb_reports(report_dir: Path) -> List[PDBReport]:
    """Load all PDB reports from directory.

    Args:
        report_dir: Directory containing PDB report JSON files.

    Returns:
        List of PDBReport objects.
    """
    reports = []
    for json_path in sorted(report_dir.glob("*.json")):
        if json_path.name == "aggregate.json":
            continue
        try:
            report = load_pdb_report(json_path)
            reports.append(report)
        except Exception as e:
            print(f"Warning: Failed to load {json_path}: {e}")
    return reports

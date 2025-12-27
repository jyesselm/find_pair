"""Diagnostic dataclasses for pair analysis."""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Tuple


@dataclass
class ExpectedHBond:
    """Expected hydrogen bond for canonical base pairing.

    Attributes:
        donor_atom: Name of donor atom (e.g., "N1", "N6").
        acceptor_atom: Name of acceptor atom (e.g., "O6", "N7").
        source: Source of expectation ("canonical" or "dssr").
    """

    donor_atom: str
    acceptor_atom: str
    source: str


@dataclass
class HBondDiagnostics:
    """Diagnostics for hydrogen bonding in a base pair.

    Attributes:
        expected_hbonds: Expected H-bonds for this sequence.
        found_hbonds: Detected H-bonds with donor/acceptor atoms, distance, context.
        missing_hbonds: Expected H-bonds not found.
        extra_hbonds: Unexpected H-bonds found.
        wrong_atoms: Dict mapping issue description to details.
        distance_issues: List of (atom, distance, issue_description) tuples.
        overloaded_acceptors: Acceptor atoms receiving >2 H-bonds.
    """

    expected_hbonds: List[ExpectedHBond] = field(default_factory=list)
    found_hbonds: List[Dict[str, Any]] = field(default_factory=list)
    missing_hbonds: List[ExpectedHBond] = field(default_factory=list)
    extra_hbonds: List[Dict[str, Any]] = field(default_factory=list)
    wrong_atoms: Dict[str, str] = field(default_factory=dict)
    distance_issues: List[Tuple[str, float, str]] = field(default_factory=list)
    overloaded_acceptors: List[str] = field(default_factory=list)


@dataclass
class GeometricDiagnostics:
    """Geometric diagnostics for base pair quality.

    Attributes:
        rmsd_cww: RMSD to cWW template.
        rmsd_best: RMSD to best-fitting template.
        best_lw: Leontis-Westhof class of best-fitting template.
        rmsd_gap: Difference (rmsd_cww - rmsd_best).
        interbase_angle: Angle between base planes in degrees.
        n1n9_distance: Distance between glycosidic N atoms (N1/N9).
        is_geometric_outlier: True if angle >15 degrees or distance abnormal.
    """

    rmsd_cww: float
    rmsd_best: float
    best_lw: str
    rmsd_gap: float
    interbase_angle: float
    n1n9_distance: float
    is_geometric_outlier: bool


@dataclass
class MissAnnotation:
    """Complete annotation for a false negative or false positive.

    Attributes:
        res_id1: Residue ID of first base (format: "chain-name-seq[ins]").
        res_id2: Residue ID of second base.
        sequence: Sequence pair (e.g., "GC").
        our_prediction: Our LW class prediction.
        dssr_class: DSSR's LW class.
        saenger: Saenger classification or "--".
        reasons: List of reason codes explaining the miss.
        hbond_diagnostics: H-bond analysis.
        geometric_diagnostics: Geometric quality metrics.
    """

    res_id1: str
    res_id2: str
    sequence: str
    our_prediction: str
    dssr_class: str
    saenger: str
    reasons: List[str] = field(default_factory=list)
    hbond_diagnostics: HBondDiagnostics = field(default_factory=HBondDiagnostics)
    geometric_diagnostics: GeometricDiagnostics = field(
        default_factory=lambda: GeometricDiagnostics(
            rmsd_cww=0.0,
            rmsd_best=0.0,
            best_lw="",
            rmsd_gap=0.0,
            interbase_angle=0.0,
            n1n9_distance=0.0,
            is_geometric_outlier=False,
        )
    )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization.

        Returns:
            Dictionary representation suitable for JSON output.
        """
        return {
            "res_id1": self.res_id1,
            "res_id2": self.res_id2,
            "sequence": self.sequence,
            "our_prediction": self.our_prediction,
            "dssr_class": self.dssr_class,
            "saenger": self.saenger,
            "reasons": self.reasons,
            "hbond_diagnostics": self._serialize_hbond_diagnostics(),
            "geometric_diagnostics": self._serialize_geometric_diagnostics(),
        }

    def _serialize_hbond_diagnostics(self) -> Dict[str, Any]:
        """Serialize H-bond diagnostics to dict."""
        return {
            "expected_hbonds": [
                {
                    "donor_atom": hb.donor_atom,
                    "acceptor_atom": hb.acceptor_atom,
                    "source": hb.source,
                }
                for hb in self.hbond_diagnostics.expected_hbonds
            ],
            "found_hbonds": self.hbond_diagnostics.found_hbonds,
            "missing_hbonds": [
                {
                    "donor_atom": hb.donor_atom,
                    "acceptor_atom": hb.acceptor_atom,
                    "source": hb.source,
                }
                for hb in self.hbond_diagnostics.missing_hbonds
            ],
            "extra_hbonds": self.hbond_diagnostics.extra_hbonds,
            "wrong_atoms": self.hbond_diagnostics.wrong_atoms,
            "distance_issues": [
                {"atom": atom, "distance": dist, "issue": issue}
                for atom, dist, issue in self.hbond_diagnostics.distance_issues
            ],
            "overloaded_acceptors": self.hbond_diagnostics.overloaded_acceptors,
        }

    def _serialize_geometric_diagnostics(self) -> Dict[str, Any]:
        """Serialize geometric diagnostics to dict."""
        return {
            "rmsd_cww": self.geometric_diagnostics.rmsd_cww,
            "rmsd_best": self.geometric_diagnostics.rmsd_best,
            "best_lw": self.geometric_diagnostics.best_lw,
            "rmsd_gap": self.geometric_diagnostics.rmsd_gap,
            "interbase_angle": self.geometric_diagnostics.interbase_angle,
            "n1n9_distance": self.geometric_diagnostics.n1n9_distance,
            "is_geometric_outlier": self.geometric_diagnostics.is_geometric_outlier,
        }


@dataclass
class PDBReport:
    """Summary report for cWW identification in a single PDB.

    Attributes:
        pdb_id: PDB identifier.
        total_canonical_cww: Total canonical cWW pairs in DSSR reference.
        true_positives: Correct cWW identifications.
        false_negatives: Canonical cWW pairs missed.
        false_positives: Non-cWW pairs incorrectly called cWW.
        fn_annotations: Detailed annotations for false negatives.
        fp_annotations: Detailed annotations for false positives.
    """

    pdb_id: str
    total_canonical_cww: int
    true_positives: int
    false_negatives: int
    false_positives: int
    fn_annotations: List[MissAnnotation] = field(default_factory=list)
    fp_annotations: List[MissAnnotation] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization.

        Returns:
            Dictionary representation suitable for JSON output.
        """
        total = self.total_canonical_cww
        tp = self.true_positives
        fp = self.false_positives

        return {
            "pdb_id": self.pdb_id,
            "summary": {
                "total_canonical_cww": total,
                "true_positives": tp,
                "false_negatives": self.false_negatives,
                "false_positives": fp,
                "sensitivity": tp / total if total > 0 else 0.0,
                "precision": tp / (tp + fp) if (tp + fp) > 0 else 0.0,
            },
            "false_negatives": [ann.to_dict() for ann in self.fn_annotations],
            "false_positives": [ann.to_dict() for ann in self.fp_annotations],
        }

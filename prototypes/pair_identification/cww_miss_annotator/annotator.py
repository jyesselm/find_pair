"""Main annotator that orchestrates cWW miss analysis.

This module provides the MissAnnotator class that combines H-bond and geometric
analysis to produce detailed annotations explaining why Watson-Crick base pairs
are misclassified (false negatives or false positives).
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from cww_miss_annotator.diagnostics import (
    MissAnnotation,
    HBondDiagnostics,
    GeometricDiagnostics,
    PDBReport,
)
from cww_miss_annotator.hbond_analyzer import HBondAnalyzer
from cww_miss_annotator.geometric_analyzer import GeometricAnalyzer
from cww_miss_annotator.loaders import (
    DSSRPair,
    SlotHBond,
    load_dssr_pairs,
    load_slot_hbonds,
)
from core.pdb_parser import parse_pdb


def get_edge_pair(lw_class: str) -> str:
    """Extract edge pair from LW class (ignoring cis/trans).

    cWW and tWW both use Watson-Watson edges, so they should be treated
    as equivalent when comparing base atoms only (sugars distinguish them).

    Args:
        lw_class: LW classification (e.g., "cWW", "tWW", "cWS")

    Returns:
        Edge pair without cis/trans (e.g., "WW", "WS", "HS")
    """
    if len(lw_class) >= 3:
        return lw_class[1:3]  # Extract edge letters (WW, WS, HS, etc.)
    return lw_class


class MissAnnotator:
    """Annotates misclassified cWW pairs with detailed diagnostics.

    This orchestrator combines H-bond analysis and geometric analysis to
    determine why a Watson-Crick base pair was misclassified. It compares
    our predictions against DSSR reference data.

    Attributes:
        pdb_dir: Directory containing PDB files.
        hbond_dir: Directory containing slot H-bond JSON files.
        dssr_dir: Directory containing DSSR reference JSON files.
        hbond_analyzer: Analyzer for H-bond pattern diagnostics.
        geometric_analyzer: Analyzer for geometric quality metrics.
    """

    REASON_CODES = [
        "no_hbonds",  # No H-bonds detected
        "missing_hbonds",  # Some expected H-bonds missing
        "wrong_atoms",  # H-bonds to wrong acceptor atoms
        "extra_hbonds",  # Unexpected H-bonds detected
        "distance_issues",  # H-bond distances out of range
        "overloaded_acceptor",  # Acceptor has >2 H-bonds
        "rmsd_prefers_other",  # Better RMSD to non-cWW template
        "geometric_outlier",  # Angle or N1N9 distance abnormal
        "non_canonical",  # DSSR Saenger is "--"
        "dssr_questionable",  # DSSR likely wrong (high RMSD + no H-bonds)
    ]

    def __init__(
        self,
        pdb_dir: Path = Path("data/pdb"),
        hbond_dir: Path = Path("data/json/slot_hbonds"),
        dssr_dir: Path = Path("data/json_dssr"),
        idealized_dir: Path = Path("basepair-idealized"),
        exemplar_dir: Path = Path("basepair-exemplars"),
    ):
        """Initialize the miss annotator.

        Args:
            pdb_dir: Directory containing PDB files. Default: data/pdb
            hbond_dir: Directory with slot H-bond JSON. Default: data/json/slot_hbonds
            dssr_dir: Directory with DSSR reference JSON. Default: data/json_dssr
            idealized_dir: Directory with idealized templates. Default: basepair-idealized
            exemplar_dir: Directory with exemplar templates. Default: basepair-exemplars
        """
        self.pdb_dir = Path(pdb_dir)
        self.hbond_dir = Path(hbond_dir)
        self.dssr_dir = Path(dssr_dir)

        self.hbond_analyzer = HBondAnalyzer()
        self.geometric_analyzer = GeometricAnalyzer(
            idealized_dir=Path(idealized_dir), exemplar_dir=Path(exemplar_dir)
        )

    def annotate_pdb(self, pdb_id: str) -> Optional[PDBReport]:
        """Annotate all cWW classification differences for a PDB.

        Args:
            pdb_id: PDB identifier (e.g., "1EHZ")

        Returns:
            PDBReport with summary counts (TP, FN, FP) and detailed annotations
            for false negatives (we miss cWW) and false positives (we incorrectly
            predict cWW). Returns None if DSSR data not available.
        """
        # Load data files
        dssr_path = self.dssr_dir / f"{pdb_id}.json"
        hbond_path = self.hbond_dir / f"{pdb_id}.json"
        pdb_path = self.pdb_dir / f"{pdb_id}.pdb"

        if not dssr_path.exists():
            return None

        # Load DSSR canonical cWW pairs
        dssr_pairs = load_dssr_pairs(dssr_path, lw_filter="cWW")

        # Load slot H-bonds
        slot_hbonds = {}
        if hbond_path.exists():
            slot_hbonds = load_slot_hbonds(hbond_path)

        # Load PDB residues for geometric analysis
        residues = {}
        if pdb_path.exists():
            try:
                residues = parse_pdb(pdb_path)
            except Exception:
                pass

        # Classify each pair and compare to DSSR
        fn_annotations = []  # DSSR says cWW, we say something else
        fp_annotations = []  # We say cWW, DSSR says something else
        true_positives = 0

        for (res_id1, res_id2), dssr_pair in dssr_pairs.items():
            # Get our H-bonds for this pair
            hbonds = slot_hbonds.get((res_id1, res_id2), [])

            # Analyze H-bonds
            sequence = dssr_pair.bp.replace("-", "")
            h_diag = self.hbond_analyzer.analyze(
                sequence=sequence,
                found_hbonds=hbonds,
                dssr_hbonds_desc=dssr_pair.hbonds_desc,
            )

            # Analyze geometry
            res1_atoms = self._extract_atom_coords(residues, res_id1)
            res2_atoms = self._extract_atom_coords(residues, res_id2)

            g_diag = self.geometric_analyzer.analyze(
                res1_atoms=res1_atoms,
                res2_atoms=res2_atoms,
                sequence=sequence,
                interbase_angle=dssr_pair.interbase_angle,
                n1n9_distance=dssr_pair.n1n9_distance,
            )

            # Determine if this is a miss
            reasons = self._categorize_miss(h_diag, g_diag, dssr_pair)

            if reasons:
                # This is a false negative (we would miss this cWW pair)
                our_prediction = g_diag.best_lw if g_diag.best_lw else "unknown"

                annotation = MissAnnotation(
                    res_id1=res_id1,
                    res_id2=res_id2,
                    sequence=sequence,
                    our_prediction=our_prediction,
                    dssr_class="cWW",
                    saenger=dssr_pair.saenger,
                    reasons=reasons,
                    hbond_diagnostics=h_diag,
                    geometric_diagnostics=g_diag,
                )
                fn_annotations.append(annotation)
            else:
                true_positives += 1

        # TODO: Check for false positives
        # (pairs where we predict cWW but DSSR doesn't)
        # This requires running our classifier on all pairs

        return PDBReport(
            pdb_id=pdb_id,
            total_canonical_cww=len(dssr_pairs),
            true_positives=true_positives,
            false_negatives=len(fn_annotations),
            false_positives=len(fp_annotations),
            fn_annotations=fn_annotations,
            fp_annotations=fp_annotations,
        )

    def _extract_atom_coords(
        self, residues: Dict, res_id: str
    ) -> Dict[str, "np.ndarray"]:
        """Extract atom coordinates from residue dictionary.

        Args:
            residues: Dict mapping res_id to Residue objects
            res_id: Residue ID to extract

        Returns:
            Dict mapping atom name to 3D coordinates (numpy arrays).
            Returns empty dict if residue not found.
        """
        if res_id not in residues:
            return {}

        return {name: atom.coords for name, atom in residues[res_id].atoms.items()}

    def _categorize_miss(
        self, h_diag: HBondDiagnostics, g_diag: GeometricDiagnostics, dssr_pair: DSSRPair
    ) -> List[str]:
        """Determine why a pair would be misclassified.

        Args:
            h_diag: H-bond diagnostics
            g_diag: Geometric diagnostics
            dssr_pair: DSSR reference data

        Returns:
            List of reason codes explaining the miss. Empty list if no miss
            (i.e., pair would be correctly classified).
        """
        reasons = []

        # H-bond based reasons
        if not h_diag.found_hbonds:
            reasons.append("no_hbonds")
        else:
            if h_diag.missing_hbonds:
                reasons.append("missing_hbonds")
            if h_diag.wrong_atoms:
                reasons.append("wrong_atoms")
            if h_diag.extra_hbonds:
                reasons.append("extra_hbonds")
            if h_diag.distance_issues:
                reasons.append("distance_issues")
            if h_diag.overloaded_acceptors:
                reasons.append("overloaded_acceptor")

        # Geometry based reasons
        # Only flag rmsd_prefers_other if the best template has DIFFERENT edges
        # cWW and tWW have same edges (WW), so treat them as equivalent
        if g_diag.rmsd_gap > 0.5:
            best_edges = get_edge_pair(g_diag.best_lw)
            target_edges = "WW"  # We're looking for Watson-Watson pairs
            if best_edges != target_edges:
                reasons.append("rmsd_prefers_other")
        # Only flag geometric_outlier if RMSD is also poor (> 0.5Å)
        # If RMSD to cWW is good, the pair is clearly cWW regardless of minor variations
        if g_diag.is_geometric_outlier and g_diag.rmsd_cww > 0.5:
            reasons.append("geometric_outlier")

        # DSSR classification based
        if dssr_pair.saenger == "--":
            reasons.append("non_canonical")

        # Flag likely DSSR errors: high RMSD + no H-bonds = probably not a real cWW
        # These shouldn't count against our accuracy
        if g_diag.rmsd_cww > 1.0 and not h_diag.found_hbonds:
            reasons.append("dssr_questionable")

        return reasons

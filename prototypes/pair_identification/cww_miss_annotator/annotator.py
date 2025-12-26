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
from cww_miss_annotator.bp_score import compute_bp_score
from core.pdb_parser import parse_pdb


def normalize_res_id_for_hbond_lookup(res_id: str) -> str:
    """Normalize DNA res_ids to match slot H-bond file format.

    Slot H-bond files use single-letter base codes (G, C, A, T, U) while
    DSSR uses DNA-specific codes (DG, DC, DA, DT). This converts DSSR
    format to the format used in slot H-bond files.

    Args:
        res_id: Residue ID in format "chain-base-seqnum" (e.g., "E-DG-2")

    Returns:
        Normalized res_id (e.g., "E-G-2")

    Examples:
        >>> normalize_res_id_for_hbond_lookup("E-DG-2")
        'E-G-2'
        >>> normalize_res_id_for_hbond_lookup("A-G-10")
        'A-G-10'
    """
    parts = res_id.split("-")
    if len(parts) >= 2:
        base = parts[1]
        # Convert DNA bases (DG, DC, DA, DT) to single-letter codes
        if base.upper() in ("DG", "DC", "DA", "DT"):
            parts[1] = base[1].upper()  # Take second letter (G, C, A, T)
            return "-".join(parts)
    return res_id


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
        "no_hbonds",  # No H-bonds detected between bases
        "poor_planarity",  # Bases not coplanar (interbase angle > 20°)
        "missing_hbonds",  # Some expected H-bonds missing
        "wrong_hbonds",  # H-bonds to wrong acceptor atoms
        "extra_hbonds",  # Unexpected H-bonds detected
        "long_hbonds",  # H-bonds detected but distances > 4.0Å
        "short_hbonds",  # H-bonds detected but distances < 2.0Å (clash)
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
            # Normalize IDs for fallback lookup (DNA prefixes -> RNA codes)
            norm_id1 = normalize_res_id_for_hbond_lookup(res_id1)
            norm_id2 = normalize_res_id_for_hbond_lookup(res_id2)

            # Get our H-bonds for this pair
            # Try original IDs first, then normalized IDs (some files use DNA prefixes, some don't)
            hbonds = slot_hbonds.get((res_id1, res_id2), [])
            if not hbonds:
                hbonds = slot_hbonds.get((norm_id1, norm_id2), [])

            # Analyze H-bonds
            sequence = dssr_pair.bp.replace("-", "")
            h_diag = self.hbond_analyzer.analyze(
                sequence=sequence,
                found_hbonds=hbonds,
                dssr_hbonds_desc=dssr_pair.hbonds_desc,
            )

            # Analyze geometry (try original IDs, then normalized)
            res1_atoms = self._extract_atom_coords(residues, res_id1)
            if not res1_atoms:
                res1_atoms = self._extract_atom_coords(residues, norm_id1)
            res2_atoms = self._extract_atom_coords(residues, res_id2)
            if not res2_atoms:
                res2_atoms = self._extract_atom_coords(residues, norm_id2)

            g_diag = self.geometric_analyzer.analyze(
                res1_atoms=res1_atoms,
                res2_atoms=res2_atoms,
                sequence=sequence,
                interbase_angle=dssr_pair.interbase_angle,
                n1n9_distance=dssr_pair.n1n9_distance,
            )

            # Determine if this is a miss
            # Pass atom coords for extended H-bond search when geometry is good
            reasons = self._categorize_miss(h_diag, g_diag, dssr_pair, res1_atoms, res2_atoms)

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
        self, h_diag: HBondDiagnostics, g_diag: GeometricDiagnostics, dssr_pair: DSSRPair,
        res1_atoms: Dict = None, res2_atoms: Dict = None,
    ) -> List[str]:
        """Determine why a pair would be misclassified.

        Args:
            h_diag: H-bond diagnostics
            g_diag: Geometric diagnostics
            dssr_pair: DSSR reference data
            res1_atoms: Optional atom coordinates for residue 1 (for extended H-bond search)
            res2_atoms: Optional atom coordinates for residue 2 (for extended H-bond search)

        Returns:
            List of reason codes explaining the miss. Empty list if no miss
            (i.e., pair would be correctly classified).
        """
        sequence = dssr_pair.bp.replace("-", "")

        # Compute BP score first - if score >= threshold, pair is good
        # found_hbonds is already a list of dicts with context, distance, alignment keys
        # Pass atom coords for extended H-bond search when geometry is good but H-bonds sparse
        bp_score, bp_components = compute_bp_score(
            sequence, g_diag.rmsd_cww, h_diag.found_hbonds,
            res1_atoms=res1_atoms, res2_atoms=res2_atoms,
            interbase_angle=g_diag.interbase_angle,
        )

        # Use lower threshold (0.60) when extended search found H-bonds
        # These pairs have good geometry but stretched H-bonds - likely real cWW
        threshold = 0.60 if bp_components.get('extended_search', False) else 0.70
        if bp_score >= threshold:
            return []

        reasons = []

        # H-bond based reasons
        if not h_diag.found_hbonds:
            reasons.append("no_hbonds")
        else:
            if h_diag.missing_hbonds:
                reasons.append("missing_hbonds")
            if h_diag.wrong_atoms:
                reasons.append("wrong_hbonds")
            if h_diag.extra_hbonds:
                reasons.append("extra_hbonds")
            # Separate long vs short distance issues
            if h_diag.distance_issues:
                has_long = any(issue[2] == "too_long" for issue in h_diag.distance_issues)
                has_short = any(issue[2] == "too_short" for issue in h_diag.distance_issues)
                if has_long:
                    reasons.append("long_hbonds")
                if has_short:
                    reasons.append("short_hbonds")
            if h_diag.overloaded_acceptors:
                reasons.append("overloaded_acceptor")

        # Geometry based reasons
        if g_diag.rmsd_gap > 0.5:
            best_edges = get_edge_pair(g_diag.best_lw)
            target_edges = "WW"
            if best_edges != target_edges:
                reasons.append("rmsd_prefers_other")

        if g_diag.is_geometric_outlier:
            reasons.append("geometric_outlier")

        # Check planarity - bases should be coplanar (angle < 20°)
        if g_diag.interbase_angle >= 20.0:
            reasons.append("poor_planarity")

        # DSSR classification based
        if dssr_pair.saenger == "--":
            reasons.append("non_canonical")

        # Flag likely DSSR errors: high RMSD + no H-bonds = probably not a real cWW
        if g_diag.rmsd_cww > 1.0 and not h_diag.found_hbonds:
            reasons.append("dssr_questionable")

        return reasons

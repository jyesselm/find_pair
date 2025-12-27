"""Main annotator that orchestrates cWW miss analysis.

This module combines H-bond analysis, geometric analysis, and scoring
to produce detailed annotations explaining why Watson-Crick base pairs
are misclassified (false negatives or false positives).
"""

import re
import numpy as np
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from prototypes.pair_identification.analysis.diagnostics import (
    ExpectedHBond,
    HBondDiagnostics,
    GeometricDiagnostics,
    MissAnnotation,
    PDBReport,
)
from prototypes.pair_identification.core.pdb_parser import parse_pdb
from prototypes.pair_identification.core.residue import Residue, Atom
from prototypes.pair_identification.template_aligner import TemplateAligner

# Expected H-bond patterns for canonical Watson-Crick pairs
CWW_EXPECTED_PATTERNS: Dict[str, List[Tuple[str, str]]] = {
    "GC": [("N1", "N3"), ("N2", "O2"), ("N4", "O6")],
    "CG": [("N4", "O6"), ("N1", "N3"), ("N2", "O2")],
    "AU": [("N6", "O4"), ("N3", "N1")],
    "UA": [("N6", "O4"), ("N3", "N1")],
}

# Normal N1-N9 distance ranges for cWW pairs (Angstroms)
CWW_N1N9_RANGE = {
    "GC": (7.5, 10.0),
    "CG": (7.5, 10.0),
    "AU": (7.5, 10.0),
    "UA": (7.5, 10.0),
    "GU": (7.0, 10.0),
    "UG": (7.0, 10.0),
}

# Threshold for interbase angle outliers (degrees)
CWW_ANGLE_THRESHOLD = 30.0

# Glycosidic nitrogen mapping
GLYCOSIDIC_NITROGEN = {
    "A": "N9",
    "G": "N9",
    "C": "N1",
    "U": "N1",
    "T": "N1",
}


@dataclass
class DSSRPair:
    """DSSR base pair record."""

    nt1: str
    nt2: str
    lw_class: str
    saenger: str
    hbonds_desc: str
    interbase_angle: float
    n1n9_distance: float
    bp: str


@dataclass
class SlotHBond:
    """Slot-based hydrogen bond record."""

    donor_res_id: str
    donor_atom: str
    acceptor_res_id: str
    acceptor_atom: str
    distance: float
    context: str
    h_slot: Optional[int]
    lp_slot: Optional[int]
    alignment: Optional[float]


def load_dssr_pairs(
    dssr_path: Path, lw_filter: str = "cWW"
) -> Dict[Tuple[str, str], DSSRPair]:
    """Load DSSR pairs filtered by Leontis-Westhof class.

    Args:
        dssr_path: Path to DSSR JSON file.
        lw_filter: LW class to filter (default "cWW").

    Returns:
        Dict keyed by (res_id1, res_id2) with DSSRPair values.
    """
    import json
    from prototypes.pair_identification.core.identifiers import dssr_to_res_id

    with open(dssr_path) as f:
        data = json.load(f)

    pairs = {}
    standard_seqs = {"GC", "CG", "AU", "UA"}

    for pair in data.get("pairs", []):
        if pair.get("LW", "") != lw_filter:
            continue

        bp = pair.get("bp", "")
        bp_clean = bp.replace("-", "").replace("+", "")
        if bp_clean not in standard_seqs:
            continue

        nt1_slot = dssr_to_res_id(pair.get("nt1", ""))
        nt2_slot = dssr_to_res_id(pair.get("nt2", ""))
        if not nt1_slot or not nt2_slot:
            continue

        dssr_pair = DSSRPair(
            nt1=nt1_slot,
            nt2=nt2_slot,
            lw_class=lw_filter,
            saenger=pair.get("Saenger", ""),
            hbonds_desc=pair.get("hbonds_desc", ""),
            interbase_angle=pair.get("interBase_angle", 0.0),
            n1n9_distance=pair.get("N1N9_dist", 0.0),
            bp=bp,
        )
        pairs[(nt1_slot, nt2_slot)] = dssr_pair

    return pairs


def load_slot_hbonds(hbond_path: Path) -> Dict[Tuple[str, str], List[SlotHBond]]:
    """Load slot-based H-bond data.

    Args:
        hbond_path: Path to slot hbonds JSON file.

    Returns:
        Dict keyed by (res_id_i, res_id_j) with list of SlotHBond values.
    """
    import json

    with open(hbond_path) as f:
        data = json.load(f)

    hbonds = {}
    for record in data:
        res_id_i = record.get("res_id_i", "")
        res_id_j = record.get("res_id_j", "")

        slot_hbonds = []
        for hbond in record.get("hbonds", []):
            slot_hbond = SlotHBond(
                donor_res_id=hbond.get("donor_res_id", res_id_i),
                donor_atom=hbond.get("donor_atom", ""),
                acceptor_res_id=hbond.get("acceptor_res_id", res_id_j),
                acceptor_atom=hbond.get("acceptor_atom", ""),
                distance=hbond.get("distance", 0.0),
                context=hbond.get("context", ""),
                h_slot=hbond.get("h_slot"),
                lp_slot=hbond.get("lp_slot"),
                alignment=hbond.get("alignment"),
            )
            slot_hbonds.append(slot_hbond)

        hbonds[(res_id_i, res_id_j)] = slot_hbonds
        hbonds[(res_id_j, res_id_i)] = slot_hbonds

    return hbonds


def parse_dssr_hbonds(hbonds_desc: str) -> List[Tuple[str, str, float]]:
    """Parse DSSR hbonds_desc string.

    Args:
        hbonds_desc: DSSR format like "O6(carbonyl)-N4(amino)[2.83]".

    Returns:
        List of (donor_atom, acceptor_atom, distance) tuples.
    """
    if not hbonds_desc:
        return []

    hbonds = []
    pattern = r"([A-Z]\d+(?:\*)?)\([^)]+\)-([A-Z]\d+(?:\*)?)\([^)]+\)\[(\d+\.\d+)\]"

    for match in re.finditer(pattern, hbonds_desc):
        atom1 = match.group(1).rstrip("*")
        atom2 = match.group(2).rstrip("*")
        distance = float(match.group(3))
        hbonds.append((atom1, atom2, distance))

    simple_pattern = r"([A-Z]\d+(?:\*)?)-([A-Z]\d+(?:\*)?)\[(\d+\.\d+)\]"
    for match in re.finditer(simple_pattern, hbonds_desc):
        atom1 = match.group(1).rstrip("*")
        atom2 = match.group(2).rstrip("*")
        distance = float(match.group(3))
        if (atom1, atom2, distance) not in hbonds:
            hbonds.append((atom1, atom2, distance))

    return hbonds


class MissAnnotator:
    """Annotates misclassified cWW pairs with detailed diagnostics.

    Combines H-bond analysis, geometric analysis, and scoring to determine
    why a Watson-Crick base pair was misclassified.

    Attributes:
        pdb_dir: Directory containing PDB files.
        hbond_dir: Directory containing slot H-bond JSON files.
        dssr_dir: Directory containing DSSR reference JSON files.
        idealized_dir: Directory with idealized templates.
        exemplar_dir: Directory with exemplar templates.
        template_aligner: TemplateAligner for geometric analysis.
    """

    REASON_CODES = [
        "no_hbonds",
        "poor_planarity",
        "missing_hbonds",
        "wrong_hbonds",
        "extra_hbonds",
        "long_hbonds",
        "short_hbonds",
        "overloaded_acceptor",
        "rmsd_prefers_other",
        "geometric_outlier",
        "non_canonical",
        "dssr_questionable",
    ]

    def __init__(
        self,
        pdb_dir: Path = Path("data/pdb"),
        hbond_dir: Path = Path("data/json/slot_hbonds"),
        dssr_dir: Path = Path("data/json_dssr"),
        idealized_dir: Path = Path("basepair-idealized"),
        exemplar_dir: Path = Path("basepair-exemplars"),
    ):
        """Initialize the miss annotator."""
        self.pdb_dir = Path(pdb_dir)
        self.hbond_dir = Path(hbond_dir)
        self.dssr_dir = Path(dssr_dir)

        try:
            self.template_aligner = TemplateAligner(
                Path(idealized_dir), Path(exemplar_dir)
            )
            self.has_templates = True
        except Exception as e:
            self.template_aligner = None
            self.has_templates = False

    def annotate_pdb(self, pdb_id: str) -> Optional[PDBReport]:
        """Annotate all cWW classification differences for a PDB.

        Args:
            pdb_id: PDB identifier (e.g., "1EHZ").

        Returns:
            PDBReport with summary counts and detailed annotations,
            or None if DSSR data not available.
        """
        dssr_path = self.dssr_dir / f"{pdb_id}.json"
        if not dssr_path.exists():
            return None

        hbond_path = self.hbond_dir / f"{pdb_id}.json"
        pdb_path = self.pdb_dir / f"{pdb_id}.pdb"

        dssr_pairs = load_dssr_pairs(dssr_path, lw_filter="cWW")
        slot_hbonds = load_slot_hbonds(hbond_path) if hbond_path.exists() else {}
        residues = parse_pdb(pdb_path) if pdb_path.exists() else {}

        fn_annotations = []
        true_positives = 0

        for (res_id1, res_id2), dssr_pair in dssr_pairs.items():
            norm_id1 = self._normalize_res_id(res_id1)
            norm_id2 = self._normalize_res_id(res_id2)

            hbonds = slot_hbonds.get((res_id1, res_id2), [])
            if not hbonds:
                hbonds = slot_hbonds.get((norm_id1, norm_id2), [])

            sequence = dssr_pair.bp.replace("-", "")
            h_diag = self._analyze_hbonds(sequence, hbonds, dssr_pair.hbonds_desc)

            res1_atoms = self._extract_atom_coords(residues, res_id1, norm_id1)
            res2_atoms = self._extract_atom_coords(residues, res_id2, norm_id2)

            g_diag = self._analyze_geometry(
                res1_atoms,
                res2_atoms,
                sequence,
                dssr_pair.interbase_angle,
                dssr_pair.n1n9_distance,
            )

            reasons = self._categorize_miss(
                h_diag, g_diag, dssr_pair, res1_atoms, res2_atoms
            )

            if reasons:
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

        return PDBReport(
            pdb_id=pdb_id,
            total_canonical_cww=len(dssr_pairs),
            true_positives=true_positives,
            false_negatives=len(fn_annotations),
            false_positives=0,
            fn_annotations=fn_annotations,
            fp_annotations=[],
        )

    def _normalize_res_id(self, res_id: str) -> str:
        """Normalize DNA res_ids to match slot H-bond format."""
        parts = res_id.split("-")
        if len(parts) >= 2:
            base = parts[1]
            if base.upper() in ("DG", "DC", "DA", "DT"):
                parts[1] = base[1].upper()
                return "-".join(parts)
        return res_id

    def _extract_atom_coords(
        self, residues: Dict, res_id: str, norm_id: str
    ) -> Dict[str, np.ndarray]:
        """Extract atom coordinates, trying both IDs."""
        if res_id in residues:
            return {name: atom.coords for name, atom in residues[res_id].atoms.items()}
        if norm_id in residues:
            return {
                name: atom.coords for name, atom in residues[norm_id].atoms.items()
            }
        return {}

    def _analyze_hbonds(
        self, sequence: str, found_hbonds: List[SlotHBond], dssr_hbonds_desc: str
    ) -> HBondDiagnostics:
        """Analyze H-bond pattern for a base pair."""
        expected = CWW_EXPECTED_PATTERNS.get(sequence, [])
        expected_hbonds = [
            ExpectedHBond(donor_atom=d, acceptor_atom=a, source="canonical")
            for d, a in expected
        ]

        found_dicts = [
            {
                "donor_res_id": hb.donor_res_id,
                "donor_atom": hb.donor_atom,
                "acceptor_res_id": hb.acceptor_res_id,
                "acceptor_atom": hb.acceptor_atom,
                "distance": hb.distance,
                "context": hb.context,
                "h_slot": hb.h_slot,
                "lp_slot": hb.lp_slot,
                "alignment": hb.alignment,
            }
            for hb in found_hbonds
        ]

        diagnostics = HBondDiagnostics(
            expected_hbonds=expected_hbonds, found_hbonds=found_dicts
        )

        missing = []
        extra = []
        for exp_hb in expected_hbonds:
            if not self._match_hbond(found_hbonds, exp_hb.donor_atom, exp_hb.acceptor_atom):
                missing.append(exp_hb)

        for hb in found_hbonds:
            if not self._is_expected(hb, expected):
                extra.append(
                    {
                        "donor_atom": hb.donor_atom,
                        "acceptor_atom": hb.acceptor_atom,
                        "distance": hb.distance,
                        "context": hb.context,
                    }
                )

        diagnostics.missing_hbonds = missing
        diagnostics.extra_hbonds = extra

        wrong_atoms = {}
        for hb in found_hbonds:
            wrong_desc = self._find_wrong_atom(hb, expected)
            if wrong_desc:
                key = f"{hb.donor_atom}->{hb.acceptor_atom}"
                wrong_atoms[key] = wrong_desc
        diagnostics.wrong_atoms = wrong_atoms

        distance_issues = []
        for hb in found_hbonds:
            if hb.distance < 2.0:
                distance_issues.append((f"{hb.donor_atom}-{hb.acceptor_atom}", hb.distance, "too_short"))
            elif hb.distance > 4.0:
                distance_issues.append((f"{hb.donor_atom}-{hb.acceptor_atom}", hb.distance, "too_long"))
        diagnostics.distance_issues = distance_issues

        acceptor_counts = Counter(hb.acceptor_atom for hb in found_hbonds)
        diagnostics.overloaded_acceptors = [
            atom for atom, count in acceptor_counts.items() if count > 2
        ]

        return diagnostics

    def _match_hbond(
        self, found_hbonds: List[SlotHBond], donor: str, acceptor: str
    ) -> Optional[SlotHBond]:
        """Find a matching H-bond in the found list."""
        for hb in found_hbonds:
            if hb.donor_atom == donor and hb.acceptor_atom == acceptor:
                return hb
            if hb.donor_atom == acceptor and hb.acceptor_atom == donor:
                return hb
        return None

    def _is_expected(
        self, hb: SlotHBond, expected: List[Tuple[str, str]]
    ) -> bool:
        """Check if an H-bond is in the expected pattern."""
        for donor, acceptor in expected:
            if hb.donor_atom == donor and hb.acceptor_atom == acceptor:
                return True
            if hb.donor_atom == acceptor and hb.acceptor_atom == donor:
                return True
        return False

    def _find_wrong_atom(
        self, hb: SlotHBond, expected: List[Tuple[str, str]]
    ) -> Optional[str]:
        """Check if donor matches but acceptor differs from expected."""
        for donor, acceptor in expected:
            if hb.donor_atom == donor and hb.acceptor_atom != acceptor:
                return f"Expected {donor}->{acceptor}, got {donor}->{hb.acceptor_atom}"
            if hb.donor_atom == acceptor and hb.acceptor_atom != donor:
                return f"Expected {acceptor}->{donor}, got {acceptor}->{hb.acceptor_atom}"
        return None

    def _analyze_geometry(
        self,
        res1_atoms: Dict[str, np.ndarray],
        res2_atoms: Dict[str, np.ndarray],
        sequence: str,
        interbase_angle: Optional[float],
        n1n9_distance: Optional[float],
    ) -> GeometricDiagnostics:
        """Analyze geometric fit to templates."""
        if not self.has_templates or not res1_atoms or not res2_atoms:
            return GeometricDiagnostics(
                rmsd_cww=999.0,
                rmsd_best=999.0,
                best_lw="unknown",
                rmsd_gap=0.0,
                interbase_angle=interbase_angle or 0.0,
                n1n9_distance=n1n9_distance or 0.0,
                is_geometric_outlier=False,
            )

        if n1n9_distance is None:
            n1n9_distance = self._compute_n1n9_distance(res1_atoms, res2_atoms, sequence)
        if interbase_angle is None:
            interbase_angle = self._compute_interbase_angle(res1_atoms, res2_atoms)

        res1 = self._create_residue("1", sequence[0], res1_atoms)
        res2 = self._create_residue("2", sequence[1], res2_atoms)

        cww_rmsd = self._align_to_lw_class(res1, res2, sequence, "cWW")

        all_rmsds = {}
        for lw_class in self.template_aligner.DEFAULT_LW_CLASSES:
            rmsd = self._align_to_lw_class(res1, res2, sequence, lw_class)
            if rmsd < float("inf"):
                all_rmsds[lw_class] = rmsd

        if not all_rmsds:
            return GeometricDiagnostics(
                rmsd_cww=cww_rmsd,
                rmsd_best=cww_rmsd,
                best_lw="cWW",
                rmsd_gap=0.0,
                interbase_angle=interbase_angle,
                n1n9_distance=n1n9_distance,
                is_geometric_outlier=False,
            )

        sorted_rmsds = sorted(all_rmsds.items(), key=lambda x: x[1])
        best_lw, best_rmsd = sorted_rmsds[0]
        rmsd_gap = cww_rmsd - best_rmsd

        n1n9_outlier = self._is_n1n9_outlier(sequence, n1n9_distance)
        angle_outlier = abs(interbase_angle) > CWW_ANGLE_THRESHOLD
        is_geometric_outlier = n1n9_outlier or angle_outlier

        return GeometricDiagnostics(
            rmsd_cww=cww_rmsd,
            rmsd_best=best_rmsd,
            best_lw=best_lw,
            rmsd_gap=rmsd_gap,
            interbase_angle=interbase_angle,
            n1n9_distance=n1n9_distance,
            is_geometric_outlier=is_geometric_outlier,
        )

    def _align_to_lw_class(
        self, res1: Residue, res2: Residue, sequence: str, lw_class: str
    ) -> float:
        """Align to template for specific LW class."""
        template_path = self.template_aligner._find_template(sequence, lw_class)
        if template_path is None:
            return float("inf")

        try:
            rmsd, _ = self.template_aligner.align_to_template(res1, res2, template_path)
            return rmsd
        except Exception:
            return float("inf")

    def _compute_n1n9_distance(
        self, res1_atoms: Dict, res2_atoms: Dict, sequence: str
    ) -> float:
        """Compute distance between glycosidic nitrogens."""
        n1_atom = GLYCOSIDIC_NITROGEN.get(sequence[0])
        n2_atom = GLYCOSIDIC_NITROGEN.get(sequence[1])

        if not n1_atom or not n2_atom:
            return 0.0
        if n1_atom not in res1_atoms or n2_atom not in res2_atoms:
            return 0.0

        return float(np.linalg.norm(res1_atoms[n1_atom] - res2_atoms[n2_atom]))

    def _compute_interbase_angle(self, res1_atoms: Dict, res2_atoms: Dict) -> float:
        """Compute angle between base planes."""
        ring_atoms = ["C2", "C4", "C5", "C6", "N1", "N3"]
        coords1 = [res1_atoms[atom] for atom in ring_atoms if atom in res1_atoms]
        coords2 = [res2_atoms[atom] for atom in ring_atoms if atom in res2_atoms]

        if len(coords1) < 3 or len(coords2) < 3:
            return 0.0

        normal1 = self._fit_plane_normal(np.array(coords1))
        normal2 = self._fit_plane_normal(np.array(coords2))

        cos_angle = np.clip(np.dot(normal1, normal2), -1.0, 1.0)
        angle = np.arccos(abs(cos_angle))
        return float(np.degrees(angle))

    @staticmethod
    def _fit_plane_normal(coords: np.ndarray) -> np.ndarray:
        """Fit plane to points and return normal vector."""
        centroid = np.mean(coords, axis=0)
        centered = coords - centroid
        _, _, Vt = np.linalg.svd(centered)
        normal = Vt[-1, :]
        return normal / np.linalg.norm(normal)

    @staticmethod
    def _is_n1n9_outlier(sequence: str, distance: float) -> bool:
        """Check if N1-N9 distance is outside expected range."""
        if sequence not in CWW_N1N9_RANGE:
            return False
        min_dist, max_dist = CWW_N1N9_RANGE[sequence]
        return distance < min_dist or distance > max_dist

    @staticmethod
    def _create_residue(res_id: str, base_type: str, atoms: Dict) -> Residue:
        """Create Residue object from atom dictionary."""
        residue = Residue(res_id=res_id, base_type=base_type)
        for atom_name, coords in atoms.items():
            residue.atoms[atom_name] = Atom(name=atom_name, coords=coords)
        return residue

    def _categorize_miss(
        self,
        h_diag: HBondDiagnostics,
        g_diag: GeometricDiagnostics,
        dssr_pair: DSSRPair,
        res1_atoms: Dict,
        res2_atoms: Dict,
    ) -> List[str]:
        """Determine why a pair would be misclassified.

        Returns:
            List of reason codes, empty if no miss.
        """
        from prototypes.pair_identification.cww_miss_annotator.bp_score import compute_bp_score

        sequence = dssr_pair.bp.replace("-", "")
        bp_score, bp_components = compute_bp_score(
            sequence,
            g_diag.rmsd_cww,
            h_diag.found_hbonds,
            res1_atoms=res1_atoms,
            res2_atoms=res2_atoms,
            interbase_angle=g_diag.interbase_angle,
        )

        threshold = 0.60 if bp_components.get("extended_search", False) else 0.70
        if bp_score >= threshold:
            return []

        reasons = []

        if not h_diag.found_hbonds:
            reasons.append("no_hbonds")
        else:
            if h_diag.missing_hbonds:
                reasons.append("missing_hbonds")
            if h_diag.wrong_atoms:
                reasons.append("wrong_hbonds")
            if h_diag.extra_hbonds:
                reasons.append("extra_hbonds")

            has_long = any(issue[2] == "too_long" for issue in h_diag.distance_issues)
            has_short = any(issue[2] == "too_short" for issue in h_diag.distance_issues)
            if has_long:
                reasons.append("long_hbonds")
            if has_short:
                reasons.append("short_hbonds")

            if h_diag.overloaded_acceptors:
                reasons.append("overloaded_acceptor")

        # Only flag rmsd_prefers_other if best match is a DIFFERENT edge pair
        # (cWW vs tWW is NOT flagged since they're both "WW" edge pair)
        if g_diag.rmsd_gap > 0.5 and g_diag.best_edge_pair != "WW":
            reasons.append("rmsd_prefers_other")

        if g_diag.is_geometric_outlier:
            reasons.append("geometric_outlier")

        if g_diag.interbase_angle >= 20.0:
            reasons.append("poor_planarity")

        if dssr_pair.saenger == "--":
            reasons.append("non_canonical")

        if g_diag.rmsd_cww > 1.0 and not h_diag.found_hbonds:
            reasons.append("dssr_questionable")

        return reasons

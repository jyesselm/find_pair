"""High-level pair visualization API."""

import sys
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

import numpy as np

# Add parent paths for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from core.pdb_parser import (
    parse_pdb,
    parse_template_pdb,
    write_pair_pdb,
    extract_pair_from_pdb,
)
from core.alignment import align_atom_dicts, kabsch_align
from core.residue import Residue
from core.identifiers import extract_sequence

from .pymol_generator import PyMOLScriptGenerator, HBondViz, TemplateViz


@dataclass
class VisualizationResult:
    """Result of visualizing a pair."""

    pdb_id: str
    res_id1: str
    res_id2: str
    sequence: str

    target_pdb: Path
    template_pdbs: Dict[str, Path]  # lw_class -> path
    pymol_script: Path

    rmsd_values: Dict[str, float]
    best_lw: str
    best_rmsd: float


class PairVisualizer:
    """Visualize base pairs with templates and H-bonds."""

    # Ring atoms for alignment
    RING_ATOMS = ["C2", "C4", "C5", "C6", "N1", "N3", "N7", "C8", "N9"]

    # LW classes to show
    DEFAULT_LW_CLASSES = ["cWW", "tWW", "cWS", "tHS"]

    # Template colors
    TEMPLATE_COLORS = {
        "cWW": "green",
        "tWW": "orange",
        "cWS": "magenta",
        "tWH": "yellow",
        "tHS": "pink",
        "cHS": "salmon",
    }

    def __init__(
        self,
        pdb_dir: Path = Path("data/pdb"),
        template_dir: Path = Path("basepair-idealized"),
        hbond_dir: Path = Path("data/json/slot_hbonds"),
        output_dir: Path = Path("viz_output"),
    ):
        self.pdb_dir = Path(pdb_dir)
        self.template_dir = Path(template_dir)
        self.hbond_dir = Path(hbond_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.pymol_gen = PyMOLScriptGenerator(self.output_dir)

    def visualize_pair(
        self,
        pdb_id: str,
        res_id1: str,
        res_id2: str,
        lw_classes: Optional[List[str]] = None,
        include_hbonds: bool = True,
    ) -> Optional[VisualizationResult]:
        """Generate complete visualization for a base pair.

        Args:
            pdb_id: PDB identifier (e.g., "1EHZ")
            res_id1: First residue ID (e.g., "A-G-1")
            res_id2: Second residue ID (e.g., "A-C-72")
            lw_classes: LW classes to show (default: cWW, tWW, cWS, tHS)
            include_hbonds: Whether to load and show H-bonds

        Returns:
            VisualizationResult with paths to generated files
        """
        if lw_classes is None:
            lw_classes = self.DEFAULT_LW_CLASSES

        # Load PDB
        pdb_path = self.pdb_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            pdb_path = self.pdb_dir / f"{pdb_id.lower()}.pdb"
        if not pdb_path.exists():
            print(f"PDB not found: {pdb_path}")
            return None

        res1, res2 = extract_pair_from_pdb(pdb_path, res_id1, res_id2)
        if res1 is None or res2 is None:
            print(f"Residues not found: {res_id1}, {res_id2}")
            return None

        sequence = res1.base_type + res2.base_type

        # Extract target atoms
        target_res1 = {n: a.coords for n, a in res1.atoms.items()}
        target_res2 = {n: a.coords for n, a in res2.atoms.items()}

        # Write target PDB
        target_name = f"target_{pdb_id}_{sequence}"
        target_pdb = self.output_dir / f"{target_name}.pdb"
        write_pair_pdb(
            target_res1, target_res2, res1.base_type, res2.base_type, target_pdb
        )

        # Align templates and write aligned PDBs
        template_pdbs = {}
        rmsd_values = {}
        template_vizs = []

        for lw_class in lw_classes:
            template_path = self._find_template(sequence, lw_class)
            if template_path is None:
                continue

            # Load and align template
            tmpl_res1, tmpl_res2 = parse_template_pdb(template_path)

            aligned_res1, aligned_res2, rmsd = self._align_template_to_target(
                tmpl_res1, tmpl_res2, target_res1, target_res2
            )

            if rmsd < float("inf"):
                # Write aligned template
                aligned_name = f"aligned_{pdb_id}_{lw_class}_{sequence}"
                aligned_path = self.output_dir / f"{aligned_name}.pdb"
                write_pair_pdb(
                    aligned_res1,
                    aligned_res2,
                    res1.base_type,
                    res2.base_type,
                    aligned_path,
                )

                template_pdbs[lw_class] = aligned_path
                rmsd_values[lw_class] = rmsd

                color = self.TEMPLATE_COLORS.get(lw_class, "gray")
                template_vizs.append(
                    TemplateViz(
                        name=lw_class,
                        pdb_path=aligned_path,
                        rmsd=rmsd,
                        color=color,
                        enabled=True,
                    )
                )

        # Find best template
        if rmsd_values:
            best_lw = min(rmsd_values, key=rmsd_values.get)
            best_rmsd = rmsd_values[best_lw]
        else:
            best_lw = "unknown"
            best_rmsd = float("inf")

        # Load H-bonds if requested
        observed_hbonds = []
        expected_hbonds = []

        if include_hbonds:
            observed_hbonds, expected_hbonds = self._load_hbonds(
                pdb_id, res_id1, res_id2, sequence
            )

        # Generate PyMOL script
        script_name = f"view_{pdb_id}_{res_id1}_{res_id2}".replace("-", "_")
        title = f"{pdb_id} {res_id1} - {res_id2} ({sequence})"

        pymol_script = self.pymol_gen.generate_pair_script(
            target_pdb=target_pdb,
            templates=template_vizs,
            observed_hbonds=observed_hbonds,
            expected_hbonds=expected_hbonds,
            output_name=script_name,
            title=title,
        )

        return VisualizationResult(
            pdb_id=pdb_id,
            res_id1=res_id1,
            res_id2=res_id2,
            sequence=sequence,
            target_pdb=target_pdb,
            template_pdbs=template_pdbs,
            pymol_script=pymol_script,
            rmsd_values=rmsd_values,
            best_lw=best_lw,
            best_rmsd=best_rmsd,
        )

    def _find_template(self, sequence: str, lw_class: str) -> Optional[Path]:
        """Find template PDB for sequence and LW class."""
        # Try idealized templates
        template_path = self.template_dir / lw_class / f"{sequence}.pdb"
        if template_path.exists():
            return template_path

        # Try reversed sequence
        rev_seq = sequence[::-1]
        template_path = self.template_dir / lw_class / f"{rev_seq}.pdb"
        if template_path.exists():
            return template_path

        return None

    def _align_template_to_target(
        self,
        tmpl_res1: Dict[str, np.ndarray],
        tmpl_res2: Dict[str, np.ndarray],
        target_res1: Dict[str, np.ndarray],
        target_res2: Dict[str, np.ndarray],
    ) -> Tuple[Dict, Dict, float]:
        """Align template atoms to target atoms."""
        # Combine atoms from both residues for alignment
        combined_tmpl = {}
        combined_target = {}

        for atom in self.RING_ATOMS:
            if atom in tmpl_res1 and atom in target_res1:
                combined_tmpl[f"1_{atom}"] = tmpl_res1[atom]
                combined_target[f"1_{atom}"] = target_res1[atom]
            if atom in tmpl_res2 and atom in target_res2:
                combined_tmpl[f"2_{atom}"] = tmpl_res2[atom]
                combined_target[f"2_{atom}"] = target_res2[atom]

        if len(combined_tmpl) < 6:
            return tmpl_res1, tmpl_res2, float("inf")

        # Align using combined atoms
        transformed, rmsd, _ = align_atom_dicts(combined_tmpl, combined_target)

        # Apply same transformation to all atoms
        # Get transformation parameters from aligned combined atoms
        common_atoms = sorted(set(combined_tmpl.keys()) & set(combined_target.keys()))
        source_coords = np.array([combined_tmpl[a] for a in common_atoms])
        target_coords = np.array([combined_target[a] for a in common_atoms])

        R, centroid_target, centroid_source = kabsch_align(source_coords, target_coords)

        # Transform all template atoms
        aligned_res1 = {}
        for atom, coords in tmpl_res1.items():
            aligned_res1[atom] = (np.array(coords) - centroid_source) @ R + centroid_target

        aligned_res2 = {}
        for atom, coords in tmpl_res2.items():
            aligned_res2[atom] = (np.array(coords) - centroid_source) @ R + centroid_target

        return aligned_res1, aligned_res2, rmsd

    def _load_hbonds(
        self,
        pdb_id: str,
        res_id1: str,
        res_id2: str,
        sequence: str,
    ) -> Tuple[List[HBondViz], List[HBondViz]]:
        """Load observed and expected H-bonds."""
        from hbond.patterns import get_cww_expected

        observed = []
        expected = []

        # Load slot H-bonds
        hbond_path = self.hbond_dir / f"{pdb_id}.json"
        if hbond_path.exists():
            with open(hbond_path) as f:
                hbond_data = json.load(f)

            # Find H-bonds for this pair
            for entry in hbond_data:
                if (
                    entry["res_id_i"] == res_id1 and entry["res_id_j"] == res_id2
                ) or (entry["res_id_i"] == res_id2 and entry["res_id_j"] == res_id1):

                    for hb in entry.get("hbonds", []):
                        if hb.get("context") == "base_base":
                            # Use explicit donor/acceptor residue IDs from JSON
                            donor_res_id = hb.get("donor_res_id", entry["res_id_i"])
                            acceptor_res_id = hb.get("acceptor_res_id", entry["res_id_j"])

                            # Map to residue numbers (1 or 2)
                            donor_res = 1 if donor_res_id == res_id1 else 2
                            acceptor_res = 1 if acceptor_res_id == res_id1 else 2

                            observed.append(
                                HBondViz(
                                    donor_res=donor_res,
                                    donor_atom=hb["donor_atom"],
                                    acceptor_res=acceptor_res,
                                    acceptor_atom=hb["acceptor_atom"],
                                    distance=hb["distance"],
                                    color="yellow",
                                )
                            )
                    break

        # Get expected H-bonds for cWW
        cww_expected = get_cww_expected(sequence)

        for donor, acceptor in cww_expected:
            # Check if this is in observed
            found = any(
                (h.donor_atom == donor and h.acceptor_atom == acceptor)
                or (h.donor_atom == acceptor and h.acceptor_atom == donor)
                for h in observed
            )

            expected.append(
                HBondViz(
                    donor_res=1,
                    donor_atom=donor,
                    acceptor_res=2,
                    acceptor_atom=acceptor,
                    distance=0.0,
                    color="green" if found else "red",
                )
            )

        return observed, expected

#!/usr/bin/env python3
"""Geometric analyzer for cWW miss diagnostics.

This module computes RMSD-based diagnostics using template alignment
to determine if a pair is a geometric outlier or fits better to a
different LW class.
"""

import sys
import numpy as np
from pathlib import Path
from typing import Dict, Optional, Tuple

# Import from parent directory core modules
sys.path.insert(0, str(Path(__file__).parent.parent))
from core.residue import Residue, Atom
from template_aligner import TemplateAligner

# Import diagnostics from current directory
from .diagnostics import GeometricDiagnostics


# Normal N1-N9 distance ranges for cWW pairs (Angstroms)
# Widened thresholds - if RMSD to cWW template is good, the pair is cWW
CWW_N1N9_RANGE = {
    "GC": (7.5, 10.0),
    "CG": (7.5, 10.0),
    "AU": (7.5, 10.0),
    "UA": (7.5, 10.0),
    "GU": (7.0, 10.0),
    "UG": (7.0, 10.0),
}

# Threshold for interbase angle outliers (degrees)
# Widened from 15 to 30 - real cWW pairs can have significant buckling
CWW_ANGLE_THRESHOLD = 30.0

# Glycosidic nitrogen mapping
GLYCOSIDIC_NITROGEN = {
    "A": "N9",
    "G": "N9",
    "C": "N1",
    "U": "N1",
    "T": "N1",
}


class GeometricAnalyzer:
    """Analyzes geometric fit to templates using RMSD."""

    def __init__(
        self,
        idealized_dir: Optional[Path] = None,
        exemplar_dir: Optional[Path] = None,
    ):
        """Initialize with template directories.

        Args:
            idealized_dir: Directory with idealized templates.
                Default: basepair-idealized
            exemplar_dir: Directory with exemplar templates.
                Default: basepair-exemplars
        """
        if idealized_dir is None:
            idealized_dir = Path("basepair-idealized")
        if exemplar_dir is None:
            exemplar_dir = Path("basepair-exemplars")

        self.idealized_dir = idealized_dir
        self.exemplar_dir = exemplar_dir

        # Initialize template aligner
        try:
            self.aligner = TemplateAligner(idealized_dir, exemplar_dir)
            self.has_templates = True
        except Exception as e:
            self.aligner = None
            self.has_templates = False
            self.template_error = str(e)

    def analyze(
        self,
        res1_atoms: Dict[str, np.ndarray],
        res2_atoms: Dict[str, np.ndarray],
        sequence: str,
        interbase_angle: Optional[float] = None,
        n1n9_distance: Optional[float] = None,
    ) -> GeometricDiagnostics:
        """Analyze geometric fit to templates.

        Args:
            res1_atoms: Dict mapping atom name to 3D coordinates for residue 1
            res2_atoms: Dict mapping atom name to 3D coordinates for residue 2
            sequence: Two-letter sequence (e.g., "GC")
            interbase_angle: From DSSR if available, otherwise computed
            n1n9_distance: From DSSR if available, otherwise computed

        Returns:
            GeometricDiagnostics with RMSD analysis
        """
        # Check if templates are available
        if not self.has_templates:
            return GeometricDiagnostics(
                rmsd_cww=999.0,
                rmsd_best=999.0,
                best_lw="unknown",
                rmsd_gap=0.0,
                interbase_angle=interbase_angle or 0.0,
                n1n9_distance=n1n9_distance or 0.0,
                is_geometric_outlier=False,
            )

        # Compute N1-N9 distance if not provided
        if n1n9_distance is None:
            n1n9_distance = self._compute_n1n9_distance(
                res1_atoms, res2_atoms, sequence
            )

        # Compute interbase angle if not provided
        if interbase_angle is None:
            interbase_angle = self._compute_interbase_angle(
                res1_atoms, res2_atoms
            )

        # Convert to Residue objects for template aligner
        res1 = self._create_residue("1", sequence[0], res1_atoms)
        res2 = self._create_residue("2", sequence[1], res2_atoms)

        # Try aligning to cWW template
        cww_rmsd = self._align_to_lw_class(res1, res2, sequence, "cWW")

        # Try aligning to all LW classes
        all_rmsds = {}
        for lw_class in self.aligner.DEFAULT_LW_CLASSES:
            rmsd = self._align_to_lw_class(res1, res2, sequence, lw_class)
            if rmsd < float('inf'):
                all_rmsds[lw_class] = rmsd

        # Find best and second best
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
        second_lw = sorted_rmsds[1][0] if len(sorted_rmsds) > 1 else None
        second_rmsd = sorted_rmsds[1][1] if len(sorted_rmsds) > 1 else None

        # Compute RMSD gap (positive = cWW worse than best)
        rmsd_gap = cww_rmsd - best_rmsd

        # Check for geometric outliers
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
        self,
        res1: Residue,
        res2: Residue,
        sequence: str,
        lw_class: str,
    ) -> float:
        """Align to template for specific LW class.

        Args:
            res1: First residue
            res2: Second residue
            sequence: Two-letter sequence
            lw_class: LW class to try

        Returns:
            RMSD to template, or inf if template not found
        """
        template_path = self.aligner._find_template(sequence, lw_class)
        if template_path is None:
            return float('inf')

        try:
            rmsd, num_atoms = self.aligner.align_to_template(
                res1, res2, template_path
            )
            return rmsd
        except Exception:
            return float('inf')

    def _compute_n1n9_distance(
        self,
        res1_atoms: Dict[str, np.ndarray],
        res2_atoms: Dict[str, np.ndarray],
        sequence: str,
    ) -> float:
        """Compute distance between glycosidic nitrogens.

        Args:
            res1_atoms: First residue atoms
            res2_atoms: Second residue atoms
            sequence: Two-letter sequence

        Returns:
            N1-N9 distance in Angstroms, or 0.0 if atoms missing
        """
        base1, base2 = sequence[0], sequence[1]

        n1_atom = GLYCOSIDIC_NITROGEN.get(base1)
        n2_atom = GLYCOSIDIC_NITROGEN.get(base2)

        if n1_atom is None or n2_atom is None:
            return 0.0

        if n1_atom not in res1_atoms or n2_atom not in res2_atoms:
            return 0.0

        coords1 = res1_atoms[n1_atom]
        coords2 = res2_atoms[n2_atom]

        return float(np.linalg.norm(coords1 - coords2))

    def _compute_interbase_angle(
        self,
        res1_atoms: Dict[str, np.ndarray],
        res2_atoms: Dict[str, np.ndarray],
    ) -> float:
        """Compute angle between base planes.

        Uses best-fit plane through ring atoms.

        Args:
            res1_atoms: First residue atoms
            res2_atoms: Second residue atoms

        Returns:
            Angle in degrees, or 0.0 if insufficient atoms
        """
        # Get ring atoms for plane fitting
        ring_atoms = ["C2", "C4", "C5", "C6", "N1", "N3"]

        # Collect coordinates
        coords1 = [
            res1_atoms[atom] for atom in ring_atoms if atom in res1_atoms
        ]
        coords2 = [
            res2_atoms[atom] for atom in ring_atoms if atom in res2_atoms
        ]

        if len(coords1) < 3 or len(coords2) < 3:
            return 0.0

        # Fit planes using SVD
        normal1 = self._fit_plane_normal(np.array(coords1))
        normal2 = self._fit_plane_normal(np.array(coords2))

        # Compute angle between normals
        cos_angle = np.dot(normal1, normal2)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle = np.arccos(abs(cos_angle))  # Use abs for undirected angle

        return float(np.degrees(angle))

    @staticmethod
    def _fit_plane_normal(coords: np.ndarray) -> np.ndarray:
        """Fit plane to points and return normal vector.

        Args:
            coords: N x 3 array of coordinates

        Returns:
            Unit normal vector (3,)
        """
        # Center the points
        centroid = np.mean(coords, axis=0)
        centered = coords - centroid

        # SVD - normal is last right singular vector
        _, _, Vt = np.linalg.svd(centered)
        normal = Vt[-1, :]

        return normal / np.linalg.norm(normal)

    @staticmethod
    def _is_n1n9_outlier(sequence: str, distance: float) -> bool:
        """Check if N1-N9 distance is outside expected range.

        Args:
            sequence: Two-letter sequence
            distance: N1-N9 distance in Angstroms

        Returns:
            True if outside expected range for cWW
        """
        if sequence not in CWW_N1N9_RANGE:
            return False

        min_dist, max_dist = CWW_N1N9_RANGE[sequence]
        return distance < min_dist or distance > max_dist

    @staticmethod
    def _create_residue(
        res_id: str,
        base_type: str,
        atoms: Dict[str, np.ndarray],
    ) -> Residue:
        """Create Residue object from atom dictionary.

        Args:
            res_id: Residue ID
            base_type: Single-letter base type
            atoms: Dict mapping atom name to coordinates

        Returns:
            Residue object for template aligner
        """
        residue = Residue(res_id=res_id, base_type=base_type)

        for atom_name, coords in atoms.items():
            residue.atoms[atom_name] = Atom(
                name=atom_name,
                coords=coords,
            )

        return residue


def main():
    """Example usage."""
    import argparse
    from core.pdb_parser import parse_pdb

    parser = argparse.ArgumentParser(
        description="Analyze geometric fit to templates"
    )
    parser.add_argument(
        "--pdb", type=str, required=True,
        help="PDB ID or path to PDB file"
    )
    parser.add_argument(
        "--res1", type=str, required=True,
        help="First residue ID (e.g., A-G-1)"
    )
    parser.add_argument(
        "--res2", type=str, required=True,
        help="Second residue ID (e.g., A-C-72)"
    )
    parser.add_argument(
        "--idealized-dir", type=Path, default=Path("basepair-idealized"),
        help="Directory with idealized templates"
    )
    parser.add_argument(
        "--exemplar-dir", type=Path, default=Path("basepair-exemplars"),
        help="Directory with exemplar templates"
    )
    parser.add_argument(
        "--pdb-dir", type=Path, default=Path("data/pdb"),
        help="Directory with PDB files"
    )

    args = parser.parse_args()

    # Load PDB
    if Path(args.pdb).exists():
        pdb_path = Path(args.pdb)
    else:
        pdb_path = args.pdb_dir / f"{args.pdb}.pdb"
        if not pdb_path.exists():
            pdb_path = args.pdb_dir / f"{args.pdb.lower()}.pdb"

    residues = parse_pdb(pdb_path)

    if args.res1 not in residues:
        print(f"Residue {args.res1} not found in PDB")
        return
    if args.res2 not in residues:
        print(f"Residue {args.res2} not found in PDB")
        return

    res1 = residues[args.res1]
    res2 = residues[args.res2]

    # Get atom dictionaries
    res1_atoms = {name: atom.coords for name, atom in res1.atoms.items()}
    res2_atoms = {name: atom.coords for name, atom in res2.atoms.items()}

    sequence = res1.base_type + res2.base_type

    # Analyze
    analyzer = GeometricAnalyzer(args.idealized_dir, args.exemplar_dir)
    diagnostics = analyzer.analyze(res1_atoms, res2_atoms, sequence)

    # Print results
    print(f"\nGeometric Analysis for {args.res1} - {args.res2} ({sequence}):")
    print(f"\n  cWW RMSD: {diagnostics.rmsd_cww:.3f} Å")
    print(f"  Best RMSD: {diagnostics.rmsd_best:.3f} Å ({diagnostics.best_lw})")
    print(f"  RMSD gap: {diagnostics.rmsd_gap:.3f} Å")
    print(f"\n  N1-N9 distance: {diagnostics.n1n9_distance:.2f} Å")
    print(f"  Interbase angle: {diagnostics.interbase_angle:.1f}°")
    print(f"\n  Geometric outlier: {diagnostics.is_geometric_outlier}")


if __name__ == "__main__":
    main()

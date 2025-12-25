#!/usr/bin/env python3
"""Template-based pair classification using Kabsch alignment.

Aligns target pairs to templates for each LW class and computes RMSD.
The best-matching template indicates the LW classification.
"""

import json
import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

# Import from core modules
from core.residue import Residue, Atom
from core.pdb_parser import parse_pdb, parse_template_pdb
from core.alignment import kabsch_align, compute_rmsd


@dataclass
class AlignmentResult:
    """Result of aligning a target pair to a template."""
    lw_class: str
    sequence: str
    rmsd: float
    num_atoms_aligned: int
    template_path: Optional[Path] = None

    @property
    def score(self) -> float:
        """Score combining RMSD and atom count. Lower is better."""
        # Penalize templates with few atoms (less reliable)
        # RMSD alone isn't enough - 6 atoms with 0.01Å RMSD is worse than
        # 15 atoms with 0.15Å RMSD
        min_atoms = 10
        if self.num_atoms_aligned < min_atoms:
            # Heavy penalty for too few atoms
            return self.rmsd + (min_atoms - self.num_atoms_aligned) * 0.5
        return self.rmsd


@dataclass
class ClassificationResult:
    """Result of classifying a pair by trying all LW templates."""
    target_res1: str
    target_res2: str
    sequence: str

    # Best match
    best_lw: str
    best_rmsd: float

    # Second best match
    second_lw: Optional[str] = None
    second_rmsd: Optional[float] = None

    # All results sorted by RMSD
    all_results: List[AlignmentResult] = field(default_factory=list)

    @property
    def confidence(self) -> float:
        """Confidence based on gap between best and second best RMSD."""
        if self.second_rmsd is None:
            return 1.0
        gap = self.second_rmsd - self.best_rmsd
        # Normalize: 0.5Å gap = high confidence
        return min(1.0, gap / 0.5)


class TemplateAligner:
    """Align target pairs to templates and compute RMSD."""

    # Ring atoms used for alignment
    RING_ATOMS = ["C2", "C4", "C5", "C6", "N1", "N3", "N7", "C8", "N9"]

    # LW classes to try
    DEFAULT_LW_CLASSES = ["cWW", "tWW", "cWH", "tWH", "cWS", "tWS",
                          "cHH", "tHH", "cHS", "tHS", "cSS", "tSS"]

    def __init__(
        self,
        idealized_dir: Path = Path("basepair-idealized"),
        exemplar_dir: Path = Path("basepair-exemplars"),
    ):
        self.idealized_dir = idealized_dir
        self.exemplar_dir = exemplar_dir
        self._template_cache: Dict[str, Tuple[Dict, Dict]] = {}

    def _load_template(self, template_path: Path) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        """Load template PDB and return atom coords for residue 1 and 2."""
        if str(template_path) in self._template_cache:
            return self._template_cache[str(template_path)]

        res1_atoms, res2_atoms = parse_template_pdb(template_path)
        self._template_cache[str(template_path)] = (res1_atoms, res2_atoms)
        return res1_atoms, res2_atoms

    def _find_template(self, sequence: str, lw_class: str) -> Optional[Path]:
        """Find template PDB file for given sequence and LW class."""
        # Try idealized first
        idealized_path = self.idealized_dir / lw_class / f"{sequence}.pdb"
        if idealized_path.exists():
            return idealized_path

        # Try with underscore format
        idealized_path = self.idealized_dir / lw_class / f"{sequence[0]}_{sequence[1].lower()}.pdb"
        if idealized_path.exists():
            return idealized_path

        # Try exemplars with various naming conventions
        for pattern in [
            f"{sequence[0]}-{sequence[1]}-{lw_class}.pdb",
            f"{sequence[0]}plus{sequence[1]}-{lw_class}.pdb",
            f"{sequence[0].lower()}-{sequence[1]}-{lw_class}.pdb",
            f"{sequence}-{lw_class}.pdb",
        ]:
            exemplar_path = self.exemplar_dir / pattern
            if exemplar_path.exists():
                return exemplar_path

        return None


    def align_to_template(
        self,
        target_res1: Residue,
        target_res2: Residue,
        template_path: Path,
    ) -> Tuple[float, int]:
        """Align target pair to template and compute RMSD.

        Args:
            target_res1: First residue of target pair
            target_res2: Second residue of target pair
            template_path: Path to template PDB

        Returns:
            rmsd: RMSD after optimal superposition
            num_atoms: Number of atoms used in alignment
        """
        template_res1, template_res2 = self._load_template(template_path)

        # Collect matching atoms
        template_points = []
        target_points = []

        for atom_name in self.RING_ATOMS:
            if atom_name in template_res1 and atom_name in target_res1.atoms:
                template_points.append(template_res1[atom_name])
                target_points.append(target_res1.atoms[atom_name].coords)
            if atom_name in template_res2 and atom_name in target_res2.atoms:
                template_points.append(template_res2[atom_name])
                target_points.append(target_res2.atoms[atom_name].coords)

        if len(template_points) < 4:
            return float('inf'), 0

        template_points = np.array(template_points)
        target_points = np.array(target_points)

        # Align template to target using core module
        R, centroid_target, centroid_template = kabsch_align(template_points, target_points)

        # Transform template points
        aligned_template = (template_points - centroid_template) @ R + centroid_target

        # Compute RMSD
        rmsd = compute_rmsd(aligned_template, target_points)

        return rmsd, len(template_points)

    def classify_pair(
        self,
        target_res1: Residue,
        target_res2: Residue,
        lw_classes: Optional[List[str]] = None,
    ) -> ClassificationResult:
        """Classify a pair by trying all LW class templates.

        Args:
            target_res1: First residue
            target_res2: Second residue
            lw_classes: LW classes to try (default: all 12)

        Returns:
            ClassificationResult with best match and confidence
        """
        if lw_classes is None:
            lw_classes = self.DEFAULT_LW_CLASSES

        # Determine sequence
        sequence = target_res1.base_type + target_res2.base_type

        results = []

        for lw in lw_classes:
            # Try forward sequence
            template_path = self._find_template(sequence, lw)
            if template_path:
                rmsd, num_atoms = self.align_to_template(
                    target_res1, target_res2, template_path
                )
                results.append(AlignmentResult(
                    lw_class=lw,
                    sequence=sequence,
                    rmsd=rmsd,
                    num_atoms_aligned=num_atoms,
                    template_path=template_path,
                ))

            # Try reversed sequence (swap residues)
            rev_sequence = sequence[1] + sequence[0]
            template_path = self._find_template(rev_sequence, lw)
            if template_path:
                rmsd, num_atoms = self.align_to_template(
                    target_res2, target_res1, template_path
                )
                results.append(AlignmentResult(
                    lw_class=lw,
                    sequence=rev_sequence + "(rev)",
                    rmsd=rmsd,
                    num_atoms_aligned=num_atoms,
                    template_path=template_path,
                ))

        # Sort by score (best first) - considers both RMSD and atom count
        results.sort(key=lambda r: r.score)

        if not results:
            return ClassificationResult(
                target_res1=target_res1.res_id,
                target_res2=target_res2.res_id,
                sequence=sequence,
                best_lw="unknown",
                best_rmsd=float('inf'),
            )

        best = results[0]
        second = results[1] if len(results) > 1 else None

        return ClassificationResult(
            target_res1=target_res1.res_id,
            target_res2=target_res2.res_id,
            sequence=sequence,
            best_lw=best.lw_class,
            best_rmsd=best.rmsd,
            second_lw=second.lw_class if second else None,
            second_rmsd=second.rmsd if second else None,
            all_results=results,
        )


# Backward compatibility: re-export parse_pdb as parse_pdb_residues
def parse_pdb_residues(pdb_path: Path) -> Dict[str, Residue]:
    """Parse PDB file and return residues keyed by res_id.

    This is a backward-compatible wrapper around core.pdb_parser.parse_pdb.
    """
    return parse_pdb(pdb_path)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Classify base pairs using template alignment"
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

    residues = parse_pdb_residues(pdb_path)

    if args.res1 not in residues:
        print(f"Residue {args.res1} not found in PDB")
        return
    if args.res2 not in residues:
        print(f"Residue {args.res2} not found in PDB")
        return

    res1 = residues[args.res1]
    res2 = residues[args.res2]

    # Classify
    aligner = TemplateAligner(args.idealized_dir, args.exemplar_dir)
    result = aligner.classify_pair(res1, res2)

    print(f"\nClassification for {result.target_res1} - {result.target_res2} ({result.sequence}):")
    print(f"\n  Best match: {result.best_lw} (RMSD = {result.best_rmsd:.3f} Å)")
    if result.second_lw:
        print(f"  Second best: {result.second_lw} (RMSD = {result.second_rmsd:.3f} Å)")
    print(f"  Confidence: {result.confidence:.2f}")

    print(f"\n  All results (sorted by score):")
    for r in result.all_results[:10]:
        print(f"    {r.lw_class:6s} {r.sequence:8s} RMSD={r.rmsd:.3f}Å score={r.score:.3f} ({r.num_atoms_aligned} atoms)")


if __name__ == "__main__":
    main()

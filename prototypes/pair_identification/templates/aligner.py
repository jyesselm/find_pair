"""Template alignment using Kabsch algorithm."""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import numpy as np

from prototypes.pair_identification.core.residue import Residue
from prototypes.pair_identification.core.alignment import kabsch_align, compute_rmsd
from prototypes.pair_identification.templates.loader import TemplateLoader


@dataclass
class AlignmentResult:
    """Result of aligning a target pair to a template.

    Attributes:
        lw_class: Leontis-Westhof class of template.
        sequence: Template sequence (may include "(rev)" suffix).
        rmsd: RMSD after optimal superposition (Angstroms).
        num_atoms_aligned: Number of atoms used in alignment.
        template_path: Path to template PDB file.
    """

    lw_class: str
    sequence: str
    rmsd: float
    num_atoms_aligned: int
    template_path: Optional[Path] = None

    @property
    def score(self) -> float:
        """Compute alignment score (lower is better).

        Penalizes templates with few atoms, since RMSD alone can be
        misleading when computed over very few atoms.

        Returns:
            Combined score incorporating RMSD and atom count penalty.
        """
        min_atoms = 10
        if self.num_atoms_aligned < min_atoms:
            atom_penalty = (min_atoms - self.num_atoms_aligned) * 0.5
            return self.rmsd + atom_penalty
        return self.rmsd


@dataclass
class ClassificationResult:
    """Result of classifying a pair by trying all LW templates.

    Attributes:
        target_res1: First residue ID.
        target_res2: Second residue ID.
        sequence: Target sequence.
        best_lw: Best matching LW class.
        best_rmsd: RMSD of best match.
        second_lw: Second best LW class.
        second_rmsd: RMSD of second best match.
        all_results: All alignment results sorted by score.
    """

    target_res1: str
    target_res2: str
    sequence: str
    best_lw: str
    best_rmsd: float
    second_lw: Optional[str] = None
    second_rmsd: Optional[float] = None
    all_results: List[AlignmentResult] = None

    def __post_init__(self):
        """Initialize all_results if None."""
        if self.all_results is None:
            self.all_results = []

    @property
    def confidence(self) -> float:
        """Compute confidence based on gap between best and second best.

        Returns:
            Confidence value in [0, 1], where 1 is highest confidence.
        """
        if self.second_rmsd is None:
            return 1.0
        gap = self.second_rmsd - self.best_rmsd
        return min(1.0, gap / 0.5)


class TemplateAligner:
    """Align target pairs to templates and compute RMSD.

    Uses Kabsch algorithm to find optimal superposition of target pair
    onto each template, then computes RMSD over ring atoms.
    """

    RING_ATOMS = ["C2", "C4", "C5", "C6", "N1", "N3", "N7", "C8", "N9"]

    DEFAULT_LW_CLASSES = [
        "cWW", "tWW", "cWH", "tWH", "cWS", "tWS",
        "cHH", "tHH", "cHS", "tHS", "cSS", "tSS"
    ]

    def __init__(
        self,
        idealized_dir: Path = Path("basepair-idealized"),
        exemplar_dir: Path = Path("basepair-exemplars"),
    ):
        """Initialize template aligner.

        Args:
            idealized_dir: Directory with idealized templates.
            exemplar_dir: Directory with exemplar templates.
        """
        self.loader = TemplateLoader(idealized_dir, exemplar_dir)

    def align_to_template(
        self,
        target_res1: Residue,
        target_res2: Residue,
        template_path: Path,
    ) -> tuple[float, int]:
        """Align target pair to template and compute RMSD.

        Args:
            target_res1: First residue of target pair.
            target_res2: Second residue of target pair.
            template_path: Path to template PDB.

        Returns:
            Tuple of (rmsd, num_atoms) where rmsd is the RMSD after
            optimal superposition and num_atoms is the number of atoms aligned.
        """
        template_res1, template_res2 = self.loader.load_template(template_path)

        template_points, target_points = self._collect_matching_atoms(
            template_res1, template_res2, target_res1, target_res2
        )

        if len(template_points) < 4:
            return float('inf'), 0

        template_points = np.array(template_points)
        target_points = np.array(target_points)

        R, centroid_target, centroid_template = kabsch_align(
            template_points, target_points
        )

        aligned_template = (
            (template_points - centroid_template) @ R + centroid_target
        )

        rmsd = compute_rmsd(aligned_template, target_points)

        return rmsd, len(template_points)

    def _collect_matching_atoms(
        self, template_res1, template_res2, target_res1, target_res2
    ) -> tuple[list, list]:
        """Collect matching ring atoms from template and target."""
        template_points = []
        target_points = []

        for atom_name in self.RING_ATOMS:
            if atom_name in template_res1 and atom_name in target_res1.atoms:
                template_points.append(template_res1[atom_name])
                target_points.append(target_res1.atoms[atom_name].coords)
            if atom_name in template_res2 and atom_name in target_res2.atoms:
                template_points.append(template_res2[atom_name])
                target_points.append(target_res2.atoms[atom_name].coords)

        return template_points, target_points

    def classify_pair(
        self,
        target_res1: Residue,
        target_res2: Residue,
        lw_classes: Optional[List[str]] = None,
    ) -> ClassificationResult:
        """Classify a pair by trying all LW class templates.

        Tries both forward and reverse sequence orientations for each LW class.

        Args:
            target_res1: First residue.
            target_res2: Second residue.
            lw_classes: LW classes to try (default: all 12).

        Returns:
            ClassificationResult with best match and confidence.
        """
        if lw_classes is None:
            lw_classes = self.DEFAULT_LW_CLASSES

        sequence = target_res1.base_type + target_res2.base_type

        results = []

        for lw in lw_classes:
            self._try_forward(
                lw, sequence, target_res1, target_res2, results
            )
            self._try_reverse(
                lw, sequence, target_res1, target_res2, results
            )

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

    def _try_forward(
        self,
        lw: str,
        sequence: str,
        res1: Residue,
        res2: Residue,
        results: List[AlignmentResult],
    ) -> None:
        """Try forward sequence orientation."""
        template_path = self.loader.find_template(lw, sequence)
        if not template_path:
            return

        rmsd, num_atoms = self.align_to_template(res1, res2, template_path)
        results.append(
            AlignmentResult(
                lw_class=lw,
                sequence=sequence,
                rmsd=rmsd,
                num_atoms_aligned=num_atoms,
                template_path=template_path,
            )
        )

    def _try_reverse(
        self,
        lw: str,
        sequence: str,
        res1: Residue,
        res2: Residue,
        results: List[AlignmentResult],
    ) -> None:
        """Try reversed sequence orientation."""
        rev_sequence = sequence[1] + sequence[0]
        template_path = self.loader.find_template(lw, rev_sequence)
        if not template_path:
            return

        rmsd, num_atoms = self.align_to_template(res2, res1, template_path)
        results.append(
            AlignmentResult(
                lw_class=lw,
                sequence=rev_sequence + "(rev)",
                rmsd=rmsd,
                num_atoms_aligned=num_atoms,
                template_path=template_path,
            )
        )

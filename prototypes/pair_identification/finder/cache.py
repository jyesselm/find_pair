"""Pair candidate caching for efficient pair finding.

This module provides fast access to pre-validated pair candidates by building
a spatial index (KDTree) and caching geometric validation results. Supports
querying by residue for greedy selection algorithms.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from scipy.spatial import KDTree

from prototypes.pair_identification.frame_loader import FrameLoader, ReferenceFrame
from prototypes.pair_identification.validation import (
    GeometricValidator,
    ValidationResult,
)


@dataclass
class AtomCoords:
    """Coordinates for key atoms of a residue.

    Attributes:
        c1_prime: C1' position (sugar carbon).
        n1_or_n9: N1 (pyrimidines) or N9 (purines) position.
        ring_atoms: Dictionary of ring atom positions.
        hbond_atoms: Dictionary of H-bond donor/acceptor positions.
    """

    c1_prime: Optional[np.ndarray] = None
    n1_or_n9: Optional[np.ndarray] = None
    ring_atoms: Dict[str, np.ndarray] = field(default_factory=dict)
    hbond_atoms: Dict[str, np.ndarray] = field(default_factory=dict)


@dataclass
class CandidateInfo:
    """Cached validation result for a pair candidate.

    Attributes:
        res_id1: First residue ID.
        res_id2: Second residue ID.
        res_name1: First residue single-letter code.
        res_name2: Second residue single-letter code.
        frame1: Reference frame for residue 1.
        frame2: Reference frame for residue 2.
        validation: Geometric validation result.
        quality_score: Overall quality score (0-1, higher is better).
        lw_class: Leontis-Westhof classification (optional).
    """

    res_id1: str
    res_id2: str
    res_name1: str
    res_name2: str
    frame1: ReferenceFrame
    frame2: ReferenceFrame
    validation: ValidationResult
    quality_score: float = 0.0
    lw_class: Optional[str] = None

    @property
    def sequence(self) -> str:
        """Get two-letter sequence code."""
        return self.res_name1 + self.res_name2


class PairCache:
    """Cache for pre-computed pair validation results.

    This class builds a spatial index on residue frame origins, finds all
    candidate pairs within distance, validates them geometrically, and
    provides efficient lookup by residue ID.

    Attributes:
        pdb_id: PDB identifier.
        json_dir: Directory containing JSON output.
        frames: Reference frames by residue ID.
        atoms: Atom coordinates by residue ID.
        candidates: All validated candidates.
        index: Spatial index for residue lookup.
    """

    def __init__(self, pdb_id: str, json_dir: Path | str):
        """Initialize pair cache.

        Args:
            pdb_id: PDB identifier (e.g., "1EHZ").
            json_dir: Directory containing JSON output.
        """
        self.pdb_id = pdb_id
        self.json_dir = Path(json_dir)
        self.frames: Dict[str, ReferenceFrame] = {}
        self.atoms: Dict[str, AtomCoords] = {}
        self.candidates: List[CandidateInfo] = []
        self._candidates_by_res: Dict[str, List[CandidateInfo]] = {}

    def build(
        self, max_distance: float = 15.0, validator: Optional[GeometricValidator] = None
    ) -> None:
        """Pre-compute all valid pairs within distance.

        Args:
            max_distance: Maximum distance between frame origins (Angstroms).
            validator: Geometric validator (default: GeometricValidator()).

        Raises:
            FileNotFoundError: If required JSON files do not exist.
            ValueError: If JSON format is invalid.
        """
        if validator is None:
            validator = GeometricValidator()

        self._load_frames()
        if not self.frames:
            raise ValueError(f"No frames loaded for {self.pdb_id}")

        self._load_atoms()
        self._find_candidates(max_distance, validator)
        self._build_index()

    def get_candidates_for(self, res_id: str) -> List[CandidateInfo]:
        """Get all valid candidates for a residue.

        Args:
            res_id: Residue ID to query.

        Returns:
            List of CandidateInfo where res_id appears.
        """
        return self._candidates_by_res.get(res_id, [])

    def get_valid_candidates(self) -> List[CandidateInfo]:
        """Get all candidates with valid geometry.

        Returns:
            List of CandidateInfo where validation.is_valid is True.
        """
        return [c for c in self.candidates if c.validation.is_valid]

    def _load_frames(self) -> None:
        """Load reference frames from ls_fitting JSON."""
        frame_loader = FrameLoader(self.json_dir)
        self.frames = frame_loader.load_frames(self.pdb_id)

    def _load_atoms(self) -> None:
        """Load atom coordinates from pdb_atoms JSON.

        Raises:
            FileNotFoundError: If pdb_atoms JSON does not exist.
            ValueError: If JSON format is invalid.
        """
        import json

        atoms_file = self.json_dir / "pdb_atoms" / f"{self.pdb_id}.json"
        if not atoms_file.exists():
            raise FileNotFoundError(f"Atoms JSON not found: {atoms_file}")

        with open(atoms_file) as f:
            data = json.load(f)

        if not isinstance(data, list) or len(data) == 0:
            raise ValueError("Expected non-empty list in atoms JSON")
        if "atoms" not in data[0]:
            raise ValueError("Missing 'atoms' key in JSON data")

        residue_atoms: Dict[str, List[dict]] = {}
        for atom in data[0]["atoms"]:
            res_id = atom["res_id"]
            if res_id not in residue_atoms:
                residue_atoms[res_id] = []
            residue_atoms[res_id].append(atom)

        for res_id, atoms in residue_atoms.items():
            self.atoms[res_id] = self._extract_coords(atoms)

    def _extract_coords(self, atoms: List[dict]) -> AtomCoords:
        """Extract atom coordinates from atom list."""
        coords = AtomCoords()
        residue_name = atoms[0]["residue_name"]
        is_purine = residue_name in ["A", "G", "DA", "DG"]

        for atom in atoms:
            atom_name = atom["atom_name"]
            xyz = np.array(atom["xyz"], dtype=np.float64)

            if atom_name == "C1'":
                coords.c1_prime = xyz
            if is_purine and atom_name == "N9":
                coords.n1_or_n9 = xyz
            elif not is_purine and atom_name == "N1":
                coords.n1_or_n9 = xyz
            if atom_name in ["C2", "C4", "C5", "C6", "C8", "N1", "N3", "N7", "N9"]:
                coords.ring_atoms[atom_name] = xyz
            if atom_name in ["N1", "N2", "N3", "N4", "N6", "N7", "O2", "O4", "O6"]:
                coords.hbond_atoms[atom_name] = xyz

        return coords

    def _find_candidates(
        self, max_distance: float, validator: GeometricValidator
    ) -> None:
        """Find and validate all pairs within distance."""
        res_ids = list(self.frames.keys())
        origins = np.array([self.frames[r].origin for r in res_ids])
        kdtree = KDTree(origins)

        for i, res1_id in enumerate(res_ids):
            neighbors = kdtree.query_ball_point(origins[i], r=max_distance)

            for j in neighbors:
                if j <= i:
                    continue

                res2_id = res_ids[j]
                candidate = self._validate_pair(
                    res1_id, res2_id, validator
                )
                if candidate is not None:
                    self.candidates.append(candidate)

    def _validate_pair(
        self, res1_id: str, res2_id: str, validator: GeometricValidator
    ) -> Optional[CandidateInfo]:
        """Validate a single pair candidate."""
        if res1_id not in self.atoms or res2_id not in self.atoms:
            return None

        n1n9_pos1 = self.atoms[res1_id].n1_or_n9
        n1n9_pos2 = self.atoms[res2_id].n1_or_n9
        if n1n9_pos1 is None or n1n9_pos2 is None:
            return None

        frame1 = self.frames[res1_id]
        frame2 = self.frames[res2_id]
        validation = validator.validate(frame1, frame2, n1n9_pos1, n1n9_pos2)

        res1_name = self._extract_residue_name(res1_id)
        res2_name = self._extract_residue_name(res2_id)

        return CandidateInfo(
            res_id1=res1_id,
            res_id2=res2_id,
            res_name1=res1_name,
            res_name2=res2_name,
            frame1=frame1,
            frame2=frame2,
            validation=validation,
        )

    def _extract_residue_name(self, res_id: str) -> str:
        """Extract single-letter residue name from res_id.

        Args:
            res_id: Residue ID (format: "chain-name-seq[ins]").

        Returns:
            Single-letter residue name (e.g., "G", "C").
        """
        parts = res_id.split("-")
        if len(parts) < 2:
            raise ValueError(f"Invalid res_id format: {res_id}")

        name = parts[1]
        if name.startswith("D") and len(name) == 2:
            return name[1]
        return name

    def _build_index(self) -> None:
        """Build lookup index for efficient querying by residue."""
        self._candidates_by_res.clear()

        for candidate in self.candidates:
            for res_id in [candidate.res_id1, candidate.res_id2]:
                if res_id not in self._candidates_by_res:
                    self._candidates_by_res[res_id] = []
                self._candidates_by_res[res_id].append(candidate)

    def get_residue_ids(self) -> Set[str]:
        """Get all residue IDs that have candidates.

        Returns:
            Set of residue IDs.
        """
        return set(self._candidates_by_res.keys())

"""Canonical Residue and Atom data classes for pair identification."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional
import numpy as np

from .constants import (
    PURINE_RING_ATOMS,
    PYRIMIDINE_RING_ATOMS,
    GLYCOSIDIC_N,
    BASE_ATOMS,
    PURINES,
    PYRIMIDINES,
    MODIFIED_BASE_MAP,
    DNA_PREFIX_MAP,
)


@dataclass
class Atom:
    """Represents a single atom with coordinates."""
    name: str
    coords: np.ndarray
    element: str = ""

    def distance_to(self, other: "Atom") -> float:
        """Compute distance to another atom."""
        return float(np.linalg.norm(self.coords - other.coords))

    def __post_init__(self):
        """Ensure coords is numpy array."""
        if not isinstance(self.coords, np.ndarray):
            self.coords = np.array(self.coords, dtype=np.float64)


@dataclass
class Residue:
    """Represents a nucleotide residue with atoms."""
    res_id: str
    base_type: str
    atoms: Dict[str, Atom] = field(default_factory=dict)

    @property
    def chain(self) -> str:
        """Extract chain from res_id."""
        parts = self.res_id.split("-")
        return parts[0] if parts else ""

    @property
    def seq_num(self) -> str:
        """Extract sequence number from res_id."""
        parts = self.res_id.split("-")
        return parts[-1] if len(parts) >= 3 else ""

    def get_glycosidic_n(self) -> Optional[Atom]:
        """Get the glycosidic nitrogen (N1 for pyrimidines, N9 for purines)."""
        n_name = GLYCOSIDIC_N.get(self.base_type.upper())
        return self.atoms.get(n_name) if n_name else None

    def get_ring_atoms(self) -> Dict[str, np.ndarray]:
        """Get ring atom coordinates for alignment."""
        if self.base_type.upper() in PURINES:
            ring_names = PURINE_RING_ATOMS
        else:
            ring_names = PYRIMIDINE_RING_ATOMS

        return {name: self.atoms[name].coords
                for name in ring_names if name in self.atoms}

    def get_base_atoms(self) -> Dict[str, Atom]:
        """Get only base atoms (exclude sugar/phosphate)."""
        base_atom_names = BASE_ATOMS.get(self.base_type.upper(), set())
        return {name: atom for name, atom in self.atoms.items()
                if name in base_atom_names}

    def get_coords_array(self, atom_names: Optional[List[str]] = None) -> np.ndarray:
        """Get coordinates as Nx3 array for specified atoms.

        Args:
            atom_names: List of atom names to include. If None, use all atoms.

        Returns:
            Nx3 array of coordinates.
        """
        if atom_names is None:
            atom_names = list(self.atoms.keys())

        coords = []
        for name in atom_names:
            if name in self.atoms:
                coords.append(self.atoms[name].coords)

        return np.array(coords) if coords else np.empty((0, 3))

    def add_atom(self, name: str, coords: np.ndarray, element: str = ""):
        """Add an atom to the residue."""
        self.atoms[name] = Atom(name=name, coords=coords, element=element)

    def has_atom(self, name: str) -> bool:
        """Check if residue has a specific atom."""
        return name in self.atoms

    def __repr__(self):
        return f"Residue({self.res_id}, {self.base_type}, {len(self.atoms)} atoms)"


def is_purine(base_type: str) -> bool:
    """Check if base type is a purine (A, G, I)."""
    return base_type.upper() in PURINES


def is_pyrimidine(base_type: str) -> bool:
    """Check if base type is a pyrimidine (C, U, T)."""
    return base_type.upper() in PYRIMIDINES


def normalize_base_type(residue_name: str) -> str:
    """Normalize residue name to single letter base type.

    Handles:
    - Standard: A, G, C, U, T, I, P
    - DNA prefix: DA, DG, DC, DT -> A, G, C, T
    - Modified bases: 2MG -> G, PSU -> U, etc.

    Args:
        residue_name: Raw residue name from PDB.

    Returns:
        Normalized single-letter base type.
    """
    name = residue_name.upper().strip()

    # Standard single letters (RNA, DNA, modified)
    standard = PURINES | PYRIMIDINES | {"I", "P"}
    if name in standard:
        return name

    # DNA prefixes
    if name in DNA_PREFIX_MAP:
        return DNA_PREFIX_MAP[name]

    # Common modified bases -> parent
    if name in MODIFIED_BASE_MAP:
        return MODIFIED_BASE_MAP[name]

    # Try first letter as fallback
    if name and name[0] in standard:
        return name[0]

    return name

"""Geometry utilities for H-bond slot prediction.

Predicts hydrogen and lone pair positions from molecular geometry.
"""

from dataclasses import dataclass, field
from typing import Dict, List, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    pass


def normalize(v: np.ndarray) -> np.ndarray:
    """Normalize a vector to unit length.

    Args:
        v: Input vector.

    Returns:
        Unit vector in same direction as v.
    """
    norm = np.linalg.norm(v)
    if norm < 1e-10:
        return v
    return v / norm


def angle_between(v1: np.ndarray, v2: np.ndarray) -> float:
    """Compute angle in degrees between two vectors.

    Args:
        v1: First vector.
        v2: Second vector.

    Returns:
        Angle in degrees [0, 180].
    """
    v1_u = normalize(v1)
    v2_u = normalize(v2)
    dot = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    return np.degrees(np.arccos(dot))


def rotate_vector(v: np.ndarray, axis: np.ndarray, angle_deg: float) -> np.ndarray:
    """Rotate vector around axis by angle using Rodrigues formula.

    Args:
        v: Vector to rotate.
        axis: Rotation axis (will be normalized).
        angle_deg: Rotation angle in degrees.

    Returns:
        Rotated vector.
    """
    angle_rad = np.radians(angle_deg)
    axis = normalize(axis)
    cos_a = np.cos(angle_rad)
    sin_a = np.sin(angle_rad)

    term1 = v * cos_a
    term2 = np.cross(axis, v) * sin_a
    term3 = axis * np.dot(axis, v) * (1 - cos_a)
    return term1 + term2 + term3


@dataclass
class HSlot:
    """Hydrogen donor slot with direction and bond tracking.

    Supports bifurcation: one H can donate to two acceptors if
    angularly separated.
    """

    direction: np.ndarray
    max_bonds: int = 2
    bond_directions: List[np.ndarray] = field(default_factory=list)

    def is_available(self) -> bool:
        """Check if slot can accept another bond."""
        return len(self.bond_directions) < self.max_bonds

    def can_add_bond(self, new_direction: np.ndarray, min_angle: float = 45.0) -> bool:
        """Check if a new bond can be added (bifurcation check).

        Args:
            new_direction: Direction vector to the new acceptor.
            min_angle: Minimum angle between bonds for bifurcation (degrees).

        Returns:
            True if bond can be added, False otherwise.
        """
        if len(self.bond_directions) == 0:
            return True
        if len(self.bond_directions) >= self.max_bonds:
            return False

        for existing_dir in self.bond_directions:
            angle = angle_between(existing_dir, new_direction)
            if angle < min_angle:
                return False
        return True

    def add_bond(self, direction: np.ndarray):
        """Record a bond using this slot.

        Args:
            direction: Direction vector to the acceptor.
        """
        self.bond_directions.append(normalize(direction))


@dataclass
class LPSlot:
    """Lone pair acceptor slot with direction and bond tracking.

    Supports bifurcation: one LP can accept from two donors if
    angularly separated.
    """

    direction: np.ndarray
    max_bonds: int = 2
    bond_directions: List[np.ndarray] = field(default_factory=list)

    def is_available(self) -> bool:
        """Check if slot can accept another bond."""
        return len(self.bond_directions) < self.max_bonds

    def can_add_bond(self, new_direction: np.ndarray, min_angle: float = 45.0) -> bool:
        """Check if a new bond can be added (bifurcation check).

        Args:
            new_direction: Direction vector from the new donor.
            min_angle: Minimum angle between bonds for bifurcation (degrees).

        Returns:
            True if bond can be added, False otherwise.
        """
        if len(self.bond_directions) == 0:
            return True
        if len(self.bond_directions) >= self.max_bonds:
            return False

        for existing_dir in self.bond_directions:
            angle = angle_between(existing_dir, new_direction)
            if angle < min_angle:
                return False
        return True

    def add_bond(self, direction: np.ndarray):
        """Record a bond using this slot.

        Args:
            direction: Direction vector from the donor.
        """
        self.bond_directions.append(normalize(direction))


# Base connectivity for standard nucleotides
# Maps (base_type, atom_name) -> list of bonded atom names
BASE_CONNECTIVITY: Dict[tuple, List[str]] = {
    # Adenine
    ("A", "N6"): ["C6"],
    ("A", "N1"): ["C2", "C6"],
    ("A", "N3"): ["C2", "C4"],
    ("A", "N7"): ["C5", "C8"],
    ("A", "C2"): ["N1", "N3"],
    ("A", "C8"): ["N7", "N9"],
    # Guanine
    ("G", "N1"): ["C2", "C6"],
    ("G", "N2"): ["C2"],
    ("G", "O6"): ["C6"],
    ("G", "N3"): ["C2", "C4"],
    ("G", "N7"): ["C5", "C8"],
    ("G", "C8"): ["N7", "N9"],
    # Cytosine
    ("C", "N4"): ["C4"],
    ("C", "N3"): ["C2", "C4"],
    ("C", "O2"): ["C2"],
    ("C", "C5"): ["C4", "C6"],
    ("C", "C6"): ["C5", "N1"],
    # Uracil/Thymine
    ("U", "N3"): ["C2", "C4"],
    ("U", "O2"): ["C2"],
    ("U", "O4"): ["C4"],
    ("U", "C5"): ["C4", "C6"],
    ("U", "C6"): ["C5", "N1"],
    ("T", "N3"): ["C2", "C4"],
    ("T", "O2"): ["C2"],
    ("T", "O4"): ["C4"],
    ("T", "C6"): ["C5", "N1"],
    # Ribose atoms (common to all nucleotides)
    ("A", "O2'"): ["C2'"],
    ("G", "O2'"): ["C2'"],
    ("C", "O2'"): ["C2'"],
    ("U", "O2'"): ["C2'"],
    ("T", "O2'"): ["C2'"],
    ("A", "O4'"): ["C1'", "C4'"],
    ("G", "O4'"): ["C1'", "C4'"],
    ("C", "O4'"): ["C1'", "C4'"],
    ("U", "O4'"): ["C1'", "C4'"],
    ("T", "O4'"): ["C1'", "C4'"],
    ("A", "O3'"): ["C3'"],
    ("G", "O3'"): ["C3'"],
    ("C", "O3'"): ["C3'"],
    ("U", "O3'"): ["C3'"],
    ("T", "O3'"): ["C3'"],
    ("A", "O5'"): ["C5'"],
    ("G", "O5'"): ["C5'"],
    ("C", "O5'"): ["C5'"],
    ("U", "O5'"): ["C5'"],
    ("T", "O5'"): ["C5'"],
    # Pseudouridine (P)
    ("P", "N1"): ["C2", "C6"],
    ("P", "N3"): ["C2", "C4"],
    ("P", "O2"): ["C2"],
    ("P", "O4"): ["C4"],
    ("P", "O2'"): ["C2'"],
    ("P", "O4'"): ["C1'", "C4'"],
    ("P", "O3'"): ["C3'"],
    ("P", "O5'"): ["C5'"],
    # Inosine (I)
    ("I", "N1"): ["C2", "C6"],
    ("I", "O6"): ["C6"],
    ("I", "N3"): ["C2", "C4"],
    ("I", "N7"): ["C5", "C8"],
    ("I", "O2'"): ["C2'"],
    ("I", "O4'"): ["C1'", "C4'"],
    ("I", "O3'"): ["C3'"],
    ("I", "O5'"): ["C5'"],
    # DNA bases (deoxyribose - no O2')
    ("DA", "N6"): ["C6"],
    ("DA", "N1"): ["C2", "C6"],
    ("DA", "N3"): ["C2", "C4"],
    ("DA", "N7"): ["C5", "C8"],
    ("DA", "O4'"): ["C1'", "C4'"],
    ("DA", "O3'"): ["C3'"],
    ("DA", "O5'"): ["C5'"],
    ("DG", "N1"): ["C2", "C6"],
    ("DG", "N2"): ["C2"],
    ("DG", "O6"): ["C6"],
    ("DG", "N3"): ["C2", "C4"],
    ("DG", "N7"): ["C5", "C8"],
    ("DG", "O4'"): ["C1'", "C4'"],
    ("DG", "O3'"): ["C3'"],
    ("DG", "O5'"): ["C5'"],
    ("DC", "N4"): ["C4"],
    ("DC", "N3"): ["C2", "C4"],
    ("DC", "O2"): ["C2"],
    ("DC", "O4'"): ["C1'", "C4'"],
    ("DC", "O3'"): ["C3'"],
    ("DC", "O5'"): ["C5'"],
    ("DT", "N3"): ["C2", "C4"],
    ("DT", "O2"): ["C2"],
    ("DT", "O4"): ["C4"],
    ("DT", "O4'"): ["C1'", "C4'"],
    ("DT", "O3'"): ["C3'"],
    ("DT", "O5'"): ["C5'"],
}


# Donor capacity: how many H atoms can donate
DONOR_CAPACITY: Dict[tuple, int] = {
    # NH2 amino groups - 2 hydrogens
    ("A", "N6"): 2,
    ("C", "N4"): 2,
    ("G", "N2"): 2,
    # Imino NH - 1 hydrogen
    ("G", "N1"): 1,
    ("U", "N3"): 1,
    ("T", "N3"): 1,
    # Ribose O2' hydroxyl - can donate 1
    ("A", "O2'"): 1,
    ("G", "O2'"): 1,
    ("C", "O2'"): 1,
    ("U", "O2'"): 1,
    ("T", "O2'"): 1,
    ("P", "O2'"): 1,
    ("I", "O2'"): 1,
    # Pseudouridine (P)
    ("P", "N1"): 1,
    ("P", "N3"): 1,
    # Inosine (I)
    ("I", "N1"): 1,
    # DNA bases
    ("DA", "N6"): 2,
    ("DG", "N1"): 1,
    ("DG", "N2"): 2,
    ("DC", "N4"): 2,
    ("DT", "N3"): 1,
}


# Acceptor capacity: how many lone pairs can accept
ACCEPTOR_CAPACITY: Dict[tuple, int] = {
    # sp2 carbonyl oxygens - 2 lone pairs
    ("G", "O6"): 2,
    ("U", "O2"): 2,
    ("U", "O4"): 2,
    ("C", "O2"): 2,
    ("T", "O2"): 2,
    ("T", "O4"): 2,
    # sp2 ring nitrogens - 1 lone pair
    ("A", "N1"): 1,
    ("A", "N3"): 1,
    ("A", "N7"): 1,
    ("G", "N3"): 1,
    ("G", "N7"): 1,
    ("C", "N3"): 1,
    # Ribose atoms
    ("A", "O2'"): 2,
    ("G", "O2'"): 2,
    ("C", "O2'"): 2,
    ("U", "O2'"): 2,
    ("T", "O2'"): 2,
    ("A", "O4'"): 1,
    ("G", "O4'"): 1,
    ("C", "O4'"): 1,
    ("U", "O4'"): 1,
    ("T", "O4'"): 1,
    # Phosphate oxygens - 3 lone pairs
    ("A", "OP1"): 3,
    ("G", "OP1"): 3,
    ("C", "OP1"): 3,
    ("U", "OP1"): 3,
    ("T", "OP1"): 3,
    ("A", "OP2"): 3,
    ("G", "OP2"): 3,
    ("C", "OP2"): 3,
    ("U", "OP2"): 3,
    ("T", "OP2"): 3,
    ("A", "O1P"): 3,
    ("G", "O1P"): 3,
    ("C", "O1P"): 3,
    ("U", "O1P"): 3,
    ("A", "O2P"): 3,
    ("G", "O2P"): 3,
    ("C", "O2P"): 3,
    ("U", "O2P"): 3,
    # Pseudouridine (P)
    ("P", "O2"): 2,
    ("P", "O4"): 2,
    ("P", "O2'"): 2,
    ("P", "O4'"): 1,
    ("P", "OP1"): 3,
    ("P", "OP2"): 3,
    # Inosine (I)
    ("I", "O6"): 2,
    ("I", "N3"): 1,
    ("I", "N7"): 1,
    ("I", "O2'"): 2,
    ("I", "O4'"): 1,
    ("I", "OP1"): 3,
    ("I", "OP2"): 3,
    # DNA bases
    ("DA", "N1"): 1,
    ("DA", "N3"): 1,
    ("DA", "N7"): 1,
    ("DA", "O4'"): 1,
    ("DG", "O6"): 2,
    ("DG", "N3"): 1,
    ("DG", "N7"): 1,
    ("DG", "O4'"): 1,
    ("DC", "O2"): 2,
    ("DC", "N3"): 1,
    ("DC", "O4'"): 1,
    ("DT", "O2"): 2,
    ("DT", "O4"): 2,
    ("DT", "O4'"): 1,
}


def compute_base_normal(atoms: Dict[str, np.ndarray]) -> np.ndarray:
    """Compute the normal vector to the base plane.

    Uses ring atoms (C2, C4, C6, N1, N3) which are present in all bases.

    Args:
        atoms: Dict mapping atom names to positions.

    Returns:
        Unit normal vector to the base plane.
    """
    ring_atoms = ["C2", "C4", "C6", "N1", "N3"]
    positions = []
    for name in ring_atoms:
        if name in atoms:
            positions.append(atoms[name])

    if len(positions) < 3:
        return np.array([0.0, 0.0, 1.0])

    v1 = positions[1] - positions[0]
    v2 = positions[2] - positions[0]
    normal = np.cross(v1, v2)
    return normalize(normal)


def predict_h_slots(
    base_type: str,
    atom_name: str,
    atoms: Dict[str, np.ndarray],
    base_normal: np.ndarray,
) -> List[HSlot]:
    """Predict hydrogen slot positions for a donor atom.

    Args:
        base_type: Single letter base code (A, G, C, U, T, etc.).
        atom_name: Name of the donor atom (e.g., 'N6', 'N1').
        atoms: Dict mapping atom names to positions.
        base_normal: Normal vector to the base plane.

    Returns:
        List of HSlot objects with predicted H directions.
    """
    key = (base_type.upper(), atom_name.strip())

    if key not in DONOR_CAPACITY:
        return []

    capacity = DONOR_CAPACITY[key]
    connectivity = BASE_CONNECTIVITY.get(key, [])

    if atom_name.strip() not in atoms:
        return []

    donor_pos = atoms[atom_name.strip()]

    antecedent_positions = []
    for ant_name in connectivity:
        if ant_name in atoms:
            antecedent_positions.append(atoms[ant_name])

    if not antecedent_positions:
        return []

    slots = []

    if capacity == 2 and len(antecedent_positions) == 1:
        # sp2 NH2: Two H atoms at 120° from C-N bond
        ant_to_donor = normalize(donor_pos - antecedent_positions[0])
        h1_dir = rotate_vector(ant_to_donor, base_normal, 120.0)
        h2_dir = rotate_vector(ant_to_donor, base_normal, -120.0)
        slots.append(HSlot(direction=h1_dir))
        slots.append(HSlot(direction=h2_dir))

    elif capacity == 1 and len(antecedent_positions) == 2:
        # sp2 imino NH: H points away from ring
        avg_ant = (antecedent_positions[0] + antecedent_positions[1]) / 2
        h_dir = normalize(donor_pos - avg_ant)
        slots.append(HSlot(direction=h_dir))

    elif capacity == 1 and len(antecedent_positions) == 1:
        # Simple case: H opposite to single antecedent
        h_dir = normalize(donor_pos - antecedent_positions[0])
        slots.append(HSlot(direction=h_dir))

    return slots


def _get_phosphate_lp_slots() -> List[LPSlot]:
    """Get LP slots for phosphate oxygens (isotropic model).

    Returns:
        Three orthogonal LPSlot objects.
    """
    return [
        LPSlot(direction=np.array([1.0, 0.0, 0.0])),
        LPSlot(direction=np.array([0.0, 1.0, 0.0])),
        LPSlot(direction=np.array([0.0, 0.0, 1.0])),
    ]


def _get_ribose_lp_slots(base_normal: np.ndarray) -> List[LPSlot]:
    """Get LP slots for ribose oxygens (sp3 with flexible geometry).

    Args:
        base_normal: Normal vector to base plane.

    Returns:
        Two perpendicular LPSlot objects.
    """
    perp1 = np.cross(base_normal, np.array([1.0, 0.0, 0.0]))
    if np.linalg.norm(perp1) < 0.1:
        perp1 = np.cross(base_normal, np.array([0.0, 1.0, 0.0]))
    perp1 = normalize(perp1)
    perp2 = np.cross(base_normal, perp1)
    return [LPSlot(direction=perp1), LPSlot(direction=perp2)]


def _get_fallback_lp_slots(capacity: int) -> List[LPSlot]:
    """Get fallback LP slots (isotropic model).

    Args:
        capacity: Number of lone pairs.

    Returns:
        List of LPSlot objects with isotropic directions.
    """
    slots = [LPSlot(direction=np.array([1.0, 0.0, 0.0]))]
    if capacity >= 2:
        slots.append(LPSlot(direction=np.array([0.0, 1.0, 0.0])))
    return slots


def _get_geometry_based_lp_slots(
    capacity: int,
    acceptor_pos: np.ndarray,
    antecedent_positions: List[np.ndarray],
    base_normal: np.ndarray,
) -> List[LPSlot]:
    """Get LP slots based on molecular geometry.

    Args:
        capacity: Number of lone pairs.
        acceptor_pos: Position of acceptor atom.
        antecedent_positions: Positions of bonded atoms.
        base_normal: Normal vector to base plane.

    Returns:
        List of LPSlot objects.
    """
    slots = []

    if capacity == 2 and len(antecedent_positions) == 1:
        # sp2 carbonyl: Two lone pairs at 120° from C=O
        ant_to_acc = normalize(acceptor_pos - antecedent_positions[0])
        lp1_dir = rotate_vector(ant_to_acc, base_normal, 120.0)
        lp2_dir = rotate_vector(ant_to_acc, base_normal, -120.0)
        slots.append(LPSlot(direction=lp1_dir))
        slots.append(LPSlot(direction=lp2_dir))

    elif capacity == 1 and len(antecedent_positions) == 2:
        # sp2 ring nitrogen: One lone pair pointing out of ring
        avg_ant = (antecedent_positions[0] + antecedent_positions[1]) / 2
        lp_dir = normalize(acceptor_pos - avg_ant)
        slots.append(LPSlot(direction=lp_dir, max_bonds=1))

    elif capacity == 1 and len(antecedent_positions) == 1:
        # Fallback
        lp_dir = normalize(acceptor_pos - antecedent_positions[0])
        slots.append(LPSlot(direction=lp_dir))

    return slots


def predict_lp_slots(
    base_type: str,
    atom_name: str,
    atoms: Dict[str, np.ndarray],
    base_normal: np.ndarray,
) -> List[LPSlot]:
    """Predict lone pair slot directions for an acceptor atom.

    Args:
        base_type: Single letter base code.
        atom_name: Name of acceptor atom (e.g., 'O6', 'N1').
        atoms: Dict mapping atom names to positions.
        base_normal: Normal vector to base plane.

    Returns:
        List of LPSlot objects with predicted LP directions.
    """
    key = (base_type.upper(), atom_name.strip())

    if key not in ACCEPTOR_CAPACITY:
        return []

    capacity = ACCEPTOR_CAPACITY[key]
    connectivity = BASE_CONNECTIVITY.get(key, [])

    if atom_name.strip() not in atoms:
        return []

    acceptor_pos = atoms[atom_name.strip()]

    # Handle special cases
    if atom_name.strip() in ("OP1", "OP2", "O1P", "O2P"):
        return _get_phosphate_lp_slots()

    if atom_name.strip() in ("O2'", "O4'", "O3'", "O5'"):
        return _get_ribose_lp_slots(base_normal)

    # Get antecedent positions
    antecedent_positions = []
    for ant_name in connectivity:
        if ant_name in atoms:
            antecedent_positions.append(atoms[ant_name])

    if not antecedent_positions:
        return _get_fallback_lp_slots(capacity)

    return _get_geometry_based_lp_slots(
        capacity, acceptor_pos, antecedent_positions, base_normal
    )

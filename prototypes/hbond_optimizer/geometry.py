"""
Geometry utilities for H-bond prediction.

Provides functions to compute hydrogen positions and lone pair directions
from known molecular geometry.
"""

import json
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass, field


# Load CIF-derived donor/acceptor data for modified residues
_CIF_DONORS: Dict[str, Dict[str, int]] = {}
_CIF_ACCEPTORS: Dict[str, Dict[str, int]] = {}
_CIF_DATA_LOADED = False


def _load_cif_data():
    """Load CIF donor/acceptor data from JSON files."""
    global _CIF_DONORS, _CIF_ACCEPTORS, _CIF_DATA_LOADED
    if _CIF_DATA_LOADED:
        return

    cif_dir = Path(__file__).parent
    donors_path = cif_dir / "cif_donors_chain.json"
    acceptors_path = cif_dir / "cif_acceptors.json"

    if donors_path.exists():
        with open(donors_path) as f:
            _CIF_DONORS = json.load(f)

    if acceptors_path.exists():
        with open(acceptors_path) as f:
            _CIF_ACCEPTORS = json.load(f)

    _CIF_DATA_LOADED = True


def get_donor_capacity(residue_code: str, atom_name: str) -> int:
    """
    Get donor capacity for an atom, checking CIF data first for modified residues.

    Args:
        residue_code: Full residue code (e.g., "1MA", "A", "G")
        atom_name: Atom name (e.g., "N6", "O2'")

    Returns:
        Donor capacity (number of H atoms that can donate)
    """
    _load_cif_data()

    # Check CIF data first (for modified residues)
    if residue_code in _CIF_DONORS:
        cif_donors = _CIF_DONORS[residue_code]
        if atom_name in cif_donors:
            return cif_donors[atom_name]
        # Atom not a donor in this residue
        return 0

    # Fall back to standard base type lookup
    # Map common modified residue prefixes to parent base
    base_type = _get_parent_base_type(residue_code)
    if base_type:
        key = (base_type, atom_name.strip())
        return DONOR_CAPACITY.get(key, 0)

    return 0


def get_acceptor_capacity(residue_code: str, atom_name: str) -> int:
    """
    Get acceptor capacity for an atom, checking CIF data first for modified residues.

    Args:
        residue_code: Full residue code (e.g., "1MA", "A", "G")
        atom_name: Atom name (e.g., "O6", "N1")

    Returns:
        Acceptor capacity (number of lone pairs that can accept)
    """
    _load_cif_data()

    # Check CIF data first (for modified residues)
    if residue_code in _CIF_ACCEPTORS:
        cif_acceptors = _CIF_ACCEPTORS[residue_code]
        if atom_name in cif_acceptors:
            return cif_acceptors[atom_name]
        # Atom not an acceptor in this residue
        return 0

    # Fall back to standard base type lookup
    base_type = _get_parent_base_type(residue_code)
    if base_type:
        key = (base_type, atom_name.strip())
        return ACCEPTOR_CAPACITY.get(key, 0)

    return 0


def _get_parent_base_type(residue_code: str) -> Optional[str]:
    """
    Get the parent base type for a residue code.

    Standard nucleotides return themselves (A, G, C, U, T).
    Modified residues are mapped to their parent base.
    """
    # Standard nucleotides
    if residue_code in ('A', 'G', 'C', 'U', 'T', 'P', 'I'):
        return residue_code

    # DNA variants
    if residue_code in ('DA', 'DG', 'DC', 'DT'):
        return residue_code[1]  # DA -> A, etc.

    # Try to infer from residue code
    # Many modified residues have the parent base as last character
    # e.g., 1MA (1-methyladenosine), 5MC (5-methylcytidine)
    if len(residue_code) >= 2:
        last_char = residue_code[-1].upper()
        if last_char in ('A', 'G', 'C', 'U', 'T'):
            return last_char

    # Check if first char is the base type
    if len(residue_code) >= 1:
        first_char = residue_code[0].upper()
        if first_char in ('A', 'G', 'C', 'U', 'T'):
            return first_char

    return None


def normalize(v: np.ndarray) -> np.ndarray:
    """Normalize a vector to unit length."""
    norm = np.linalg.norm(v)
    if norm < 1e-10:
        return v
    return v / norm


def angle_between(v1: np.ndarray, v2: np.ndarray) -> float:
    """Return angle in degrees between vectors v1 and v2."""
    v1_u = normalize(v1)
    v2_u = normalize(v2)
    dot = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    return np.degrees(np.arccos(dot))


def rotate_vector(v: np.ndarray, axis: np.ndarray, angle_deg: float) -> np.ndarray:
    """Rotate vector v around axis by angle_deg degrees (Rodrigues formula)."""
    angle_rad = np.radians(angle_deg)
    axis = normalize(axis)
    cos_a = np.cos(angle_rad)
    sin_a = np.sin(angle_rad)
    return v * cos_a + np.cross(axis, v) * sin_a + axis * np.dot(axis, v) * (1 - cos_a)


@dataclass
class HSlot:
    """A hydrogen slot on a donor atom.

    Supports bifurcation: one H can donate to two acceptors if angularly separated.
    """
    direction: np.ndarray  # Unit vector pointing from donor toward where H is
    used: bool = False
    # Track directions to acceptors for bifurcation check
    bond_directions: List[np.ndarray] = field(default_factory=list)
    max_bonds: int = 2  # Allow bifurcation (1 H to 2 acceptors)

    def is_available(self) -> bool:
        """Check if slot can accept another bond."""
        return len(self.bond_directions) < self.max_bonds

    def can_add_bond(self, new_direction: np.ndarray, min_angle: float = 60.0) -> bool:
        """Check if a new bond can be added (bifurcation check).

        Args:
            new_direction: Direction vector to the new acceptor
            min_angle: Minimum angle between bonds for bifurcation (degrees)
        """
        if len(self.bond_directions) == 0:
            return True
        if len(self.bond_directions) >= self.max_bonds:
            return False
        # Check angular separation from existing bonds
        for existing_dir in self.bond_directions:
            angle = angle_between(existing_dir, new_direction)
            if angle < min_angle:
                return False  # Too close, would be stacking
        return True

    def add_bond(self, direction: np.ndarray):
        """Record a bond using this slot."""
        self.bond_directions.append(normalize(direction))
        self.used = True


@dataclass
class LPSlot:
    """A lone pair slot on an acceptor atom.

    Supports bifurcation: one LP can accept from two donors if angularly separated.
    """
    direction: np.ndarray  # Unit vector pointing from acceptor in LP direction
    used: bool = False
    # Track directions from donors for bifurcation check
    bond_directions: List[np.ndarray] = field(default_factory=list)
    max_bonds: int = 2  # Allow bifurcation (2 donors to 1 LP)

    def is_available(self) -> bool:
        """Check if slot can accept another bond."""
        return len(self.bond_directions) < self.max_bonds

    def can_add_bond(self, new_direction: np.ndarray, min_angle: float = 60.0) -> bool:
        """Check if a new bond can be added (bifurcation check).

        Args:
            new_direction: Direction vector from the new donor
            min_angle: Minimum angle between bonds for bifurcation (degrees)
        """
        if len(self.bond_directions) == 0:
            return True
        if len(self.bond_directions) >= self.max_bonds:
            return False
        # Check angular separation from existing bonds
        for existing_dir in self.bond_directions:
            angle = angle_between(existing_dir, new_direction)
            if angle < min_angle:
                return False  # Too close, would be stacking
        return True

    def add_bond(self, direction: np.ndarray):
        """Record a bond using this slot."""
        self.bond_directions.append(normalize(direction))
        self.used = True


# Base connectivity for standard nucleotides
# Maps (base_type, atom_name) -> list of bonded atom names
BASE_CONNECTIVITY = {
    # Adenine
    ('A', 'N6'): ['C6'],           # Amino nitrogen, bonded to C6
    ('A', 'N1'): ['C2', 'C6'],     # Ring N, bonded to C2 and C6
    ('A', 'N3'): ['C2', 'C4'],     # Ring N
    ('A', 'N7'): ['C5', 'C8'],     # Ring N
    ('A', 'C2'): ['N1', 'N3'],     # For C-H
    ('A', 'C8'): ['N7', 'N9'],     # For C-H

    # Guanine
    ('G', 'N1'): ['C2', 'C6'],     # Imino N-H
    ('G', 'N2'): ['C2'],           # Amino nitrogen
    ('G', 'O6'): ['C6'],           # Carbonyl oxygen
    ('G', 'N3'): ['C2', 'C4'],     # Ring N acceptor
    ('G', 'N7'): ['C5', 'C8'],     # Ring N acceptor
    ('G', 'C8'): ['N7', 'N9'],     # For C-H

    # Cytosine
    ('C', 'N4'): ['C4'],           # Amino nitrogen
    ('C', 'N3'): ['C2', 'C4'],     # Ring N acceptor
    ('C', 'O2'): ['C2'],           # Carbonyl oxygen
    ('C', 'C5'): ['C4', 'C6'],     # For C-H
    ('C', 'C6'): ['C5', 'N1'],     # For C-H

    # Uracil/Thymine
    ('U', 'N3'): ['C2', 'C4'],     # Imino N-H
    ('U', 'O2'): ['C2'],           # Carbonyl oxygen
    ('U', 'O4'): ['C4'],           # Carbonyl oxygen
    ('U', 'C5'): ['C4', 'C6'],     # For C-H
    ('U', 'C6'): ['C5', 'N1'],     # For C-H

    ('T', 'N3'): ['C2', 'C4'],
    ('T', 'O2'): ['C2'],
    ('T', 'O4'): ['C4'],
    ('T', 'C6'): ['C5', 'N1'],

    # Ribose atoms (common to all nucleotides)
    # O2' hydroxyl bonded to C2'
    ('A', "O2'"): ["C2'"],
    ('G', "O2'"): ["C2'"],
    ('C', "O2'"): ["C2'"],
    ('U', "O2'"): ["C2'"],
    ('T', "O2'"): ["C2'"],
    # O4' ring oxygen bonded to C1' and C4'
    ('A', "O4'"): ["C1'", "C4'"],
    ('G', "O4'"): ["C1'", "C4'"],
    ('C', "O4'"): ["C1'", "C4'"],
    ('U', "O4'"): ["C1'", "C4'"],
    ('T', "O4'"): ["C1'", "C4'"],
    # O3' bonded to C3' (also to phosphate)
    ('A', "O3'"): ["C3'"],
    ('G', "O3'"): ["C3'"],
    ('C', "O3'"): ["C3'"],
    ('U', "O3'"): ["C3'"],
    ('T', "O3'"): ["C3'"],
    # O5' bonded to C5' (also to phosphate)
    ('A', "O5'"): ["C5'"],
    ('G', "O5'"): ["C5'"],
    ('C', "O5'"): ["C5'"],
    ('U', "O5'"): ["C5'"],
    ('T', "O5'"): ["C5'"],

    # Pseudouridine (P) - C-glycosidic bond at C5, not N1
    # Base atoms
    ('P', 'N1'): ['C2', 'C6'],     # N1 is free (not glycosidic), has H
    ('P', 'N3'): ['C2', 'C4'],     # Imino N-H (same as uridine)
    ('P', 'O2'): ['C2'],           # Carbonyl oxygen
    ('P', 'O4'): ['C4'],           # Carbonyl oxygen
    # Ribose atoms
    ('P', "O2'"): ["C2'"],
    ('P', "O4'"): ["C1'", "C4'"],
    ('P', "O3'"): ["C3'"],
    ('P', "O5'"): ["C5'"],

    # Inosine (I) - like guanine but no N2 amino
    # Base atoms
    ('I', 'N1'): ['C2', 'C6'],     # Imino N-H
    ('I', 'O6'): ['C6'],           # Carbonyl oxygen
    ('I', 'N3'): ['C2', 'C4'],     # Ring N acceptor
    ('I', 'N7'): ['C5', 'C8'],     # Ring N acceptor
    # Ribose atoms
    ('I', "O2'"): ["C2'"],
    ('I', "O4'"): ["C1'", "C4'"],
    ('I', "O3'"): ["C3'"],
    ('I', "O5'"): ["C5'"],
}

# Amino acid connectivity for H-bond geometry prediction
# Maps (residue_code, atom_name) -> list of bonded atom names
PROTEIN_CONNECTIVITY = {
    # Backbone (common to all amino acids)
    # N (backbone amide) bonded to CA
    ('ALA', 'N'): ['CA'],
    ('ARG', 'N'): ['CA'],
    ('ASN', 'N'): ['CA'],
    ('ASP', 'N'): ['CA'],
    ('CYS', 'N'): ['CA'],
    ('GLN', 'N'): ['CA'],
    ('GLU', 'N'): ['CA'],
    ('GLY', 'N'): ['CA'],
    ('HIS', 'N'): ['CA'],
    ('ILE', 'N'): ['CA'],
    ('LEU', 'N'): ['CA'],
    ('LYS', 'N'): ['CA'],
    ('MET', 'N'): ['CA'],
    ('PHE', 'N'): ['CA'],
    ('PRO', 'N'): ['CA'],  # Proline has no H on N (imino)
    ('SER', 'N'): ['CA'],
    ('THR', 'N'): ['CA'],
    ('TRP', 'N'): ['CA'],
    ('TYR', 'N'): ['CA'],
    ('VAL', 'N'): ['CA'],

    # O (backbone carbonyl) bonded to C
    ('ALA', 'O'): ['C'],
    ('ARG', 'O'): ['C'],
    ('ASN', 'O'): ['C'],
    ('ASP', 'O'): ['C'],
    ('CYS', 'O'): ['C'],
    ('GLN', 'O'): ['C'],
    ('GLU', 'O'): ['C'],
    ('GLY', 'O'): ['C'],
    ('HIS', 'O'): ['C'],
    ('ILE', 'O'): ['C'],
    ('LEU', 'O'): ['C'],
    ('LYS', 'O'): ['C'],
    ('MET', 'O'): ['C'],
    ('PHE', 'O'): ['C'],
    ('PRO', 'O'): ['C'],
    ('SER', 'O'): ['C'],
    ('THR', 'O'): ['C'],
    ('TRP', 'O'): ['C'],
    ('TYR', 'O'): ['C'],
    ('VAL', 'O'): ['C'],

    # ARG - guanidinium group
    ('ARG', 'NE'): ['CD', 'CZ'],      # NE bonded to CD and CZ
    ('ARG', 'NH1'): ['CZ'],           # NH2 amino
    ('ARG', 'NH2'): ['CZ'],           # NH2 amino

    # ASN - amide side chain
    ('ASN', 'ND2'): ['CG'],           # NH2 amino
    ('ASN', 'OD1'): ['CG'],           # Carbonyl

    # ASP - carboxylate
    ('ASP', 'OD1'): ['CG'],
    ('ASP', 'OD2'): ['CG'],

    # GLN - amide side chain
    ('GLN', 'NE2'): ['CD'],           # NH2 amino
    ('GLN', 'OE1'): ['CD'],           # Carbonyl

    # GLU - carboxylate
    ('GLU', 'OE1'): ['CD'],
    ('GLU', 'OE2'): ['CD'],

    # HIS - imidazole
    ('HIS', 'ND1'): ['CG', 'CE1'],    # Can be donor or acceptor depending on tautomer
    ('HIS', 'NE2'): ['CD2', 'CE1'],   # Can be donor or acceptor depending on tautomer

    # LYS - amino
    ('LYS', 'NZ'): ['CE'],            # NH3+ amino

    # SER - hydroxyl
    ('SER', 'OG'): ['CB'],

    # THR - hydroxyl
    ('THR', 'OG1'): ['CB'],

    # TRP - indole NH
    ('TRP', 'NE1'): ['CD1', 'CE2'],

    # TYR - hydroxyl
    ('TYR', 'OH'): ['CZ'],
}

# Donor capacity: how many H atoms can donate
# Only N and O donors (no C-H)
DONOR_CAPACITY = {
    # NH2 amino groups - 2 hydrogens
    ('A', 'N6'): 2,
    ('C', 'N4'): 2,
    ('G', 'N2'): 2,

    # Imino NH - 1 hydrogen
    ('G', 'N1'): 1,
    ('U', 'N3'): 1,
    ('T', 'N3'): 1,

    # Ribose O2' hydroxyl - can donate 1
    ('A', "O2'"): 1,
    ('G', "O2'"): 1,
    ('C', "O2'"): 1,
    ('U', "O2'"): 1,
    ('T', "O2'"): 1,
    ('P', "O2'"): 1,  # Pseudouridine
    ('I', "O2'"): 1,  # Inosine

    # Ribose O3' hydroxyl (3' end) - can donate 1
    ('A', "O3'"): 1,
    ('G', "O3'"): 1,
    ('C', "O3'"): 1,
    ('U', "O3'"): 1,
    ('T', "O3'"): 1,
    ('P', "O3'"): 1,
    ('I', "O3'"): 1,

    # Ribose O5' hydroxyl (5' end) - can donate 1
    ('A', "O5'"): 1,
    ('G', "O5'"): 1,
    ('C', "O5'"): 1,
    ('U', "O5'"): 1,
    ('T', "O5'"): 1,
    ('P', "O5'"): 1,
    ('I', "O5'"): 1,

    # Pseudouridine (P) - C-glycosidic bond, so N1 is free to donate
    ('P', 'N1'): 1,   # N1 is now a donor (not glycosidic)
    ('P', 'N3'): 1,   # Same as uridine

    # Inosine (I) - like guanine but without N2 amino
    ('I', 'N1'): 1,   # Imino NH
}

# Acceptor capacity: how many lone pairs can accept
ACCEPTOR_CAPACITY = {
    # sp2 carbonyl oxygens - 2 lone pairs
    ('G', 'O6'): 2,
    ('U', 'O2'): 2,
    ('U', 'O4'): 2,
    ('C', 'O2'): 2,
    ('T', 'O2'): 2,
    ('T', 'O4'): 2,

    # sp2 ring nitrogens - 1 lone pair (in plane)
    ('A', 'N1'): 1,
    ('A', 'N3'): 1,
    ('A', 'N7'): 1,
    ('G', 'N3'): 1,
    ('G', 'N7'): 1,
    ('C', 'N3'): 1,

    # Ribose atoms (sp3 oxygens) - 2 lone pairs each
    ('A', "O2'"): 2,
    ('G', "O2'"): 2,
    ('C', "O2'"): 2,
    ('U', "O2'"): 2,
    ('T', "O2'"): 2,
    ('A', "O4'"): 1,
    ('G', "O4'"): 1,
    ('C', "O4'"): 1,
    ('U', "O4'"): 1,
    ('T', "O4'"): 1,
    # O3' hydroxyl - 2 lone pairs
    ('A', "O3'"): 2,
    ('G', "O3'"): 2,
    ('C', "O3'"): 2,
    ('U', "O3'"): 2,
    ('T', "O3'"): 2,
    ('P', "O3'"): 2,
    ('I', "O3'"): 2,
    # O5' hydroxyl - 2 lone pairs
    ('A', "O5'"): 2,
    ('G', "O5'"): 2,
    ('C', "O5'"): 2,
    ('U', "O5'"): 2,
    ('T', "O5'"): 2,
    ('P', "O5'"): 2,
    ('I', "O5'"): 2,

    # Phosphate oxygens - 3 lone pairs each (tetrahedral)
    ('A', 'OP1'): 3,
    ('G', 'OP1'): 3,
    ('C', 'OP1'): 3,
    ('U', 'OP1'): 3,
    ('T', 'OP1'): 3,
    ('A', 'OP2'): 3,
    ('G', 'OP2'): 3,
    ('C', 'OP2'): 3,
    ('U', 'OP2'): 3,
    ('T', 'OP2'): 3,
    # Legacy names
    ('A', 'O1P'): 3,
    ('G', 'O1P'): 3,
    ('C', 'O1P'): 3,
    ('U', 'O1P'): 3,
    ('A', 'O2P'): 3,
    ('G', 'O2P'): 3,
    ('C', 'O2P'): 3,
    ('U', 'O2P'): 3,

    # Pseudouridine (P) - similar to uridine
    ('P', 'O2'): 2,   # Carbonyl oxygen
    ('P', 'O4'): 2,   # Carbonyl oxygen
    ('P', "O2'"): 2,  # Ribose
    ('P', "O4'"): 1,  # Ribose ring
    ('P', 'OP1'): 3,  # Phosphate
    ('P', 'OP2'): 3,

    # Inosine (I) - similar to guanine
    ('I', 'O6'): 2,   # Carbonyl oxygen
    ('I', 'N3'): 1,   # Ring N acceptor
    ('I', 'N7'): 1,   # Ring N acceptor
    ('I', "O2'"): 2,  # Ribose
    ('I', "O4'"): 1,  # Ribose ring
    ('I', 'OP1'): 3,  # Phosphate
    ('I', 'OP2'): 3,
}


def compute_base_normal(atoms: dict) -> np.ndarray:
    """
    Compute the normal vector to the base plane.

    Uses C2, C4, C6 atoms which are present in all bases.
    """
    # Try to use ring atoms
    ring_atoms = ['C2', 'C4', 'C6', 'N1', 'N3']
    positions = []
    for name in ring_atoms:
        if name in atoms:
            positions.append(atoms[name])

    if len(positions) < 3:
        return np.array([0.0, 0.0, 1.0])  # Default

    # Fit plane using first 3 atoms
    v1 = positions[1] - positions[0]
    v2 = positions[2] - positions[0]
    normal = np.cross(v1, v2)
    return normalize(normal)


def predict_h_slots(base_type: str, atom_name: str,
                    atoms: dict, base_normal: np.ndarray) -> List[HSlot]:
    """
    Predict hydrogen slot positions for a donor atom.

    Args:
        base_type: Single letter base code (A, G, C, U, T) OR amino acid 3-letter code
        atom_name: Name of the donor atom (e.g., 'N6', 'N1', 'N', 'NZ')
        atoms: Dict mapping atom names to positions (np.ndarray)
        base_normal: Normal vector to the base plane (may be None for proteins)

    Returns:
        List of HSlot objects with predicted H directions
    """
    base_type_upper = base_type.upper()
    atom_name_stripped = atom_name.strip()
    key = (base_type_upper, atom_name_stripped)

    # Check nucleotide connectivity first, then protein
    if key in DONOR_CAPACITY:
        capacity = DONOR_CAPACITY[key]
        connectivity = BASE_CONNECTIVITY.get(key, [])
    elif key in PROTEIN_CONNECTIVITY:
        # Use CIF-based capacity for amino acids
        capacity = get_donor_capacity(base_type_upper, atom_name_stripped)
        connectivity = PROTEIN_CONNECTIVITY.get(key, [])
    else:
        # Try CIF data for modified residues
        capacity = get_donor_capacity(base_type_upper, atom_name_stripped)
        connectivity = BASE_CONNECTIVITY.get(key, []) or PROTEIN_CONNECTIVITY.get(key, [])
        if capacity == 0:
            return []

    if capacity == 0:
        return []

    if atom_name_stripped not in atoms:
        return []

    donor_pos = atoms[atom_name_stripped]

    # Get antecedent atom positions
    antecedent_positions = []
    for ant_name in connectivity:
        if ant_name in atoms:
            antecedent_positions.append(atoms[ant_name])

    if not antecedent_positions:
        return []

    slots = []

    # Determine if we have a valid base normal for rotation
    has_normal = base_normal is not None and np.linalg.norm(base_normal) > 0.1

    if capacity == 2 and len(antecedent_positions) == 1:
        # sp2 NH2: Two H atoms at 120° from the C-N bond
        ant_to_donor = normalize(donor_pos - antecedent_positions[0])

        if has_normal:
            # Rotate ±120° around base normal
            h1_dir = rotate_vector(ant_to_donor, base_normal, 120.0)
            h2_dir = rotate_vector(ant_to_donor, base_normal, -120.0)
        else:
            # For proteins without base normal, use arbitrary perpendicular axis
            perp = np.cross(ant_to_donor, np.array([0.0, 0.0, 1.0]))
            if np.linalg.norm(perp) < 0.1:
                perp = np.cross(ant_to_donor, np.array([0.0, 1.0, 0.0]))
            perp = normalize(perp)
            h1_dir = rotate_vector(ant_to_donor, perp, 120.0)
            h2_dir = rotate_vector(ant_to_donor, perp, -120.0)

        slots.append(HSlot(direction=h1_dir))
        slots.append(HSlot(direction=h2_dir))

    elif capacity == 3 and len(antecedent_positions) == 1:
        # sp3 NH3+ (like LYS NZ): Three H atoms in tetrahedral arrangement
        ant_to_donor = normalize(donor_pos - antecedent_positions[0])

        # Create 3 H directions at ~109° from the C-N bond
        # Use arbitrary perpendicular as rotation axis
        perp = np.cross(ant_to_donor, np.array([0.0, 0.0, 1.0]))
        if np.linalg.norm(perp) < 0.1:
            perp = np.cross(ant_to_donor, np.array([0.0, 1.0, 0.0]))
        perp = normalize(perp)

        # Tilt away from ant_to_donor by ~70° (180-109)
        h_base = rotate_vector(ant_to_donor, perp, 70.0)

        # Rotate around ant_to_donor axis
        h1_dir = h_base
        h2_dir = rotate_vector(h_base, ant_to_donor, 120.0)
        h3_dir = rotate_vector(h_base, ant_to_donor, -120.0)

        slots.append(HSlot(direction=h1_dir))
        slots.append(HSlot(direction=h2_dir))
        slots.append(HSlot(direction=h3_dir))

    elif capacity == 1 and len(antecedent_positions) == 2:
        # sp2 imino NH or ring C-H: H points away from ring
        avg_ant = (antecedent_positions[0] + antecedent_positions[1]) / 2
        h_dir = normalize(donor_pos - avg_ant)
        slots.append(HSlot(direction=h_dir))

    elif capacity == 1 and len(antecedent_positions) == 1:
        # Simple case: H opposite to single antecedent
        h_dir = normalize(donor_pos - antecedent_positions[0])
        slots.append(HSlot(direction=h_dir))

    return slots


def predict_lp_slots(base_type: str, atom_name: str,
                     atoms: dict, base_normal: np.ndarray) -> List[LPSlot]:
    """
    Predict lone pair slot directions for an acceptor atom.

    Args:
        base_type: Single letter base code OR amino acid 3-letter code
        atom_name: Name of acceptor atom (e.g., 'O6', 'N1', 'O', 'OD1')
        atoms: Dict mapping atom names to positions
        base_normal: Normal vector to base plane (may be None for proteins)

    Returns:
        List of LPSlot objects with predicted LP directions
    """
    base_type_upper = base_type.upper()
    atom_name_stripped = atom_name.strip()
    key = (base_type_upper, atom_name_stripped)

    # Check nucleotide capacity first, then protein
    if key in ACCEPTOR_CAPACITY:
        capacity = ACCEPTOR_CAPACITY[key]
        connectivity = BASE_CONNECTIVITY.get(key, [])
    elif key in PROTEIN_CONNECTIVITY:
        # Use CIF-based capacity for amino acids
        capacity = get_acceptor_capacity(base_type_upper, atom_name_stripped)
        connectivity = PROTEIN_CONNECTIVITY.get(key, [])
    else:
        # Try CIF data for modified residues
        capacity = get_acceptor_capacity(base_type_upper, atom_name_stripped)
        connectivity = BASE_CONNECTIVITY.get(key, []) or PROTEIN_CONNECTIVITY.get(key, [])
        if capacity == 0:
            return []

    if capacity == 0:
        return []

    if atom_name_stripped not in atoms:
        return []

    acceptor_pos = atoms[atom_name_stripped]

    antecedent_positions = []
    for ant_name in connectivity:
        if ant_name in atoms:
            antecedent_positions.append(atoms[ant_name])

    slots = []

    # Determine if we have a valid base normal
    has_normal = base_normal is not None and np.linalg.norm(base_normal) > 0.1

    # Handle phosphate oxygens specially - use isotropic LP model
    if atom_name_stripped in ('OP1', 'OP2', 'O1P', 'O2P'):
        # Phosphate oxygens have 3 lone pairs in rough tetrahedral arrangement
        # Use 3 orthogonal directions as approximation
        slots.append(LPSlot(direction=np.array([1.0, 0.0, 0.0])))
        slots.append(LPSlot(direction=np.array([0.0, 1.0, 0.0])))
        slots.append(LPSlot(direction=np.array([0.0, 0.0, 1.0])))
        return slots

    # Handle ribose O2' and O4' - sp3 with flexible geometry
    if atom_name_stripped in ("O2'", "O4'", "O3'", "O5'"):
        if has_normal:
            # Use 2 orthogonal directions perpendicular to base normal
            perp1 = np.cross(base_normal, np.array([1.0, 0.0, 0.0]))
            if np.linalg.norm(perp1) < 0.1:
                perp1 = np.cross(base_normal, np.array([0.0, 1.0, 0.0]))
            perp1 = normalize(perp1)
            perp2 = np.cross(base_normal, perp1)
            slots.append(LPSlot(direction=perp1))
            slots.append(LPSlot(direction=perp2))
        else:
            # Fallback: isotropic
            slots.append(LPSlot(direction=np.array([1.0, 0.0, 0.0])))
            slots.append(LPSlot(direction=np.array([0.0, 1.0, 0.0])))
        return slots

    if not antecedent_positions:
        # Fallback: isotropic
        slots.append(LPSlot(direction=np.array([1.0, 0.0, 0.0])))
        if capacity >= 2:
            slots.append(LPSlot(direction=np.array([0.0, 1.0, 0.0])))
        return slots

    if capacity == 2 and len(antecedent_positions) == 1:
        # sp2 carbonyl oxygen: Two lone pairs at 120° from C=O
        ant_to_acc = normalize(acceptor_pos - antecedent_positions[0])

        if has_normal:
            lp1_dir = rotate_vector(ant_to_acc, base_normal, 120.0)
            lp2_dir = rotate_vector(ant_to_acc, base_normal, -120.0)
        else:
            # For proteins without base normal, use arbitrary perpendicular axis
            perp = np.cross(ant_to_acc, np.array([0.0, 0.0, 1.0]))
            if np.linalg.norm(perp) < 0.1:
                perp = np.cross(ant_to_acc, np.array([0.0, 1.0, 0.0]))
            perp = normalize(perp)
            lp1_dir = rotate_vector(ant_to_acc, perp, 120.0)
            lp2_dir = rotate_vector(ant_to_acc, perp, -120.0)

        slots.append(LPSlot(direction=lp1_dir))
        slots.append(LPSlot(direction=lp2_dir))

    elif capacity == 1 and len(antecedent_positions) == 2:
        # sp2 ring nitrogen: One lone pair pointing out of ring
        avg_ant = (antecedent_positions[0] + antecedent_positions[1]) / 2
        lp_dir = normalize(acceptor_pos - avg_ant)
        slots.append(LPSlot(direction=lp_dir))

    elif capacity == 1 and len(antecedent_positions) == 1:
        # Fallback
        lp_dir = normalize(acceptor_pos - antecedent_positions[0])
        slots.append(LPSlot(direction=lp_dir))

    return slots


def score_hbond_alignment(donor_pos: np.ndarray, acceptor_pos: np.ndarray,
                          h_slots: List[HSlot], lp_slots: List[LPSlot],
                          ignore_used: bool = False
                          ) -> Tuple[int, int, float]:
    """
    Score how well an H-bond aligns with H and LP slots.

    Args:
        donor_pos: Position of donor atom
        acceptor_pos: Position of acceptor atom
        h_slots: Available H slots on donor
        lp_slots: Available LP slots on acceptor
        ignore_used: If True, consider all slots regardless of used status.
                     Useful for diagnostic purposes to compute "theoretical" alignment.

    Returns:
        (best_h_idx, best_lp_idx, alignment_score)
        Score is sum of two dot products, max = 2.0 for perfect alignment
    """
    if not h_slots or not lp_slots:
        return 0, 0, 0.0

    donor_to_acceptor = normalize(acceptor_pos - donor_pos)
    acceptor_to_donor = -donor_to_acceptor

    # Find best matching H slot
    best_h_score = -2.0
    best_h_idx = 0
    for i, h_slot in enumerate(h_slots):
        if not ignore_used and h_slot.used:
            continue
        alignment = np.dot(h_slot.direction, donor_to_acceptor)
        if alignment > best_h_score:
            best_h_score = alignment
            best_h_idx = i

    # Find best matching LP slot
    best_lp_score = -2.0
    best_lp_idx = 0
    for i, lp_slot in enumerate(lp_slots):
        if not ignore_used and lp_slot.used:
            continue
        alignment = np.dot(lp_slot.direction, acceptor_to_donor)
        if alignment > best_lp_score:
            best_lp_score = alignment
            best_lp_idx = i

    total_score = best_h_score + best_lp_score
    return best_h_idx, best_lp_idx, total_score

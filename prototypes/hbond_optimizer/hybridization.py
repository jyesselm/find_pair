"""
Hybridization detection and H-bond inference from molecular geometry.

This module provides functions to detect hybridization state from bond angles
and infer H-bond donor/acceptor capability without explicit hydrogen atoms.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from enum import Enum


class Hybridization(Enum):
    """Hybridization states."""
    SP = "sp"      # Linear, 180 degrees
    SP2 = "sp2"    # Trigonal planar, 120 degrees
    SP3 = "sp3"    # Tetrahedral, 109.5 degrees
    UNKNOWN = "unknown"


@dataclass
class AtomInfo:
    """Information about an atom for H-bond inference."""
    name: str
    element: str
    position: np.ndarray
    hybridization: Hybridization = Hybridization.UNKNOWN
    num_bonds: int = 0
    bonded_atoms: List[str] = None

    def __post_init__(self):
        if self.bonded_atoms is None:
            self.bonded_atoms = []


def get_element(atom_name: str) -> str:
    """Extract element from atom name (e.g., 'N6' -> 'N', "O2'" -> 'O')."""
    name = atom_name.strip()
    if not name:
        return ""

    # Handle common cases
    first_char = name[0].upper()
    if first_char in ('N', 'O', 'C', 'S', 'P', 'H'):
        return first_char

    # Two-letter elements
    if len(name) >= 2:
        two_char = name[:2].upper()
        if two_char in ('CL', 'BR', 'FE', 'ZN', 'MG', 'CA', 'MN', 'CO', 'NI', 'CU'):
            return two_char

    return first_char


def angle_between_vectors(v1: np.ndarray, v2: np.ndarray) -> float:
    """Calculate angle in degrees between two vectors."""
    v1_norm = np.linalg.norm(v1)
    v2_norm = np.linalg.norm(v2)
    if v1_norm < 1e-10 or v2_norm < 1e-10:
        return 0.0

    cos_angle = np.clip(np.dot(v1, v2) / (v1_norm * v2_norm), -1.0, 1.0)
    return np.degrees(np.arccos(cos_angle))


def detect_hybridization_from_angles(angles: List[float]) -> Hybridization:
    """
    Detect hybridization from bond angles.

    Args:
        angles: List of bond angles in degrees

    Returns:
        Detected hybridization state
    """
    if not angles:
        return Hybridization.UNKNOWN

    avg_angle = np.mean(angles)

    # sp: ~180 degrees (linear)
    if avg_angle > 160:
        return Hybridization.SP

    # sp2: ~120 degrees (trigonal planar)
    if 110 <= avg_angle <= 140:
        return Hybridization.SP2

    # sp3: ~109.5 degrees (tetrahedral)
    if 100 <= avg_angle < 115:
        return Hybridization.SP3

    # Borderline cases - use angle ranges
    if avg_angle >= 115:
        return Hybridization.SP2

    return Hybridization.SP3


def find_bonded_atoms(
    atom_name: str,
    atom_pos: np.ndarray,
    all_atoms: Dict[str, np.ndarray],
    max_bond_distance: float = 1.9
) -> List[Tuple[str, np.ndarray]]:
    """
    Find atoms bonded to the given atom based on distance.

    Args:
        atom_name: Name of the central atom
        atom_pos: Position of the central atom
        all_atoms: Dict of all atom names to positions
        max_bond_distance: Maximum distance for a covalent bond (Angstroms)

    Returns:
        List of (atom_name, position) tuples for bonded atoms
    """
    bonded = []
    for other_name, other_pos in all_atoms.items():
        if other_name == atom_name:
            continue

        dist = np.linalg.norm(other_pos - atom_pos)
        if dist <= max_bond_distance:
            bonded.append((other_name, other_pos))

    return bonded


def calculate_bond_angles(
    center_pos: np.ndarray,
    bonded_positions: List[np.ndarray]
) -> List[float]:
    """
    Calculate all bond angles at a central atom.

    Args:
        center_pos: Position of the central atom
        bonded_positions: Positions of atoms bonded to center

    Returns:
        List of angles between pairs of bonds
    """
    angles = []
    n = len(bonded_positions)

    if n < 2:
        return angles

    for i in range(n):
        for j in range(i + 1, n):
            v1 = bonded_positions[i] - center_pos
            v2 = bonded_positions[j] - center_pos
            angle = angle_between_vectors(v1, v2)
            angles.append(angle)

    return angles


def detect_atom_hybridization(
    atom_name: str,
    atoms: Dict[str, np.ndarray],
    max_bond_distance: float = 1.9
) -> Tuple[Hybridization, int, List[str]]:
    """
    Detect hybridization of an atom from its geometry.

    Args:
        atom_name: Name of the atom to analyze
        atoms: Dict of all atom names to positions
        max_bond_distance: Maximum covalent bond distance

    Returns:
        (hybridization, num_bonds, list of bonded atom names)
    """
    if atom_name not in atoms:
        return Hybridization.UNKNOWN, 0, []

    atom_pos = atoms[atom_name]
    bonded = find_bonded_atoms(atom_name, atom_pos, atoms, max_bond_distance)

    if not bonded:
        return Hybridization.UNKNOWN, 0, []

    bonded_names = [name for name, _ in bonded]
    bonded_positions = [pos for _, pos in bonded]
    num_bonds = len(bonded)

    # Calculate bond angles
    angles = calculate_bond_angles(atom_pos, bonded_positions)

    if not angles:
        # Single bond - can't determine hybridization from geometry alone
        # Use element-based heuristics
        element = get_element(atom_name)
        if element == 'O':
            # Single-bonded oxygen is likely sp3 (hydroxyl or ether)
            return Hybridization.SP3, num_bonds, bonded_names
        elif element == 'N':
            # Could be sp2 or sp3 depending on context
            return Hybridization.SP3, num_bonds, bonded_names
        return Hybridization.UNKNOWN, num_bonds, bonded_names

    hybridization = detect_hybridization_from_angles(angles)
    return hybridization, num_bonds, bonded_names


def infer_donor_capacity(
    element: str,
    hybridization: Hybridization,
    num_bonds: int,
    bonded_elements: List[str]
) -> int:
    """
    Infer H-bond donor capacity from atom properties.

    An atom can donate if it has H atoms attached.
    We infer the number of H from the element, hybridization, and bond count.

    Args:
        element: Element symbol (N, O, etc.)
        hybridization: Detected hybridization state
        num_bonds: Number of heavy atom bonds detected
        bonded_elements: Elements of bonded atoms

    Returns:
        Number of H-bonds the atom can donate (0 if not a donor)
    """
    if element not in ('N', 'O', 'S'):
        return 0  # Only N, O, S can be donors in typical H-bonds

    if element == 'N':
        if hybridization == Hybridization.SP3:
            # sp3 nitrogen: NH3 = 3H, NH2R = 2H, NHR2 = 1H, NR3 = 0H
            expected_bonds = 4  # N typically has 4 bonds total (including H)
            implicit_h = max(0, expected_bonds - num_bonds - 1)  # -1 for lone pair consideration
            # Actually for sp3 N: 3 bonds + 1 LP
            implicit_h = max(0, 3 - num_bonds)
            return implicit_h

        elif hybridization == Hybridization.SP2:
            # sp2 nitrogen (amide, imine, aromatic)
            # Check if it's in a ring (aromatic) - these typically don't have H
            if num_bonds == 2:
                # Could be NH (imino) or =N- (imine without H)
                # If one bond is to C and looks like C=N, likely no H
                # If it's like amide N-C(=O), it has 1 H
                # Conservative: assume 1 H unless clearly 3 bonds
                return 1
            elif num_bonds == 3:
                # Likely no H (e.g., tertiary amide or aromatic N)
                return 0
            else:
                return max(0, 2 - num_bonds)

        elif hybridization == Hybridization.SP:
            # Linear nitrogen (nitrile, isocyanide) - no H
            return 0

    elif element == 'O':
        if hybridization == Hybridization.SP3:
            # sp3 oxygen: OH2 = 2H, OH = 1H, O-R2 = 0H (ether)
            # Oxygen has 2 bonds + 2 lone pairs
            if num_bonds == 1:
                # Hydroxyl - has 1 H
                return 1
            else:
                # Ether or other - no H to donate
                return 0

        elif hybridization == Hybridization.SP2:
            # sp2 oxygen (carbonyl, carboxyl)
            # Carbonyl O=C has no H
            # Carboxylic OH might, but that's sp3
            return 0

    elif element == 'S':
        if hybridization == Hybridization.SP3 and num_bonds == 1:
            # Thiol - has 1 H
            return 1
        return 0

    return 0


def infer_acceptor_capacity(
    element: str,
    hybridization: Hybridization,
    num_bonds: int
) -> int:
    """
    Infer H-bond acceptor capacity from atom properties.

    An atom can accept H-bonds through its lone pairs.

    Args:
        element: Element symbol
        hybridization: Detected hybridization state
        num_bonds: Number of heavy atom bonds

    Returns:
        Number of H-bonds the atom can accept (lone pair count for H-bonding)
    """
    if element not in ('N', 'O', 'S', 'F', 'CL'):
        return 0

    if element == 'O':
        if hybridization == Hybridization.SP3:
            # sp3 oxygen has 2 lone pairs (water, alcohol, ether)
            return 2
        elif hybridization == Hybridization.SP2:
            # sp2 oxygen (carbonyl) has 2 lone pairs in the plane
            return 2
        else:
            return 2  # Default for oxygen

    elif element == 'N':
        if hybridization == Hybridization.SP3:
            # sp3 nitrogen has 1 lone pair (amine)
            return 1
        elif hybridization == Hybridization.SP2:
            # sp2 nitrogen in ring or imine has 1 lone pair
            if num_bonds == 2:
                # Ring N with 2 bonds - 1 LP in plane
                return 1
            elif num_bonds == 3:
                # sp2 N with 3 bonds - like amide N, still has some LP character
                # but typically considered poor acceptor
                return 0
        elif hybridization == Hybridization.SP:
            # sp nitrogen (nitrile) - 1 LP
            return 1
        return 1  # Default for nitrogen

    elif element == 'S':
        # Sulfur has 2 lone pairs
        return 2

    elif element in ('F', 'CL'):
        # Halogens have 3 lone pairs but typically only 1 used for H-bonding
        return 1

    return 0


def analyze_residue_hbond_properties(
    residue_name: str,
    atoms: Dict[str, np.ndarray],
    max_bond_distance: float = 1.9
) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Analyze a residue to determine H-bond donor and acceptor atoms.

    Args:
        residue_name: Name of the residue (for logging)
        atoms: Dict mapping atom names to positions
        max_bond_distance: Maximum covalent bond distance

    Returns:
        (donors, acceptors) - dicts mapping atom names to capacities
    """
    donors = {}
    acceptors = {}

    for atom_name, atom_pos in atoms.items():
        element = get_element(atom_name)

        # Skip non-heteroatoms for H-bonding
        if element not in ('N', 'O', 'S', 'F'):
            continue

        # Detect hybridization
        hybridization, num_bonds, bonded_names = detect_atom_hybridization(
            atom_name, atoms, max_bond_distance
        )

        # Get elements of bonded atoms
        bonded_elements = [get_element(name) for name in bonded_names]

        # Infer donor capacity
        donor_cap = infer_donor_capacity(element, hybridization, num_bonds, bonded_elements)
        if donor_cap > 0:
            donors[atom_name] = donor_cap

        # Infer acceptor capacity
        acceptor_cap = infer_acceptor_capacity(element, hybridization, num_bonds)
        if acceptor_cap > 0:
            acceptors[atom_name] = acceptor_cap

    return donors, acceptors


# Special cases for amino acids backbone
BACKBONE_DONORS = {'N': 1}  # Backbone NH (except proline)
BACKBONE_ACCEPTORS = {'O': 2}  # Backbone carbonyl

# Common functional group patterns
AMINO_PATTERN = {  # -NH2
    'donor': 2,  # 2 hydrogens
    'acceptor': 1  # 1 lone pair (though usually poor)
}

IMINO_PATTERN = {  # >NH
    'donor': 1,
    'acceptor': 1
}

CARBONYL_PATTERN = {  # C=O
    'donor': 0,
    'acceptor': 2
}

HYDROXYL_PATTERN = {  # -OH
    'donor': 1,
    'acceptor': 2
}

CARBOXYL_PATTERN = {  # -COO-
    'donor': 0,  # Deprotonated typically
    'acceptor': 2  # Per oxygen
}


def detect_base_type(atoms: Dict[str, np.ndarray]) -> Optional[str]:
    """
    Detect the base type from atoms present.

    Returns:
        'A', 'G', 'C', 'U', 'T', or None if unknown
    """
    atom_names = set(a.strip() for a in atoms.keys())

    # Purines have N9, N7, N3, N1
    is_purine = 'N9' in atom_names and 'N7' in atom_names

    # Check for specific atoms to distinguish bases
    if is_purine:
        if 'N6' in atom_names:
            return 'A'  # Adenine has amino N6
        elif 'O6' in atom_names and 'N2' in atom_names:
            return 'G'  # Guanine has O6 and amino N2
        elif 'O6' in atom_names:
            return 'I'  # Inosine has O6 but no N2

    else:
        # Pyrimidines
        if 'N4' in atom_names:
            return 'C'  # Cytosine has amino N4
        elif 'O4' in atom_names:
            if 'C7' in atom_names or 'C5M' in atom_names:
                return 'T'  # Thymine has C7/methyl
            else:
                return 'U'  # Uracil

    return None


def classify_nitrogen(
    n_atom: str,
    atoms: Dict[str, np.ndarray],
    max_bond_dist: float = 1.9
) -> Tuple[int, int]:
    """
    Classify a nitrogen atom and return (donor_capacity, acceptor_capacity).

    Uses geometry, atom naming, AND base type detection to distinguish between:
    - Primary amine (NH2): donor=2, acceptor=0
    - Imino (NH): donor=1, acceptor=0
    - Ring nitrogen: donor=0, acceptor=1
    - Glycosidic nitrogen (N9): donor=0, acceptor=0
    """
    if n_atom not in atoms:
        return 0, 0

    n_pos = atoms[n_atom]
    bonded = find_bonded_atoms(n_atom, n_pos, atoms, max_bond_dist)
    num_bonds = len(bonded)
    stripped_name = n_atom.strip()

    if num_bonds == 0:
        return 0, 0

    # Detect base type for context
    base_type = detect_base_type(atoms)

    # N9 is glycosidic in purines - not a donor or acceptor
    if stripped_name == 'N9':
        return 0, 0

    # Amino nitrogens - named N2, N4, N6 in bases
    if stripped_name in ('N2', 'N4', 'N6'):
        return 2, 0  # Amino group -NH2

    # Handle specific ring nitrogens based on base type
    if stripped_name == 'N1':
        if base_type == 'G':
            # Guanine N1 is imino (has H)
            return 1, 0
        elif base_type in ('A', 'C', 'U', 'T', 'I'):
            # Adenine N1 is ring acceptor (pyrimidine N1 is glycosidic-like)
            if base_type == 'A':
                return 0, 1  # Acceptor
            else:
                return 0, 0  # Glycosidic, not involved in base pairing

    if stripped_name == 'N3':
        if base_type in ('U', 'T'):
            # Uracil/Thymine N3 is imino (has H)
            return 1, 0
        elif base_type in ('A', 'G', 'C', 'I'):
            # Adenine/Guanine/Cytosine N3 is ring acceptor
            return 0, 1

    if stripped_name == 'N7':
        # N7 in purines is always acceptor
        return 0, 1

    # For amino acid backbone N
    if stripped_name == 'N':
        # Backbone amide - typically 1 H (except proline)
        return 1, 0

    # Generic classification for unknown nitrogens
    if num_bonds >= 2:
        bonded_positions = [pos for _, pos in bonded]
        angles = calculate_bond_angles(n_pos, bonded_positions)
        avg_angle = np.mean(angles) if angles else 109.5

        if avg_angle > 115:  # sp2-like
            if num_bonds == 2:
                # sp2 N with 2 bonds: check for carbonyl neighbors
                has_carbonyl_neighbor = False
                for name, _ in bonded:
                    if get_element(name) == 'C':
                        c_pos = atoms.get(name)
                        if c_pos is not None:
                            c_bonded = find_bonded_atoms(name, c_pos, atoms, max_bond_dist)
                            for c_neighbor, c_neighbor_pos in c_bonded:
                                if get_element(c_neighbor) == 'O':
                                    co_dist = np.linalg.norm(c_neighbor_pos - c_pos)
                                    if co_dist < 1.35:
                                        has_carbonyl_neighbor = True
                                        break

                if has_carbonyl_neighbor:
                    return 1, 0  # Imino adjacent to carbonyl
                else:
                    return 0, 1  # Ring nitrogen

            elif num_bonds == 3:
                return 0, 1  # Tertiary N

        else:  # sp3-like
            if num_bonds == 1:
                return 2, 0  # Primary amine
            elif num_bonds == 2:
                return 1, 0  # Secondary amine
            elif num_bonds == 3:
                return 0, 1  # Tertiary amine

    elif num_bonds == 1:
        return 2, 0  # Primary amine

    return 0, 0


def classify_oxygen(
    o_atom: str,
    atoms: Dict[str, np.ndarray],
    max_bond_dist: float = 1.9
) -> Tuple[int, int]:
    """
    Classify an oxygen atom and return (donor_capacity, acceptor_capacity).

    Uses geometry AND atom naming to distinguish between:
    - Carbonyl (C=O): donor=0, acceptor=2 (O2, O4, O6 in bases)
    - Hydroxyl (-OH): donor=1, acceptor=2 (O2', O3', O5')
    - Ether (-O-): donor=0, acceptor=1-2 (O4' ring)
    - Phosphate: donor=0, acceptor=2-3 (OP1, OP2)
    """
    if o_atom not in atoms:
        return 0, 0

    o_pos = atoms[o_atom]
    bonded = find_bonded_atoms(o_atom, o_pos, atoms, max_bond_dist)
    num_bonds = len(bonded)
    stripped_name = o_atom.strip()

    if num_bonds == 0:
        return 0, 0

    # Use atom naming conventions for nucleic acids
    # O2, O4, O6 are carbonyl oxygens in bases - no H
    if stripped_name in ('O2', 'O4', 'O6'):
        return 0, 2  # Carbonyl - acceptor only

    # O2', O3', O5' are hydroxyl oxygens - have H
    if stripped_name in ("O2'", "O3'", "O5'"):
        return 1, 2  # Hydroxyl - donor + acceptor

    # O4' is ring oxygen in sugar - no H, less accessible
    if stripped_name in ("O4'",):
        return 0, 1  # Ether - acceptor only, reduced capacity

    # OP1, OP2, O1P, O2P, O3P are phosphate oxygens
    if stripped_name in ('OP1', 'OP2', 'OP3', 'O1P', 'O2P', 'O3P'):
        return 0, 2  # Phosphate - acceptor only

    # For unknown oxygens, use geometry
    bonded_elements = [get_element(name) for name, _ in bonded]

    if num_bonds == 1:
        bond_name, bond_pos = bonded[0]
        bond_element = get_element(bond_name)
        bond_dist = np.linalg.norm(bond_pos - o_pos)

        if bond_element == 'P':
            # Phosphate oxygen - acceptor only
            return 0, 2

        elif bond_element == 'C':
            # Check C-O bond length to distinguish carbonyl from hydroxyl
            if bond_dist < 1.30:
                # Double bond (carbonyl) - very short
                return 0, 2

            # Check if the carbon is sp2 (has other double bonds)
            c_pos = atoms.get(bond_name)
            if c_pos is not None:
                c_bonded = find_bonded_atoms(bond_name, c_pos, atoms, max_bond_dist)

                # Count oxygens and nitrogens bonded to this carbon
                o_count = sum(1 for name, _ in c_bonded if get_element(name) == 'O')
                n_count = sum(1 for name, _ in c_bonded if get_element(name) == 'N')

                # If carbon has 2 oxygens (carboxyl) or is in aromatic ring
                if o_count >= 2:
                    # Find the other oxygen
                    for other_name, other_pos in c_bonded:
                        if get_element(other_name) == 'O' and other_name != o_atom:
                            other_dist = np.linalg.norm(other_pos - c_pos)
                            if bond_dist < other_dist:
                                # This is the shorter bond - carbonyl
                                return 0, 2
                            else:
                                # This is the longer bond - could be OH or O-
                                return 0, 2  # Usually deprotonated at pH 7

                # If carbon is bonded to N (amide/urea), this O is carbonyl
                if n_count >= 1:
                    return 0, 2

                # Check number of carbons bonded to this carbon
                c_count = sum(1 for name, _ in c_bonded if get_element(name) == 'C')
                if c_count >= 2 and bond_dist < 1.35:
                    # Likely carbonyl in aromatic ring
                    return 0, 2

            # Default for single-bonded O to C with long bond: hydroxyl
            if bond_dist >= 1.35:
                return 1, 2

            return 0, 2  # Conservative: assume carbonyl

        elif bond_element == 'S':
            return 0, 2  # Sulfoxide/sulfone

        elif bond_element == 'N':
            return 0, 2  # N-oxide

        # Default for unknown single-bonded O: assume hydroxyl
        return 1, 2

    elif num_bonds == 2:
        # Ether or ring oxygen
        return 0, 1  # Reduced acceptor capacity for ethers

    return 0, 0


if __name__ == "__main__":
    # Test with a simple example
    print("Hybridization Detection Module")
    print("=" * 40)

    # Example: Adenine N6 amino group
    test_atoms = {
        'C6': np.array([0.0, 0.0, 0.0]),
        'N6': np.array([1.35, 0.0, 0.0]),  # Amino N attached to C6
        'N1': np.array([-0.7, 1.2, 0.0]),
        'C5': np.array([-0.7, -1.2, 0.0]),
    }

    donors, acceptors = analyze_residue_hbond_properties("ADE", test_atoms)

    print(f"Detected donors: {donors}")
    print(f"Detected acceptors: {acceptors}")

    # Test nitrogen classification
    print("\nNitrogen classification for N6:")
    d, a = classify_nitrogen('N6', test_atoms)
    print(f"  Donor capacity: {d}")
    print(f"  Acceptor capacity: {a}")

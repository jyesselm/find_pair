#!/usr/bin/env python3
"""
Generate idealized base pair structures with optimal H-bond geometry.

Takes exemplar base pairs and creates idealized versions where:
1. Both bases are perfectly planar (z=0)
2. Hydrogen bonds have optimal distances (~2.9 Å)

Usage:
    python tools/generate_idealized_basepairs.py
    python tools/generate_idealized_basepairs.py --lw-class cWW -v
"""

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import minimize

# Add tools directory to path
sys.path.insert(0, str(Path(__file__).parent))

from bp_catalog_parser import load_database, ExemplarEntry
from bp_catalog_extractor import BaseTemplates, encode_sequence_for_filename


@dataclass
class AttachedAtomInfo:
    """Info about attached atom for angle calculation."""
    attached: str       # Attached atom name
    ideal_angle: float  # Ideal angle in degrees
    atom_type: str      # Type (amino, carbonyl, ring_N, ring_NH)


@dataclass
class HBondDef:
    """Single H-bond definition."""
    donor_base: int      # 1 or 2
    donor_atom: str      # Atom name on donor base
    acceptor_atom: str   # Atom name on acceptor base
    ideal_dist: float    # Ideal distance in Angstroms
    ideal_donor_angle: float = 120.0    # Ideal X-D...A angle
    ideal_acceptor_angle: float = 120.0  # Ideal D...A-Y angle


@dataclass
class HBondPattern:
    """H-bond pattern for a specific LW class + sequence."""
    lw_class: str
    sequence: str
    hbonds: List[HBondDef]


@dataclass
class HBondGeometry:
    """Measured H-bond geometry."""
    donor_atom: str
    acceptor_atom: str
    distance: float
    donor_angle: Optional[float] = None   # X-D...A angle
    acceptor_angle: Optional[float] = None  # D...A-Y angle
    ideal_dist: float = 2.9
    ideal_donor_angle: float = 120.0
    ideal_acceptor_angle: float = 120.0


class HBondPatternDB:
    """Database of H-bond patterns."""

    def __init__(self, pattern_file: Path):
        with open(pattern_file) as f:
            self.data = json.load(f)
        self.patterns = self.data.get("patterns", {})
        self.donor_attached = self._load_attached_atoms("donor_attached_atoms")
        self.acceptor_attached = self._load_attached_atoms("acceptor_attached_atoms")

    def _load_attached_atoms(self, key: str) -> Dict[str, AttachedAtomInfo]:
        """Load attached atom info from JSON."""
        result = {}
        data = self.data.get(key, {})
        for atom_name, info in data.items():
            if atom_name.startswith("_"):
                continue
            result[atom_name] = AttachedAtomInfo(
                attached=info.get("attached", ""),
                ideal_angle=info.get("ideal_angle", 120.0),
                atom_type=info.get("type", "unknown")
            )
        return result

    def get_donor_attached(self, donor_atom: str) -> Optional[AttachedAtomInfo]:
        """Get attached atom info for a donor."""
        return self.donor_attached.get(donor_atom)

    def get_acceptor_attached(self, acceptor_atom: str) -> Optional[AttachedAtomInfo]:
        """Get attached atom info for an acceptor."""
        return self.acceptor_attached.get(acceptor_atom)

    def get_pattern(self, lw_class: str, sequence: str) -> Optional[HBondPattern]:
        """Get H-bond pattern for a specific LW class and sequence."""
        # Normalize sequence to uppercase for lookup
        seq_upper = sequence.upper()

        if lw_class not in self.patterns:
            return None

        lw_patterns = self.patterns[lw_class]
        if seq_upper not in lw_patterns:
            return None

        hbonds = []
        for hb_data in lw_patterns[seq_upper]:
            donor_atom = hb_data["donor_atom"]
            acceptor_atom = hb_data["acceptor_atom"]

            # Get ideal angles from attached atom info
            donor_info = self.get_donor_attached(donor_atom)
            acceptor_info = self.get_acceptor_attached(acceptor_atom)

            ideal_donor_angle = donor_info.ideal_angle if donor_info else 120.0
            ideal_acceptor_angle = acceptor_info.ideal_angle if acceptor_info else 120.0

            hbonds.append(HBondDef(
                donor_base=hb_data["donor_base"],
                donor_atom=donor_atom,
                acceptor_atom=acceptor_atom,
                ideal_dist=hb_data.get("ideal_dist", 2.9),
                ideal_donor_angle=ideal_donor_angle,
                ideal_acceptor_angle=ideal_acceptor_angle
            ))

        return HBondPattern(lw_class=lw_class, sequence=seq_upper, hbonds=hbonds)


def get_base_atoms(templates: BaseTemplates, base_type: str) -> Dict[str, np.ndarray]:
    """Get atom coordinates for a base from templates (already at z=0)."""
    base = base_type.upper()

    # Map modified bases to parent
    if base not in templates.templates:
        parent_map = {
            "1MA": "A", "MIA": "A", "M1A": "A", "6MA": "A", "M6A": "A",
            "1MG": "G", "2MG": "G", "7MG": "G", "M2G": "G", "OMG": "G",
            "5MC": "C", "OMC": "C", "5CM": "C",
            "PSU": "U", "H2U": "U", "5MU": "U", "OMU": "U",
        }
        base = parent_map.get(base, "A")

    return templates.templates.get(base, {}).copy()


def calculate_angle(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """
    Calculate angle at p2 formed by p1-p2-p3 in degrees.

    Returns angle between vectors (p1-p2) and (p3-p2).
    """
    v1 = p1 - p2
    v2 = p3 - p2

    v1_norm = np.linalg.norm(v1)
    v2_norm = np.linalg.norm(v2)

    if v1_norm < 1e-6 or v2_norm < 1e-6:
        return 0.0

    cos_angle = np.dot(v1, v2) / (v1_norm * v2_norm)
    cos_angle = np.clip(cos_angle, -1.0, 1.0)

    return math.degrees(math.acos(cos_angle))


# Mapping of donor atoms to their attached atoms for angle calculation
DONOR_ATTACHED_ATOMS = {
    "N1": "C2",  # Ring NH (G)
    "N2": "C2",  # Amino (G)
    "N3": "C2",  # Ring NH (U, T)
    "N4": "C4",  # Amino (C)
    "N6": "C6",  # Amino (A)
}

# Mapping of acceptor atoms to their attached atoms for angle calculation
ACCEPTOR_ATTACHED_ATOMS = {
    "O2": "C2",  # Carbonyl
    "O4": "C4",  # Carbonyl
    "O6": "C6",  # Carbonyl
    "N1": "C6",  # Ring N (A) - different attached than when donor
    "N3": "C4",  # Ring N (C, G when acceptor)
    "N7": "C5",  # Ring N (purines)
}


def calculate_hbond_angles(
    donor_coords: Dict[str, np.ndarray],
    acceptor_coords: Dict[str, np.ndarray],
    donor_atom: str,
    acceptor_atom: str
) -> Tuple[Optional[float], Optional[float]]:
    """
    Calculate H-bond angles.

    Returns:
        (donor_angle, acceptor_angle) where:
        - donor_angle: X-D...A angle (angle at donor heavy atom)
        - acceptor_angle: D...A-Y angle (angle at acceptor heavy atom)
        Returns None for either if atoms not found.
    """
    if donor_atom not in donor_coords or acceptor_atom not in acceptor_coords:
        return None, None

    d_coord = donor_coords[donor_atom]
    a_coord = acceptor_coords[acceptor_atom]

    donor_angle = None
    acceptor_angle = None

    # Calculate donor angle (X-D...A)
    attached_donor = DONOR_ATTACHED_ATOMS.get(donor_atom)
    if attached_donor and attached_donor in donor_coords:
        x_coord = donor_coords[attached_donor]
        donor_angle = calculate_angle(x_coord, d_coord, a_coord)

    # Calculate acceptor angle (D...A-Y)
    attached_acceptor = ACCEPTOR_ATTACHED_ATOMS.get(acceptor_atom)
    if attached_acceptor and attached_acceptor in acceptor_coords:
        y_coord = acceptor_coords[attached_acceptor]
        acceptor_angle = calculate_angle(d_coord, a_coord, y_coord)

    return donor_angle, acceptor_angle


def transform_coords(
    coords: Dict[str, np.ndarray],
    x: float,
    y: float,
    theta: float,
    flip_y: bool = False
) -> Dict[str, np.ndarray]:
    """
    Transform coordinates by rotation about z-axis and translation in xy-plane.

    Args:
        coords: Dict of atom_name -> [x, y, z] coordinates
        x, y: Translation in xy-plane
        theta: Rotation angle about z-axis (radians)
        flip_y: If True, mirror coordinates in y-axis (x -> -x) before rotation

    Returns:
        Transformed coordinates (all still in z=0 plane)
    """
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)

    result = {}
    for atom_name, coord in coords.items():
        cx, cy = coord[0], coord[1]

        # Optional y-axis flip (mirror)
        if flip_y:
            cx = -cx

        # Rotate about z-axis
        new_x = cx * cos_t - cy * sin_t
        new_y = cx * sin_t + cy * cos_t
        new_z = 0.0  # Force z=0 for perfect planarity

        # Translate
        result[atom_name] = np.array([new_x + x, new_y + y, new_z])

    return result


def calculate_hbond_distances(
    base1_coords: Dict[str, np.ndarray],
    base2_coords: Dict[str, np.ndarray],
    pattern: HBondPattern
) -> List[Tuple[float, float]]:
    """
    Calculate actual vs ideal H-bond distances.

    Returns list of (actual_distance, ideal_distance) tuples.
    """
    distances = []

    for hb in pattern.hbonds:
        if hb.donor_base == 1:
            donor_coords = base1_coords
            acceptor_coords = base2_coords
        else:
            donor_coords = base2_coords
            acceptor_coords = base1_coords

        if hb.donor_atom not in donor_coords or hb.acceptor_atom not in acceptor_coords:
            continue

        d_coord = donor_coords[hb.donor_atom]
        a_coord = acceptor_coords[hb.acceptor_atom]

        dist = np.linalg.norm(d_coord - a_coord)
        distances.append((dist, hb.ideal_dist))

    return distances


def calculate_hbond_geometry(
    base1_coords: Dict[str, np.ndarray],
    base2_coords: Dict[str, np.ndarray],
    pattern: HBondPattern
) -> List[HBondGeometry]:
    """
    Calculate full H-bond geometry including distances and angles.

    Returns list of HBondGeometry objects.
    """
    geometries = []

    for hb in pattern.hbonds:
        if hb.donor_base == 1:
            donor_coords = base1_coords
            acceptor_coords = base2_coords
        else:
            donor_coords = base2_coords
            acceptor_coords = base1_coords

        if hb.donor_atom not in donor_coords or hb.acceptor_atom not in acceptor_coords:
            continue

        d_coord = donor_coords[hb.donor_atom]
        a_coord = acceptor_coords[hb.acceptor_atom]

        dist = np.linalg.norm(d_coord - a_coord)

        # Calculate angles
        donor_angle, acceptor_angle = calculate_hbond_angles(
            donor_coords, acceptor_coords, hb.donor_atom, hb.acceptor_atom
        )

        geometries.append(HBondGeometry(
            donor_atom=hb.donor_atom,
            acceptor_atom=hb.acceptor_atom,
            distance=dist,
            donor_angle=donor_angle,
            acceptor_angle=acceptor_angle,
            ideal_dist=hb.ideal_dist,
            ideal_donor_angle=hb.ideal_donor_angle,
            ideal_acceptor_angle=hb.ideal_acceptor_angle
        ))

    return geometries


def calculate_hbond_score(
    params: np.ndarray,
    base1_coords: Dict[str, np.ndarray],
    base2_template: Dict[str, np.ndarray],
    pattern: HBondPattern,
    is_cis: bool = True,
    flip_y: bool = False,
    angle_weight: float = 0.5
) -> float:
    """
    Objective function: H-bond distances + angles + clash avoidance + geometry constraints.

    Args:
        params: [x, y, theta] - position/rotation of base2
        base1_coords: Fixed base1 coordinates
        base2_template: Template coordinates for base2 (will be transformed)
        pattern: H-bond pattern to optimize
        is_cis: True for cis pairs (same side), False for trans
        flip_y: Whether to mirror base2 in y-axis
        angle_weight: Weight for angle deviation penalty (default 0.5)

    Returns:
        Weighted score (lower is better)
    """
    x, y, theta = params

    # Transform base2
    base2_coords = transform_coords(base2_template, x, y, theta, flip_y=flip_y)

    score = 0.0

    # 1. H-bond distance term - penalize deviation from ideal more strongly
    hbond_atoms = set()
    max_hbond_dev = 0.0

    for hb in pattern.hbonds:
        if hb.donor_base == 1:
            donor_coords = base1_coords
            acceptor_coords = base2_coords
            hbond_atoms.add((hb.donor_atom, 1))
            hbond_atoms.add((hb.acceptor_atom, 2))
        else:
            donor_coords = base2_coords
            acceptor_coords = base1_coords
            hbond_atoms.add((hb.donor_atom, 2))
            hbond_atoms.add((hb.acceptor_atom, 1))

        if hb.donor_atom not in donor_coords or hb.acceptor_atom not in acceptor_coords:
            score += 100.0
            continue

        d_coord = donor_coords[hb.donor_atom]
        a_coord = acceptor_coords[hb.acceptor_atom]

        dist = np.linalg.norm(d_coord - a_coord)
        dev = abs(dist - hb.ideal_dist)

        # Penalize both too short and too long, but especially too short (clashing)
        if dist < hb.ideal_dist:
            score += 20.0 * dev ** 2  # Extra penalty for too short
        else:
            score += 10.0 * dev ** 2

        max_hbond_dev = max(max_hbond_dev, dev)

        # 2. H-bond angle terms
        donor_angle, acceptor_angle = calculate_hbond_angles(
            donor_coords, acceptor_coords, hb.donor_atom, hb.acceptor_atom
        )

        if donor_angle is not None:
            # Ideal donor angle is typically 120° (X-D...A for sp2 donor)
            # Penalize deviation from ideal
            angle_dev = abs(donor_angle - hb.ideal_donor_angle)
            score += angle_weight * (angle_dev / 10.0) ** 2  # Normalized by 10°

        if acceptor_angle is not None:
            # Ideal acceptor angle depends on type (120° for carbonyl, 125° for ring N)
            angle_dev = abs(acceptor_angle - hb.ideal_acceptor_angle)
            score += angle_weight * (angle_dev / 10.0) ** 2

    # Extra penalty for having any single H-bond very far from ideal
    if max_hbond_dev > 0.5:
        score += 20.0 * (max_hbond_dev - 0.5) ** 2

    # 3. Clash avoidance - penalize non-H-bond atoms that are too close
    min_dist = 2.5  # Minimum allowed distance for non-bonded atoms
    for name1, coord1 in base1_coords.items():
        if name1 == "C1'":
            continue  # Skip C1' for clash check
        for name2, coord2 in base2_coords.items():
            if name2 == "C1'":
                continue
            # Skip H-bond atom pairs
            if (name1, 1) in hbond_atoms and (name2, 2) in hbond_atoms:
                continue

            dist = np.linalg.norm(coord1 - coord2)
            if dist < min_dist:
                score += 50.0 * (min_dist - dist) ** 2

    # 4. C1'-C1' distance constraint (typical WC pair: ~10.5 Å)
    # Use lower weight so it doesn't prevent good H-bond geometry
    if "C1'" in base1_coords and "C1'" in base2_coords:
        c1_dist = np.linalg.norm(base1_coords["C1'"] - base2_coords["C1'"])
        ideal_c1_dist = 10.5
        # Only penalize if too close (allow larger distances for planar geometry)
        if c1_dist < ideal_c1_dist:
            score += 0.5 * (c1_dist - ideal_c1_dist) ** 2

    # 5. Planarity bonus - bases should stay in z=0 plane (already enforced by transform)

    return score


def get_initial_guess_for_basepair(
    lw_class: str,
    base1_type: str,
    base2_type: str
) -> Tuple[float, float, float]:
    """
    Get initial guess for base2 position based on LW class and base types.

    For Watson-Crick pairs: base2 is ~9 Å away and rotated 180°.
    The exact positioning depends on the edge interaction.
    """
    # Determine if cis or trans
    is_cis = lw_class.startswith('c')

    # Determine edge combination from LW class
    # Format: [c/t][edge1][edge2] where edges are W, H, or S
    edges = lw_class[1:3] if len(lw_class) >= 3 else "WW"

    # Base translation and rotation based on edge combination
    # These are approximate starting points for optimization

    if edges == "WW":
        # Watson-Watson: bases face each other along Watson edge
        # Typical C1'-C1' distance ~10.5 Å for WC pairs
        if is_cis:
            # cWW: Standard Watson-Crick geometry
            # Base2 is translated in -x direction and rotated 180°
            return (-9.0, -2.0, math.pi)
        else:
            # tWW: Trans Watson-Watson (reverse WC)
            return (-9.0, 2.0, 0.0)

    elif edges == "WH":
        # Watson-Hoogsteen
        if is_cis:
            return (-6.0, -5.0, math.pi * 0.7)
        else:
            return (-6.0, 5.0, -math.pi * 0.3)

    elif edges == "WS":
        # Watson-Sugar
        if is_cis:
            return (-8.0, -3.0, math.pi * 0.8)
        else:
            return (-8.0, 3.0, math.pi * 0.2)

    elif edges == "HH":
        # Hoogsteen-Hoogsteen
        if is_cis:
            return (-5.0, -6.0, math.pi * 0.5)
        else:
            return (-5.0, 6.0, -math.pi * 0.5)

    elif edges == "HS":
        # Hoogsteen-Sugar
        if is_cis:
            return (-6.0, -4.0, math.pi * 0.6)
        else:
            return (-6.0, 4.0, math.pi * 0.4)

    elif edges == "SS":
        # Sugar-Sugar
        if is_cis:
            return (-7.0, -5.0, math.pi * 0.7)
        else:
            return (-7.0, 5.0, math.pi * 0.3)

    else:
        # Default: standard WC-like geometry
        return (-9.0, -2.0, math.pi)


def calculate_wc_initial_position(
    base1_coords: Dict[str, np.ndarray],
    base2_template: Dict[str, np.ndarray],
    pattern: HBondPattern
) -> Tuple[float, float, float]:
    """
    Calculate a good initial position for Watson-Crick-like pairs based on H-bond atoms.

    Places base2 such that the first H-bond donor/acceptor are at ~3 Å distance.
    """
    if not pattern.hbonds:
        return (-9.0, -2.0, math.pi)

    # Use first H-bond to determine approximate positioning
    hb = pattern.hbonds[0]

    if hb.donor_base == 1:
        donor_atom = hb.donor_atom
        acceptor_atom = hb.acceptor_atom
        if donor_atom not in base1_coords or acceptor_atom not in base2_template:
            return (-9.0, -2.0, math.pi)

        # Target: place acceptor_atom of base2 near donor_atom of base1
        donor_pos = base1_coords[donor_atom]
        acceptor_template_pos = base2_template[acceptor_atom]

        # Translate base2 so acceptor is near donor (with ~3 Å separation along x)
        dx = donor_pos[0] - acceptor_template_pos[0] - 3.0
        dy = donor_pos[1] - acceptor_template_pos[1]

    else:
        donor_atom = hb.donor_atom
        acceptor_atom = hb.acceptor_atom
        if donor_atom not in base2_template or acceptor_atom not in base1_coords:
            return (-9.0, -2.0, math.pi)

        acceptor_pos = base1_coords[acceptor_atom]
        donor_template_pos = base2_template[donor_atom]

        dx = acceptor_pos[0] - donor_template_pos[0] - 3.0
        dy = acceptor_pos[1] - donor_template_pos[1]

    return (dx, dy, math.pi)


def optimize_base2_position(
    base1_coords: Dict[str, np.ndarray],
    base2_template: Dict[str, np.ndarray],
    pattern: HBondPattern,
    initial_guess: Tuple[float, float, float],
    is_cis: bool = True
) -> Tuple[float, float, float, float, bool]:
    """
    Optimize base2 position to minimize H-bond distance deviations.

    Uses grid search over multiple starting positions to avoid local minima.
    Also tries with and without y-axis flip.

    Returns:
        (x, y, theta, final_score, flip_y)
    """
    # Bounds: x,y in reasonable range, theta full rotation
    bounds = [
        (-15, 5),    # x
        (-10, 10),   # y
        (-math.pi, 2 * math.pi)  # theta
    ]

    best_result = None
    best_score = float('inf')
    best_flip = False

    # Get intelligent starting point based on H-bond geometry
    wc_init = calculate_wc_initial_position(base1_coords, base2_template, pattern)

    # Grid search over starting positions - reduced for efficiency
    x_starts = [wc_init[0], -8, -4]
    y_starts = [wc_init[1], -2, 2]
    theta_starts = [wc_init[2], 0, math.pi, -math.pi/2]

    # Try both with and without flip
    for flip_y in [False, True]:
        for x_start in x_starts:
            for y_start in y_starts:
                for theta_start in theta_starts:
                    x0 = np.array([x_start, y_start, theta_start])

                    try:
                        result = minimize(
                            calculate_hbond_score,
                            x0,
                            args=(base1_coords, base2_template, pattern, is_cis, flip_y),
                            method='L-BFGS-B',
                            bounds=bounds,
                            options={'maxiter': 1000, 'ftol': 1e-8}
                        )

                        if result.fun < best_score:
                            best_score = result.fun
                            best_result = result
                            best_flip = flip_y
                    except Exception:
                        continue

    if best_result is None:
        return initial_guess[0], initial_guess[1], initial_guess[2], float('inf'), False

    return best_result.x[0], best_result.x[1], best_result.x[2], best_result.fun, best_flip


def write_idealized_pdb(
    output_path: Path,
    base1_coords: Dict[str, np.ndarray],
    base2_coords: Dict[str, np.ndarray],
    entry: ExemplarEntry,
    pattern: HBondPattern,
    hbond_geometries: List[HBondGeometry]
) -> None:
    """Write idealized base pair to PDB file with H-bond geometry info."""
    with open(output_path, 'w') as f:
        f.write(f"REMARK   1 Idealized base pair: {entry.bp_key}\n")
        f.write(f"REMARK   2 LW classification: {entry.lw_class}\n")
        f.write(f"REMARK   3 Sequence: {entry.sequence}\n")
        f.write(f"REMARK   4 Based on exemplar: {entry.pdb_id}\n")
        f.write(f"REMARK   5 All atoms in z=0 plane (perfectly planar)\n")
        f.write(f"REMARK   6 H-bonds optimized (dist Å, X-D...A°, D...A-Y°):\n")

        for geom in hbond_geometries:
            donor_ang_str = f"{geom.donor_angle:.1f}" if geom.donor_angle else "N/A"
            acc_ang_str = f"{geom.acceptor_angle:.1f}" if geom.acceptor_angle else "N/A"
            f.write(f"REMARK   7   {geom.donor_atom}-{geom.acceptor_atom}: "
                   f"d={geom.distance:.2f}({geom.ideal_dist:.2f}) "
                   f"ang={donor_ang_str}({geom.ideal_donor_angle:.0f})/"
                   f"{acc_ang_str}({geom.ideal_acceptor_angle:.0f})\n")

        # Write base1 atoms
        atom_num = 1
        base1_type = entry.sequence[0].upper()
        for atom_name, coord in sorted(base1_coords.items()):
            line = f"ATOM  {atom_num:5d}  {atom_name:<3s} {base1_type:>3s} A   1    "
            line += f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
            line += f"  1.00  0.00           {atom_name[0]:>2s}\n"
            f.write(line)
            atom_num += 1

        # Write base2 atoms
        base2_type = entry.sequence[1].upper() if len(entry.sequence) > 1 else base1_type
        for atom_name, coord in sorted(base2_coords.items()):
            line = f"ATOM  {atom_num:5d}  {atom_name:<3s} {base2_type:>3s} B   2    "
            line += f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
            line += f"  1.00  0.00           {atom_name[0]:>2s}\n"
            f.write(line)
            atom_num += 1

        f.write("END\n")


def extract_base_atoms_from_pdb(pdb_path: Path) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    Extract base atoms from exemplar PDB.

    Uses order of residue appearance in the file (first residue = base1).
    Returns full 3D coordinates (not projected to z=0).

    Returns (base1_coords, base2_coords).
    """
    base_atom_names = {
        "N1", "N2", "N3", "N4", "N6", "N7", "N9",
        "C2", "C4", "C5", "C6", "C8",
        "O2", "O4", "O6",
        "C1'"  # Include C1' for reference
    }

    base1_coords = {}
    base2_coords = {}
    first_residue_key = None

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            name = line[12:16].strip()

            if name not in base_atom_names:
                continue

            # Get residue key (chain + resnum)
            chain_id = line[21] if len(line) > 21 else 'A'
            res_num = line[22:26].strip()
            residue_key = f"{chain_id}_{res_num}"

            # Track first residue seen
            if first_residue_key is None:
                first_residue_key = residue_key

            # Parse coordinates
            parts = line[30:].split()
            if len(parts) < 3:
                continue
            try:
                x = float(parts[0])
                y = float(parts[1])
                z = float(parts[2])
            except (ValueError, IndexError):
                continue

            # First residue in file = base1, second = base2
            if residue_key == first_residue_key:
                base1_coords[name] = np.array([x, y, z])
            else:
                base2_coords[name] = np.array([x, y, z])

    return base1_coords, base2_coords


class IdealizedBasePairGenerator:
    """Main generator for idealized base pairs."""

    def __init__(
        self,
        template_dir: Path,
        pattern_file: Path,
        exemplar_dir: Path,
        output_dir: Path
    ):
        self.templates = BaseTemplates(template_dir)
        self.pattern_db = HBondPatternDB(pattern_file)
        self.exemplar_dir = exemplar_dir
        self.output_dir = output_dir

    def generate(self, entry: ExemplarEntry, verbose: bool = False) -> Optional[Path]:
        """Generate idealized structure for a single entry."""
        # Get H-bond pattern
        pattern = self.pattern_db.get_pattern(entry.lw_class, entry.sequence)
        if pattern is None:
            if verbose:
                print(f"  No H-bond pattern for {entry.lw_class}/{entry.sequence}")
            return None

        safe_seq = encode_sequence_for_filename(entry.sequence)
        exemplar_path = self.exemplar_dir / entry.lw_class / f"{safe_seq}.pdb"

        # Always use template-based optimization for best H-bond geometry
        # (Exemplar projection to z=0 doesn't preserve H-bond distances since
        # WC pairs have z-separation between bases, not xy-separation)
        base1_type = entry.sequence[0].upper()
        base2_type = entry.sequence[1].upper() if len(entry.sequence) > 1 else base1_type

        base1_coords = get_base_atoms(self.templates, base1_type)
        base2_template = get_base_atoms(self.templates, base2_type)

        if not base1_coords or not base2_template:
            if verbose:
                print(f"  Missing template for {entry.sequence}")
            return None

        # Get initial guess based on LW class geometry
        initial_guess = get_initial_guess_for_basepair(
            entry.lw_class, base1_type, base2_type
        )

        # Determine if cis or trans
        is_cis = entry.lw_class.startswith('c')

        # Optimize base2 position (also tries with/without y-flip)
        x, y, theta, score, flip_y = optimize_base2_position(
            base1_coords, base2_template, pattern, initial_guess, is_cis
        )

        # Transform base2
        base2_coords = transform_coords(base2_template, x, y, theta, flip_y=flip_y)

        # Calculate final H-bond geometry (distances and angles)
        hbond_geometries = calculate_hbond_geometry(base1_coords, base2_coords, pattern)

        # Create output directory
        lw_dir = self.output_dir / entry.lw_class
        lw_dir.mkdir(parents=True, exist_ok=True)

        # Write output
        output_file = lw_dir / f"{safe_seq}.pdb"
        write_idealized_pdb(output_file, base1_coords, base2_coords, entry, pattern, hbond_geometries)

        if verbose:
            avg_dist_dev = np.mean([abs(g.distance - g.ideal_dist) for g in hbond_geometries]) if hbond_geometries else 0
            avg_donor_ang_dev = np.mean([abs(g.donor_angle - g.ideal_donor_angle)
                                          for g in hbond_geometries if g.donor_angle]) if hbond_geometries else 0
            avg_acc_ang_dev = np.mean([abs(g.acceptor_angle - g.ideal_acceptor_angle)
                                        for g in hbond_geometries if g.acceptor_angle]) if hbond_geometries else 0
            print(f"  {entry.lw_class}/{entry.sequence}: score={score:.4f}, "
                  f"dist_dev={avg_dist_dev:.3f}Å, ang_dev={avg_donor_ang_dev:.1f}°/{avg_acc_ang_dev:.1f}°")

        return output_file


def compare_idealized_to_exemplar(
    idealized_path: Path,
    exemplar_path: Path,
    pattern: HBondPattern,
    verbose: bool = False
) -> Optional[Dict]:
    """
    Compare H-bond geometry between idealized and exemplar PDBs.

    Returns comparison dict or None if comparison failed.
    """
    if not idealized_path.exists() or not exemplar_path.exists():
        return None

    # Extract coordinates from both PDBs
    ideal_base1, ideal_base2 = extract_base_atoms_from_pdb(idealized_path)
    exem_base1, exem_base2 = extract_base_atoms_from_pdb(exemplar_path)

    if not ideal_base1 or not ideal_base2 or not exem_base1 or not exem_base2:
        return None

    # Calculate geometry for both
    ideal_geom = calculate_hbond_geometry(ideal_base1, ideal_base2, pattern)
    exem_geom = calculate_hbond_geometry(exem_base1, exem_base2, pattern)

    if not ideal_geom or not exem_geom:
        return None

    comparison = {
        'lw_class': pattern.lw_class,
        'sequence': pattern.sequence,
        'idealized_path': str(idealized_path),
        'exemplar_path': str(exemplar_path),
        'hbonds': []
    }

    for i, (ig, eg) in enumerate(zip(ideal_geom, exem_geom)):
        hb_compare = {
            'donor_atom': ig.donor_atom,
            'acceptor_atom': ig.acceptor_atom,
            'idealized': {
                'distance': round(ig.distance, 3),
                'donor_angle': round(ig.donor_angle, 1) if ig.donor_angle else None,
                'acceptor_angle': round(ig.acceptor_angle, 1) if ig.acceptor_angle else None,
            },
            'exemplar': {
                'distance': round(eg.distance, 3),
                'donor_angle': round(eg.donor_angle, 1) if eg.donor_angle else None,
                'acceptor_angle': round(eg.acceptor_angle, 1) if eg.acceptor_angle else None,
            },
            'ideal_values': {
                'distance': ig.ideal_dist,
                'donor_angle': ig.ideal_donor_angle,
                'acceptor_angle': ig.ideal_acceptor_angle,
            },
            'delta': {
                'distance': round(ig.distance - eg.distance, 3),
                'donor_angle': round(ig.donor_angle - eg.donor_angle, 1)
                              if ig.donor_angle and eg.donor_angle else None,
                'acceptor_angle': round(ig.acceptor_angle - eg.acceptor_angle, 1)
                                  if ig.acceptor_angle and eg.acceptor_angle else None,
            }
        }
        comparison['hbonds'].append(hb_compare)

    # Calculate summary statistics
    dist_devs_ideal = [abs(g.distance - g.ideal_dist) for g in ideal_geom]
    dist_devs_exem = [abs(g.distance - g.ideal_dist) for g in exem_geom]

    ang_devs_ideal = []
    ang_devs_exem = []
    for g in ideal_geom:
        if g.donor_angle:
            ang_devs_ideal.append(abs(g.donor_angle - g.ideal_donor_angle))
        if g.acceptor_angle:
            ang_devs_ideal.append(abs(g.acceptor_angle - g.ideal_acceptor_angle))
    for g in exem_geom:
        if g.donor_angle:
            ang_devs_exem.append(abs(g.donor_angle - g.ideal_donor_angle))
        if g.acceptor_angle:
            ang_devs_exem.append(abs(g.acceptor_angle - g.ideal_acceptor_angle))

    comparison['summary'] = {
        'idealized': {
            'avg_dist_dev': round(np.mean(dist_devs_ideal), 3) if dist_devs_ideal else 0,
            'max_dist_dev': round(max(dist_devs_ideal), 3) if dist_devs_ideal else 0,
            'avg_angle_dev': round(np.mean(ang_devs_ideal), 1) if ang_devs_ideal else 0,
            'max_angle_dev': round(max(ang_devs_ideal), 1) if ang_devs_ideal else 0,
        },
        'exemplar': {
            'avg_dist_dev': round(np.mean(dist_devs_exem), 3) if dist_devs_exem else 0,
            'max_dist_dev': round(max(dist_devs_exem), 3) if dist_devs_exem else 0,
            'avg_angle_dev': round(np.mean(ang_devs_exem), 1) if ang_devs_exem else 0,
            'max_angle_dev': round(max(ang_devs_exem), 1) if ang_devs_exem else 0,
        }
    }

    return comparison


def run_comparison(
    pattern_db: HBondPatternDB,
    idealized_dir: Path,
    exemplar_dir: Path,
    entries: Dict[str, List[ExemplarEntry]],
    lw_class_filter: Optional[str] = None,
    verbose: bool = False
) -> List[Dict]:
    """Run comparison between idealized and exemplar base pairs."""
    comparisons = []

    for lw_class, entry_list in entries.items():
        if lw_class_filter and lw_class != lw_class_filter:
            continue

        # Get unique sequences
        seen_seqs = set()
        for entry in entry_list:
            seq_upper = entry.sequence.upper()
            if seq_upper in seen_seqs:
                continue
            seen_seqs.add(seq_upper)

            pattern = pattern_db.get_pattern(lw_class, seq_upper)
            if pattern is None:
                continue

            safe_seq = encode_sequence_for_filename(entry.sequence)
            idealized_path = idealized_dir / lw_class / f"{safe_seq}.pdb"
            exemplar_path = exemplar_dir / lw_class / f"{safe_seq}.pdb"

            result = compare_idealized_to_exemplar(
                idealized_path, exemplar_path, pattern, verbose
            )

            if result:
                comparisons.append(result)
                if verbose:
                    summary = result['summary']
                    print(f"  {lw_class}/{seq_upper}: "
                          f"ideal dist_dev={summary['idealized']['avg_dist_dev']:.3f}Å "
                          f"ang_dev={summary['idealized']['avg_angle_dev']:.1f}° | "
                          f"exemplar dist_dev={summary['exemplar']['avg_dist_dev']:.3f}Å "
                          f"ang_dev={summary['exemplar']['avg_angle_dev']:.1f}°")

    return comparisons


def main():
    parser = argparse.ArgumentParser(
        description="Generate idealized base pair structures and compare to exemplars"
    )
    parser.add_argument(
        "--template-dir",
        default="resources/templates",
        help="Directory containing base templates"
    )
    parser.add_argument(
        "--pattern-file",
        default="tools/hbond_patterns.json",
        help="H-bond pattern definitions file"
    )
    parser.add_argument(
        "--exemplar-dir",
        default="basepair-catalog-exemplars",
        help="Directory containing exemplar PDBs"
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="basepair-idealized",
        help="Output directory for idealized PDBs"
    )
    parser.add_argument(
        "--database", "-d",
        default="data/bp_database.txt",
        help="BP database file"
    )
    parser.add_argument(
        "--lw-class",
        help="Only process specific LW class"
    )
    parser.add_argument(
        "--skip-modeled",
        action="store_true",
        help="Skip modeled entries"
    )
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Compare idealized structures to exemplars (requires both to exist)"
    )
    parser.add_argument(
        "--compare-output",
        default="basepair-comparison.json",
        help="Output file for comparison results (default: basepair-comparison.json)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )

    args = parser.parse_args()

    project_root = Path(__file__).parent.parent
    template_dir = project_root / args.template_dir
    pattern_file = project_root / args.pattern_file
    exemplar_dir = project_root / args.exemplar_dir
    output_dir = project_root / args.output_dir
    db_path = project_root / args.database

    # Load database
    print(f"Loading database from {db_path}...")
    entries = load_database(db_path)

    # Filter
    if args.lw_class:
        if args.lw_class not in entries:
            print(f"Error: LW class '{args.lw_class}' not found")
            return 1
        entries = {args.lw_class: entries[args.lw_class]}

    # Collect unique entries (one per sequence per LW class)
    unique_entries: Dict[str, ExemplarEntry] = {}
    for lw_class, entry_list in entries.items():
        for entry in entry_list:
            if args.skip_modeled and entry.is_modeled:
                continue
            key = f"{lw_class}/{entry.sequence.upper()}"
            if key not in unique_entries:
                unique_entries[key] = entry

    print(f"Processing {len(unique_entries)} unique base pair types...")

    # Create generator
    generator = IdealizedBasePairGenerator(
        template_dir=template_dir,
        pattern_file=pattern_file,
        exemplar_dir=exemplar_dir,
        output_dir=output_dir
    )

    # Generate
    output_dir.mkdir(parents=True, exist_ok=True)
    generated = 0
    failed = 0
    summary = []

    for key, entry in sorted(unique_entries.items()):
        result = generator.generate(entry, verbose=args.verbose)
        if result:
            generated += 1
            summary.append({
                'lw_class': entry.lw_class,
                'sequence': entry.sequence,
                'output_file': str(result.relative_to(output_dir))
            })
        else:
            failed += 1

    # Write summary
    summary_file = output_dir / "summary.json"
    with open(summary_file, 'w') as f:
        json.dump({
            'generated': generated,
            'failed': failed,
            'idealized_pairs': summary
        }, f, indent=2)

    print(f"\nGeneration complete!")
    print(f"  Generated: {generated}")
    print(f"  Failed: {failed} (no H-bond pattern)")
    print(f"  Output: {output_dir}")

    # Run comparison mode if requested
    if args.compare:
        print(f"\nRunning comparison with exemplars...")
        pattern_db = HBondPatternDB(pattern_file)

        comparisons = run_comparison(
            pattern_db=pattern_db,
            idealized_dir=output_dir,
            exemplar_dir=exemplar_dir,
            entries=entries,
            lw_class_filter=args.lw_class,
            verbose=args.verbose
        )

        # Calculate overall statistics
        if comparisons:
            all_ideal_dist_devs = []
            all_exem_dist_devs = []
            all_ideal_ang_devs = []
            all_exem_ang_devs = []

            for comp in comparisons:
                all_ideal_dist_devs.append(comp['summary']['idealized']['avg_dist_dev'])
                all_exem_dist_devs.append(comp['summary']['exemplar']['avg_dist_dev'])
                all_ideal_ang_devs.append(comp['summary']['idealized']['avg_angle_dev'])
                all_exem_ang_devs.append(comp['summary']['exemplar']['avg_angle_dev'])

            overall_summary = {
                'total_comparisons': len(comparisons),
                'idealized': {
                    'mean_dist_dev': round(np.mean(all_ideal_dist_devs), 3),
                    'mean_angle_dev': round(np.mean(all_ideal_ang_devs), 1),
                },
                'exemplar': {
                    'mean_dist_dev': round(np.mean(all_exem_dist_devs), 3),
                    'mean_angle_dev': round(np.mean(all_exem_ang_devs), 1),
                }
            }

            # Write comparison output
            compare_output_path = project_root / args.compare_output
            with open(compare_output_path, 'w') as f:
                json.dump({
                    'overall_summary': overall_summary,
                    'comparisons': comparisons
                }, f, indent=2)

            print(f"\nComparison complete!")
            print(f"  Compared: {len(comparisons)} base pair types")
            print(f"  Idealized: mean dist_dev={overall_summary['idealized']['mean_dist_dev']:.3f}Å, "
                  f"mean ang_dev={overall_summary['idealized']['mean_angle_dev']:.1f}°")
            print(f"  Exemplar:  mean dist_dev={overall_summary['exemplar']['mean_dist_dev']:.3f}Å, "
                  f"mean ang_dev={overall_summary['exemplar']['mean_angle_dev']:.1f}°")
            print(f"  Output: {compare_output_path}")
        else:
            print("  No comparisons could be made (missing idealized or exemplar files)")

    return 0


if __name__ == "__main__":
    sys.exit(main())

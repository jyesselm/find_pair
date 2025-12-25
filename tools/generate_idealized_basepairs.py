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
class HBondDef:
    """Single H-bond definition."""
    donor_base: int      # 1 or 2
    donor_atom: str      # Atom name on donor base
    acceptor_atom: str   # Atom name on acceptor base
    ideal_dist: float    # Ideal distance in Angstroms


@dataclass
class HBondPattern:
    """H-bond pattern for a specific LW class + sequence."""
    lw_class: str
    sequence: str
    hbonds: List[HBondDef]


class HBondPatternDB:
    """Database of H-bond patterns."""

    def __init__(self, pattern_file: Path):
        with open(pattern_file) as f:
            self.data = json.load(f)
        self.patterns = self.data.get("patterns", {})

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
            hbonds.append(HBondDef(
                donor_base=hb_data["donor_base"],
                donor_atom=hb_data["donor_atom"],
                acceptor_atom=hb_data["acceptor_atom"],
                ideal_dist=hb_data.get("ideal_dist", 2.9)
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


def calculate_hbond_score(
    params: np.ndarray,
    base1_coords: Dict[str, np.ndarray],
    base2_template: Dict[str, np.ndarray],
    pattern: HBondPattern,
    is_cis: bool = True,
    flip_y: bool = False
) -> float:
    """
    Objective function: H-bond distances + clash avoidance + geometry constraints.

    Args:
        params: [x, y, theta] - position/rotation of base2
        base1_coords: Fixed base1 coordinates
        base2_template: Template coordinates for base2 (will be transformed)
        pattern: H-bond pattern to optimize
        is_cis: True for cis pairs (same side), False for trans
        flip_y: Whether to mirror base2 in y-axis

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

    # Extra penalty for having any single H-bond very far from ideal
    if max_hbond_dev > 0.5:
        score += 20.0 * (max_hbond_dev - 0.5) ** 2

    # 2. Clash avoidance - penalize non-H-bond atoms that are too close
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

    # 3. C1'-C1' distance constraint (typical WC pair: ~10.5 Å)
    # Use lower weight so it doesn't prevent good H-bond geometry
    if "C1'" in base1_coords and "C1'" in base2_coords:
        c1_dist = np.linalg.norm(base1_coords["C1'"] - base2_coords["C1'"])
        ideal_c1_dist = 10.5
        # Only penalize if too close (allow larger distances for planar geometry)
        if c1_dist < ideal_c1_dist:
            score += 0.5 * (c1_dist - ideal_c1_dist) ** 2

    # 4. Planarity bonus - bases should stay in z=0 plane (already enforced by transform)

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
    hbond_distances: List[Tuple[float, float]]
) -> None:
    """Write idealized base pair to PDB file."""
    with open(output_path, 'w') as f:
        f.write(f"REMARK   1 Idealized base pair: {entry.bp_key}\n")
        f.write(f"REMARK   2 LW classification: {entry.lw_class}\n")
        f.write(f"REMARK   3 Sequence: {entry.sequence}\n")
        f.write(f"REMARK   4 Based on exemplar: {entry.pdb_id}\n")
        f.write(f"REMARK   5 All atoms in z=0 plane (perfectly planar)\n")
        f.write(f"REMARK   6 H-bonds optimized:\n")

        for i, (actual, ideal) in enumerate(hbond_distances):
            hb = pattern.hbonds[i] if i < len(pattern.hbonds) else None
            if hb:
                f.write(f"REMARK   7   {hb.donor_atom}-{hb.acceptor_atom}: "
                       f"{actual:.2f} A (ideal: {ideal:.2f})\n")

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

        # Calculate final H-bond distances
        hbond_distances = calculate_hbond_distances(base1_coords, base2_coords, pattern)

        # Create output directory
        lw_dir = self.output_dir / entry.lw_class
        lw_dir.mkdir(parents=True, exist_ok=True)

        # Write output
        output_file = lw_dir / f"{safe_seq}.pdb"
        write_idealized_pdb(output_file, base1_coords, base2_coords, entry, pattern, hbond_distances)

        if verbose:
            avg_dev = np.mean([abs(a - i) for a, i in hbond_distances]) if hbond_distances else 0
            print(f"  {entry.lw_class}/{entry.sequence}: score={score:.4f}, "
                  f"avg_dev={avg_dev:.3f} A, {len(hbond_distances)} H-bonds")

        return output_file


def main():
    parser = argparse.ArgumentParser(
        description="Generate idealized base pair structures"
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

    print(f"\nComplete!")
    print(f"  Generated: {generated}")
    print(f"  Failed: {failed} (no H-bond pattern)")
    print(f"  Output: {output_dir}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

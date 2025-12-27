#!/usr/bin/env python3
"""Generate PyMOL visualization of template alignments."""

import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple
import sys

sys.path.insert(0, str(Path(__file__).parent))
from template_aligner import parse_pdb_residues, Residue


def load_template_atoms(template_path: Path) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], str, str]:
    """Load template PDB and return atom coords + residue names."""
    res1_atoms = {}
    res2_atoms = {}
    res1_name = ""
    res2_name = ""

    with open(template_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            res_seq = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            if res_seq == 1:
                res1_atoms[atom_name] = np.array([x, y, z])
                res1_name = res_name
            else:
                res2_atoms[atom_name] = np.array([x, y, z])
                res2_name = res_name

    return res1_atoms, res2_atoms, res1_name, res2_name


def kabsch_align(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute optimal rotation and translation to align P onto Q."""
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    H = P_centered.T @ Q_centered
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    t = centroid_Q - R @ centroid_P
    return R, t


def write_pdb(atoms: Dict[str, np.ndarray], res_name: str, chain: str, res_num: int, output_path: Path):
    """Write atoms to PDB format."""
    with open(output_path, 'w') as f:
        atom_num = 1
        for atom_name, coords in atoms.items():
            f.write(f"ATOM  {atom_num:5d}  {atom_name:<3s} {res_name:3s} {chain}{res_num:4d}    "
                    f"{coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}  1.00  0.00\n")
            atom_num += 1
        f.write("END\n")


def write_pair_pdb(res1_atoms: Dict[str, np.ndarray], res2_atoms: Dict[str, np.ndarray],
                   res1_name: str, res2_name: str, chain: str, output_path: Path):
    """Write a base pair to PDB format."""
    with open(output_path, 'w') as f:
        atom_num = 1
        for atom_name, coords in res1_atoms.items():
            f.write(f"ATOM  {atom_num:5d}  {atom_name:<3s} {res1_name:3s} {chain}   1    "
                    f"{coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}  1.00  0.00\n")
            atom_num += 1
        for atom_name, coords in res2_atoms.items():
            f.write(f"ATOM  {atom_num:5d}  {atom_name:<3s} {res2_name:3s} {chain}   2    "
                    f"{coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}  1.00  0.00\n")
            atom_num += 1
        f.write("END\n")


def align_template_to_target(
    target_res1: Residue,
    target_res2: Residue,
    template_path: Path,
    ring_atoms: List[str]
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], str, str, float]:
    """Align template to target and return transformed coordinates."""
    template_res1, template_res2, res1_name, res2_name = load_template_atoms(template_path)

    # Collect matching ring atoms for alignment
    template_points = []
    target_points = []

    for atom_name in ring_atoms:
        if atom_name in template_res1 and atom_name in target_res1.atoms:
            template_points.append(template_res1[atom_name])
            target_points.append(target_res1.atoms[atom_name].coords)
        if atom_name in template_res2 and atom_name in target_res2.atoms:
            template_points.append(template_res2[atom_name])
            target_points.append(target_res2.atoms[atom_name].coords)

    if len(template_points) < 4:
        return None, None, None, None, float('inf')

    template_points = np.array(template_points)
    target_points = np.array(target_points)

    # Compute alignment (align template TO target)
    # kabsch_align returns R, t such that R @ template + t ≈ target
    centroid_template = np.mean(template_points, axis=0)
    centroid_target = np.mean(target_points, axis=0)

    T_centered = template_points - centroid_template
    P_centered = target_points - centroid_target

    H = T_centered.T @ P_centered
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # Transform ALL template atoms (not just ring atoms)
    aligned_res1 = {}
    for atom_name, coords in template_res1.items():
        aligned_res1[atom_name] = R @ (coords - centroid_template) + centroid_target

    aligned_res2 = {}
    for atom_name, coords in template_res2.items():
        aligned_res2[atom_name] = R @ (coords - centroid_template) + centroid_target

    # Compute RMSD on ring atoms
    aligned_ring = np.array([R @ (p - centroid_template) + centroid_target for p in template_points])
    rmsd = float(np.sqrt(np.mean(np.sum((aligned_ring - target_points)**2, axis=1))))

    return aligned_res1, aligned_res2, res1_name, res2_name, rmsd


def generate_visualization(
    pdb_id: str,
    res_id1: str,
    res_id2: str,
    pdb_dir: Path,
    idealized_dir: Path,
    output_dir: Path,
):
    """Generate PyMOL visualization files."""

    pdb_path = pdb_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"

    if not pdb_path.exists():
        print(f"PDB file not found: {pdb_path}")
        return

    residues = parse_pdb_residues(pdb_path)

    if res_id1 not in residues or res_id2 not in residues:
        print(f"Residue not found: {res_id1} or {res_id2}")
        return

    res1 = residues[res_id1]
    res2 = residues[res_id2]
    seq = res1.base_type + res2.base_type

    output_dir.mkdir(parents=True, exist_ok=True)

    ring_atoms = ["C2", "C4", "C5", "C6", "N1", "N3", "N7", "C8", "N9"]

    # Extract target pair atoms
    target_res1_atoms = {name: atom.coords for name, atom in res1.atoms.items()}
    target_res2_atoms = {name: atom.coords for name, atom in res2.atoms.items()}

    # Write target pair
    target_pdb = output_dir / f"target_{pdb_id}_{seq}.pdb"
    write_pair_pdb(target_res1_atoms, target_res2_atoms,
                   res1.base_type, res2.base_type, "T", target_pdb)
    print(f"Wrote target: {target_pdb}")

    # Align templates
    templates = {
        "cWW": idealized_dir / "cWW" / f"{seq}.pdb",
        "tWW": idealized_dir / "tWW" / f"{seq}.pdb",
        "cWS": idealized_dir / "cWS" / f"{seq}.pdb",
    }

    aligned_files = {"target": target_pdb}
    rmsd_values = {}

    for lw_class, template_path in templates.items():
        if not template_path.exists():
            # Try alternate naming
            template_path = idealized_dir / lw_class / f"{seq[0]}_{seq[1].lower()}.pdb"

        if not template_path.exists():
            print(f"No {lw_class} template for {seq}")
            continue

        aligned_res1, aligned_res2, res1_name, res2_name, rmsd = align_template_to_target(
            res1, res2, template_path, ring_atoms
        )

        if aligned_res1 is None:
            continue

        # Write aligned template (include pdb_id to avoid overwrites)
        aligned_pdb = output_dir / f"aligned_{pdb_id}_{lw_class}_{seq}.pdb"
        write_pair_pdb(aligned_res1, aligned_res2, res1_name, res2_name,
                       lw_class[0].upper(), aligned_pdb)
        print(f"Wrote {lw_class} aligned (RMSD={rmsd:.3f}Å): {aligned_pdb}")

        aligned_files[lw_class] = aligned_pdb
        rmsd_values[lw_class] = rmsd

    # Generate PyMOL script
    pml_script = output_dir / f"view_{pdb_id}_{seq}.pml"
    with open(pml_script, 'w') as f:
        f.write(f"# PyMOL visualization of {pdb_id} {res_id1}-{res_id2} ({seq})\n")
        f.write(f"# Template alignment comparison\n\n")

        # Load structures
        f.write(f"load {target_pdb.absolute()}, target\n")

        colors = {"cWW": "green", "tWW": "orange", "cWS": "magenta"}

        for lw_class, pdb_file in aligned_files.items():
            if lw_class == "target":
                continue
            f.write(f"load {pdb_file.absolute()}, {lw_class}\n")

        f.write("\n# Style settings\n")
        f.write("hide everything\n")
        f.write("show sticks, all\n")
        f.write("set stick_radius, 0.15\n")

        # Color target
        f.write("\n# Color target (cyan)\n")
        f.write("color cyan, target\n")

        # Color templates
        f.write("\n# Color templates\n")
        for lw_class, color in colors.items():
            if lw_class in aligned_files:
                f.write(f"color {color}, {lw_class}\n")

        # Add labels
        f.write("\n# Labels\n")
        f.write("set label_size, 20\n")

        # Show H-bond atoms
        f.write("\n# Highlight key atoms\n")
        f.write("select hbond_atoms, name N1+N2+N3+N4+N6+N7+O2+O4+O6\n")
        f.write("show spheres, hbond_atoms and target\n")
        f.write("set sphere_scale, 0.2, hbond_atoms\n")

        # Set view
        f.write("\n# Center and zoom\n")
        f.write("center target\n")
        f.write("zoom target, 5\n")

        # Add title
        f.write("\n# Add title\n")
        rmsd_str = ", ".join([f"{lw}={rmsd_values.get(lw, 0):.2f}Å" for lw in ["cWW", "tWW", "cWS"] if lw in rmsd_values])
        f.write(f'set_name target, "Target {pdb_id} {seq}"\n')

        f.write("\n# Legend\n")
        f.write("# Target = cyan\n")
        for lw_class, color in colors.items():
            if lw_class in rmsd_values:
                f.write(f"# {lw_class} = {color} (RMSD={rmsd_values[lw_class]:.3f}Å)\n")

        f.write("\n# Toggle visibility\n")
        f.write("# Use: disable cWW / enable cWW to toggle templates\n")
        f.write("# All templates enabled by default\n")

    print(f"\nWrote PyMOL script: {pml_script}")
    print(f"\nTo visualize, run: pymol {pml_script}")

    # Also create a summary
    print(f"\n{'='*50}")
    print(f"RMSD Summary for {pdb_id} {res_id1}-{res_id2} ({seq}):")
    print(f"{'='*50}")
    for lw_class in ["cWW", "tWW", "cWS"]:
        if lw_class in rmsd_values:
            marker = " <-- BEST" if rmsd_values[lw_class] == min(rmsd_values.values()) else ""
            print(f"  {lw_class}: {rmsd_values[lw_class]:.4f}Å{marker}")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Generate alignment visualization")
    parser.add_argument("--pdb", type=str, required=True, help="PDB ID")
    parser.add_argument("--res1", type=str, required=True, help="First residue ID")
    parser.add_argument("--res2", type=str, required=True, help="Second residue ID")
    parser.add_argument("--pdb-dir", type=Path, default=Path("data/pdb"))
    parser.add_argument("--idealized-dir", type=Path, default=Path("basepair-idealized"))
    parser.add_argument("--output-dir", type=Path, default=Path("viz_output"))

    args = parser.parse_args()

    generate_visualization(
        args.pdb, args.res1, args.res2,
        args.pdb_dir, args.idealized_dir, args.output_dir
    )


if __name__ == "__main__":
    main()

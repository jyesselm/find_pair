"""Unified PDB parsing utilities.

Consolidates 4 different PDB parsing implementations across the codebase
into a single, canonical interface.
"""

from pathlib import Path
from typing import Dict, Optional, Tuple
import numpy as np

from .residue import Residue, normalize_base_type


def parse_pdb(path: Path) -> Dict[str, Residue]:
    """Parse PDB file into Residue objects keyed by res_id.

    Uses canonical res_id format: "chain-base-seq[ins]" (e.g., "A-G-1", "A-C-10A").

    Args:
        path: Path to PDB file

    Returns:
        Dict mapping res_id to Residue objects
    """
    residues: Dict[str, Residue] = {}

    with open(path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21:22].strip() or "A"
            res_seq = line[22:26].strip()
            ins_code = line[26:27].strip()

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            element = line[76:78].strip() if len(line) > 76 else ""

            # Build canonical res_id
            seq_with_ins = res_seq + ins_code if ins_code else res_seq
            base_type = normalize_base_type(res_name)
            res_id = f"{chain_id}-{base_type}-{seq_with_ins}"

            # Create or update residue
            if res_id not in residues:
                residues[res_id] = Residue(res_id=res_id, base_type=base_type)

            residues[res_id].add_atom(
                name=atom_name, coords=np.array([x, y, z]), element=element
            )

    return residues


def parse_template_pdb(
    path: Path,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Parse a template PDB file with exactly 2 residues.

    Templates have residue numbers 1 and 2. Returns raw coordinate
    dictionaries for efficient alignment operations.

    Args:
        path: Path to template PDB file

    Returns:
        Tuple of (res1_atoms, res2_atoms) where each is a dict
        mapping atom name to coordinates
    """
    res1_atoms: Dict[str, np.ndarray] = {}
    res2_atoms: Dict[str, np.ndarray] = {}

    with open(path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            atom_name = line[12:16].strip()
            res_seq = int(line[22:26])

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            coords = np.array([x, y, z])

            if res_seq == 1:
                res1_atoms[atom_name] = coords
            elif res_seq == 2:
                res2_atoms[atom_name] = coords

    return res1_atoms, res2_atoms


def write_pair_pdb(
    res1_atoms: Dict[str, np.ndarray],
    res2_atoms: Dict[str, np.ndarray],
    res1_name: str,
    res2_name: str,
    output_path: Path,
    chain: str = "A",
) -> None:
    """Write a base pair to PDB format.

    Args:
        res1_atoms: Dict mapping atom name to coords for residue 1
        res2_atoms: Dict mapping atom name to coords for residue 2
        res1_name: Residue name for residue 1 (e.g., "G")
        res2_name: Residue name for residue 2 (e.g., "C")
        output_path: Path to write PDB file
        chain: Chain identifier
    """
    atom_num = 1

    with open(output_path, "w") as f:
        # Write residue 1
        for atom_name in sorted(res1_atoms.keys()):
            coords = res1_atoms[atom_name]
            line = _format_atom_line(atom_num, atom_name, res1_name, chain, 1, coords)
            f.write(line + "\n")
            atom_num += 1

        # Write residue 2
        for atom_name in sorted(res2_atoms.keys()):
            coords = res2_atoms[atom_name]
            line = _format_atom_line(atom_num, atom_name, res2_name, chain, 2, coords)
            f.write(line + "\n")
            atom_num += 1

        f.write("END\n")


def _format_atom_line(
    atom_num: int,
    atom_name: str,
    res_name: str,
    chain: str,
    res_seq: int,
    coords: np.ndarray,
) -> str:
    """Format a PDB ATOM line following PDB format specification."""
    # Atom name formatting (columns 13-16)
    # Atoms with 1-3 char names are left-justified in columns 14-16
    # 4-char names fill columns 13-16
    if len(atom_name) < 4:
        atom_name_fmt = f" {atom_name:<3}"
    else:
        atom_name_fmt = atom_name[:4]

    # Element from atom name (usually first character)
    element = atom_name[0] if atom_name else ""

    return (
        f"ATOM  {atom_num:5d} {atom_name_fmt}{res_name:>3} {chain}{res_seq:4d}    "
        f"{coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}"
        f"  1.00  0.00          {element:>2}"
    )


def extract_pair_from_pdb(
    pdb_path: Path,
    res_id1: str,
    res_id2: str,
) -> Tuple[Optional[Residue], Optional[Residue]]:
    """Extract two specific residues from a PDB file.

    Args:
        pdb_path: Path to PDB file
        res_id1: First residue ID (e.g., "A-G-1")
        res_id2: Second residue ID (e.g., "A-C-72")

    Returns:
        Tuple of (residue1, residue2), either may be None if not found
    """
    residues = parse_pdb(pdb_path)
    return residues.get(res_id1), residues.get(res_id2)

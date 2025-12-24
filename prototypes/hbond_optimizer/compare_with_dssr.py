#!/usr/bin/env python3
"""
Compare H-bond optimizer output with DSSR output.

DSSR H-bond format: "atom1(type)-atom2(type)[distance]"
Example: "N6(amino)-O4(carbonyl)[2.86]"

Usage:
    python compare_with_dssr.py <pdb_file> [dssr_output]

If dssr_output not provided, runs x3dna-dssr on the PDB.
"""

import sys
import re
import subprocess
import tempfile
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass

sys.path.insert(0, str(Path(__file__).parent))

try:
    from .optimizer import HBondOptimizer, Residue, HBond
    from .geometry import DONOR_CAPACITY, ACCEPTOR_CAPACITY
except ImportError:
    from optimizer import HBondOptimizer, Residue, HBond
    from geometry import DONOR_CAPACITY, ACCEPTOR_CAPACITY


@dataclass
class DSSRHBond:
    """An H-bond from DSSR output."""
    res1: str          # e.g., "A.G1"
    res2: str          # e.g., "B.C2"
    donor_atom: str    # e.g., "N6"
    acceptor_atom: str # e.g., "O4"
    distance: float


def parse_dssr_output(dssr_text: str, exclude_questionable: bool = True) -> List[DSSRHBond]:
    """
    Parse DSSR --get-hbond output to extract H-bonds.

    Format:
       11   485  #1     p    2.776 O:N O2@A.DC1 N2@B.DG24

    Fields: atom1_idx, atom2_idx, #N, type, distance, donor_type:acc_type, atom1@res1, atom2@res2

    Types:
        p = pair (base pair H-bond)
        o = other (non-pair, e.g., ribose-ribose)
        x = questionable (e.g., OP1-OP2)

    Args:
        exclude_questionable: If True, exclude 'x' type H-bonds
    """
    hbonds = []

    # Pattern for H-bond lines from --get-hbond
    # Format: idx1 idx2 #N type dist D:A atom1@res1 atom2@res2
    hbond_pattern = re.compile(
        r'^\s*\d+\s+\d+\s+#\d+\s+(\w+)\s+(\d+\.\d+)\s+\w:\w\s+(\w+)@(\S+)\s+(\w+)@(\S+)',
        re.MULTILINE
    )

    for match in hbond_pattern.finditer(dssr_text):
        hb_type = match.group(1)  # 'p', 'o', or 'x'
        dist = float(match.group(2))
        atom1 = match.group(3)  # e.g., "O2"
        res1 = match.group(4)   # e.g., "A.DC1"
        atom2 = match.group(5)  # e.g., "N2"
        res2 = match.group(6)   # e.g., "B.DG24"

        # Skip questionable H-bonds if requested
        if exclude_questionable and hb_type == 'x':
            continue

        # Normalize residue IDs (remove 'D' prefix for DNA)
        res1 = normalize_res_id(res1)
        res2 = normalize_res_id(res2)

        hbonds.append(DSSRHBond(
            res1=res1,
            res2=res2,
            donor_atom=atom1,
            acceptor_atom=atom2,
            distance=dist
        ))

    return hbonds


def normalize_res_id(res_id: str) -> str:
    """
    Normalize DSSR residue ID to standard format.

    E.g., "A.DC1" -> "A.C1", "B.DG24" -> "B.G24"
    """
    # Handle DNA prefix (D before base letter)
    parts = res_id.split('.')
    if len(parts) == 2:
        chain = parts[0]
        res = parts[1]
        # Remove 'D' prefix for DNA bases
        if len(res) >= 2 and res[0] == 'D' and res[1] in 'ACGTU':
            res = res[1:]
        return f"{chain}.{res}"
    return res_id


def parse_dssr_base_pairs(dssr_text: str) -> List[Tuple[str, str, str]]:
    """
    Parse base pairs from DSSR output.

    Returns list of (res1, res2, pair_type) tuples.
    Example: ("A.G1", "B.C2", "G-C")
    """
    pairs = []

    # Look for "List of N base pair" section
    in_bp_section = False

    for line in dssr_text.split('\n'):
        if 'List of' in line and 'base pair' in line:
            in_bp_section = True
            continue
        if in_bp_section and line.startswith('***'):
            in_bp_section = False
            continue

        if in_bp_section:
            # Format: "   1 A.G1             B.C2             G-C WC..."
            match = re.match(r'\s+\d+\s+(\S+)\s+(\S+)\s+(\S+-\S+)', line)
            if match:
                pairs.append((match.group(1), match.group(2), match.group(3)))

    return pairs


DSSR_PATH = "/Users/jyesselman2/Library/CloudStorage/Dropbox/2_code/python/py_dssr/pydssr/resources/dssr/osx/x3dna-dssr"


def run_dssr(pdb_path: Path) -> str:
    """Run DSSR on a PDB file and return output."""
    dssr_cmd = DSSR_PATH if Path(DSSR_PATH).exists() else 'x3dna-dssr'

    try:
        result = subprocess.run(
            [dssr_cmd, '-i=' + str(pdb_path), '--non-pair', '--get-hbond'],
            capture_output=True,
            text=True,
            timeout=120
        )
        if result.returncode != 0:
            print(f"DSSR stderr: {result.stderr}")
        return result.stdout
    except FileNotFoundError:
        print(f"Error: DSSR not found at {dssr_cmd}")
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print("Error: DSSR timed out")
        sys.exit(1)


def parse_pdb_residues(pdb_path: Path) -> Dict[str, Residue]:
    """Parse PDB and return residues keyed by DSSR-style ID (e.g., "A.G1")."""
    residues = {}

    RES_MAP = {
        'ADE': 'A', 'A': 'A', 'DA': 'A', 'RA': 'A',
        'GUA': 'G', 'G': 'G', 'DG': 'G', 'RG': 'G',
        'CYT': 'C', 'C': 'C', 'DC': 'C', 'RC': 'C',
        'URA': 'U', 'U': 'U', 'RU': 'U',
        'THY': 'T', 'T': 'T', 'DT': 'T',
    }

    # Atoms we track (N and O only, no C-H)
    TRACKED_ATOMS = {
        # Base N atoms
        'N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9',
        # Base O atoms
        'O2', 'O4', 'O6',
        # Ribose O atoms
        "O2'", "O4'", "O3'", "O5'",
        # Phosphate O atoms
        'OP1', 'OP2', 'O1P', 'O2P',
        # Ring carbons for geometry
        'C2', 'C4', 'C5', 'C6', 'C8',
    }

    current_res_key = None

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue

            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21]
            res_seq = line[22:26].strip()
            ins_code = line[26].strip()

            base_type = RES_MAP.get(res_name.upper())
            if not base_type:
                continue

            # DSSR-style residue ID: "chain.base#"
            dssr_id = f"{chain}.{base_type}{res_seq}{ins_code}".rstrip()

            if dssr_id not in residues:
                residues[dssr_id] = Residue(
                    res_id=dssr_id,
                    base_type=base_type,
                    atoms={}
                )

            if atom_name in TRACKED_ATOMS:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                residues[dssr_id].atoms[atom_name] = np.array([x, y, z])

    return residues


def compare_with_dssr(optimizer: HBondOptimizer,
                      residues: Dict[str, Residue],
                      dssr_hbonds: List[DSSRHBond]):
    """Compare optimizer output with DSSR H-bonds."""

    print("\n" + "="*80)
    print("COMPARISON WITH DSSR")
    print("="*80)

    # Group DSSR H-bonds by residue pair
    from collections import defaultdict
    dssr_by_pair = defaultdict(list)
    for hb in dssr_hbonds:
        pair_key = tuple(sorted([hb.res1, hb.res2]))
        dssr_by_pair[pair_key].append(hb)

    total_dssr = 0
    total_opt = 0
    total_match = 0

    # Process each unique residue pair from DSSR
    for pair_key in sorted(dssr_by_pair.keys()):
        res1_id, res2_id = pair_key

        if res1_id not in residues or res2_id not in residues:
            continue

        res1 = residues[res1_id]
        res2 = residues[res2_id]

        # Get DSSR H-bonds for this pair
        dssr_for_pair = dssr_by_pair[pair_key]

        # Run optimizer
        optimizer.add_residue(res1)
        optimizer.add_residue(res2)
        opt_hbonds = optimizer.optimize_pair(res1.res_id, res2.res_id)

        # Convert to comparable sets (normalize atom pairs)
        dssr_set = set()
        for hb in dssr_for_pair:
            # Normalize: sort atoms to handle direction differences
            atoms = tuple(sorted([hb.donor_atom, hb.acceptor_atom]))
            dssr_set.add(atoms)

        opt_set = set()
        for hb in opt_hbonds:
            atoms = tuple(sorted([hb.donor_atom, hb.acceptor_atom]))
            opt_set.add(atoms)

        matches = dssr_set & opt_set
        dssr_only = dssr_set - opt_set
        opt_only = opt_set - dssr_set

        total_dssr += len(dssr_set)
        total_opt += len(opt_set)
        total_match += len(matches)

        if dssr_only or opt_only:  # Only show mismatches
            print(f"\n{res1_id} <-> {res2_id}")
            print(f"  DSSR ({len(dssr_set)}):      ", end="")
            for hb in dssr_for_pair:
                print(f"{hb.donor_atom}-{hb.acceptor_atom}({hb.distance:.2f}) ", end="")
            print()

            print(f"  Optimizer ({len(opt_set)}): ", end="")
            for hb in opt_hbonds:
                print(f"{hb.donor_atom}-{hb.acceptor_atom}({hb.distance:.2f}) ", end="")
            print()

            if matches:
                print(f"  ✓ Matches: {matches}")
            if opt_only:
                print(f"  + Optimizer only: {opt_only}")
            if dssr_only:
                print(f"  - DSSR only: {dssr_only}")

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"DSSR H-bonds:      {total_dssr}")
    print(f"Optimizer:         {total_opt}")
    print(f"Matching:          {total_match}")
    if total_dssr > 0:
        print(f"Recall:            {total_match/total_dssr*100:.1f}%")
    if total_opt > 0:
        print(f"Precision:         {total_match/total_opt*100:.1f}%")


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Compare H-bond optimizer with DSSR')
    parser.add_argument('pdb', help='PDB file or ID')
    parser.add_argument('--dssr', help='Pre-computed DSSR output file')
    parser.add_argument('--max-dist', type=float, default=4.0, help='Max distance (default: 4.0)')
    parser.add_argument('--min-align', type=float, default=0.3, help='Min alignment (default: 0.3)')
    args = parser.parse_args()

    pdb_path = Path(args.pdb)
    if not pdb_path.exists():
        # Try data/pdb directory
        base_dir = Path(__file__).parent.parent.parent
        pdb_path = base_dir / "data" / "pdb" / f"{args.pdb}.pdb"

    if not pdb_path.exists():
        print(f"PDB file not found: {pdb_path}")
        sys.exit(1)

    # Get DSSR output
    if args.dssr:
        with open(args.dssr) as f:
            dssr_text = f.read()
    else:
        print(f"Running DSSR on {pdb_path}...")
        dssr_text = run_dssr(pdb_path)

    print(f"PDB: {pdb_path}")

    # Parse
    residues = parse_pdb_residues(pdb_path)
    print(f"Loaded {len(residues)} residues")

    dssr_hbonds = parse_dssr_output(dssr_text)
    print(f"DSSR found {len(dssr_hbonds)} H-bonds")

    # Compare with configurable thresholds
    print(f"Using max_distance={args.max_dist}, min_alignment={args.min_align}")

    optimizer = HBondOptimizer(max_distance=args.max_dist, min_alignment=args.min_align)
    compare_with_dssr(optimizer, residues, dssr_hbonds)


if __name__ == '__main__':
    main()

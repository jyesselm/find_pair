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
    res1: str          # e.g., "A.G1" (normalized for matching)
    res2: str          # e.g., "B.C2" (normalized for matching)
    donor_atom: str    # e.g., "N6"
    acceptor_atom: str # e.g., "O4"
    distance: float
    res1_original: str = ""  # Original DSSR ID, e.g., "A.7MG1"
    res2_original: str = ""  # Original DSSR ID, e.g., "B.5MC2"


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

        # Store original IDs before normalization
        res1_orig = res1
        res2_orig = res2

        # Normalize residue IDs for matching with our parsed residues
        res1 = normalize_res_id(res1)
        res2 = normalize_res_id(res2)

        hbonds.append(DSSRHBond(
            res1=res1,
            res2=res2,
            donor_atom=atom1,
            acceptor_atom=atom2,
            distance=dist,
            res1_original=res1_orig,
            res2_original=res2_orig
        ))

    return hbonds


def normalize_res_id(res_id: str) -> str:
    """
    Normalize DSSR residue ID to standard format.

    E.g., "A.DC1" -> "A.C1", "B.DG24" -> "B.G24", "R.PSU655" -> "R.P655", "A.A2M5" -> "A.A5"

    Handles:
    - DNA prefixes (DC, DG, DA, DT, DU)
    - Modified residue names with digits (PSU, 5MC, 7MG, A2M, H2U, etc.)
    - Insertion codes with caret (A.A22^L -> A.A22L)
    - Slash notation for alternate chains (A.U31/6 -> A.U31)
    """
    # Import registry for modified residue lookup
    try:
        from .modified_registry import get_registry
    except ImportError:
        from modified_registry import get_registry

    registry = get_registry()

    parts = res_id.split('.')
    if len(parts) != 2:
        return res_id

    chain = parts[0]
    res = parts[1]

    # Handle insertion code caret: A22^L -> A22L
    res = res.replace('^', '')

    # Handle slash notation: U31/6 -> U31
    if '/' in res:
        res = res.split('/')[0]

    # Standard bases (quick check)
    STANDARD_BASES = {'A', 'G', 'C', 'U', 'T'}
    DNA_PREFIXES = {'DA', 'DG', 'DC', 'DT', 'DU'}

    # PRIORITY 1: Check DNA prefix FIRST (DA, DG, DC, DT, DU)
    # DNA prefixes are always exactly 2 chars
    if len(res) >= 3 and res[:2].upper() in DNA_PREFIXES:
        candidate_num = res[2:]
        if candidate_num and (candidate_num[0].isdigit() or
                              (candidate_num[0] == '-' and len(candidate_num) > 1 and candidate_num[1].isdigit())):
            return f"{chain}.{res[1].upper()}{candidate_num}"

    # PRIORITY 2: Check standard 1-character bases (A, G, C, U, T)
    # BUT only if the remainder is purely numeric (no letters that could be part of modified name)
    # This handles "A231" as "A" + "231" but NOT "A2M5"
    if len(res) >= 2 and res[0].upper() in STANDARD_BASES:
        candidate_num = res[1:]
        # Only match if the rest is purely numeric (or starts with minus for negative numbers)
        # This prevents "A2M5" from matching as "A" + "2M5"
        if candidate_num:
            first_char = candidate_num[0]
            # Check if it's a pure number or negative number
            if first_char.isdigit():
                # Make sure there are no letters after the first digit (except insertion codes at end)
                # This allows "A231" and "A231A" but not "A2M5"
                rest = candidate_num[1:]
                # Check if rest contains only digits and possibly an insertion code letter at the end
                if all(c.isdigit() for c in rest) or (rest and all(c.isdigit() for c in rest[:-1]) and rest[-1].isalpha()):
                    return f"{chain}.{res[0].upper()}{candidate_num}"
            elif first_char == '-' and len(candidate_num) > 1 and candidate_num[1].isdigit():
                return f"{chain}.{res[0].upper()}{candidate_num}"

    # PRIORITY 3: Check for 3-char modified residues (A2M, 5MC, G7M, H2U, PSU, etc.)
    # These have the format: 3-char-code + residue-number
    if len(res) >= 4:  # Need at least 3 chars for code + 1 for number
        candidate_name = res[:3]
        candidate_num = res[3:]
        if candidate_num and (candidate_num[0].isdigit() or
                              (candidate_num[0] == '-' and len(candidate_num) > 1 and candidate_num[1].isdigit())):
            parent = registry.get_parent_base(candidate_name.upper())
            if parent:
                return f"{chain}.{parent}{candidate_num}"

    # PRIORITY 4 (Fallback): Try progressively longer prefixes for unusual modified residues
    # Handles cases not covered above (4+ char modified residue codes like IRES, 5BRU, etc.)
    for i in range(len(res), 0, -1):
        candidate_name = res[:i]
        candidate_num = res[i:]

        # Must have a number at the end (can start with minus for negative residue numbers)
        if not candidate_num:
            continue
        # Valid residue numbers: "1", "10", "-1", "-10", "1A" (with insertion code)
        if not (candidate_num[0].isdigit() or
                (candidate_num[0] == '-' and len(candidate_num) > 1 and candidate_num[1].isdigit())):
            continue

        # Check modified residue registry
        name_upper = candidate_name.upper()
        parent = registry.get_parent_base(name_upper)
        if parent:
            return f"{chain}.{parent}{candidate_num}"

    # Fallback: return as-is
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


def parse_pdb_residues(pdb_path: Path, include_protein: bool = True) -> Dict[str, Residue]:
    """Parse PDB and return residues keyed by DSSR-style ID (e.g., "A.G1").

    Uses ModifiedResidueRegistry to handle modified nucleotides like 5MC, H2U, PSU.

    Args:
        pdb_path: Path to PDB file
        include_protein: If True, also parse amino acid residues for RNA-protein H-bonds
    """
    # Import registry for modified residue lookup
    try:
        from .modified_registry import get_registry
    except ImportError:
        from modified_registry import get_registry

    registry = get_registry()
    residues = {}

    # Standard nucleotide residue mapping (quick lookup)
    NUC_RES_MAP = {
        'ADE': 'A', 'A': 'A', 'DA': 'A', 'RA': 'A',
        'GUA': 'G', 'G': 'G', 'DG': 'G', 'RG': 'G',
        'CYT': 'C', 'C': 'C', 'DC': 'C', 'RC': 'C',
        'URA': 'U', 'U': 'U', 'RU': 'U',
        'THY': 'T', 'T': 'T', 'DT': 'T',
    }

    # Standard amino acids (3-letter codes)
    AMINO_ACIDS = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    }

    # Atoms we track for nucleotides
    NUC_TRACKED_ATOMS = {
        # Base N atoms
        'N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9',
        # Base O atoms
        'O2', 'O4', 'O6',
        # Ribose O atoms
        "O2'", "O4'", "O3'", "O5'",
        # Ribose C atoms (needed for H-bond geometry)
        "C1'", "C2'", "C3'", "C4'", "C5'",
        # Phosphate O atoms
        'OP1', 'OP2', 'O1P', 'O2P',
        # Ring carbons for geometry
        'C2', 'C4', 'C5', 'C6', 'C8',
    }

    # Atoms we track for amino acids
    PROTEIN_TRACKED_ATOMS = {
        # Backbone atoms
        'N', 'CA', 'C', 'O', 'OXT',
        # Side chain donors/acceptors by residue type
        # ARG: guanidinium
        'NE', 'NH1', 'NH2', 'CZ',
        # ASN: amide
        'ND2', 'OD1', 'CG',
        # ASP: carboxylate
        'OD1', 'OD2',
        # CYS: thiol (can be donor in some cases)
        'SG',
        # GLN: amide
        'NE2', 'OE1', 'CD',
        # GLU: carboxylate
        'OE1', 'OE2',
        # HIS: imidazole
        'ND1', 'NE2', 'CE1', 'CD2', 'CG',
        # LYS: amino
        'NZ', 'CE',
        # SER: hydroxyl
        'OG', 'CB',
        # THR: hydroxyl
        'OG1', 'CB', 'CG2',
        # TRP: indole NH
        'NE1', 'CD1', 'CE2',
        # TYR: hydroxyl
        'OH', 'CZ', 'CE1', 'CE2',
    }

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue

            # Handle alternate conformers - only take ' ' or 'A' (major conformer)
            alt_loc = line[16] if len(line) > 16 else ' '
            if alt_loc not in (' ', 'A'):
                continue

            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21]
            res_seq = line[22:26].strip()
            ins_code = line[26].strip()

            res_name_upper = res_name.upper()

            # Check if it's a nucleotide
            base_type = NUC_RES_MAP.get(res_name_upper)
            is_nucleotide = True

            # If not found, try the modified residue registry
            if not base_type:
                base_type = registry.get_parent_base(res_name_upper)

            # Check if it's an amino acid
            if not base_type and include_protein and res_name_upper in AMINO_ACIDS:
                base_type = res_name_upper  # Use 3-letter code as base_type for proteins
                is_nucleotide = False

            if not base_type:
                continue

            # Build residue ID
            # For nucleotides: "chain.base#" (e.g., "A.G1")
            # For proteins: "chain.RES#" (e.g., "A.ALA1")
            dssr_id = f"{chain}.{base_type}{res_seq}{ins_code}".rstrip()

            if dssr_id not in residues:
                residues[dssr_id] = Residue(
                    res_id=dssr_id,
                    base_type=base_type,
                    atoms={},
                    residue_code=res_name_upper  # Store original 3-letter code
                )

            # Track appropriate atoms based on residue type
            tracked = NUC_TRACKED_ATOMS if is_nucleotide else PROTEIN_TRACKED_ATOMS
            if atom_name in tracked:
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

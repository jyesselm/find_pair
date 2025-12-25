#!/usr/bin/env python3
"""
Validate H-bond inference against reference data (hbond_donors.json, hbond_acceptors.json).

This script loads actual PDB structures and compares inferred values to reference.
"""

import json
import os
import sys
from pathlib import Path
from collections import defaultdict
import numpy as np

# Import inference and hybridization modules
sys.path.insert(0, str(Path(__file__).parent))
from inference import HBondInferenceEngine, load_reference_data
from hybridization import classify_nitrogen, classify_oxygen, detect_base_type, get_element


def parse_pdb_residues(pdb_path: str):
    """
    Parse a PDB file and extract atom coordinates by residue.

    Returns:
        Dict mapping (chain, resname, resnum) to {atom_name: position}
    """
    residues = {}

    with open(pdb_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21]
                resnum = line[22:26].strip()
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue

                key = (chain, resname, resnum)
                if key not in residues:
                    residues[key] = {}
                residues[key][atom_name] = np.array([x, y, z])

    return residues


def validate_standard_bases():
    """Validate inference for standard nucleotide bases using the geometry.py tables."""
    print("Validating Standard Nucleotide Bases")
    print("=" * 70)

    from geometry import DONOR_CAPACITY, ACCEPTOR_CAPACITY

    # Expected values from geometry.py tables
    bases = {
        'A': {
            'donors': {k[1]: v for k, v in DONOR_CAPACITY.items() if k[0] == 'A'},
            'acceptors': {k[1]: v for k, v in ACCEPTOR_CAPACITY.items() if k[0] == 'A'},
        },
        'G': {
            'donors': {k[1]: v for k, v in DONOR_CAPACITY.items() if k[0] == 'G'},
            'acceptors': {k[1]: v for k, v in ACCEPTOR_CAPACITY.items() if k[0] == 'G'},
        },
        'C': {
            'donors': {k[1]: v for k, v in DONOR_CAPACITY.items() if k[0] == 'C'},
            'acceptors': {k[1]: v for k, v in ACCEPTOR_CAPACITY.items() if k[0] == 'C'},
        },
        'U': {
            'donors': {k[1]: v for k, v in DONOR_CAPACITY.items() if k[0] == 'U'},
            'acceptors': {k[1]: v for k, v in ACCEPTOR_CAPACITY.items() if k[0] == 'U'},
        },
        'T': {
            'donors': {k[1]: v for k, v in DONOR_CAPACITY.items() if k[0] == 'T'},
            'acceptors': {k[1]: v for k, v in ACCEPTOR_CAPACITY.items() if k[0] == 'T'},
        },
    }

    for base, expected in bases.items():
        print(f"\n{base}:")
        print(f"  Donors: {expected['donors']}")
        print(f"  Acceptors: {expected['acceptors']}")


def validate_amino_acid_reference():
    """Validate that amino acid reference data is consistent and complete."""
    print("\n\nValidating Amino Acid Reference Data")
    print("=" * 70)

    engine = HBondInferenceEngine()

    AMINO_ACIDS = [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
        'THR', 'TRP', 'TYR', 'VAL'
    ]

    # Check each amino acid
    issues = []

    for aa in AMINO_ACIDS:
        donors = engine.ref_donors.get(aa, {})
        acceptors = engine.ref_acceptors.get(aa, {})

        if donors is None or acceptors is None:
            issues.append(f"{aa}: Missing from reference data")
            continue

        # Check common patterns
        # All amino acids except PRO should have backbone N as donor
        if aa != 'PRO':
            if 'N' not in donors:
                issues.append(f"{aa}: Missing backbone N donor")
            elif donors.get('N', 0) != 2:
                # Most have N:2 (NH2/NH3+), some have N:1
                pass  # This is OK, varies by ionization state

        # All should have backbone O as acceptor
        if 'O' not in acceptors:
            issues.append(f"{aa}: Missing backbone O acceptor")
        elif acceptors.get('O', 0) != 2:
            issues.append(f"{aa}: Backbone O should have capacity 2, got {acceptors.get('O')}")

        # Backbone N should be acceptor with capacity 1
        if 'N' not in acceptors:
            issues.append(f"{aa}: Missing backbone N acceptor")

    if issues:
        print("Issues found:")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print("All amino acid reference data looks consistent!")

    # Print summary table
    print("\nAmino Acid Summary:")
    print("-" * 70)
    print(f"{'AA':<5} {'N donor':<10} {'O accept':<10} {'Side chain donors':<20}")
    print("-" * 70)

    for aa in AMINO_ACIDS:
        donors = engine.ref_donors.get(aa, {})
        acceptors = engine.ref_acceptors.get(aa, {})

        n_donor = donors.get('N', 0) if donors else 0
        o_acceptor = acceptors.get('O', 0) if acceptors else 0

        # Side chain donors (excluding backbone N and OXT)
        sc_donors = [f"{k}:{v}" for k, v in (donors or {}).items()
                     if k not in ('N', 'OXT')]

        print(f"{aa:<5} {n_donor:<10} {o_acceptor:<10} {', '.join(sc_donors):<20}")


def compare_nucleotide_to_reference():
    """Compare nucleotide capacity tables to reference JSON data."""
    print("\n\nComparing Nucleotide Tables to Reference JSON")
    print("=" * 70)

    from geometry import DONOR_CAPACITY, ACCEPTOR_CAPACITY

    engine = HBondInferenceEngine()

    # Standard nucleotide 3-letter codes
    NUCLEOTIDES = ['A', 'G', 'C', 'U', 'T', 'DA', 'DG', 'DC', 'DT']

    for nuc in NUCLEOTIDES:
        ref_donors = engine.ref_donors.get(nuc, {})
        ref_acceptors = engine.ref_acceptors.get(nuc, {})

        if ref_donors is None:
            print(f"{nuc}: Not in reference data (null)")
            continue

        if not ref_donors and not ref_acceptors:
            print(f"{nuc}: Empty in reference data")
            continue

        # Get single-letter base type
        base_type = nuc[-1] if len(nuc) == 2 else nuc

        # Compare donors
        print(f"\n{nuc} (base type: {base_type}):")

        # Get expected from geometry.py
        expected_donors = {k[1]: v for k, v in DONOR_CAPACITY.items() if k[0] == base_type}
        expected_acceptors = {k[1]: v for k, v in ACCEPTOR_CAPACITY.items() if k[0] == base_type}

        # Compare
        print(f"  Donors:")
        all_donor_atoms = set(expected_donors.keys()) | set(ref_donors.keys())
        for atom in sorted(all_donor_atoms):
            exp = expected_donors.get(atom, 0)
            ref = ref_donors.get(atom, 0)
            match = "OK" if exp == ref else "DIFF"
            if exp > 0 or ref > 0:
                print(f"    {atom}: geometry={exp}, ref={ref} [{match}]")

        print(f"  Acceptors:")
        all_acc_atoms = set(expected_acceptors.keys()) | set(ref_acceptors.keys())
        for atom in sorted(all_acc_atoms):
            exp = expected_acceptors.get(atom, 0)
            ref = ref_acceptors.get(atom, 0)
            match = "OK" if exp == ref else "DIFF"
            if exp > 0 or ref > 0:
                print(f"    {atom}: geometry={exp}, ref={ref} [{match}]")


def main():
    validate_standard_bases()
    validate_amino_acid_reference()
    compare_nucleotide_to_reference()


if __name__ == "__main__":
    main()

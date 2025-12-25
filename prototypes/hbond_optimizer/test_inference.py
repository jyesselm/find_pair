#!/usr/bin/env python3
"""
Test inference against existing capacity tables in geometry.py.
"""

import numpy as np
from pathlib import Path
import sys

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from geometry import DONOR_CAPACITY, ACCEPTOR_CAPACITY, BASE_CONNECTIVITY
from hybridization import (
    classify_nitrogen, classify_oxygen, detect_atom_hybridization,
    get_element, Hybridization
)


def create_test_atoms_adenine():
    """Create idealized adenine atom positions."""
    # Approximate planar adenine coordinates
    return {
        'N9': np.array([0.0, 0.0, 0.0]),
        'C8': np.array([1.2, 0.5, 0.0]),
        'N7': np.array([1.2, 1.8, 0.0]),
        'C5': np.array([0.0, 2.2, 0.0]),
        'C6': np.array([-0.8, 3.4, 0.0]),
        'N6': np.array([-0.4, 4.7, 0.0]),  # Amino group
        'N1': np.array([-2.1, 3.2, 0.0]),
        'C2': np.array([-2.6, 1.9, 0.0]),
        'N3': np.array([-1.8, 0.8, 0.0]),
        'C4': np.array([-0.5, 1.0, 0.0]),
        # Ribose
        "C1'": np.array([0.0, -1.4, 0.0]),
        "O4'": np.array([1.3, -1.8, 0.0]),
        "C4'": np.array([1.4, -3.2, 0.0]),
        "C2'": np.array([-0.8, -2.6, 0.0]),
        "O2'": np.array([-2.1, -2.8, 0.0]),
        "C3'": np.array([0.0, -3.8, 0.0]),
        "O3'": np.array([-0.3, -5.1, 0.0]),
        "C5'": np.array([2.7, -3.8, 0.0]),
        "O5'": np.array([3.0, -5.1, 0.0]),
    }


def create_test_atoms_guanine():
    """Create idealized guanine atom positions."""
    return {
        'N9': np.array([0.0, 0.0, 0.0]),
        'C8': np.array([1.2, 0.5, 0.0]),
        'N7': np.array([1.2, 1.8, 0.0]),
        'C5': np.array([0.0, 2.2, 0.0]),
        'C6': np.array([-0.8, 3.4, 0.0]),
        'O6': np.array([-0.4, 4.7, 0.0]),  # Carbonyl
        'N1': np.array([-2.1, 3.2, 0.0]),  # Imino NH
        'C2': np.array([-2.6, 1.9, 0.0]),
        'N2': np.array([-3.9, 1.7, 0.0]),  # Amino group
        'N3': np.array([-1.8, 0.8, 0.0]),
        'C4': np.array([-0.5, 1.0, 0.0]),
        # Ribose (simplified)
        "C1'": np.array([0.0, -1.4, 0.0]),
        "O4'": np.array([1.3, -1.8, 0.0]),
        "C2'": np.array([-0.8, -2.6, 0.0]),
        "O2'": np.array([-2.1, -2.8, 0.0]),
        "C3'": np.array([0.0, -3.8, 0.0]),
        "O3'": np.array([-0.3, -5.1, 0.0]),
    }


def create_test_atoms_cytosine():
    """Create idealized cytosine atom positions."""
    return {
        'N1': np.array([0.0, 0.0, 0.0]),
        'C2': np.array([1.2, 0.4, 0.0]),
        'O2': np.array([2.2, -0.3, 0.0]),  # Carbonyl
        'N3': np.array([1.2, 1.8, 0.0]),
        'C4': np.array([0.0, 2.5, 0.0]),
        'N4': np.array([0.0, 3.9, 0.0]),  # Amino group
        'C5': np.array([-1.2, 1.8, 0.0]),
        'C6': np.array([-1.2, 0.4, 0.0]),
        # Ribose
        "C1'": np.array([0.0, -1.4, 0.0]),
        "O4'": np.array([1.3, -1.8, 0.0]),
        "C2'": np.array([-0.8, -2.6, 0.0]),
        "O2'": np.array([-2.1, -2.8, 0.0]),
        "C3'": np.array([0.0, -3.8, 0.0]),
        "O3'": np.array([-0.3, -5.1, 0.0]),
    }


def create_test_atoms_uracil():
    """Create idealized uracil atom positions."""
    return {
        'N1': np.array([0.0, 0.0, 0.0]),
        'C2': np.array([1.2, 0.4, 0.0]),
        'O2': np.array([2.2, -0.3, 0.0]),  # Carbonyl
        'N3': np.array([1.2, 1.8, 0.0]),  # Imino NH
        'C4': np.array([0.0, 2.5, 0.0]),
        'O4': np.array([0.0, 3.9, 0.0]),  # Carbonyl
        'C5': np.array([-1.2, 1.8, 0.0]),
        'C6': np.array([-1.2, 0.4, 0.0]),
        # Ribose
        "C1'": np.array([0.0, -1.4, 0.0]),
        "O4'": np.array([1.3, -1.8, 0.0]),
        "C2'": np.array([-0.8, -2.6, 0.0]),
        "O2'": np.array([-2.1, -2.8, 0.0]),
        "C3'": np.array([0.0, -3.8, 0.0]),
        "O3'": np.array([-0.3, -5.1, 0.0]),
    }


def test_base_atoms():
    """Test inference on base atoms against DONOR_CAPACITY and ACCEPTOR_CAPACITY."""
    print("Testing Base Atom Inference")
    print("=" * 70)

    test_cases = [
        ('A', create_test_atoms_adenine()),
        ('G', create_test_atoms_guanine()),
        ('C', create_test_atoms_cytosine()),
        ('U', create_test_atoms_uracil()),
    ]

    total_donor_matches = 0
    total_donor_mismatches = 0
    total_acceptor_matches = 0
    total_acceptor_mismatches = 0

    for base_type, atoms in test_cases:
        print(f"\n{base_type}:")
        print("-" * 50)

        # Test key atoms
        for atom_name in atoms.keys():
            element = get_element(atom_name)
            if element not in ('N', 'O'):
                continue

            # Get expected values
            expected_donor = DONOR_CAPACITY.get((base_type, atom_name.strip()), 0)
            expected_acceptor = ACCEPTOR_CAPACITY.get((base_type, atom_name.strip()), 0)

            # Infer from geometry
            if element == 'N':
                inferred_donor, inferred_acceptor = classify_nitrogen(atom_name, atoms)
            else:
                inferred_donor, inferred_acceptor = classify_oxygen(atom_name, atoms)

            # Compare
            donor_match = "OK" if inferred_donor == expected_donor else "MISMATCH"
            acceptor_match = "OK" if inferred_acceptor == expected_acceptor else "MISMATCH"

            if donor_match == "OK" and expected_donor > 0:
                total_donor_matches += 1
            elif donor_match == "MISMATCH":
                total_donor_mismatches += 1

            if acceptor_match == "OK" and expected_acceptor > 0:
                total_acceptor_matches += 1
            elif acceptor_match == "MISMATCH":
                total_acceptor_mismatches += 1

            # Only print if there's something interesting
            if expected_donor > 0 or expected_acceptor > 0 or inferred_donor > 0 or inferred_acceptor > 0:
                status = ""
                if donor_match == "MISMATCH" or acceptor_match == "MISMATCH":
                    status = " <-- FIX"
                print(f"  {atom_name:>4}: donor={inferred_donor} (exp={expected_donor}) {donor_match:>8}, "
                      f"acceptor={inferred_acceptor} (exp={expected_acceptor}) {acceptor_match:>8}{status}")

    print("\n" + "=" * 70)
    print(f"Summary:")
    print(f"  Donor matches:     {total_donor_matches}")
    print(f"  Donor mismatches:  {total_donor_mismatches}")
    print(f"  Acceptor matches:  {total_acceptor_matches}")
    print(f"  Acceptor mismatches: {total_acceptor_mismatches}")


def test_specific_atoms():
    """Test specific atoms that are in DONOR_CAPACITY/ACCEPTOR_CAPACITY."""
    print("\n\nTesting Specific Key Atoms")
    print("=" * 70)

    # Test cases: (base, atom, expected_donor, expected_acceptor, description)
    key_atoms = [
        # Amino groups (should be donor=2, acceptor=1)
        ('A', 'N6', 2, 0, 'Adenine amino'),
        ('C', 'N4', 2, 0, 'Cytosine amino'),
        ('G', 'N2', 2, 0, 'Guanine amino'),

        # Imino groups (should be donor=1, acceptor=0 or 1)
        ('G', 'N1', 1, 0, 'Guanine imino'),
        ('U', 'N3', 1, 0, 'Uracil imino'),

        # Ring nitrogens (acceptors only)
        ('A', 'N1', 0, 1, 'Adenine N1'),
        ('A', 'N3', 0, 1, 'Adenine N3'),
        ('A', 'N7', 0, 1, 'Adenine N7'),
        ('G', 'N3', 0, 1, 'Guanine N3'),
        ('G', 'N7', 0, 1, 'Guanine N7'),
        ('C', 'N3', 0, 1, 'Cytosine N3'),

        # Carbonyl oxygens (acceptor=2)
        ('G', 'O6', 0, 2, 'Guanine O6'),
        ('U', 'O2', 0, 2, 'Uracil O2'),
        ('U', 'O4', 0, 2, 'Uracil O4'),
        ('C', 'O2', 0, 2, 'Cytosine O2'),

        # Ribose O2' (donor=1, acceptor=2)
        ('A', "O2'", 1, 2, "Adenine O2'"),
        ('G', "O2'", 1, 2, "Guanine O2'"),
    ]

    atoms_A = create_test_atoms_adenine()
    atoms_G = create_test_atoms_guanine()
    atoms_C = create_test_atoms_cytosine()
    atoms_U = create_test_atoms_uracil()

    base_atoms = {'A': atoms_A, 'G': atoms_G, 'C': atoms_C, 'U': atoms_U}

    passed = 0
    failed = 0

    for base, atom, exp_donor, exp_acceptor, desc in key_atoms:
        atoms = base_atoms[base]
        if atom not in atoms:
            print(f"  {desc}: SKIP (atom not in test coordinates)")
            continue

        element = get_element(atom)
        if element == 'N':
            inf_donor, inf_acceptor = classify_nitrogen(atom, atoms)
        else:
            inf_donor, inf_acceptor = classify_oxygen(atom, atoms)

        donor_ok = inf_donor == exp_donor
        acceptor_ok = inf_acceptor == exp_acceptor

        if donor_ok and acceptor_ok:
            status = "PASS"
            passed += 1
        else:
            status = "FAIL"
            failed += 1

        print(f"  {desc:25s}: donor={inf_donor}(exp={exp_donor}), "
              f"acceptor={inf_acceptor}(exp={exp_acceptor}) [{status}]")

    print(f"\nResults: {passed} passed, {failed} failed")


def analyze_geometry_issues():
    """Analyze why certain atoms might be misclassified."""
    print("\n\nAnalyzing Geometry-Based Classification")
    print("=" * 70)

    atoms_A = create_test_atoms_adenine()
    atoms_G = create_test_atoms_guanine()

    # Check N6 in adenine (amino - should have 1 heavy atom bond)
    print("\nAdenine N6 (amino group):")
    hyb, n_bonds, bonded = detect_atom_hybridization('N6', atoms_A)
    print(f"  Hybridization: {hyb}")
    print(f"  Num bonds to heavy atoms: {n_bonds}")
    print(f"  Bonded to: {bonded}")
    print(f"  Expected: 1 bond (to C6), sp2, donor=2")

    # Check N1 in adenine (ring N - should have 2 heavy atom bonds)
    print("\nAdenine N1 (ring nitrogen):")
    hyb, n_bonds, bonded = detect_atom_hybridization('N1', atoms_A)
    print(f"  Hybridization: {hyb}")
    print(f"  Num bonds to heavy atoms: {n_bonds}")
    print(f"  Bonded to: {bonded}")
    print(f"  Expected: 2 bonds (C2, C6), sp2, acceptor=1, no H")

    # Check N1 in guanine (imino - should have 2 heavy atom bonds + 1H)
    print("\nGuanine N1 (imino group):")
    hyb, n_bonds, bonded = detect_atom_hybridization('N1', atoms_G)
    print(f"  Hybridization: {hyb}")
    print(f"  Num bonds to heavy atoms: {n_bonds}")
    print(f"  Bonded to: {bonded}")
    print(f"  Expected: 2 bonds (C2, C6), sp2, donor=1 (has H)")

    # Check O6 in guanine (carbonyl)
    print("\nGuanine O6 (carbonyl):")
    hyb, n_bonds, bonded = detect_atom_hybridization('O6', atoms_G)
    print(f"  Hybridization: {hyb}")
    print(f"  Num bonds to heavy atoms: {n_bonds}")
    print(f"  Bonded to: {bonded}")
    print(f"  Expected: 1 bond (C6), sp2, acceptor=2, no H")


if __name__ == "__main__":
    test_base_atoms()
    test_specific_atoms()
    analyze_geometry_issues()

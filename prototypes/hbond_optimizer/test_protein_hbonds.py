#!/usr/bin/env python3
"""Test RNA-protein H-bond support."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from prototypes.hbond_optimizer import parse_pdb_residues
from prototypes.hbond_optimizer.optimizer import HBondOptimizer
from prototypes.hbond_optimizer.geometry import (
    predict_h_slots, predict_lp_slots, get_donor_capacity, get_acceptor_capacity
)


def test_amino_acid_parsing():
    """Test that amino acids are parsed from PDB files."""
    # Find a PDB with both RNA and protein
    project_root = Path(__file__).parent.parent.parent
    pdb_dir = project_root / "data" / "pdb"

    # Try PDBs that have both RNA and protein
    test_pdbs = ['1A34', '1A9N', '1AQ3', '1ASY']

    for pdb_id in test_pdbs:
        pdb_path = pdb_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            continue

        print(f"\n{'='*60}")
        print(f"Testing {pdb_id}")
        print('='*60)

        residues = parse_pdb_residues(pdb_path, include_protein=True)

        # Count residue types
        nucleotides = []
        amino_acids = []

        for res_id, res in residues.items():
            if res.base_type in ['A', 'G', 'C', 'U', 'T', 'P', 'I']:
                nucleotides.append(res_id)
            else:
                amino_acids.append(res_id)

        print(f"  Nucleotides: {len(nucleotides)}")
        print(f"  Amino acids: {len(amino_acids)}")

        if amino_acids:
            print(f"\n  Sample amino acids: {amino_acids[:5]}")

            # Test H-slot and LP-slot prediction for first amino acid
            for aa_id in amino_acids[:3]:
                aa = residues[aa_id]
                print(f"\n  Residue: {aa_id} ({aa.base_type})")
                print(f"    Atoms: {list(aa.atoms.keys())}")

                # Test backbone N (donor)
                if 'N' in aa.atoms and 'CA' in aa.atoms:
                    capacity = get_donor_capacity(aa.base_type, 'N')
                    print(f"    N donor capacity: {capacity}")

                    h_slots = predict_h_slots(aa.base_type, 'N', aa.atoms, None)
                    print(f"    N H-slots: {len(h_slots)}")
                    for i, slot in enumerate(h_slots):
                        print(f"      H slot {i}: dir={slot.direction}")

                # Test backbone O (acceptor)
                if 'O' in aa.atoms and 'C' in aa.atoms:
                    capacity = get_acceptor_capacity(aa.base_type, 'O')
                    print(f"    O acceptor capacity: {capacity}")

                    lp_slots = predict_lp_slots(aa.base_type, 'O', aa.atoms, None)
                    print(f"    O LP-slots: {len(lp_slots)}")
                    for i, slot in enumerate(lp_slots):
                        print(f"      LP slot {i}: dir={slot.direction}")

        if amino_acids and nucleotides:
            print(f"\n  Testing RNA-protein H-bond detection...")

            optimizer = HBondOptimizer(max_distance=4.0, min_alignment=0.3)
            for res in residues.values():
                optimizer.add_residue(res)

            # First find close RNA-protein pairs by checking distances
            import numpy as np
            close_pairs = []
            for nuc_id in nucleotides:
                nuc = residues[nuc_id]
                for aa_id in amino_acids:
                    aa = residues[aa_id]
                    # Check minimum distance between any atoms
                    for nuc_atom, nuc_pos in nuc.atoms.items():
                        for aa_atom, aa_pos in aa.atoms.items():
                            dist = np.linalg.norm(aa_pos - nuc_pos)
                            if dist < 4.0:
                                close_pairs.append((nuc_id, aa_id, nuc_atom, aa_atom, dist))

            print(f"    Found {len(close_pairs)} close atom pairs (< 4.0 Å)")

            # Try to find H-bonds between close RNA-protein pairs
            found_rna_protein = 0
            checked_pairs = set()
            for nuc_id, aa_id, _, _, _ in close_pairs[:50]:
                pair_key = (nuc_id, aa_id)
                if pair_key in checked_pairs:
                    continue
                checked_pairs.add(pair_key)

                hbonds = optimizer.optimize_pair(nuc_id, aa_id)
                if hbonds:
                    found_rna_protein += len(hbonds)
                    print(f"    Found {len(hbonds)} H-bond(s) between {nuc_id} and {aa_id}")
                    for hb in hbonds:
                        print(f"      {hb.donor_res_id}.{hb.donor_atom} -> {hb.acceptor_res_id}.{hb.acceptor_atom} dist={hb.distance:.2f} align={hb.alignment_score:.2f}")

            print(f"\n  Total RNA-protein H-bonds found: {found_rna_protein}")

        break  # Just test first available PDB


def test_amino_acid_capacities():
    """Test that CIF-based capacities work for amino acids."""
    print("\n" + "="*60)
    print("Testing amino acid donor/acceptor capacities")
    print("="*60)

    test_cases = [
        # (residue, atom, expected_donor, expected_acceptor)
        ('ALA', 'N', 2, 0),   # Backbone amide - donor
        ('ALA', 'O', 0, 2),   # Backbone carbonyl - acceptor
        ('ARG', 'NH1', 2, 0), # Guanidinium - donor
        ('ARG', 'NH2', 2, 0), # Guanidinium - donor
        ('ASN', 'ND2', 2, 0), # Amide NH2 - donor
        ('ASN', 'OD1', 0, 2), # Amide carbonyl - acceptor
        ('LYS', 'NZ', 3, 0),  # Amino NH3+ - donor
        ('SER', 'OG', 1, 2),  # Hydroxyl - both
        ('TYR', 'OH', 1, 2),  # Hydroxyl - both
    ]

    for res, atom, exp_donor, exp_acceptor in test_cases:
        donor_cap = get_donor_capacity(res, atom)
        acc_cap = get_acceptor_capacity(res, atom)

        donor_ok = "OK" if donor_cap == exp_donor else f"FAIL (got {donor_cap})"
        acc_ok = "OK" if acc_cap == exp_acceptor else f"FAIL (got {acc_cap})"

        print(f"  {res}.{atom}: donor={donor_cap} [{donor_ok}], acceptor={acc_cap} [{acc_ok}]")


if __name__ == '__main__':
    test_amino_acid_capacities()
    test_amino_acid_parsing()

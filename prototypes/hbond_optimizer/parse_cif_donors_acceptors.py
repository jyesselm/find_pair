#!/usr/bin/env python3
"""
Parse CIF files with explicit hydrogens to determine true donor/acceptor capacities.

These CIF files have explicit H atoms, so we can accurately determine:
- Donors: N or O atoms with H atoms bonded to them
- Acceptors: N or O atoms (capacity based on hybridization)

Terminal atoms that would be different in chain context:
- OP3 (with H) - terminal phosphate, wouldn't exist in chain
- O3' (with H) - in chain, connects to next phosphate (no H)
- OP2 may have H in monomer but not always in chain
"""

import re
import json
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Set, Optional
from collections import defaultdict


@dataclass
class Atom:
    """Atom from CIF file."""
    name: str
    element: str
    leaving_atom: bool = False


@dataclass
class Bond:
    """Bond from CIF file."""
    atom1: str
    atom2: str
    order: str  # SING, DOUB, TRIP


@dataclass
class ResidueInfo:
    """Donor/acceptor info for a residue."""
    residue_code: str
    residue_name: str
    donors: Dict[str, int] = field(default_factory=dict)
    acceptors: Dict[str, int] = field(default_factory=dict)
    terminal_donors: Set[str] = field(default_factory=set)


def parse_cif_file(cif_path: Path) -> Tuple[str, str, List[Atom], List[Bond]]:
    """Parse a CIF file to extract atoms and bonds."""
    with open(cif_path) as f:
        lines = f.readlines()

    # Extract residue code and name
    residue_code = cif_path.stem
    residue_name = ""

    for line in lines:
        if line.startswith('_chem_comp.id'):
            parts = line.split()
            if len(parts) >= 2:
                residue_code = parts[1]
        elif line.startswith('_chem_comp.name'):
            match = re.search(r'"([^"]+)"', line)
            if match:
                residue_name = match.group(1)

    atoms = []
    bonds = []

    # Parse atom data - lines that start with residue_code followed by atom name
    # Format: CODE ATOM_ID ALT_ID ELEMENT CHARGE ALIGN AROM LEAVING ...
    in_atom_section = False
    in_bond_section = False

    for line in lines:
        stripped = line.strip()

        # Detect section starts
        if '_chem_comp_atom.comp_id' in line:
            in_atom_section = True
            in_bond_section = False
            continue
        elif '_chem_comp_bond.comp_id' in line:
            in_atom_section = False
            in_bond_section = True
            continue
        elif stripped.startswith('loop_') or stripped.startswith('#'):
            if not '_chem_comp_atom' in stripped and not '_chem_comp_bond' in stripped:
                in_atom_section = False
                in_bond_section = False
            continue
        elif stripped.startswith('_'):
            continue  # Column header

        # Parse atom data
        if in_atom_section and stripped.startswith(residue_code + ' '):
            parts = stripped.split()
            if len(parts) >= 8:
                # Handle quoted atom names like "O5'" or "C1'"
                atom_name = parts[1].strip('"')
                element = parts[3]
                leaving = parts[7] == 'Y' if len(parts) > 7 else False
                atoms.append(Atom(name=atom_name, element=element, leaving_atom=leaving))

        # Parse bond data
        if in_bond_section and stripped.startswith(residue_code + ' '):
            parts = stripped.split()
            if len(parts) >= 4:
                # Remove quotes from atom names
                atom1 = parts[1].strip('"')
                atom2 = parts[2].strip('"')
                order = parts[3]
                bonds.append(Bond(atom1=atom1, atom2=atom2, order=order))

    return residue_code, residue_name, atoms, bonds


def analyze_donors_acceptors(atoms: List[Atom], bonds: List[Bond]) -> Tuple[Dict, Dict, Set]:
    """Analyze atoms and bonds to determine donors and acceptors."""
    atom_by_name = {a.name: a for a in atoms}

    # Build adjacency
    neighbors = defaultdict(list)
    bond_orders = {}
    for bond in bonds:
        neighbors[bond.atom1].append(bond.atom2)
        neighbors[bond.atom2].append(bond.atom1)
        bond_orders[(bond.atom1, bond.atom2)] = bond.order
        bond_orders[(bond.atom2, bond.atom1)] = bond.order

    donors = {}
    acceptors = {}
    terminal_donors = set()

    for atom in atoms:
        if atom.element not in ('N', 'O'):
            continue

        atom_name = atom.name

        # Count H atoms bonded to this atom
        h_count = sum(1 for n in neighbors[atom_name]
                     if n in atom_by_name and atom_by_name[n].element == 'H')

        # Count heavy atom neighbors
        heavy_neighbors = [n for n in neighbors[atom_name]
                          if n in atom_by_name and atom_by_name[n].element != 'H']

        # Check for double bonds
        has_double_bond = any(bond_orders.get((atom_name, n), 'SING') == 'DOUB'
                             for n in heavy_neighbors)

        # Determine donor capacity (number of H atoms)
        if h_count > 0:
            donors[atom_name] = h_count

            # Check if terminal (wouldn't have H in chain context)
            if atom.leaving_atom:
                terminal_donors.add(atom_name)
            # OP3 is 5' terminal phosphate oxygen - doesn't exist in chain
            if atom_name == 'OP3':
                terminal_donors.add(atom_name)
            # OP2 typically anionic in chain (no H)
            if atom_name == 'OP2':
                terminal_donors.add(atom_name)
            # O3' with H is terminal (in chain, connects to next phosphate)
            if atom_name == "O3'":
                terminal_donors.add(atom_name)

        # Determine acceptor capacity
        if atom.element == 'O':
            if len(heavy_neighbors) == 1:
                if has_double_bond:
                    # Carbonyl oxygen - 2 lone pairs
                    acceptors[atom_name] = 2
                else:
                    # Hydroxyl or terminal O - 2 lone pairs
                    acceptors[atom_name] = 2
            elif len(heavy_neighbors) == 2:
                # Bridging oxygen (O4', O5' in chain) - less accessible
                acceptors[atom_name] = 1

        elif atom.element == 'N':
            if h_count >= 2:
                # Amino group -NH2: 1 lone pair (though usually only donor)
                # For bases, amino N typically doesn't act as acceptor
                acceptors[atom_name] = 0
            elif h_count == 1:
                # Imino =NH: typically only donor, weak acceptor
                acceptors[atom_name] = 0
            elif len(heavy_neighbors) == 2 and h_count == 0:
                # Ring nitrogen with lone pair (like A.N1, A.N3, A.N7)
                acceptors[atom_name] = 1
            elif len(heavy_neighbors) == 3:
                # Tertiary N (like N9 connecting sugar to base)
                acceptors[atom_name] = 0

    return donors, acceptors, terminal_donors


def process_all_cifs(cif_dir: Path) -> Dict[str, ResidueInfo]:
    """Process all CIF files in directory."""
    results = {}
    cif_files = list(cif_dir.glob("*.cif"))
    print(f"Processing {len(cif_files)} CIF files...")

    for cif_path in sorted(cif_files):
        try:
            code, name, atoms, bonds = parse_cif_file(cif_path)
            donors, acceptors, terminal = analyze_donors_acceptors(atoms, bonds)

            results[code] = ResidueInfo(
                residue_code=code,
                residue_name=name,
                donors=donors,
                acceptors=acceptors,
                terminal_donors=terminal
            )
        except Exception as e:
            print(f"Error processing {cif_path.name}: {e}")

    return results


def generate_json_output(results: Dict[str, ResidueInfo], output_dir: Path):
    """Generate JSON files for donors and acceptors."""
    donors_all = {}  # All donors including terminal
    donors_chain = {}  # Donors in chain context (excluding terminal)
    acceptors_all = {}

    for code, info in results.items():
        if info.donors:
            donors_all[code] = dict(info.donors)

            # Chain context - exclude terminal donors
            chain_donors = {k: v for k, v in info.donors.items()
                           if k not in info.terminal_donors}
            if chain_donors:
                donors_chain[code] = chain_donors

        if info.acceptors:
            acceptors_all[code] = dict(info.acceptors)

    # Save all versions
    with open(output_dir / "cif_donors_all.json", 'w') as f:
        json.dump(donors_all, f, indent=2, sort_keys=True)

    with open(output_dir / "cif_donors_chain.json", 'w') as f:
        json.dump(donors_chain, f, indent=2, sort_keys=True)

    with open(output_dir / "cif_acceptors.json", 'w') as f:
        json.dump(acceptors_all, f, indent=2, sort_keys=True)

    print(f"\nSaved to {output_dir}:")
    print(f"  cif_donors_all.json - {len(donors_all)} residues")
    print(f"  cif_donors_chain.json - {len(donors_chain)} residues (chain context)")
    print(f"  cif_acceptors.json - {len(acceptors_all)} residues")

    return donors_all, donors_chain, acceptors_all


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Parse CIF files for donor/acceptor info')
    parser.add_argument('cif_dir', type=Path, help='Directory containing CIF files')
    parser.add_argument('--output-dir', type=Path, default=None,
                       help='Output directory for JSON files')
    args = parser.parse_args()

    results = process_all_cifs(args.cif_dir)
    print(f"Processed {len(results)} residues")

    # Show summary for standard bases
    print("\n" + "="*80)
    print("STANDARD BASE SUMMARY (from CIF with explicit H)")
    print("="*80)

    for base in ['A', 'G', 'C', 'U', 'T', 'DA', 'DG', 'DC', 'DT']:
        if base in results:
            res = results[base]
            print(f"\n{base}:")
            print(f"  Donors: {dict(res.donors)}")
            print(f"  Acceptors: {dict(res.acceptors)}")
            if res.terminal_donors:
                print(f"  Terminal (exclude in chain): {res.terminal_donors}")

    # Generate output
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        generate_json_output(results, args.output_dir)
    else:
        # Default output to same dir as script
        output_dir = Path(__file__).parent
        generate_json_output(results, output_dir)


if __name__ == '__main__':
    main()

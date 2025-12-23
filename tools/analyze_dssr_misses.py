#!/usr/bin/env python3
"""
Analyze why we miss H-bonds that DSSR detects.

For each DSSR H-bond we don't detect, explains:
- Was it a base-base, sugar, or backbone interaction?
- Was it detected but lost a conflict?
- Was it outside our detection range?
- What bond might have won the conflict?

Usage:
    python tools/analyze_dssr_misses.py 1EHZ
    python tools/analyze_dssr_misses.py 1EHZ --verbose
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass


def get_atom_category(atom_name: str) -> str:
    """Categorize atom as base, sugar, or backbone."""
    name = atom_name.strip().upper()
    # Sugar atoms have prime
    if "'" in name:
        return 'sugar'
    # Backbone phosphate
    if name in ['P', 'OP1', 'OP2', 'OP3', 'O1P', 'O2P', 'O3P']:
        return 'backbone'
    # Base atoms
    return 'base'


def parse_dssr_atom_id(atom_id: str) -> Tuple[str, str, str]:
    """Parse DSSR atom ID into (atom_name, chain, resnum)."""
    if '@' not in atom_id:
        return atom_id, '', ''
    atom, dssr_res = atom_id.split('@', 1)
    if '.' not in dssr_res:
        return atom, '', dssr_res
    chain, rest = dssr_res.split('.', 1)
    # Extract resnum from end
    resnum = ''
    for i in range(len(rest) - 1, -1, -1):
        if rest[i].isdigit():
            continue
        else:
            resnum = rest[i+1:]
            break
    if not resnum:
        # All digits or starts with digit
        for i, c in enumerate(rest):
            if c.isdigit():
                resnum = rest[i:]
                break
    return atom, chain, resnum


def to_dssr_res_id(res_id: str) -> str:
    """Convert our res_id to DSSR format for matching."""
    parts = res_id.split('-')
    if len(parts) < 3:
        return res_id
    chain = parts[0]
    name = '-'.join(parts[1:-1])
    num_part = parts[-1]
    # Handle insertion code
    num = ''
    insertion = ''
    for i, c in enumerate(num_part):
        if c.isdigit() or c == '-':
            num += c
        else:
            insertion = num_part[i:]
            break
    if not num:
        num = num_part
    result = f"{chain}.{name}{num}"
    if insertion:
        result += f"^{insertion}"
    return result


def normalize_atom_name(name: str) -> str:
    """Normalize atom names for matching."""
    name = name.strip()
    atom_map = {'O1P': 'OP1', 'O2P': 'OP2', 'O3P': 'OP3'}
    return atom_map.get(name, name)


def make_key(atom1_id: str, atom2_id: str) -> Tuple[str, str]:
    """Create sorted key for matching."""
    return tuple(sorted([atom1_id, atom2_id]))


@dataclass
class MissAnalysis:
    """Analysis of why we missed a DSSR H-bond."""
    dssr_hb: Dict
    reason: str
    details: str
    category: str  # base-base, base-sugar, etc.
    competing_bond: Dict = None  # Bond that might have won


def load_dssr_hbonds(pdb_id: str, project_root: Path) -> List[Dict]:
    """Load DSSR H-bonds, excluding questionable."""
    dssr_file = project_root / "data" / "json_dssr" / f"{pdb_id}.json"
    if not dssr_file.exists():
        raise FileNotFoundError(f"DSSR file not found: {dssr_file}")

    with open(dssr_file) as f:
        data = json.load(f)

    hbonds = []
    for hb in data.get('hbonds', []):
        if hb.get('donAcc_type') != 'questionable':
            hbonds.append(hb)
    return hbonds


def load_modern_hbonds_all(pdb_id: str, project_root: Path) -> Tuple[List[Dict], List[Dict]]:
    """Load modern H-bonds - both valid and invalid for analysis."""
    modern_file = project_root / "data" / "json" / "hbond_list" / f"{pdb_id}.json"
    if not modern_file.exists():
        raise FileNotFoundError(f"Modern file not found: {modern_file}")

    with open(modern_file) as f:
        data = json.load(f)

    valid = []
    invalid = []

    for group in data:
        res_id_i = group.get('res_id_i', '')
        res_id_j = group.get('res_id_j', '')
        dssr_res_i = to_dssr_res_id(res_id_i)
        dssr_res_j = to_dssr_res_id(res_id_j)

        for hb in group.get('hbonds', []):
            hb_type = hb.get('type', ' ')
            donor = normalize_atom_name(hb.get('donor_atom', ''))
            acceptor = normalize_atom_name(hb.get('acceptor_atom', ''))

            record = {
                'atom1_id': f"{donor}@{dssr_res_i}",
                'atom2_id': f"{acceptor}@{dssr_res_j}",
                'distance': hb.get('distance', 0.0),
                'type': hb_type,
            }

            if hb_type.strip():  # Valid
                valid.append(record)
            else:  # Invalid
                invalid.append(record)

    return valid, invalid


def load_base_pairs(pdb_id: str, project_root: Path) -> Set[Tuple[str, str]]:
    """Load set of base-paired residue IDs."""
    bp_file = project_root / "data" / "json" / "base_pair" / f"{pdb_id}.json"
    if not bp_file.exists():
        return set()

    with open(bp_file) as f:
        data = json.load(f)

    pairs = set()
    for bp in data:
        res_i = to_dssr_res_id(bp.get('res_id_i', ''))
        res_j = to_dssr_res_id(bp.get('res_id_j', ''))
        pairs.add((res_i, res_j))
        pairs.add((res_j, res_i))  # Both directions

    return pairs


def analyze_misses(pdb_id: str, project_root: Path, verbose: bool = False) -> List[MissAnalysis]:
    """Analyze why we miss DSSR H-bonds."""

    dssr_hbonds = load_dssr_hbonds(pdb_id, project_root)
    valid_modern, invalid_modern = load_modern_hbonds_all(pdb_id, project_root)
    base_pairs = load_base_pairs(pdb_id, project_root)

    # Build lookup sets
    valid_keys = {make_key(h['atom1_id'], h['atom2_id']) for h in valid_modern}
    invalid_lookup = {}
    for h in invalid_modern:
        key = make_key(h['atom1_id'], h['atom2_id'])
        invalid_lookup[key] = h

    # Also build lookup by residue pair to find competing bonds
    valid_by_respair = {}
    for h in valid_modern:
        atom1, _, res1 = parse_dssr_atom_id(h['atom1_id'])
        atom2, _, res2 = parse_dssr_atom_id(h['atom2_id'])
        respair = tuple(sorted([res1, res2]))
        if respair not in valid_by_respair:
            valid_by_respair[respair] = []
        valid_by_respair[respair].append(h)

    analyses = []

    for dssr_hb in dssr_hbonds:
        atom1_id = dssr_hb['atom1_id']
        atom2_id = dssr_hb['atom2_id']
        key = make_key(atom1_id, atom2_id)

        # Get atom categories
        atom1, _, res1 = parse_dssr_atom_id(atom1_id)
        atom2, _, res2 = parse_dssr_atom_id(atom2_id)
        cat1 = get_atom_category(atom1)
        cat2 = get_atom_category(atom2)
        category = '-'.join(sorted([cat1, cat2]))

        # Check if we have it as valid
        if key in valid_keys:
            continue  # We detected it, not a miss

        # Check if we detected but marked invalid
        if key in invalid_lookup:
            invalid_hb = invalid_lookup[key]
            # Find what might have won the conflict
            respair = tuple(sorted([res1, res2]))
            competing = valid_by_respair.get(respair, [])

            if competing:
                # Found a competing bond that won
                winner = min(competing, key=lambda x: x['distance'])
                analyses.append(MissAnalysis(
                    dssr_hb=dssr_hb,
                    reason="CONFLICT_LOST",
                    details=f"Lost to {winner['atom1_id']}--{winner['atom2_id']} (d={winner['distance']:.2f}Å)",
                    category=category,
                    competing_bond=winner,
                ))
            else:
                analyses.append(MissAnalysis(
                    dssr_hb=dssr_hb,
                    reason="INVALID_CLASSIFICATION",
                    details=f"Detected but classified as invalid (d={invalid_hb['distance']:.2f}Å)",
                    category=category,
                ))
        else:
            # Not detected at all - check why
            # Build DSSR-format residue IDs for lookup
            dssr_res1 = atom1_id.split('@')[1] if '@' in atom1_id else ''
            dssr_res2 = atom2_id.split('@')[1] if '@' in atom2_id else ''

            is_base_pair = (dssr_res1, dssr_res2) in base_pairs

            if not is_base_pair:
                analyses.append(MissAnalysis(
                    dssr_hb=dssr_hb,
                    reason="NOT_BASE_PAIRED",
                    details=f"Residues are not base-paired (we only check H-bonds between base pairs)",
                    category=category,
                ))
            elif category != 'base-base':
                analyses.append(MissAnalysis(
                    dssr_hb=dssr_hb,
                    reason="NON_BASE_ATOMS",
                    details=f"Base pair but {category} atoms (we detect these but may filter)",
                    category=category,
                ))
            else:
                dist = dssr_hb.get('distance', 0.0)
                if dist > 4.0:
                    analyses.append(MissAnalysis(
                        dssr_hb=dssr_hb,
                        reason="TOO_FAR",
                        details=f"Distance {dist:.2f}Å > 4.0Å cutoff",
                        category=category,
                    ))
                else:
                    analyses.append(MissAnalysis(
                        dssr_hb=dssr_hb,
                        reason="UNKNOWN",
                        details=f"Not detected, distance={dist:.2f}Å - needs investigation",
                        category=category,
                    ))

    return analyses


def print_analysis(analyses: List[MissAnalysis], verbose: bool = False):
    """Print analysis summary."""

    # Group by reason
    by_reason = {}
    for a in analyses:
        if a.reason not in by_reason:
            by_reason[a.reason] = []
        by_reason[a.reason].append(a)

    print(f"\n{'='*70}")
    print(f"DSSR MISS ANALYSIS")
    print(f"{'='*70}")
    print(f"Total misses: {len(analyses)}")
    print()

    # Summary by reason
    print("By reason:")
    for reason, items in sorted(by_reason.items(), key=lambda x: -len(x[1])):
        print(f"  {reason}: {len(items)}")
    print()

    # Summary by category
    by_category = {}
    for a in analyses:
        if a.category not in by_category:
            by_category[a.category] = []
        by_category[a.category].append(a)

    print("By atom category:")
    for cat, items in sorted(by_category.items(), key=lambda x: -len(x[1])):
        print(f"  {cat}: {len(items)}")
    print()

    # Details
    if verbose:
        for reason, items in sorted(by_reason.items(), key=lambda x: -len(x[1])):
            print(f"\n--- {reason} ({len(items)}) ---")
            for a in items[:10]:  # First 10
                hb = a.dssr_hb
                print(f"  {hb['atom1_id']} -- {hb['atom2_id']} (d={hb['distance']:.2f}Å, {hb.get('donAcc_type', '')})")
                print(f"    → {a.details}")
            if len(items) > 10:
                print(f"  ... and {len(items) - 10} more")


def main():
    parser = argparse.ArgumentParser(description="Analyze why we miss DSSR H-bonds")
    parser.add_argument("pdb_id", help="PDB ID to analyze")
    parser.add_argument("-v", "--verbose", action="store_true", help="Show detailed output")

    args = parser.parse_args()
    project_root = Path(__file__).parent.parent
    pdb_id = args.pdb_id.upper()

    try:
        analyses = analyze_misses(pdb_id, project_root, args.verbose)
        print_analysis(analyses, args.verbose)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())

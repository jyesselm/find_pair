#!/usr/bin/env python3
"""Validate cWW pairs using template overlay and H-bond patterns.

This combines:
1. Template RMSD fitting - how well does the pair geometry match idealized cWW?
2. H-bond pattern matching - does the pair have expected WC H-bonds?
"""

import json
import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
from collections import defaultdict


# Expected H-bond patterns for Watson-Crick pairs
# Format: (donor_res, donor_atom, acceptor_res, acceptor_atom)
# res is 1 or 2 (first or second residue in pair)
WC_HBOND_PATTERNS = {
    "GC": [
        (1, "O6", 2, "N4"),   # G:O6 accepts from C:N4
        (1, "N1", 2, "N3"),   # G:N1 donates to C:N3
        (1, "N2", 2, "O2"),   # G:N2 donates to C:O2
    ],
    "CG": [
        (1, "N4", 2, "O6"),   # C:N4 donates to G:O6
        (1, "N3", 2, "N1"),   # C:N3 accepts from G:N1
        (1, "O2", 2, "N2"),   # C:O2 accepts from G:N2
    ],
    "AU": [
        (1, "N6", 2, "O4"),   # A:N6 donates to U:O4
        (1, "N1", 2, "N3"),   # A:N1 accepts from U:N3
    ],
    "UA": [
        (1, "O4", 2, "N6"),   # U:O4 accepts from A:N6
        (1, "N3", 2, "N1"),   # U:N3 donates to A:N1
    ],
    "GU": [  # Wobble pair
        (1, "O6", 2, "N3"),   # G:O6 accepts from U:N3
        (1, "N1", 2, "O2"),   # G:N1 donates to U:O2
    ],
    "UG": [
        (1, "N3", 2, "O6"),   # U:N3 donates to G:O6
        (1, "O2", 2, "N1"),   # U:O2 accepts from G:N1
    ],
}

# Mapping of modified bases to standard bases for H-bond pattern lookup
MODIFIED_BASE_MAP = {
    # Guanine modifications
    "2MG": "G", "M2G": "G", "7MG": "G", "OMG": "G", "YG": "G",
    "1MG": "G", "G7M": "G", "QUO": "G",
    # Cytosine modifications
    "5MC": "C", "OMC": "C", "S4C": "C",
    # Adenine modifications
    "1MA": "A", "MIA": "A", "6MA": "A", "A2M": "A",
    # Uracil modifications
    "PSU": "U", "5MU": "U", "H2U": "U", "4SU": "U", "DHU": "U",
    "OMU": "U", "SSU": "U", "S4U": "U",
    # Thymine (DNA)
    "DT": "U",
}

def normalize_bp_type(bp_type: str) -> str:
    """Normalize base pair type, mapping modified bases to standard."""
    if len(bp_type) != 2:
        return bp_type

    b1 = bp_type[0].upper()
    b2 = bp_type[1].upper()

    # Map modified bases to standard
    b1 = MODIFIED_BASE_MAP.get(b1, b1) if len(b1) > 1 else MODIFIED_BASE_MAP.get(b1.upper(), b1.upper())
    b2 = MODIFIED_BASE_MAP.get(b2, b2) if len(b2) > 1 else MODIFIED_BASE_MAP.get(b2.upper(), b2.upper())

    # Handle single lowercase letters (modified notation)
    if b1.islower():
        b1 = MODIFIED_BASE_MAP.get(bp_type[0].upper(), bp_type[0].upper())
    if b2.islower():
        b2 = MODIFIED_BASE_MAP.get(bp_type[1].upper(), bp_type[1].upper())

    # Single letter bases
    if len(bp_type[0]) == 1 and bp_type[0].islower():
        b1 = bp_type[0].upper()  # g -> G, c -> C, etc.
    if len(bp_type[1]) == 1 and bp_type[1].islower():
        b2 = bp_type[1].upper()

    return b1 + b2

# N1/N9 reference atoms for distance calculation
N1N9_ATOMS = {
    "A": "N9", "G": "N9", "C": "N1", "U": "N1",
    "DA": "N9", "DG": "N9", "DC": "N1", "DT": "N1",
}

# Expected N1N9 distances for cWW pairs (from DSSR statistics)
N1N9_RANGES = {
    "GC": (8.5, 9.5),
    "CG": (8.5, 9.5),
    "AU": (8.5, 9.5),
    "UA": (8.5, 9.5),
    "GU": (8.5, 9.5),
    "UG": (8.5, 9.5),
}


@dataclass
class Atom:
    """An atom with coordinates."""
    name: str
    coords: np.ndarray


@dataclass
class Residue:
    """A nucleotide residue."""
    res_id: str
    base_type: str
    atoms: Dict[str, Atom] = field(default_factory=dict)

    @property
    def n1n9_atom(self) -> Optional[Atom]:
        """Get N1 or N9 atom depending on base type."""
        if self.base_type in ("A", "G"):
            return self.atoms.get("N9")
        else:
            return self.atoms.get("N1")


@dataclass
class HBond:
    """A hydrogen bond."""
    donor_atom: str
    acceptor_atom: str
    distance: float


@dataclass
class ValidationResult:
    """Result of validating a cWW pair."""
    res_id1: str
    res_id2: str
    bp_type: str

    # Template fitting
    template_rmsd: Optional[float] = None

    # N1N9 distance
    n1n9_dist: Optional[float] = None
    n1n9_in_range: bool = False

    # H-bond validation
    expected_hbonds: int = 0
    found_hbonds: int = 0
    hbond_details: List[dict] = field(default_factory=list)

    # Overall score
    is_valid_cww: bool = False
    quality_score: float = 0.0

    def __repr__(self):
        rmsd_str = f"{self.template_rmsd:.3f}" if self.template_rmsd else "N/A"
        n1n9_str = f"{self.n1n9_dist:.2f}" if self.n1n9_dist else "N/A"
        return (f"ValidationResult({self.res_id1}-{self.res_id2} {self.bp_type}: "
                f"N1N9={n1n9_str}Å, "
                f"hbonds={self.found_hbonds}/{self.expected_hbonds}, "
                f"valid={self.is_valid_cww})")


def parse_pdb_atoms(pdb_path: Path) -> Dict[str, Residue]:
    """Parse PDB file and return residues with atoms."""
    residues = {}

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM") and not line.startswith("HETATM"):
                continue

            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21]
            res_seq = line[22:26].strip()
            ins_code = line[26].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            res_id = f"{chain}-{res_name}-{res_seq}{ins_code}".rstrip()

            if res_id not in residues:
                # Get single-letter base type
                base_type = res_name[-1] if len(res_name) <= 3 else res_name[0]
                residues[res_id] = Residue(res_id=res_id, base_type=base_type)

            residues[res_id].atoms[atom_name] = Atom(
                name=atom_name,
                coords=np.array([x, y, z])
            )

    return residues


def load_modern_hbonds(hbond_path: Path) -> Dict[Tuple[str, str], List[HBond]]:
    """Load H-bonds from modern JSON output."""
    if not hbond_path.exists():
        return {}

    with open(hbond_path) as f:
        data = json.load(f)

    result = {}
    for entry in data:
        res_id1 = entry.get("res_id_i", "")
        res_id2 = entry.get("res_id_j", "")

        hbonds = []
        for h in entry.get("hbonds", []):
            if h.get("context") != "base_base":
                continue
            hbonds.append(HBond(
                donor_atom=h.get("donor_atom", ""),
                acceptor_atom=h.get("acceptor_atom", ""),
                distance=h.get("distance", 0.0),
            ))

        if hbonds:
            key = tuple(sorted([res_id1, res_id2]))
            result[key] = hbonds

    return result


def load_modern_pairs(pair_path: Path) -> List[dict]:
    """Load base pairs from modern JSON output."""
    if not pair_path.exists():
        return []

    with open(pair_path) as f:
        return json.load(f)


def calculate_n1n9_distance(res1: Residue, res2: Residue) -> Optional[float]:
    """Calculate N1-N9 or N1-N1 distance between two bases."""
    atom1 = res1.n1n9_atom
    atom2 = res2.n1n9_atom

    if atom1 is None or atom2 is None:
        return None

    return float(np.linalg.norm(atom2.coords - atom1.coords))


def check_hbond_pattern(
    bp_type: str,
    hbonds: List[HBond],
    distance_tolerance: float = 0.5,
) -> Tuple[int, int, List[dict]]:
    """Check if H-bonds match expected WC pattern.

    Returns: (expected_count, found_count, details)
    """
    # Normalize bp_type to handle modified bases
    normalized_type = normalize_bp_type(bp_type)
    pattern = WC_HBOND_PATTERNS.get(normalized_type, [])
    if not pattern:
        return 0, 0, []

    expected = len(pattern)
    found = 0
    details = []

    for donor_res, donor_atom, acceptor_res, acceptor_atom in pattern:
        # Look for matching H-bond in either direction
        matched = False
        matched_dist = None

        for hb in hbonds:
            # Check forward match
            if hb.donor_atom == donor_atom and hb.acceptor_atom == acceptor_atom:
                matched = True
                matched_dist = hb.distance
                break
            # Check reverse (in case donor/acceptor labeling differs)
            if hb.donor_atom == acceptor_atom and hb.acceptor_atom == donor_atom:
                matched = True
                matched_dist = hb.distance
                break

        details.append({
            "expected": f"{donor_atom}-{acceptor_atom}",
            "found": matched,
            "distance": matched_dist,
        })

        if matched:
            found += 1

    return expected, found, details


def validate_cww_pair(
    res1: Residue,
    res2: Residue,
    hbonds: List[HBond],
    bp_type: str,
) -> ValidationResult:
    """Validate a potential cWW pair using geometry and H-bonds."""
    result = ValidationResult(
        res_id1=res1.res_id,
        res_id2=res2.res_id,
        bp_type=bp_type,
    )

    # Calculate N1N9 distance
    result.n1n9_dist = calculate_n1n9_distance(res1, res2)
    if result.n1n9_dist is not None:
        n1n9_range = N1N9_RANGES.get(bp_type, (8.0, 10.0))
        result.n1n9_in_range = n1n9_range[0] <= result.n1n9_dist <= n1n9_range[1]

    # Check H-bond pattern
    expected, found, details = check_hbond_pattern(bp_type, hbonds)
    result.expected_hbonds = expected
    result.found_hbonds = found
    result.hbond_details = details

    # Calculate quality score
    hbond_score = found / expected if expected > 0 else 0
    n1n9_score = 1.0 if result.n1n9_in_range else 0.5
    result.quality_score = 0.6 * hbond_score + 0.4 * n1n9_score

    # Determine if valid cWW
    # Require at least 2 H-bonds and reasonable N1N9 distance
    result.is_valid_cww = (
        found >= 2 and
        result.n1n9_dist is not None and
        7.5 <= result.n1n9_dist <= 10.5
    )

    return result


def validate_pdb_pairs(
    pdb_id: str,
    pdb_dir: Path,
    pair_dir: Path,
    hbond_dir: Path,
    verbose: bool = False,
) -> List[ValidationResult]:
    """Validate all pairs in a PDB file."""
    # Load PDB structure
    pdb_path = pdb_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"
    if not pdb_path.exists():
        return []

    residues = parse_pdb_atoms(pdb_path)

    # Load pairs and H-bonds
    pairs = load_modern_pairs(pair_dir / f"{pdb_id}.json")
    hbonds = load_modern_hbonds(hbond_dir / f"{pdb_id}.json")

    results = []

    for pair in pairs:
        res_id1 = pair.get("res_id_i", "")
        res_id2 = pair.get("res_id_j", "")
        bp_type = pair.get("bp_type", "")

        if res_id1 not in residues or res_id2 not in residues:
            continue

        res1 = residues[res_id1]
        res2 = residues[res_id2]

        # Get H-bonds for this pair
        pair_key = tuple(sorted([res_id1, res_id2]))
        pair_hbonds = hbonds.get(pair_key, [])

        result = validate_cww_pair(res1, res2, pair_hbonds, bp_type)
        results.append(result)

        if verbose:
            status = "✓" if result.is_valid_cww else "✗"
            print(f"  {status} {result}")

    return results


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate cWW pairs using templates and H-bonds"
    )
    parser.add_argument(
        "--pdb-dir", type=Path,
        default=Path("data/pdb"),
        help="Directory with PDB files"
    )
    parser.add_argument(
        "--pair-dir", type=Path,
        default=Path("data/json/base_pair"),
        help="Directory with base_pair JSON files"
    )
    parser.add_argument(
        "--hbond-dir", type=Path,
        default=Path("data/json/all_hbond_list"),
        help="Directory with H-bond JSON files"
    )
    parser.add_argument(
        "--pdb", type=str, required=True,
        help="PDB ID to analyze"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Verbose output"
    )

    args = parser.parse_args()

    print(f"Validating cWW pairs for {args.pdb}")

    results = validate_pdb_pairs(
        args.pdb,
        args.pdb_dir,
        args.pair_dir,
        args.hbond_dir,
        verbose=args.verbose
    )

    if not results:
        print("No pairs found")
        return

    # Summary
    valid = sum(1 for r in results if r.is_valid_cww)
    total = len(results)

    print(f"\nSummary: {valid}/{total} pairs validated as cWW")

    # Group by bp_type
    by_type = defaultdict(list)
    for r in results:
        by_type[r.bp_type].append(r)

    print("\nBy base pair type:")
    for bp_type, type_results in sorted(by_type.items()):
        valid_count = sum(1 for r in type_results if r.is_valid_cww)
        avg_hbonds = np.mean([r.found_hbonds for r in type_results])
        print(f"  {bp_type}: {valid_count}/{len(type_results)} valid, "
              f"avg {avg_hbonds:.1f} H-bonds")


if __name__ == "__main__":
    main()

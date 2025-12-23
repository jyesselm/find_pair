"""
DSSR JSON comparison utilities.

Compares hydrogen bonds between our modern output and DSSR output.
"""

import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field


def to_dssr_res_id(res_id: str) -> str:
    """Convert our res_id format to DSSR format.

    Examples:
        "A-G-103" -> "A.G103"
        "A-C-10A" -> "A.C10^A"
        "B-2MG-25" -> "B.2MG25"
    """
    parts = res_id.split('-')
    if len(parts) < 3:
        return res_id  # Invalid format

    chain = parts[0]
    name = '-'.join(parts[1:-1])  # Handle names with dashes like "2MG"
    num_part = parts[-1]

    # Check for insertion code (letter at end)
    num = ""
    insertion = ""
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


def normalize_atom_name(atom_name: str) -> str:
    """Normalize atom name to DSSR convention.

    DSSR uses: OP1, OP2 (phosphate oxygens)
    PDB uses:  O1P, O2P (older naming)
    """
    name = atom_name.strip()
    # Map old PDB names to DSSR names
    atom_map = {
        'O1P': 'OP1',
        'O2P': 'OP2',
        'O3P': 'OP3',
        "O5'": "O5'",
        "O4'": "O4'",
        "O3'": "O3'",
        "O2'": "O2'",
        "C5'": "C5'",
        "C4'": "C4'",
        "C3'": "C3'",
        "C2'": "C2'",
        "C1'": "C1'",
    }
    return atom_map.get(name, name)


def to_dssr_atom_id(atom_name: str, res_id: str) -> str:
    """Create DSSR atom ID from atom name and our res_id.

    Example: "N7", "A-G-103" -> "N7@A.G103"
    """
    normalized = normalize_atom_name(atom_name)
    return f"{normalized}@{to_dssr_res_id(res_id)}"


def parse_dssr_atom_id(atom_id: str) -> Tuple[str, str]:
    """Parse DSSR atom ID into (atom_name, dssr_res_id).

    Example: "N7@A.G103" -> ("N7", "A.G103")
    """
    if '@' not in atom_id:
        return atom_id, ""
    atom_name, dssr_res_id = atom_id.split('@', 1)
    return atom_name, dssr_res_id


@dataclass
class HBondMatch:
    """Result of matching a single hydrogen bond."""
    dssr_atom1_id: str
    dssr_atom2_id: str
    dssr_distance: float
    modern_donor: Optional[str] = None
    modern_acceptor: Optional[str] = None
    modern_distance: Optional[float] = None
    matched: bool = False
    distance_diff: Optional[float] = None


@dataclass
class DSSRComparisonResult:
    """Result of comparing DSSR and modern hbond outputs."""
    pdb_id: str
    dssr_count: int = 0
    modern_count: int = 0
    matched_count: int = 0
    dssr_only: List[Dict] = field(default_factory=list)
    modern_only: List[Dict] = field(default_factory=list)
    distance_mismatches: List[Dict] = field(default_factory=list)

    @property
    def match_rate(self) -> float:
        if self.dssr_count == 0:
            return 0.0
        return self.matched_count / self.dssr_count


def flatten_modern_hbonds(modern_data: List[Dict]) -> List[Dict]:
    """Flatten modern grouped hbonds to individual records.

    Converts from grouped format (by residue pair) to flat list.
    Each record will have: donor_atom_id, acceptor_atom_id, distance
    """
    flat = []
    for group in modern_data:
        res_id_i = group.get('res_id_i', '')
        res_id_j = group.get('res_id_j', '')

        for hb in group.get('hbonds', []):
            donor_atom = hb.get('donor_atom', '').strip()
            acceptor_atom = hb.get('acceptor_atom', '').strip()
            distance = hb.get('distance', 0.0)

            # Create DSSR-format atom IDs
            # donor comes from res_i, acceptor from res_j (based on our convention)
            donor_atom_id = to_dssr_atom_id(donor_atom, res_id_i)
            acceptor_atom_id = to_dssr_atom_id(acceptor_atom, res_id_j)

            flat.append({
                'donor_atom_id': donor_atom_id,
                'acceptor_atom_id': acceptor_atom_id,
                'distance': distance,
                'res_id_i': res_id_i,
                'res_id_j': res_id_j,
                'donor_atom': donor_atom,
                'acceptor_atom': acceptor_atom,
            })

    return flat


def compare_hbonds(
    dssr_data: Dict,
    modern_data: List[Dict],
    distance_tolerance: float = 0.1
) -> DSSRComparisonResult:
    """Compare DSSR hbonds with our modern hbond output.

    Args:
        dssr_data: Full DSSR JSON (contains 'hbonds' key)
        modern_data: Our hbond_list JSON (list of grouped hbonds)
        distance_tolerance: Max allowed distance difference

    Returns:
        DSSRComparisonResult with match statistics
    """
    result = DSSRComparisonResult(pdb_id="")

    dssr_hbonds = dssr_data.get('hbonds', [])
    result.dssr_count = len(dssr_hbonds)

    # Flatten modern hbonds
    modern_flat = flatten_modern_hbonds(modern_data)
    result.modern_count = len(modern_flat)

    # Build lookup for modern hbonds by atom pair (both orderings)
    modern_by_atoms = {}
    for hb in modern_flat:
        key1 = (hb['donor_atom_id'], hb['acceptor_atom_id'])
        key2 = (hb['acceptor_atom_id'], hb['donor_atom_id'])
        modern_by_atoms[key1] = hb
        modern_by_atoms[key2] = hb

    matched_modern = set()

    # Match DSSR hbonds to modern
    for dssr_hb in dssr_hbonds:
        atom1_id = dssr_hb.get('atom1_id', '')
        atom2_id = dssr_hb.get('atom2_id', '')
        dssr_dist = dssr_hb.get('distance', 0.0)

        key = (atom1_id, atom2_id)
        modern_hb = modern_by_atoms.get(key)

        if modern_hb:
            matched_modern.add((modern_hb['donor_atom_id'], modern_hb['acceptor_atom_id']))
            result.matched_count += 1

            # Check distance
            dist_diff = abs(dssr_dist - modern_hb['distance'])
            if dist_diff > distance_tolerance:
                result.distance_mismatches.append({
                    'atom1_id': atom1_id,
                    'atom2_id': atom2_id,
                    'dssr_distance': dssr_dist,
                    'modern_distance': modern_hb['distance'],
                    'diff': dist_diff,
                })
        else:
            result.dssr_only.append({
                'atom1_id': atom1_id,
                'atom2_id': atom2_id,
                'distance': dssr_dist,
                'donAcc_type': dssr_hb.get('donAcc_type', ''),
            })

    # Find modern-only hbonds
    for hb in modern_flat:
        key = (hb['donor_atom_id'], hb['acceptor_atom_id'])
        if key not in matched_modern:
            result.modern_only.append({
                'donor_atom_id': hb['donor_atom_id'],
                'acceptor_atom_id': hb['acceptor_atom_id'],
                'distance': hb['distance'],
            })

    return result


def compare_hbonds_with_pairs(
    dssr_data: Dict,
    modern_data: List[Dict],
    distance_tolerance: float = 0.1
) -> DSSRComparisonResult:
    """Compare modern hbonds against DSSR hbonds + pairs.

    This combines DSSR's hbonds array with H-bonds extracted from pairs.
    """
    result = DSSRComparisonResult(pdb_id="")

    # Get DSSR hbonds from both sources
    dssr_hbonds = dssr_data.get('hbonds', [])
    pair_hbonds = extract_pair_hbonds(dssr_data)

    # Combine (pairs hbonds take precedence as they're base-pair specific)
    all_dssr = []
    seen = set()
    for hb in pair_hbonds:
        key = (hb['atom1_id'], hb['atom2_id'])
        if key not in seen:
            seen.add(key)
            seen.add((hb['atom2_id'], hb['atom1_id']))
            all_dssr.append(hb)
    for hb in dssr_hbonds:
        key = (hb['atom1_id'], hb['atom2_id'])
        if key not in seen:
            seen.add(key)
            seen.add((hb['atom2_id'], hb['atom1_id']))
            all_dssr.append(hb)

    result.dssr_count = len(all_dssr)

    # Flatten modern hbonds
    modern_flat = flatten_modern_hbonds(modern_data)
    result.modern_count = len(modern_flat)

    # Build lookup for modern hbonds
    modern_by_atoms = {}
    for hb in modern_flat:
        key1 = (hb['donor_atom_id'], hb['acceptor_atom_id'])
        key2 = (hb['acceptor_atom_id'], hb['donor_atom_id'])
        modern_by_atoms[key1] = hb
        modern_by_atoms[key2] = hb

    matched_modern = set()

    for dssr_hb in all_dssr:
        atom1_id = dssr_hb.get('atom1_id', '')
        atom2_id = dssr_hb.get('atom2_id', '')
        dssr_dist = dssr_hb.get('distance', 0.0)

        key = (atom1_id, atom2_id)
        modern_hb = modern_by_atoms.get(key)

        if modern_hb:
            matched_modern.add((modern_hb['donor_atom_id'], modern_hb['acceptor_atom_id']))
            result.matched_count += 1

            dist_diff = abs(dssr_dist - modern_hb['distance'])
            if dist_diff > distance_tolerance:
                result.distance_mismatches.append({
                    'atom1_id': atom1_id,
                    'atom2_id': atom2_id,
                    'dssr_distance': dssr_dist,
                    'modern_distance': modern_hb['distance'],
                    'diff': dist_diff,
                })
        else:
            result.dssr_only.append({
                'atom1_id': atom1_id,
                'atom2_id': atom2_id,
                'distance': dssr_dist,
                'donAcc_type': dssr_hb.get('donAcc_type', dssr_hb.get('source', '')),
            })

    for hb in modern_flat:
        key = (hb['donor_atom_id'], hb['acceptor_atom_id'])
        if key not in matched_modern:
            result.modern_only.append({
                'donor_atom_id': hb['donor_atom_id'],
                'acceptor_atom_id': hb['acceptor_atom_id'],
                'distance': hb['distance'],
            })

    return result


def compare_pdb(pdb_id: str, project_root: Path, include_pairs: bool = True, use_all_hbonds: bool = False) -> DSSRComparisonResult:
    """Compare DSSR and modern hbonds for a single PDB.

    Args:
        pdb_id: PDB ID (e.g., "1GID")
        project_root: Project root directory
        include_pairs: If True, also match against DSSR pairs section
        use_all_hbonds: If True, use all_hbond_list instead of hbond_list

    Returns:
        DSSRComparisonResult
    """
    dssr_file = project_root / "data" / "json_dssr" / f"{pdb_id}.json"

    # Choose which modern output to use
    if use_all_hbonds:
        modern_file = project_root / "data" / "json" / "all_hbond_list" / f"{pdb_id}.json"
    else:
        modern_file = project_root / "data" / "json" / "hbond_list" / f"{pdb_id}.json"

    if not dssr_file.exists():
        raise FileNotFoundError(f"DSSR JSON not found: {dssr_file}")
    if not modern_file.exists():
        raise FileNotFoundError(f"Modern {'all_hbond_list' if use_all_hbonds else 'hbond_list'} not found: {modern_file}")

    with open(dssr_file) as f:
        dssr_data = json.load(f)
    with open(modern_file) as f:
        modern_data = json.load(f)

    if include_pairs:
        result = compare_hbonds_with_pairs(dssr_data, modern_data)
    else:
        result = compare_hbonds(dssr_data, modern_data)
    result.pdb_id = pdb_id
    return result


def extract_pair_hbonds(dssr_data: Dict) -> List[Dict]:
    """Extract H-bonds from DSSR pairs section.

    DSSR stores base pair H-bonds in 'pairs' array under 'hbonds_desc' field.
    Format: "O6(carbonyl)-N4(amino)[2.81],N1(imino)-N3[2.92]"
    """
    hbonds = []
    for pair in dssr_data.get('pairs', []):
        nt1 = pair.get('nt1', '')  # e.g., "A.G103"
        nt2 = pair.get('nt2', '')  # e.g., "A.C217"
        hbonds_desc = pair.get('hbonds_desc', '')

        if not hbonds_desc:
            continue

        # Parse hbonds_desc: "O6(carbonyl)-N4(amino)[2.81],N1(imino)-N3[2.92]"
        for hb_str in hbonds_desc.split(','):
            hb_str = hb_str.strip()
            if not hb_str:
                continue

            # Extract atom names and distance
            # Format: ATOM1(type)-ATOM2(type)[distance] or ATOM1-ATOM2[distance]
            import re
            match = re.match(r'(\w+)(?:\([^)]+\))?-(\w+)(?:\([^)]+\))?\[([0-9.]+)\]', hb_str)
            if match:
                atom1, atom2, dist = match.groups()
                hbonds.append({
                    'atom1_id': f"{atom1}@{nt1}",
                    'atom2_id': f"{atom2}@{nt2}",
                    'distance': float(dist),
                    'source': 'pairs',
                })

    return hbonds


def print_comparison_summary(result: DSSRComparisonResult, verbose: bool = False):
    """Print comparison summary."""
    print(f"\n{'='*60}")
    print(f"DSSR vs Modern H-bond Comparison: {result.pdb_id}")
    print(f"{'='*60}")
    print(f"DSSR total:      {result.dssr_count} (hbonds + pairs)")
    print(f"Modern total:    {result.modern_count}")
    print(f"Matched:         {result.matched_count} ({result.match_rate:.1%} of DSSR)")
    print(f"DSSR only:       {len(result.dssr_only)} (backbone/tertiary)")
    print(f"Modern only:     {len(result.modern_only)} (additional weak bonds)")
    print(f"Distance diffs:  {len(result.distance_mismatches)}")

    if verbose:
        if result.dssr_only:
            print(f"\nDSSR-only hbonds (first 10):")
            for hb in result.dssr_only[:10]:
                print(f"  {hb['atom1_id']} -- {hb['atom2_id']} ({hb['distance']:.3f}Å) [{hb['donAcc_type']}]")

        if result.modern_only:
            print(f"\nModern-only hbonds (first 10):")
            for hb in result.modern_only[:10]:
                print(f"  {hb['donor_atom_id']} -- {hb['acceptor_atom_id']} ({hb['distance']:.3f}Å)")

        if result.distance_mismatches:
            print(f"\nDistance mismatches (first 10):")
            for hb in result.distance_mismatches[:10]:
                print(f"  {hb['atom1_id']} -- {hb['atom2_id']}: DSSR={hb['dssr_distance']:.3f} vs Modern={hb['modern_distance']:.3f} (diff={hb['diff']:.3f})")


if __name__ == "__main__":
    import sys

    project_root = Path(__file__).parent.parent
    pdb_id = sys.argv[1] if len(sys.argv) > 1 else "1GID"
    verbose = "-v" in sys.argv or "--verbose" in sys.argv
    use_all_hbonds = "--all" in sys.argv or "-a" in sys.argv

    try:
        result = compare_pdb(pdb_id, project_root, use_all_hbonds=use_all_hbonds)
        if use_all_hbonds:
            print("(Using all_hbond_list - all H-bonds mode)")
        print_comparison_summary(result, verbose)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

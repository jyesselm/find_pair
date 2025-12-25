#!/usr/bin/env python3
"""
DSSR H-Bond Comparison PyMOL Visualizer.

Compares H-bond detection between our code and DSSR, generates PyMOL script
for visual inspection of differences.

Supports two modes:
1. C++ mode (default): Compare modern C++ output with DSSR
2. Prototype mode (--use-prototype): Use Python H-bond optimizer prototype
   with detailed miss reason analysis

Usage:
    python tools/dssr_pymol_visualizer.py 1GID
    python tools/dssr_pymol_visualizer.py 1GID --use-prototype
    python tools/dssr_pymol_visualizer.py 1GID --output 1GID_hbonds.pml --verbose
    pymol data/pdb/1GID.pdb 1GID_hbonds.pml
"""

import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

# Add project root to path for prototype imports
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def normalize_atom_name(atom_name: str) -> str:
    """Normalize atom name (strip whitespace, map old PDB names)."""
    name = atom_name.strip()
    atom_map = {'O1P': 'OP1', 'O2P': 'OP2', 'O3P': 'OP3'}
    return atom_map.get(name, name)


def parse_dssr_atom_id(atom_id: str) -> Tuple[str, str, str, str]:
    """Parse DSSR atom ID into (atom_name, chain, resname, resnum).

    Example: "N7@A.G103" -> ("N7", "A", "G", "103")
    Example: "O2'@B.C25^A" -> ("O2'", "B", "C", "25A")
    Example: "OP2@A.5MC49" -> ("OP2", "A", "5MC", "49")
    """
    if '@' not in atom_id:
        return atom_id, "", "", ""

    atom_name, dssr_res_id = atom_id.split('@', 1)

    # Parse chain.resnameresnum format
    if '.' not in dssr_res_id:
        return atom_name, "", "", dssr_res_id

    chain, rest = dssr_res_id.split('.', 1)

    # Handle insertion code (^A format)
    insertion = ""
    if '^' in rest:
        rest, insertion = rest.split('^', 1)

    # Extract resname and resnum - find last contiguous digit sequence
    # This handles cases like "5MC49" -> ("5MC", "49")
    resname = ""
    resnum = ""

    # Find the last digit sequence
    last_digit_start = -1
    for i in range(len(rest) - 1, -1, -1):
        if rest[i].isdigit():
            if last_digit_start == -1:
                last_digit_start = i
        else:
            if last_digit_start != -1:
                # Found start of digit sequence
                resname = rest[:i+1]
                resnum = rest[i+1:]
                break
    else:
        # All digits or no digits found
        if last_digit_start != -1:
            # Find where digits start from the beginning
            for i, c in enumerate(rest):
                if c.isdigit():
                    resname = rest[:i] if i > 0 else ""
                    resnum = rest[i:]
                    break
        else:
            resname = rest

    if insertion:
        resnum += insertion

    return atom_name, chain, resname, resnum


def to_pymol_selection(pdb_id: str, atom_id: str) -> str:
    """Convert DSSR atom ID to PyMOL selection string.

    Example: "N7@A.G103" -> "/1GID//A/103/N7"
    """
    atom_name, chain, resname, resnum = parse_dssr_atom_id(atom_id)
    return f"/{pdb_id}//{chain}/{resnum}/{atom_name}"


def to_dssr_res_id(res_id: str) -> str:
    """Convert our res_id format to DSSR format.

    Examples:
        "A-G-103" -> "A.G103"
        "A-C-10A" -> "A.C10^A"
    """
    parts = res_id.split('-')
    if len(parts) < 3:
        return res_id

    chain = parts[0]
    name = '-'.join(parts[1:-1])
    num_part = parts[-1]

    # Check for insertion code
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


def to_dssr_atom_id(atom_name: str, res_id: str) -> str:
    """Create DSSR atom ID from atom name and our res_id."""
    normalized = normalize_atom_name(atom_name)
    return f"{normalized}@{to_dssr_res_id(res_id)}"


@dataclass
class HBond:
    """Normalized H-bond record."""
    atom1_id: str  # DSSR format: "N7@A.G103"
    atom2_id: str
    distance: float
    source: str  # 'dssr', 'modern', 'dssr_pairs'
    don_acc_type: str = ""  # 'standard', 'acceptable', 'questionable'


def load_dssr_hbonds(pdb_id: str, project_root: Path,
                     exclude_questionable: bool = True) -> List[HBond]:
    """Load H-bonds from DSSR JSON, optionally excluding questionable."""
    dssr_file = project_root / "data" / "json_dssr" / f"{pdb_id}.json"
    if not dssr_file.exists():
        raise FileNotFoundError(f"DSSR JSON not found: {dssr_file}")

    with open(dssr_file) as f:
        dssr_data = json.load(f)

    hbonds = []

    # Get H-bonds from 'hbonds' array
    for hb in dssr_data.get('hbonds', []):
        don_acc_type = hb.get('donAcc_type', '')
        if exclude_questionable and don_acc_type == 'questionable':
            continue

        hbonds.append(HBond(
            atom1_id=hb.get('atom1_id', ''),
            atom2_id=hb.get('atom2_id', ''),
            distance=hb.get('distance', 0.0),
            source='dssr',
            don_acc_type=don_acc_type,
        ))

    # Also get H-bonds from pairs section
    for pair in dssr_data.get('pairs', []):
        nt1 = pair.get('nt1', '')
        nt2 = pair.get('nt2', '')
        hbonds_desc = pair.get('hbonds_desc', '')

        if not hbonds_desc:
            continue

        for hb_str in hbonds_desc.split(','):
            hb_str = hb_str.strip()
            if not hb_str:
                continue

            match = re.match(r"(\w+[']?)(?:\([^)]+\))?-(\w+[']?)(?:\([^)]+\))?\[([0-9.]+)\]", hb_str)
            if match:
                atom1, atom2, dist = match.groups()
                hbonds.append(HBond(
                    atom1_id=f"{atom1}@{nt1}",
                    atom2_id=f"{atom2}@{nt2}",
                    distance=float(dist),
                    source='dssr_pairs',
                    don_acc_type='pairs',
                ))

    return hbonds


def load_modern_hbonds(pdb_id: str, project_root: Path) -> List[HBond]:
    """Load H-bonds from modern hbond_list JSON.

    Handles both old format (type='-'/'*'/' ') and new format (classification='standard'/etc).
    For old format: includes STANDARD (type='-') and NON_STANDARD (type='*').
    For new format: includes all H-bonds with classification field.
    """
    # Try all_hbond_list first, fall back to hbond_list
    modern_file = project_root / "data" / "json" / "all_hbond_list" / f"{pdb_id}.json"
    if not modern_file.exists():
        modern_file = project_root / "data" / "json" / "hbond_list" / f"{pdb_id}.json"
    if not modern_file.exists():
        raise FileNotFoundError(f"Modern H-bond list not found for {pdb_id}")

    with open(modern_file) as f:
        modern_data = json.load(f)

    hbonds = []
    for group in modern_data:
        res_id_i = group.get('res_id_i', '')
        res_id_j = group.get('res_id_j', '')

        for hb in group.get('hbonds', []):
            # Check for new format (classification field)
            classification = hb.get('classification', '')
            # Check for old format (type field)
            hb_type = hb.get('type', '')

            # Determine if this is a valid H-bond
            if classification:
                # New format: include all with classification
                don_acc_type = classification
            elif hb_type.strip():
                # Old format: only include '-' or '*'
                don_acc_type = 'standard' if hb_type == '-' else 'non_standard'
            else:
                # Skip invalid (no classification and empty type)
                continue

            donor_atom = hb.get('donor_atom', '').strip()
            acceptor_atom = hb.get('acceptor_atom', '').strip()
            distance = hb.get('distance', 0.0)

            donor_atom_id = to_dssr_atom_id(donor_atom, res_id_i)
            acceptor_atom_id = to_dssr_atom_id(acceptor_atom, res_id_j)

            hbonds.append(HBond(
                atom1_id=donor_atom_id,
                atom2_id=acceptor_atom_id,
                distance=distance,
                source='modern',
                don_acc_type=don_acc_type,
            ))

    return hbonds


def make_key(atom1_id: str, atom2_id: str) -> Tuple[str, str]:
    """Create normalized key for H-bond lookup (sorted order)."""
    return tuple(sorted([atom1_id, atom2_id]))


def match_hbonds(dssr_hbonds: List[HBond], modern_hbonds: List[HBond]) -> Tuple[List[HBond], List[HBond], List[Tuple[HBond, HBond]]]:
    """Match H-bonds between DSSR and modern.

    Returns:
        (dssr_only, modern_only, matched_pairs)
    """
    # Build lookup for DSSR (deduplicated)
    dssr_by_key = {}
    for hb in dssr_hbonds:
        key = make_key(hb.atom1_id, hb.atom2_id)
        if key not in dssr_by_key:
            dssr_by_key[key] = hb

    # Build lookup for modern (deduplicated)
    modern_by_key = {}
    for hb in modern_hbonds:
        key = make_key(hb.atom1_id, hb.atom2_id)
        if key not in modern_by_key:
            modern_by_key[key] = hb

    dssr_only = []
    modern_only = []
    matched = []

    # Find DSSR-only and matched
    for key, dssr_hb in dssr_by_key.items():
        if key in modern_by_key:
            matched.append((dssr_hb, modern_by_key[key]))
        else:
            dssr_only.append(dssr_hb)

    # Find modern-only
    for key, modern_hb in modern_by_key.items():
        if key not in dssr_by_key:
            modern_only.append(modern_hb)

    return dssr_only, modern_only, matched


def generate_pymol_script(pdb_id: str,
                          dssr_only: List[HBond],
                          modern_only: List[HBond],
                          matched: List[Tuple[HBond, HBond]],
                          pdb_path: str,
                          show_matched: bool = False) -> str:
    """Generate PyMOL script for visualization."""
    lines = []
    lines.append(f"# DSSR vs Modern H-bond Comparison: {pdb_id}")
    lines.append(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"#")
    lines.append(f"# Usage: pymol {pdb_path} {pdb_id}_hbonds.pml")
    lines.append(f"#")
    lines.append(f"# Summary:")
    lines.append(f"#   DSSR-only (RED):    {len(dssr_only)} - H-bonds we miss")
    lines.append(f"#   Modern-only (BLUE): {len(modern_only)} - Extra detections")
    lines.append(f"#   Matched (GREEN):    {len(matched)} - Both agree")
    lines.append(f"#")
    lines.append(f"# Color scheme:")
    lines.append(f"#   RED   = DSSR-only (we miss these)")
    lines.append(f"#   BLUE  = Modern-only (we detect, DSSR doesn't)")
    lines.append(f"#   GREEN = Matched (both agree)")
    lines.append("")

    # Setup display
    lines.append("# Setup display")
    lines.append("hide all")
    lines.append("show sticks")
    lines.append(f"color white, {pdb_id}")
    lines.append("set dash_gap, 0.2")
    lines.append("set dash_radius, 0.1")
    lines.append("set label_size, 12")
    lines.append("set label_color, white")
    lines.append("")

    # DSSR-only H-bonds (RED)
    if dssr_only:
        lines.append(f"# DSSR-ONLY H-bonds (RED) - {len(dssr_only)} bonds we miss")
        for i, hb in enumerate(dssr_only, 1):
            sel1 = to_pymol_selection(pdb_id, hb.atom1_id)
            sel2 = to_pymol_selection(pdb_id, hb.atom2_id)
            atom1, _, _, res1 = parse_dssr_atom_id(hb.atom1_id)
            atom2, _, _, res2 = parse_dssr_atom_id(hb.atom2_id)
            lines.append(f"# {hb.atom1_id} -- {hb.atom2_id} (dist={hb.distance:.2f}, type={hb.don_acc_type})")
            lines.append(f"distance dssr_{i}, {sel1}, {sel2}")
            lines.append(f"color red, dssr_{i}")
        lines.append("")

    # Modern-only H-bonds (BLUE)
    if modern_only:
        lines.append(f"# MODERN-ONLY H-bonds (BLUE) - {len(modern_only)} extra detections")
        for i, hb in enumerate(modern_only, 1):
            sel1 = to_pymol_selection(pdb_id, hb.atom1_id)
            sel2 = to_pymol_selection(pdb_id, hb.atom2_id)
            lines.append(f"# {hb.atom1_id} -- {hb.atom2_id} (dist={hb.distance:.2f})")
            lines.append(f"distance modern_{i}, {sel1}, {sel2}")
            lines.append(f"color blue, modern_{i}")
        lines.append("")

    # Matched H-bonds (GREEN) - optional
    if show_matched and matched:
        lines.append(f"# MATCHED H-bonds (GREEN) - {len(matched)} bonds both detect")
        for i, (dssr_hb, modern_hb) in enumerate(matched, 1):
            sel1 = to_pymol_selection(pdb_id, dssr_hb.atom1_id)
            sel2 = to_pymol_selection(pdb_id, dssr_hb.atom2_id)
            lines.append(f"# {dssr_hb.atom1_id} -- {dssr_hb.atom2_id} (dssr={dssr_hb.distance:.2f}, modern={modern_hb.distance:.2f})")
            lines.append(f"distance match_{i}, {sel1}, {sel2}")
            lines.append(f"color green, match_{i}")
        lines.append("")

    # Create groups
    lines.append("# Create groups for easy toggling")
    if dssr_only:
        lines.append("group dssr_only, dssr_*")
    if modern_only:
        lines.append("group modern_only, modern_*")
    if show_matched and matched:
        lines.append("group matched, match_*")
    lines.append("")

    # Zoom to fit
    lines.append("# Zoom to fit")
    lines.append(f"zoom {pdb_id}")
    lines.append("")

    # Usage tips
    lines.append("# Usage tips:")
    lines.append("#   hide dssr_only     - hide DSSR-only bonds")
    lines.append("#   show dssr_only     - show DSSR-only bonds")
    lines.append("#   disable modern_only - disable modern-only group")
    lines.append("#   enable matched     - enable matched group")

    return "\n".join(lines)


def print_summary(dssr_only: List[HBond], modern_only: List[HBond],
                  matched: List[Tuple[HBond, HBond]], verbose: bool = False):
    """Print comparison summary to console."""
    total_dssr = len(dssr_only) + len(matched)
    total_modern = len(modern_only) + len(matched)

    print(f"\n{'='*60}")
    print(f"DSSR vs Modern H-bond Comparison")
    print(f"{'='*60}")
    print(f"DSSR total:      {total_dssr} (excluding questionable)")
    print(f"Modern total:    {total_modern}")
    print(f"Matched:         {len(matched)} ({len(matched)/total_dssr*100:.1f}% of DSSR)")
    print(f"DSSR only:       {len(dssr_only)} (we miss these)")
    print(f"Modern only:     {len(modern_only)} (extra detections)")

    if verbose:
        if dssr_only:
            print(f"\nDSSR-only (first 10):")
            for hb in dssr_only[:10]:
                print(f"  {hb.atom1_id} -- {hb.atom2_id} ({hb.distance:.2f}Å) [{hb.don_acc_type}]")
            if len(dssr_only) > 10:
                print(f"  ... and {len(dssr_only) - 10} more")

        if modern_only:
            print(f"\nModern-only (first 10):")
            for hb in modern_only[:10]:
                print(f"  {hb.atom1_id} -- {hb.atom2_id} ({hb.distance:.2f}Å)")
            if len(modern_only) > 10:
                print(f"  ... and {len(modern_only) - 10} more")


# Miss reason to color mapping for prototype mode
REASON_COLORS = {
    'DONOR_OVERSATURATED': 'red',
    'ACCEPTOR_OVERSATURATED': 'red',
    'POOR_ALIGNMENT': 'orange',
    'DISTANCE_TOO_FAR': 'yellow',
    'BIFURCATION_ANGLE_SMALL': 'purple',
    'NOT_VALID_DONOR_ACCEPTOR': 'gray',
    'NOT_BASE_PAIRED': 'white',
    'CONFLICT_LOST': 'salmon',
    'UNKNOWN': 'pink',
}


def generate_prototype_pymol_script(pdb_id: str,
                                     analyses,
                                     matched_count: int,
                                     pdb_path: str) -> str:
    """Generate PyMOL script with miss-reason color coding from prototype analysis."""
    lines = []
    lines.append(f"# DSSR Miss Analysis (Prototype): {pdb_id}")
    lines.append(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"#")
    lines.append(f"# Usage: pymol {pdb_path} {pdb_id}_hbonds.pml")
    lines.append(f"#")
    lines.append(f"# Summary:")
    lines.append(f"#   Matched: {matched_count}")
    lines.append(f"#   Missed: {len(analyses)}")
    lines.append(f"#")
    lines.append(f"# Color scheme by miss reason:")
    lines.append(f"#   RED    = Oversaturated (donor or acceptor)")
    lines.append(f"#   ORANGE = Poor alignment score")
    lines.append(f"#   YELLOW = Distance too far")
    lines.append(f"#   PURPLE = Bifurcation angle too small")
    lines.append(f"#   SALMON = Lost to competing bond")
    lines.append(f"#   WHITE  = Not base-paired / not in residue set")
    lines.append(f"#   GRAY   = Not valid donor/acceptor")
    lines.append("")

    # Setup display
    lines.append("# Setup display")
    lines.append("hide all")
    lines.append("show sticks")
    lines.append(f"color white, {pdb_id}")
    lines.append("set dash_gap, 0.2")
    lines.append("set dash_radius, 0.1")
    lines.append("set label_size, 12")
    lines.append("set label_color, white")
    lines.append("")

    # Group by reason
    by_reason = {}
    for a in analyses:
        reason = a.reason.value
        if reason not in by_reason:
            by_reason[reason] = []
        by_reason[reason].append(a)

    # Generate H-bond visualizations grouped by reason
    for reason, items in sorted(by_reason.items(), key=lambda x: -len(x[1])):
        color = REASON_COLORS.get(reason, 'pink')
        lines.append(f"# {reason} ({len(items)}) - {color.upper()}")

        for i, a in enumerate(items, 1):
            hb = a.dssr_hb
            atom1_id = hb.get('atom1_id', '')
            atom2_id = hb.get('atom2_id', '')

            sel1 = to_pymol_selection(pdb_id, atom1_id)
            sel2 = to_pymol_selection(pdb_id, atom2_id)

            obj_name = f"{reason.lower()}_{i}"
            lines.append(f"# {atom1_id} -- {atom2_id} (d={hb.get('distance', 0):.2f}A)")
            lines.append(f"# -> {a.details}")
            lines.append(f"distance {obj_name}, {sel1}, {sel2}")
            lines.append(f"color {color}, {obj_name}")

        lines.append("")

    # Create groups
    lines.append("# Create groups by reason for easy toggling")
    for reason in by_reason.keys():
        lines.append(f"group {reason.lower()}_group, {reason.lower()}_*")
    lines.append("")

    # Zoom
    lines.append(f"zoom {pdb_id}")
    lines.append("")

    # Usage tips
    lines.append("# Usage tips:")
    lines.append("#   disable <reason>_group  - hide that reason")
    lines.append("#   enable <reason>_group   - show that reason")
    lines.append("#   Example: disable not_base_paired_group")

    return "\n".join(lines)


def run_prototype_analysis(pdb_id: str, project_root: Path, verbose: bool = False):
    """Run prototype-based miss analysis and return results."""
    from tools.analyze_dssr_misses import PrototypeMissAnalyzer

    analyzer = PrototypeMissAnalyzer(pdb_id, project_root)
    analyses = analyzer.analyze_all_misses()
    matched = len(analyzer.dssr_hbonds) - len(analyses)

    return analyses, matched, len(analyzer.dssr_hbonds)


def run_prototype_comparison(pdb_id: str, project_root: Path) -> Tuple[List[HBond], List[HBond], List[Tuple[HBond, HBond]]]:
    """Run prototype optimizer and compare with DSSR. Returns (dssr_only, extras, matched)."""
    import sys
    sys.path.insert(0, str(project_root / 'prototypes' / 'hbond_optimizer'))
    from benchmark import load_dssr_json
    from compare_with_dssr import parse_pdb_residues
    from optimizer import HBondOptimizer

    pdb_path = project_root / 'data' / 'pdb' / f'{pdb_id}.pdb'
    dssr_dir = project_root / 'data' / 'json_dssr'

    # Parse residues
    residues = parse_pdb_residues(pdb_path)

    # Create optimizer
    optimizer = HBondOptimizer(max_distance=4.0, min_alignment=0.3)
    for res in residues.values():
        optimizer.add_residue(res)

    # Load DSSR H-bonds
    dssr_hbonds_raw = load_dssr_json(pdb_id, dssr_dir)

    # Convert to HBond format and build lookup
    # Use frozenset of atoms for order-independent matching
    dssr_atom_pairs = {}  # frozenset((atom1, atom2)) -> hb
    dssr_by_pair = {}
    dssr_hbonds = []
    for hb in dssr_hbonds_raw:
        atom1_id = f"{hb.donor_atom}@{hb.res1}"
        atom2_id = f"{hb.acceptor_atom}@{hb.res2}"
        dssr_hbonds.append(HBond(
            atom1_id=atom1_id,
            atom2_id=atom2_id,
            distance=hb.distance,
            source='dssr',
            don_acc_type='standard'
        ))
        # Use frozenset for order-independent matching
        atom_key = frozenset([atom1_id, atom2_id])
        dssr_atom_pairs[atom_key] = hb
        pair_key = (hb.res1, hb.res2)
        if pair_key not in dssr_by_pair:
            dssr_by_pair[pair_key] = []
        dssr_by_pair[pair_key].append(hb)

    # Run optimizer on DSSR pairs
    matched = []
    extras = []
    dssr_matched_keys = set()

    for pair_key in dssr_by_pair.keys():
        res1_id, res2_id = pair_key
        hbonds = optimizer.optimize_pair(res1_id, res2_id)

        for hb in hbonds:
            atom1_id = f"{hb.donor_atom}@{hb.donor_res_id}"
            atom2_id = f"{hb.acceptor_atom}@{hb.acceptor_res_id}"
            atom_key = frozenset([atom1_id, atom2_id])

            if atom_key in dssr_atom_pairs:
                dssr_matched_keys.add(atom_key)
                matched.append((
                    HBond(atom1_id=atom1_id, atom2_id=atom2_id, distance=hb.distance, source='dssr', don_acc_type=''),
                    HBond(atom1_id=atom1_id, atom2_id=atom2_id, distance=hb.distance, source='prototype', don_acc_type='')
                ))
            else:
                extras.append(HBond(
                    atom1_id=atom1_id,
                    atom2_id=atom2_id,
                    distance=hb.distance,
                    source='prototype',
                    don_acc_type=f'align={hb.alignment_score:.2f}'
                ))

    # Find DSSR-only (misses)
    dssr_only = []
    for hb in dssr_hbonds_raw:
        atom1_id = f"{hb.donor_atom}@{hb.res1}"
        atom2_id = f"{hb.acceptor_atom}@{hb.res2}"
        atom_key = frozenset([atom1_id, atom2_id])
        if atom_key not in dssr_matched_keys:
            dssr_only.append(HBond(
                atom1_id=atom1_id,
                atom2_id=atom2_id,
                distance=hb.distance,
                source='dssr',
                don_acc_type='missed'
            ))

    return dssr_only, extras, matched


def main():
    parser = argparse.ArgumentParser(
        description="Compare DSSR and modern H-bond detection, generate PyMOL visualization"
    )
    parser.add_argument("pdb_id", help="PDB ID to analyze (e.g., 1GID)")
    parser.add_argument("-o", "--output", help="Output PyMOL script file (default: <pdb_id>_hbonds.pml)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Show detailed output")
    parser.add_argument("--show-matched", action="store_true", help="Include matched H-bonds in PyMOL script")
    parser.add_argument("--include-questionable", action="store_true", help="Include DSSR questionable H-bonds")
    parser.add_argument("--use-prototype", action="store_true",
                        help="Use Python H-bond optimizer prototype with miss-reason analysis")

    args = parser.parse_args()

    project_root = Path(__file__).parent.parent
    pdb_id = args.pdb_id.upper()
    pdb_path = f"data/pdb/{pdb_id}.pdb"

    if args.use_prototype:
        # Use prototype for comparison (shows misses and extras)
        try:
            print(f"Running prototype comparison for {pdb_id}...")
            dssr_only, modern_only, matched = run_prototype_comparison(pdb_id, project_root)

            total_dssr = len(dssr_only) + len(matched)
            total_modern = len(modern_only) + len(matched)

            print(f"\n{'='*60}")
            print(f"Prototype vs DSSR H-bond Comparison")
            print(f"{'='*60}")
            print(f"DSSR total:      {total_dssr}")
            print(f"Prototype total: {total_modern}")
            print(f"Matched:         {len(matched)} ({len(matched)/total_dssr*100:.1f}% of DSSR)")
            print(f"DSSR only:       {len(dssr_only)} (we miss these)")
            print(f"Prototype only:  {len(modern_only)} (extras)")

            if args.verbose:
                if dssr_only:
                    print(f"\nDSSR-only (misses):")
                    for hb in dssr_only[:10]:
                        print(f"  {hb.atom1_id} -- {hb.atom2_id} ({hb.distance:.2f}Å)")
                    if len(dssr_only) > 10:
                        print(f"  ... and {len(dssr_only) - 10} more")

                if modern_only:
                    print(f"\nPrototype-only (extras):")
                    for hb in modern_only[:10]:
                        print(f"  {hb.atom1_id} -- {hb.atom2_id} ({hb.distance:.2f}Å) [{hb.don_acc_type}]")
                    if len(modern_only) > 10:
                        print(f"  ... and {len(modern_only) - 10} more")

            # Generate script (reuse existing function)
            script = generate_pymol_script(
                pdb_id, dssr_only, modern_only, matched,
                pdb_path, show_matched=args.show_matched
            )

        except Exception as e:
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()
            return 1
    else:
        # Use C++ mode (original behavior)
        try:
            dssr_hbonds = load_dssr_hbonds(
                pdb_id, project_root,
                exclude_questionable=not args.include_questionable
            )
            modern_hbonds = load_modern_hbonds(pdb_id, project_root)
        except FileNotFoundError as e:
            print(f"Error: {e}")
            return 1

        # Match H-bonds
        dssr_only, modern_only, matched = match_hbonds(dssr_hbonds, modern_hbonds)

        # Print summary
        print_summary(dssr_only, modern_only, matched, args.verbose)

        # Generate PyMOL script
        script = generate_pymol_script(
            pdb_id, dssr_only, modern_only, matched,
            pdb_path, show_matched=args.show_matched
        )

    # Write output
    output_file = args.output or f"{pdb_id}_hbonds.pml"
    with open(output_file, 'w') as f:
        f.write(script)

    print(f"\nPyMOL script written to: {output_file}")
    print(f"Usage: pymol {pdb_path} {output_file}")

    return 0


if __name__ == "__main__":
    exit(main())

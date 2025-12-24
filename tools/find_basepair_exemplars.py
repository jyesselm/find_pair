#!/usr/bin/env python3
"""
Find the best exemplar base pairs for each Leontis-Westhof classification type.

Scans all DSSR JSON files, computes ideal geometry statistics per bp-LW type,
scores each pair by deviation from ideal, and extracts the best exemplars as PDB files.

All exemplars are transformed to a common reference frame:
- Origin at [0, 0, 0]
- Identity orientation (x=[1,0,0], y=[0,1,0], z=[0,0,1])

Usage:
    python tools/find_basepair_exemplars.py
    python tools/find_basepair_exemplars.py --min-samples 5 --top-n 3
    python tools/find_basepair_exemplars.py --scan-only  # Just print statistics
"""

import argparse
import json
import math
import os
import re
import statistics
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple


@dataclass
class BasePairFrame:
    """Reference frame for a base pair."""
    origin: Tuple[float, float, float]
    x_axis: Tuple[float, float, float]
    y_axis: Tuple[float, float, float]
    z_axis: Tuple[float, float, float]


@dataclass
class BasePairRecord:
    """A single base pair record from DSSR JSON."""
    pdb_id: str
    bp_key: str           # e.g., "G-C-cWW"
    bp_type: str          # e.g., "G-C"
    lw_class: str         # e.g., "cWW"
    nt1_id: str           # DSSR residue ID (e.g., "A.G1")
    nt2_id: str           # DSSR residue ID (e.g., "A.C72")
    n1n9_dist: float      # N1-N9 distance
    interbase_angle: float  # Angle between bases
    planarity: float      # Planarity deviation
    hbonds_num: int       # Number of H-bonds
    hbonds_desc: str = "" # H-bond description
    frame: Optional[BasePairFrame] = None  # Reference frame for alignment


@dataclass
class IdealGeometry:
    """Ideal geometry statistics for a bp-LW type."""
    bp_key: str
    ideal_n1n9: float
    stdev_n1n9: float
    ideal_angle: float
    stdev_angle: float
    ideal_planarity: float
    stdev_planarity: float
    expected_hbonds: int
    sample_size: int


def transform_point(point: Tuple[float, float, float],
                    frame: BasePairFrame) -> Tuple[float, float, float]:
    """Transform a point from world coordinates to the base pair frame.

    This aligns the base pair so its reference frame becomes identity at origin.
    The transformation is: p' = R^T * (p - origin)
    where R is the rotation matrix [x_axis, y_axis, z_axis] as columns.
    """
    # Translate to origin
    px = point[0] - frame.origin[0]
    py = point[1] - frame.origin[1]
    pz = point[2] - frame.origin[2]

    # Rotate by inverse (transpose) of frame rotation matrix
    # R^T * p = [x_axis · p, y_axis · p, z_axis · p]
    new_x = frame.x_axis[0] * px + frame.x_axis[1] * py + frame.x_axis[2] * pz
    new_y = frame.y_axis[0] * px + frame.y_axis[1] * py + frame.y_axis[2] * pz
    new_z = frame.z_axis[0] * px + frame.z_axis[1] * py + frame.z_axis[2] * pz

    return (new_x, new_y, new_z)


def parse_pdb_coords(line: str) -> Tuple[float, float, float]:
    """Parse x, y, z coordinates from a PDB ATOM/HETATM line."""
    x = float(line[30:38].strip())
    y = float(line[38:46].strip())
    z = float(line[46:54].strip())
    return (x, y, z)


def format_pdb_coords(x: float, y: float, z: float) -> str:
    """Format coordinates for PDB file (8.3f format for each)."""
    return f"{x:8.3f}{y:8.3f}{z:8.3f}"


def transform_pdb_line(line: str, frame: BasePairFrame) -> str:
    """Transform coordinates in a PDB ATOM/HETATM line to the base pair frame."""
    coords = parse_pdb_coords(line)
    new_coords = transform_point(coords, frame)

    # Replace coordinates in line (columns 31-54, 0-indexed: 30-53)
    new_line = line[:30] + format_pdb_coords(*new_coords) + line[54:]
    return new_line


def parse_dssr_nt_id(nt_id: str) -> Tuple[str, str, int, str]:
    """Parse DSSR nucleotide ID into (chain, resname, resnum, icode).

    Examples:
        "A.G1" -> ("A", "G", 1, "")
        "B.C25^A" -> ("B", "C", 25, "A")
        "A.5MC49" -> ("A", "5MC", 49, "")
    """
    if '.' not in nt_id:
        return "", "", 0, ""

    chain, rest = nt_id.split('.', 1)

    # Handle insertion code (^A format)
    insertion = ""
    if '^' in rest:
        rest, insertion = rest.split('^', 1)

    # Extract resname and resnum - find last contiguous digit sequence
    resname = ""
    resnum = 0

    # Find where digits start from the end
    digit_start = len(rest)
    for i in range(len(rest) - 1, -1, -1):
        if rest[i].isdigit():
            digit_start = i
        else:
            break

    if digit_start < len(rest):
        resname = rest[:digit_start]
        try:
            # Handle negative residue numbers
            if resname.endswith('-'):
                resname = resname[:-1]
                resnum = -int(rest[digit_start:])
            else:
                resnum = int(rest[digit_start:])
        except ValueError:
            resnum = 0
    else:
        resname = rest

    return chain, resname, resnum, insertion


def load_pairs_from_json(json_path: Path) -> List[BasePairRecord]:
    """Load base pair records from a single DSSR JSON file."""
    try:
        with open(json_path) as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError):
        return []

    pdb_id = json_path.stem
    pairs = []

    for pair in data.get('pairs', []):
        lw = pair.get('LW', '')

        # Skip pairs without valid LW classification
        if not lw or lw == '--' or lw == '.':
            continue

        bp_type = pair.get('bp', '')
        if not bp_type:
            continue

        # Get geometry metrics (handle None values)
        n1n9_dist = pair.get('N1N9_dist')
        interbase_angle = pair.get('interBase_angle')
        planarity = pair.get('planarity')
        hbonds_num = pair.get('hbonds_num', 0)
        hbonds_desc = pair.get('hbonds_desc', '')

        # Skip pairs with missing/invalid geometry
        if n1n9_dist is None or interbase_angle is None or planarity is None:
            continue
        if n1n9_dist <= 0 or interbase_angle <= 0:
            continue

        bp_key = f"{bp_type}-{lw}"

        # Parse frame data for alignment
        frame_data = pair.get('frame', {})
        frame = None
        if frame_data:
            origin = frame_data.get('origin')
            x_axis = frame_data.get('x_axis')
            y_axis = frame_data.get('y_axis')
            z_axis = frame_data.get('z_axis')
            if origin and x_axis and y_axis and z_axis:
                frame = BasePairFrame(
                    origin=tuple(origin),
                    x_axis=tuple(x_axis),
                    y_axis=tuple(y_axis),
                    z_axis=tuple(z_axis),
                )

        pairs.append(BasePairRecord(
            pdb_id=pdb_id,
            bp_key=bp_key,
            bp_type=bp_type,
            lw_class=lw,
            nt1_id=pair.get('nt1', ''),
            nt2_id=pair.get('nt2', ''),
            n1n9_dist=n1n9_dist,
            interbase_angle=interbase_angle,
            planarity=planarity,
            hbonds_num=hbonds_num,
            hbonds_desc=hbonds_desc,
            frame=frame,
        ))

    return pairs


def load_all_pairs(json_dir: Path, workers: int = None) -> List[BasePairRecord]:
    """Load all base pair records from all DSSR JSON files."""
    json_files = list(json_dir.glob("*.json"))

    if not json_files:
        raise FileNotFoundError(f"No JSON files found in {json_dir}")

    all_pairs = []

    if workers is None:
        workers = os.cpu_count() or 4

    print(f"Loading pairs from {len(json_files)} DSSR JSON files using {workers} workers...")

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(load_pairs_from_json, f): f for f in json_files}

        for i, future in enumerate(as_completed(futures), 1):
            try:
                pairs = future.result()
                all_pairs.extend(pairs)
            except Exception as e:
                print(f"  Error processing {futures[future]}: {e}")

            if i % 500 == 0 or i == len(json_files):
                print(f"  Processed {i}/{len(json_files)} files, {len(all_pairs)} pairs so far")

    return all_pairs


def compute_ideal_geometries(pairs: List[BasePairRecord],
                              min_samples: int = 10) -> Dict[str, IdealGeometry]:
    """Compute ideal geometry statistics for each bp-LW type."""
    # Group pairs by bp_key
    grouped = defaultdict(list)
    for pair in pairs:
        grouped[pair.bp_key].append(pair)

    ideals = {}

    for bp_key, group in grouped.items():
        if len(group) < min_samples:
            continue

        n1n9_vals = [p.n1n9_dist for p in group]
        angle_vals = [p.interbase_angle for p in group]
        planarity_vals = [p.planarity for p in group]
        hbonds_vals = [p.hbonds_num for p in group]

        # Calculate mean and stdev
        ideal_n1n9 = statistics.mean(n1n9_vals)
        stdev_n1n9 = statistics.stdev(n1n9_vals) if len(n1n9_vals) > 1 else 0.1

        ideal_angle = statistics.mean(angle_vals)
        stdev_angle = statistics.stdev(angle_vals) if len(angle_vals) > 1 else 1.0

        ideal_planarity = statistics.mean(planarity_vals)
        stdev_planarity = statistics.stdev(planarity_vals) if len(planarity_vals) > 1 else 0.1

        # Mode for hbonds
        hbond_counts = Counter(hbonds_vals)
        expected_hbonds = hbond_counts.most_common(1)[0][0]

        ideals[bp_key] = IdealGeometry(
            bp_key=bp_key,
            ideal_n1n9=ideal_n1n9,
            stdev_n1n9=stdev_n1n9,
            ideal_angle=ideal_angle,
            stdev_angle=stdev_angle,
            ideal_planarity=ideal_planarity,
            stdev_planarity=stdev_planarity,
            expected_hbonds=expected_hbonds,
            sample_size=len(group),
        )

    return ideals


def score_basepair(pair: BasePairRecord, ideal: IdealGeometry) -> float:
    """Score a base pair by deviation from ideal geometry (lower = better)."""
    # Normalized deviations
    n1n9_dev = abs(pair.n1n9_dist - ideal.ideal_n1n9) / max(ideal.stdev_n1n9, 0.01)
    angle_dev = abs(pair.interbase_angle - ideal.ideal_angle) / max(ideal.stdev_angle, 0.1)
    planarity_dev = pair.planarity / max(ideal.stdev_planarity, 0.01)
    hbond_penalty = max(0, ideal.expected_hbonds - pair.hbonds_num)

    # Weighted sum (weights can be tuned)
    score = (1.0 * n1n9_dev +
             1.0 * angle_dev +
             1.5 * planarity_dev +
             2.0 * hbond_penalty)

    return score


def find_best_exemplars(pairs: List[BasePairRecord],
                        ideals: Dict[str, IdealGeometry],
                        top_n: int = 1) -> Dict[str, List[Tuple[BasePairRecord, float]]]:
    """Find the best exemplars for each bp-LW type."""
    # Group pairs by bp_key
    grouped = defaultdict(list)
    for pair in pairs:
        if pair.bp_key in ideals:
            ideal = ideals[pair.bp_key]
            score = score_basepair(pair, ideal)
            grouped[pair.bp_key].append((pair, score))

    # Sort by score and take top N
    best = {}
    for bp_key, scored_pairs in grouped.items():
        scored_pairs.sort(key=lambda x: x[1])
        best[bp_key] = scored_pairs[:top_n]

    return best


def extract_residue_atoms(pdb_path: Path, chain: str, resnum: int,
                          icode: str = "") -> List[str]:
    """Extract ATOM/HETATM lines for a specific residue from PDB file."""
    atoms = []

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue

            # PDB format: columns are fixed width
            # Chain ID: column 22 (0-indexed: 21)
            # Residue number: columns 23-26 (0-indexed: 22-25)
            # Insertion code: column 27 (0-indexed: 26)

            pdb_chain = line[21] if len(line) > 21 else ''

            try:
                pdb_resnum = int(line[22:26].strip())
            except ValueError:
                continue

            pdb_icode = line[26].strip() if len(line) > 26 else ''

            # Check for alternate conformations - take first (blank or 'A')
            altloc = line[16] if len(line) > 16 else ' '
            if altloc not in (' ', 'A'):
                continue

            if pdb_chain == chain and pdb_resnum == resnum and pdb_icode == icode:
                atoms.append(line)

    return atoms


def write_exemplar_pdb(output_path: Path, atoms1: List[str], atoms2: List[str],
                       pair: BasePairRecord) -> None:
    """Write exemplar base pair to PDB file with header.

    Coordinates are transformed to align the base pair reference frame
    to the identity frame at origin.
    """
    with open(output_path, 'w') as f:
        f.write(f"REMARK   1 Base pair exemplar: {pair.bp_key}\n")
        f.write(f"REMARK   2 Source PDB: {pair.pdb_id}\n")
        f.write(f"REMARK   3 Residue 1: {pair.nt1_id}\n")
        f.write(f"REMARK   4 Residue 2: {pair.nt2_id}\n")
        f.write(f"REMARK   5 LW classification: {pair.lw_class}\n")
        f.write(f"REMARK   6 N1N9 distance: {pair.n1n9_dist:.3f}\n")
        f.write(f"REMARK   7 Interbase angle: {pair.interbase_angle:.3f}\n")
        f.write(f"REMARK   8 Planarity: {pair.planarity:.3f}\n")
        f.write(f"REMARK   9 H-bonds: {pair.hbonds_num} ({pair.hbonds_desc})\n")
        f.write(f"REMARK  10 Aligned to origin with identity frame\n")

        # Write atoms with renumbered atom serial and transformed coordinates
        atom_num = 1
        for line in atoms1 + atoms2:
            # Transform coordinates if we have frame data
            if pair.frame:
                line = transform_pdb_line(line, pair.frame)

            # Renumber atom serial (columns 7-11, 0-indexed: 6-10)
            new_line = line[:6] + f"{atom_num:5d}" + line[11:]
            f.write(new_line)
            atom_num += 1

        f.write("END\n")


def sanitize_filename(bp_key: str) -> str:
    """Sanitize bp_key for use as filename."""
    # Replace problematic characters
    return bp_key.replace('+', 'plus').replace('/', '-').replace('*', 'star')


def main():
    parser = argparse.ArgumentParser(
        description="Find best exemplar base pairs for each LW classification type"
    )
    parser.add_argument("--json-dir", default="data/json_dssr",
                        help="Directory containing DSSR JSON files")
    parser.add_argument("--pdb-dir", default="data/pdb",
                        help="Directory containing PDB files")
    parser.add_argument("--output-dir", default="basepair-exemplars",
                        help="Output directory for exemplar PDBs")
    parser.add_argument("--top-n", type=int, default=1,
                        help="Number of top exemplars to extract per type")
    parser.add_argument("--min-samples", type=int, default=10,
                        help="Minimum samples required to compute ideal geometry")
    parser.add_argument("--workers", "-w", type=int, default=None,
                        help="Number of parallel workers (default: CPU count)")
    parser.add_argument("--scan-only", action="store_true",
                        help="Only scan and print statistics, don't extract PDBs")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    project_root = Path(__file__).parent.parent
    json_dir = project_root / args.json_dir
    pdb_dir = project_root / args.pdb_dir
    output_dir = project_root / args.output_dir

    # Load all pairs
    all_pairs = load_all_pairs(json_dir, args.workers)
    print(f"\nLoaded {len(all_pairs)} base pairs with valid LW classification")

    # Compute ideal geometries
    print(f"\nComputing ideal geometries (min samples: {args.min_samples})...")
    ideals = compute_ideal_geometries(all_pairs, args.min_samples)
    print(f"Found {len(ideals)} bp-LW types with sufficient samples")

    # Print statistics by LW class
    lw_counts = defaultdict(int)
    lw_types = defaultdict(set)
    for bp_key, ideal in ideals.items():
        parts = bp_key.rsplit('-', 1)
        if len(parts) == 2:
            lw = parts[1]
            lw_counts[lw] += ideal.sample_size
            lw_types[lw].add(bp_key)

    print("\nStatistics by LW classification:")
    print(f"{'LW':<8} {'Types':<8} {'Total Pairs':<12}")
    print("-" * 30)
    for lw in sorted(lw_counts.keys()):
        print(f"{lw:<8} {len(lw_types[lw]):<8} {lw_counts[lw]:<12}")

    if args.verbose:
        print("\nDetailed ideal geometries:")
        print(f"{'Type':<15} {'N':<8} {'N1N9':<12} {'Angle':<12} {'Planarity':<12} {'HBonds':<8}")
        print("-" * 70)
        for bp_key in sorted(ideals.keys()):
            ideal = ideals[bp_key]
            print(f"{bp_key:<15} {ideal.sample_size:<8} "
                  f"{ideal.ideal_n1n9:.2f}±{ideal.stdev_n1n9:.2f}  "
                  f"{ideal.ideal_angle:.1f}±{ideal.stdev_angle:.1f}  "
                  f"{ideal.ideal_planarity:.2f}±{ideal.stdev_planarity:.2f}  "
                  f"{ideal.expected_hbonds:<8}")

    if args.scan_only:
        print("\n--scan-only specified, skipping PDB extraction")
        return 0

    # Find best exemplars
    print(f"\nFinding top {args.top_n} exemplars per type...")
    best_exemplars = find_best_exemplars(all_pairs, ideals, args.top_n)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract and write exemplar PDBs
    print(f"\nExtracting exemplar PDBs to {output_dir}...")

    summary = []
    extracted = 0
    skipped = 0

    for bp_key in sorted(best_exemplars.keys()):
        for rank, (pair, score) in enumerate(best_exemplars[bp_key], 1):
            pdb_path = pdb_dir / f"{pair.pdb_id}.pdb"

            if not pdb_path.exists():
                if args.verbose:
                    print(f"  Skipping {bp_key}: PDB {pair.pdb_id} not found")
                skipped += 1
                continue

            # Parse residue IDs
            chain1, resname1, resnum1, icode1 = parse_dssr_nt_id(pair.nt1_id)
            chain2, resname2, resnum2, icode2 = parse_dssr_nt_id(pair.nt2_id)

            # Extract atoms
            atoms1 = extract_residue_atoms(pdb_path, chain1, resnum1, icode1)
            atoms2 = extract_residue_atoms(pdb_path, chain2, resnum2, icode2)

            if not atoms1 or not atoms2:
                if args.verbose:
                    print(f"  Skipping {bp_key}: Could not extract atoms from {pair.pdb_id}")
                skipped += 1
                continue

            # Skip if no frame data for alignment
            if not pair.frame:
                if args.verbose:
                    print(f"  Skipping {bp_key}: No frame data for alignment")
                skipped += 1
                continue

            # Write exemplar PDB
            safe_name = sanitize_filename(bp_key)
            if args.top_n > 1:
                output_file = output_dir / f"{safe_name}_rank{rank}.pdb"
            else:
                output_file = output_dir / f"{safe_name}.pdb"

            write_exemplar_pdb(output_file, atoms1, atoms2, pair)
            extracted += 1

            # Add to summary
            summary.append({
                'bp_key': bp_key,
                'rank': rank,
                'score': score,
                'pdb_id': pair.pdb_id,
                'nt1_id': pair.nt1_id,
                'nt2_id': pair.nt2_id,
                'n1n9_dist': pair.n1n9_dist,
                'interbase_angle': pair.interbase_angle,
                'planarity': pair.planarity,
                'hbonds_num': pair.hbonds_num,
                'hbonds_desc': pair.hbonds_desc,
                'ideal': {
                    'n1n9': ideals[bp_key].ideal_n1n9,
                    'angle': ideals[bp_key].ideal_angle,
                    'planarity': ideals[bp_key].ideal_planarity,
                    'hbonds': ideals[bp_key].expected_hbonds,
                    'sample_size': ideals[bp_key].sample_size,
                },
                'output_file': str(output_file.name),
            })

    # Write summary JSON
    summary_file = output_dir / "summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\nComplete!")
    print(f"  Extracted: {extracted} exemplar PDBs")
    print(f"  Skipped: {skipped} (missing PDB or atoms)")
    print(f"  Summary: {summary_file}")

    return 0


if __name__ == "__main__":
    exit(main())

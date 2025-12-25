#!/usr/bin/env python3
"""
Multithreaded benchmark for H-bond optimizer vs DSSR.

Runs comparison on large PDB sets and documents all differences.
"""

import json
import sys
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from prototypes.hbond_optimizer import (
    parse_pdb_residues, run_dssr, parse_dssr_output,
    HBondOptimizer, Residue, DSSRHBond
)
from prototypes.hbond_optimizer.compare_with_dssr import normalize_res_id
from prototypes.hbond_optimizer.geometry import (
    DONOR_CAPACITY, ACCEPTOR_CAPACITY, score_hbond_alignment, normalize,
    get_donor_capacity, get_acceptor_capacity
)


@dataclass
class MissedHBond:
    """A missed H-bond with diagnostic info."""
    pdb_id: str
    res1: str           # Normalized ID for matching
    res2: str           # Normalized ID for matching
    res1_original: str  # Original DSSR ID (e.g., "A.7MG1")
    res2_original: str  # Original DSSR ID (e.g., "B.5MC2")
    donor_atom: str
    acceptor_atom: str
    distance: float
    reason: str
    details: str = ""


@dataclass
class PDBResult:
    """Results for a single PDB."""
    pdb_id: str
    dssr_count: int
    matched_count: int
    optimizer_count: int
    missed: List[MissedHBond] = field(default_factory=list)
    extra: List[Tuple[str, str, str, str, float]] = field(default_factory=list)
    error: Optional[str] = None


def extract_residue_code(dssr_res_id: str) -> str:
    """Extract residue code from DSSR ID (e.g., 'A.7MG1' -> '7MG', 'A.G1' -> 'G')."""
    import re
    parts = dssr_res_id.split('.')
    if len(parts) != 2:
        return ""
    res = parts[1]
    # Find where the number starts
    match = re.match(r'^([A-Za-z0-9]+?)(-?\d+)', res)
    if match:
        return match.group(1)
    return res


def analyze_miss(pdb_id: str, hb: DSSRHBond, residues: Dict[str, Residue],
                 optimizer: HBondOptimizer) -> MissedHBond:
    """Analyze why an H-bond was missed."""
    res1 = residues.get(hb.res1)
    res2 = residues.get(hb.res2)

    # Use original IDs for reporting
    res1_orig = hb.res1_original if hb.res1_original else hb.res1
    res2_orig = hb.res2_original if hb.res2_original else hb.res2

    # Check if residues exist
    if not res1:
        return MissedHBond(
            pdb_id=pdb_id, res1=hb.res1, res2=hb.res2,
            res1_original=res1_orig, res2_original=res2_orig,
            donor_atom=hb.donor_atom, acceptor_atom=hb.acceptor_atom,
            distance=hb.distance, reason="MISSING_RESIDUE",
            details=f"Residue {res1_orig} not found (normalized: {hb.res1})"
        )
    if not res2:
        return MissedHBond(
            pdb_id=pdb_id, res1=hb.res1, res2=hb.res2,
            res1_original=res1_orig, res2_original=res2_orig,
            donor_atom=hb.donor_atom, acceptor_atom=hb.acceptor_atom,
            distance=hb.distance, reason="MISSING_RESIDUE",
            details=f"Residue {res2_orig} not found (normalized: {hb.res2})"
        )

    # Check if atoms exist
    if hb.donor_atom not in res1.atoms and hb.donor_atom not in res2.atoms:
        return MissedHBond(
            pdb_id=pdb_id, res1=hb.res1, res2=hb.res2,
            res1_original=res1_orig, res2_original=res2_orig,
            donor_atom=hb.donor_atom, acceptor_atom=hb.acceptor_atom,
            distance=hb.distance, reason="MISSING_ATOM",
            details=f"Atom {hb.donor_atom} not found in {res1_orig}"
        )
    if hb.acceptor_atom not in res1.atoms and hb.acceptor_atom not in res2.atoms:
        return MissedHBond(
            pdb_id=pdb_id, res1=hb.res1, res2=hb.res2,
            res1_original=res1_orig, res2_original=res2_orig,
            donor_atom=hb.donor_atom, acceptor_atom=hb.acceptor_atom,
            distance=hb.distance, reason="MISSING_ATOM",
            details=f"Atom {hb.acceptor_atom} not found in {res2_orig}"
        )

    # Determine donor/acceptor assignment
    # DSSR format: atom1@res1 atom2@res2, but doesn't specify which is donor
    # Try both orientations using CIF-aware capacity functions
    donor_res, acceptor_res = None, None
    donor_atom, acceptor_atom = hb.donor_atom, hb.acceptor_atom

    # Get residue codes for CIF lookup (e.g., "1MA", "5MC")
    res1_code = res1.residue_code if res1.residue_code else res1.base_type.upper()
    res2_code = res2.residue_code if res2.residue_code else res2.base_type.upper()

    # Check if res1.donor_atom -> res2.acceptor_atom
    if get_donor_capacity(res1_code, hb.donor_atom) > 0 and get_acceptor_capacity(res2_code, hb.acceptor_atom) > 0:
        donor_res, acceptor_res = res1, res2
        donor_atom, acceptor_atom = hb.donor_atom, hb.acceptor_atom
    # Try res2.donor_atom -> res1.acceptor_atom
    elif get_donor_capacity(res2_code, hb.donor_atom) > 0 and get_acceptor_capacity(res1_code, hb.acceptor_atom) > 0:
        donor_res, acceptor_res = res2, res1
        donor_atom, acceptor_atom = hb.donor_atom, hb.acceptor_atom
    # Try swapped atoms: res1.acceptor_atom -> res2.donor_atom
    elif get_donor_capacity(res1_code, hb.acceptor_atom) > 0 and get_acceptor_capacity(res2_code, hb.donor_atom) > 0:
        donor_res, acceptor_res = res1, res2
        donor_atom, acceptor_atom = hb.acceptor_atom, hb.donor_atom
    # Try swapped: res2.acceptor_atom -> res1.donor_atom
    elif get_donor_capacity(res2_code, hb.acceptor_atom) > 0 and get_acceptor_capacity(res1_code, hb.donor_atom) > 0:
        donor_res, acceptor_res = res2, res1
        donor_atom, acceptor_atom = hb.acceptor_atom, hb.donor_atom

    if not donor_res or not acceptor_res:
        # Not a valid donor/acceptor pair in our tables
        # Include original residue codes for debugging
        res1_code = extract_residue_code(res1_orig)
        res2_code = extract_residue_code(res2_orig)
        return MissedHBond(
            pdb_id=pdb_id, res1=hb.res1, res2=hb.res2,
            res1_original=res1_orig, res2_original=res2_orig,
            donor_atom=hb.donor_atom, acceptor_atom=hb.acceptor_atom,
            distance=hb.distance, reason="NOT_DONOR_ACCEPTOR",
            details=f"{res1_code}.{hb.donor_atom}/{res2_code}.{hb.acceptor_atom} not in tables (parent: {res1.base_type}/{res2.base_type})"
        )

    # Check distance
    if hb.distance > optimizer.max_distance:
        return MissedHBond(
            pdb_id=pdb_id, res1=hb.res1, res2=hb.res2,
            res1_original=res1_orig, res2_original=res2_orig,
            donor_atom=hb.donor_atom, acceptor_atom=hb.acceptor_atom,
            distance=hb.distance, reason="DISTANCE_TOO_FAR",
            details=f"{hb.distance:.2f} > {optimizer.max_distance}"
        )

    # Check slots
    h_slots = donor_res.get_h_slots(donor_atom)
    lp_slots = acceptor_res.get_lp_slots(acceptor_atom)

    if not h_slots:
        donor_code = extract_residue_code(res1_orig if donor_res.res_id == hb.res1 else res2_orig)
        return MissedHBond(
            pdb_id=pdb_id, res1=hb.res1, res2=hb.res2,
            res1_original=res1_orig, res2_original=res2_orig,
            donor_atom=hb.donor_atom, acceptor_atom=hb.acceptor_atom,
            distance=hb.distance, reason="NO_H_SLOTS",
            details=f"{donor_code}.{donor_atom} has 0 H slots (parent: {donor_res.base_type})"
        )

    if not lp_slots:
        acc_code = extract_residue_code(res2_orig if acceptor_res.res_id == hb.res2 else res1_orig)
        return MissedHBond(
            pdb_id=pdb_id, res1=hb.res1, res2=hb.res2,
            res1_original=res1_orig, res2_original=res2_orig,
            donor_atom=hb.donor_atom, acceptor_atom=hb.acceptor_atom,
            distance=hb.distance, reason="NO_LP_SLOTS",
            details=f"{acc_code}.{acceptor_atom} has 0 LP slots (parent: {acceptor_res.base_type})"
        )

    # Check alignment - use ignore_used=True to get "theoretical" alignment
    # This tells us whether geometry would allow the bond if slots were free
    if donor_atom in donor_res.atoms and acceptor_atom in acceptor_res.atoms:
        donor_pos = donor_res.atoms[donor_atom]
        acceptor_pos = acceptor_res.atoms[acceptor_atom]

        try:
            # Compute theoretical alignment (ignoring slot usage)
            h_idx, lp_idx, theoretical_score = score_hbond_alignment(
                donor_pos, acceptor_pos, h_slots, lp_slots, ignore_used=True
            )

            # Check if any slots are actually used (i.e., taken by other bonds)
            any_h_used = any(s.used for s in h_slots)
            any_lp_used = any(s.used for s in lp_slots)

            if theoretical_score < optimizer.min_alignment:
                # Geometry itself is bad
                return MissedHBond(
                    pdb_id=pdb_id, res1=hb.res1, res2=hb.res2,
                    res1_original=res1_orig, res2_original=res2_orig,
                    donor_atom=hb.donor_atom, acceptor_atom=hb.acceptor_atom,
                    distance=hb.distance, reason="POOR_ALIGNMENT",
                    details=f"alignment={theoretical_score:.3f} < {optimizer.min_alignment}"
                )
            elif any_h_used or any_lp_used:
                # Geometry is OK but slots were taken by other bonds
                return MissedHBond(
                    pdb_id=pdb_id, res1=hb.res1, res2=hb.res2,
                    res1_original=res1_orig, res2_original=res2_orig,
                    donor_atom=hb.donor_atom, acceptor_atom=hb.acceptor_atom,
                    distance=hb.distance, reason="SLOT_CONFLICT",
                    details=f"alignment={theoretical_score:.3f} OK but slots taken (H:{any_h_used}, LP:{any_lp_used})"
                )
        except Exception as e:
            return MissedHBond(
                pdb_id=pdb_id, res1=hb.res1, res2=hb.res2,
                res1_original=res1_orig, res2_original=res2_orig,
                donor_atom=hb.donor_atom, acceptor_atom=hb.acceptor_atom,
                distance=hb.distance, reason="ALIGNMENT_ERROR",
                details=str(e)
            )

    # If we got here, it might be a slot saturation issue or conflict
    return MissedHBond(
        pdb_id=pdb_id, res1=hb.res1, res2=hb.res2,
        res1_original=res1_orig, res2_original=res2_orig,
        donor_atom=hb.donor_atom, acceptor_atom=hb.acceptor_atom,
        distance=hb.distance, reason="SLOT_CONFLICT",
        details="Slots saturated or lost in greedy selection"
    )


def load_dssr_json(pdb_id: str, dssr_json_dir: Path) -> List[DSSRHBond]:
    """Load H-bonds from pre-computed DSSR JSON file."""
    import re

    json_path = dssr_json_dir / f"{pdb_id}.json"
    if not json_path.exists():
        return []

    with open(json_path) as f:
        data = json.load(f)

    hbonds = []

    # Extract H-bonds from pairs
    for pair in data.get('pairs', []):
        hbonds_desc = pair.get('hbonds_desc', '')
        nt1 = pair.get('nt1', '')
        nt2 = pair.get('nt2', '')

        if not hbonds_desc or not nt1 or not nt2:
            continue

        # Parse hbonds_desc: "O2(carbonyl)-N2(amino)[2.77],N3-N1(imino)[2.93]"
        for hb_str in hbonds_desc.split(','):
            # Pattern: atom1(type)-atom2(type)[distance] or atom1-atom2[distance]
            match = re.match(r'([A-Za-z0-9\']+)(?:\([^)]*\))?-([A-Za-z0-9\']+)(?:\([^)]*\))?\[([0-9.]+)\]', hb_str)
            if match:
                atom1, atom2, dist = match.groups()
                # Normalize residue IDs to match PDB parsing format
                nt1_norm = normalize_res_id(nt1)
                nt2_norm = normalize_res_id(nt2)
                # The first atom is typically on nt1, second on nt2
                hbonds.append(DSSRHBond(
                    res1=nt1_norm,
                    res2=nt2_norm,
                    donor_atom=atom1,
                    acceptor_atom=atom2,
                    distance=float(dist),
                    res1_original=nt1,
                    res2_original=nt2,
                ))

    return hbonds


def process_pdb(args: Tuple[str, Path, float, float, float, Path]) -> PDBResult:
    """Process a single PDB file. Designed for multiprocessing."""
    pdb_id, pdb_path, max_distance, min_alignment, short_dist_thresh, dssr_json_dir = args

    try:
        residues = parse_pdb_residues(pdb_path)

        # Load pre-computed DSSR JSON
        dssr_hbonds = load_dssr_json(pdb_id, dssr_json_dir)

        if not dssr_hbonds:
            return PDBResult(pdb_id=pdb_id, dssr_count=0, matched_count=0, optimizer_count=0)

        optimizer = HBondOptimizer(max_distance=max_distance, min_alignment=min_alignment,
                                   short_distance_threshold=short_dist_thresh)
        for res in residues.values():
            optimizer.add_residue(res)

        # Group DSSR hbonds by residue pair
        dssr_by_pair = defaultdict(list)
        for hb in dssr_hbonds:
            pair_key = tuple(sorted([hb.res1, hb.res2]))
            dssr_by_pair[pair_key].append(hb)

        total_dssr = 0
        total_match = 0
        total_opt = 0
        missed_list = []
        extra_list = []

        for pair_key, dssr_hbs in dssr_by_pair.items():
            res1_id, res2_id = pair_key
            if res1_id not in residues or res2_id not in residues:
                # Count as missed due to missing residue
                for hb in dssr_hbs:
                    res1_orig = hb.res1_original if hb.res1_original else hb.res1
                    res2_orig = hb.res2_original if hb.res2_original else hb.res2
                    missing_id = res1_orig if res1_id not in residues else res2_orig
                    missed_list.append(MissedHBond(
                        pdb_id=pdb_id, res1=hb.res1, res2=hb.res2,
                        res1_original=res1_orig, res2_original=res2_orig,
                        donor_atom=hb.donor_atom, acceptor_atom=hb.acceptor_atom,
                        distance=hb.distance, reason="MISSING_RESIDUE",
                        details=f"Residue {missing_id} not parsed"
                    ))
                total_dssr += len(dssr_hbs)
                continue

            opt_hbonds = optimizer.optimize_pair(res1_id, res2_id)

            dssr_set = {tuple(sorted([hb.donor_atom, hb.acceptor_atom])): hb for hb in dssr_hbs}
            opt_set = {tuple(sorted([hb.donor_atom, hb.acceptor_atom])): hb for hb in opt_hbonds}

            dssr_atoms = set(dssr_set.keys())
            opt_atoms = set(opt_set.keys())

            matches = dssr_atoms & opt_atoms
            dssr_only = dssr_atoms - opt_atoms
            opt_only = opt_atoms - dssr_atoms

            total_dssr += len(dssr_atoms)
            total_match += len(matches)
            total_opt += len(opt_atoms)

            # Analyze missed H-bonds
            for atoms in dssr_only:
                hb = dssr_set[atoms]
                miss = analyze_miss(pdb_id, hb, residues, optimizer)
                missed_list.append(miss)

            # Record extra H-bonds
            for atoms in opt_only:
                hb = opt_set[atoms]
                extra_list.append((res1_id, res2_id, hb.donor_atom, hb.acceptor_atom, hb.distance))

        return PDBResult(
            pdb_id=pdb_id,
            dssr_count=total_dssr,
            matched_count=total_match,
            optimizer_count=total_opt,
            missed=missed_list,
            extra=extra_list
        )

    except Exception as e:
        return PDBResult(pdb_id=pdb_id, dssr_count=0, matched_count=0, optimizer_count=0, error=str(e))


def run_benchmark(pdb_ids: List[str], pdb_dir: Path, dssr_json_dir: Path,
                  max_workers: int = None,
                  max_distance: float = 4.0, min_alignment: float = 0.3,
                  short_distance_threshold: float = 3.2) -> Dict:
    """Run benchmark on multiple PDBs with multiprocessing."""
    if max_workers is None:
        max_workers = mp.cpu_count()

    # Prepare arguments
    args_list = []
    for pdb_id in pdb_ids:
        pdb_path = pdb_dir / f"{pdb_id}.pdb"
        if pdb_path.exists():
            args_list.append((pdb_id, pdb_path, max_distance, min_alignment, short_distance_threshold, dssr_json_dir))

    print(f"Running benchmark on {len(args_list)} PDBs with {max_workers} workers...")

    results = []
    completed = 0

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_pdb, args): args[0] for args in args_list}

        for future in as_completed(futures):
            pdb_id = futures[future]
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                results.append(PDBResult(pdb_id=pdb_id, dssr_count=0, matched_count=0,
                                        optimizer_count=0, error=str(e)))

            completed += 1
            if completed % 100 == 0:
                print(f"  Completed {completed}/{len(args_list)} PDBs...")

    return analyze_results(results)


def analyze_results(results: List[PDBResult]) -> Dict:
    """Analyze benchmark results and categorize misses."""
    total_dssr = 0
    total_match = 0
    total_opt = 0
    total_errors = 0

    # Categorize misses
    misses_by_reason = defaultdict(list)
    misses_by_atom_pair = defaultdict(list)
    misses_by_base_type = defaultdict(list)

    all_missed = []
    all_extra = []

    for result in results:
        if result.error:
            total_errors += 1
            continue

        total_dssr += result.dssr_count
        total_match += result.matched_count
        total_opt += result.optimizer_count

        for miss in result.missed:
            all_missed.append(miss)
            misses_by_reason[miss.reason].append(miss)

            atom_pair = tuple(sorted([miss.donor_atom, miss.acceptor_atom]))
            misses_by_atom_pair[atom_pair].append(miss)

        all_extra.extend(result.extra)

    recall = total_match / total_dssr * 100 if total_dssr > 0 else 0
    precision = total_match / total_opt * 100 if total_opt > 0 else 0

    return {
        "total_pdbs": len(results),
        "total_errors": total_errors,
        "total_dssr": total_dssr,
        "total_matched": total_match,
        "total_optimizer": total_opt,
        "recall": recall,
        "precision": precision,
        "misses_by_reason": {k: len(v) for k, v in misses_by_reason.items()},
        "misses_by_atom_pair": {f"{k[0]}-{k[1]}": len(v) for k, v in
                                sorted(misses_by_atom_pair.items(), key=lambda x: -len(x[1]))[:20]},
        "all_missed": all_missed,
        "all_extra": all_extra,
        "detailed_misses": misses_by_reason,
    }


def print_report(analysis: Dict):
    """Print a detailed report."""
    print("\n" + "="*80)
    print("H-BOND OPTIMIZER BENCHMARK REPORT")
    print("="*80)

    print(f"\nPDBs processed: {analysis['total_pdbs']}")
    print(f"Errors: {analysis['total_errors']}")

    print(f"\nH-Bond Statistics:")
    print(f"  DSSR H-bonds:    {analysis['total_dssr']}")
    print(f"  Matched:         {analysis['total_matched']}")
    print(f"  Optimizer total: {analysis['total_optimizer']}")
    print(f"  Missed:          {analysis['total_dssr'] - analysis['total_matched']}")
    print(f"  Extra:           {len(analysis['all_extra'])}")

    print(f"\nRecall:    {analysis['recall']:.2f}%")
    print(f"Precision: {analysis['precision']:.2f}%")

    print("\n" + "-"*80)
    print("MISSES BY REASON:")
    print("-"*80)
    for reason, count in sorted(analysis['misses_by_reason'].items(), key=lambda x: -x[1]):
        pct = count / (analysis['total_dssr'] - analysis['total_matched']) * 100 if analysis['total_dssr'] > analysis['total_matched'] else 0
        print(f"  {reason:25s}: {count:5d} ({pct:5.1f}%)")

    print("\n" + "-"*80)
    print("TOP 20 MISSED ATOM PAIRS:")
    print("-"*80)
    for pair, count in list(analysis['misses_by_atom_pair'].items())[:20]:
        print(f"  {pair:15s}: {count:5d}")

    # Detailed breakdown by reason
    print("\n" + "-"*80)
    print("DETAILED MISS ANALYSIS:")
    print("-"*80)

    detailed = analysis['detailed_misses']

    # NOT_DONOR_ACCEPTOR - likely need to add to tables
    if 'NOT_DONOR_ACCEPTOR' in detailed:
        print("\nNOT_DONOR_ACCEPTOR (need to add to geometry tables):")
        # Group by base_type.atom
        missing_entries = defaultdict(int)
        for miss in detailed['NOT_DONOR_ACCEPTOR']:
            missing_entries[miss.details] += 1
        for entry, count in sorted(missing_entries.items(), key=lambda x: -x[1])[:15]:
            print(f"  {entry}: {count}")

    # MISSING_RESIDUE - likely modified residues
    if 'MISSING_RESIDUE' in detailed:
        print("\nMISSING_RESIDUE (unrecognized residues - use original DSSR IDs):")
        missing_res = defaultdict(int)
        for miss in detailed['MISSING_RESIDUE']:
            # Use original DSSR residue IDs to show the actual modified residue codes
            for res in [miss.res1_original, miss.res2_original]:
                if res:
                    missing_res[res] += 1
        for res, count in sorted(missing_res.items(), key=lambda x: -x[1])[:20]:
            # Extract residue code for clarity
            code = extract_residue_code(res)
            print(f"  {res} ({code}): {count}")

    # NO_H_SLOTS / NO_LP_SLOTS
    for reason in ['NO_H_SLOTS', 'NO_LP_SLOTS']:
        if reason in detailed:
            print(f"\n{reason}:")
            slot_issues = defaultdict(int)
            for miss in detailed[reason]:
                slot_issues[miss.details] += 1
            for issue, count in sorted(slot_issues.items(), key=lambda x: -x[1])[:10]:
                print(f"  {issue}: {count}")

    # POOR_ALIGNMENT
    if 'POOR_ALIGNMENT' in detailed:
        print(f"\nPOOR_ALIGNMENT ({len(detailed['POOR_ALIGNMENT'])} total):")
        # Show distribution of alignment scores
        scores = []
        for miss in detailed['POOR_ALIGNMENT']:
            if 'alignment=' in miss.details:
                try:
                    score = float(miss.details.split('alignment=')[1].split()[0])
                    scores.append(score)
                except:
                    pass
        if scores:
            import statistics
            print(f"  Min alignment: {min(scores):.3f}")
            print(f"  Max alignment: {max(scores):.3f}")
            print(f"  Mean alignment: {statistics.mean(scores):.3f}")
            print(f"  Median alignment: {statistics.median(scores):.3f}")


def save_detailed_report(analysis: Dict, output_path: Path):
    """Save detailed report to JSON."""
    # Convert MissedHBond objects to dicts
    report = {
        "summary": {
            "total_pdbs": analysis['total_pdbs'],
            "total_errors": analysis['total_errors'],
            "total_dssr": analysis['total_dssr'],
            "total_matched": analysis['total_matched'],
            "total_optimizer": analysis['total_optimizer'],
            "recall": analysis['recall'],
            "precision": analysis['precision'],
        },
        "misses_by_reason": analysis['misses_by_reason'],
        "misses_by_atom_pair": analysis['misses_by_atom_pair'],
        "all_missed": [
            {
                "pdb_id": m.pdb_id,
                "res1": m.res1,
                "res2": m.res2,
                "res1_original": m.res1_original,
                "res2_original": m.res2_original,
                "donor_atom": m.donor_atom,
                "acceptor_atom": m.acceptor_atom,
                "distance": m.distance,
                "reason": m.reason,
                "details": m.details,
            }
            for m in analysis['all_missed']
        ],
    }

    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)

    print(f"\nDetailed report saved to: {output_path}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark H-bond optimizer vs DSSR')
    parser.add_argument('--test-set', default='100',
                       help='Test set: 100, fast, or comma-separated PDB IDs')
    parser.add_argument('--workers', type=int, default=None,
                       help=f'Number of workers (default: {mp.cpu_count()})')
    parser.add_argument('--max-dist', type=float, default=4.0,
                       help='Max distance threshold')
    parser.add_argument('--min-align', type=float, default=0.3,
                       help='Min alignment threshold')
    parser.add_argument('--short-dist', type=float, default=3.5,
                       help='Short distance threshold (skip alignment check below this)')
    parser.add_argument('--output', type=str, default=None,
                       help='Output JSON file for detailed report')
    args = parser.parse_args()

    # Find project root
    script_dir = Path(__file__).parent
    project_root = script_dir.parent.parent
    pdb_dir = project_root / "data" / "pdb"
    dssr_json_dir = project_root / "data" / "json_dssr"

    # Load test set
    if args.test_set == 'fast':
        with open(project_root / "data" / "valid_pdbs_fast.json") as f:
            data = json.load(f)
            pdb_ids = data['valid_pdbs_with_atoms_and_frames']
    elif args.test_set.isdigit():
        n = int(args.test_set)
        with open(project_root / "data" / "valid_pdbs_fast.json") as f:
            data = json.load(f)
            pdb_ids = data['valid_pdbs_with_atoms_and_frames'][:n]
    else:
        pdb_ids = [p.strip() for p in args.test_set.split(',')]

    print(f"Test set: {len(pdb_ids)} PDBs")
    print(f"Parameters: max_distance={args.max_dist}, min_alignment={args.min_align}, short_distance={args.short_dist}")

    # Run benchmark
    analysis = run_benchmark(
        pdb_ids, pdb_dir, dssr_json_dir,
        max_workers=args.workers,
        max_distance=args.max_dist,
        min_alignment=args.min_align,
        short_distance_threshold=args.short_dist
    )

    # Print report
    print_report(analysis)

    # Save detailed report
    if args.output:
        save_detailed_report(analysis, Path(args.output))
    else:
        default_output = project_root / "data" / "hbond_benchmark_report.json"
        save_detailed_report(analysis, default_output)


if __name__ == '__main__':
    main()

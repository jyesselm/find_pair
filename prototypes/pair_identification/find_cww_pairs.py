#!/usr/bin/env python3
"""Find cWW base pairs independently using distance filters and validation.

This script finds cWW pairs without relying on pre-computed pair lists:
1. Initial filter: C1'-C1' distance < 15Å
2. N1N9 distance filter: 8.0-10.0Å for cWW candidates
3. Template RMSD validation
4. H-bond pattern validation
5. Compare against DSSR ground truth
"""

import json
import sys
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

from core.pdb_parser import parse_pdb
from core.residue import Residue, is_purine, is_pyrimidine, normalize_base_type
from template_aligner import TemplateAligner
from cww_validator import WC_HBOND_PATTERNS, normalize_bp_type

# Import modified residue registry
sys.path.insert(0, str(Path(__file__).parent.parent / "hbond_optimizer"))
from modified_registry import get_parent_base


# N1N9 distance ranges for cWW pairs (from DSSR statistics analysis)
# Analysis showed expanding to (7.0, 11.0) adds too many false positives
N1N9_RANGE_CWW = (8.0, 10.0)

# Canonical Watson-Crick sequences (including DNA T equivalents)
CANONICAL_WC_SEQUENCES = {
    "GC", "CG", "AU", "UA", "GU", "UG",  # RNA
    "AT", "TA", "GT", "TG",  # DNA (T = U)
}


def normalize_res_id_base(res_id: str) -> str:
    """Normalize the base name in a res_id using modified registry.

    Converts modified bases to their parent: 5MC -> C, PSU -> U, etc.
    Also handles DNA: DA -> A, DG -> G, etc.

    Args:
        res_id: Residue ID in format "chain-base-num"

    Returns:
        Normalized res_id with parent base name
    """
    parts = res_id.split("-")
    if len(parts) < 2:
        return res_id

    base = parts[1]

    # Try modified registry first
    parent = get_parent_base(base)
    if parent:
        # Special cases: P (pseudouridine) -> U for H-bond purposes
        if parent == "P":
            parent = "U"
        elif parent == "I":  # Inosine -> G-like
            parent = "G"
        parts[1] = parent
        return "-".join(parts)

    # Handle DNA prefixes
    if base in ("DA", "DG", "DC", "DT"):
        parts[1] = base[1]
        return "-".join(parts)

    # Already standard
    return res_id

# Max C1'-C1' distance for initial filtering
MAX_C1_DISTANCE = 15.0


@dataclass
class CandidatePair:
    """A potential cWW pair candidate."""
    res_id1: str
    res_id2: str
    sequence: str
    c1_distance: float
    n1n9_distance: Optional[float] = None
    rmsd_cww: Optional[float] = None
    hbond_count: int = 0
    is_valid: bool = False
    score: float = 0.0


@dataclass
class FinderResult:
    """Result of pair finding for a single PDB."""
    pdb_id: str
    total_residues: int
    initial_candidates: int
    n1n9_filtered: int
    final_pairs: int
    pairs: List[CandidatePair] = field(default_factory=list)

    # Comparison with DSSR
    dssr_pairs: int = 0
    true_positives: int = 0
    false_positives: int = 0
    false_negatives: int = 0

    @property
    def precision(self) -> float:
        if self.true_positives + self.false_positives == 0:
            return 0.0
        return self.true_positives / (self.true_positives + self.false_positives)

    @property
    def recall(self) -> float:
        if self.dssr_pairs == 0:
            return 0.0
        return self.true_positives / self.dssr_pairs


def get_c1_atom(residue: Residue) -> Optional[np.ndarray]:
    """Get C1' atom coordinates from residue."""
    if "C1'" in residue.atoms:
        return residue.atoms["C1'"].coords
    return None


def get_n1n9_atom(residue: Residue) -> Optional[np.ndarray]:
    """Get N1 or N9 atom depending on base type."""
    atom = residue.get_glycosidic_n()
    return atom.coords if atom else None


def compute_base_normal(residue: Residue) -> Optional[np.ndarray]:
    """Compute base plane normal vector."""
    # Use 3 ring atoms to define plane
    if is_purine(residue.base_type):
        atoms_needed = ["N9", "C4", "C8"]
    else:
        atoms_needed = ["N1", "C2", "C4"]

    coords = []
    for name in atoms_needed:
        if name in residue.atoms:
            coords.append(residue.atoms[name].coords)

    if len(coords) < 3:
        return None

    v1 = coords[1] - coords[0]
    v2 = coords[2] - coords[0]
    normal = np.cross(v1, v2)
    norm = np.linalg.norm(normal)
    return normal / norm if norm > 0 else None


def compute_interbase_angle(res1: Residue, res2: Residue) -> float:
    """Compute angle between base planes in degrees."""
    n1 = compute_base_normal(res1)
    n2 = compute_base_normal(res2)

    if n1 is None or n2 is None:
        return 90.0  # Assume bad if can't compute

    dot = abs(np.dot(n1, n2))
    dot = min(1.0, max(-1.0, dot))
    return np.degrees(np.arccos(dot))


def find_initial_candidates(
    residues: Dict[str, Residue],
    max_distance: float = MAX_C1_DISTANCE,
) -> List[Tuple[str, str, float]]:
    """Find all residue pairs within max_distance using C1' atoms.

    Returns list of (res_id1, res_id2, distance) tuples.
    """
    candidates = []
    res_list = list(residues.items())

    for i, (res_id1, res1) in enumerate(res_list):
        c1_1 = get_c1_atom(res1)
        if c1_1 is None:
            continue

        # Only check pairs after this one (avoid duplicates)
        for res_id2, res2 in res_list[i+1:]:
            c1_2 = get_c1_atom(res2)
            if c1_2 is None:
                continue

            dist = np.linalg.norm(c1_2 - c1_1)
            if dist <= max_distance:
                candidates.append((res_id1, res_id2, dist))

    return candidates


def filter_by_n1n9(
    residues: Dict[str, Residue],
    candidates: List[Tuple[str, str, float]],
    n1n9_range: Tuple[float, float] = N1N9_RANGE_CWW,
    canonical_only: bool = True,
) -> List[CandidatePair]:
    """Filter candidates by N1N9 distance for cWW pairs.

    Args:
        residues: Dict of residues by res_id
        candidates: List of (res_id1, res_id2, c1_dist) tuples
        n1n9_range: Min/max N1N9 distance for cWW
        canonical_only: If True, only keep canonical WC sequences (GC, AU, GU, etc.)
    """
    filtered = []

    for res_id1, res_id2, c1_dist in candidates:
        res1 = residues[res_id1]
        res2 = residues[res_id2]

        # Get N1/N9 atoms
        n1 = get_n1n9_atom(res1)
        n2 = get_n1n9_atom(res2)

        if n1 is None or n2 is None:
            continue

        n1n9_dist = np.linalg.norm(n2 - n1)

        # Check if in cWW range
        if n1n9_range[0] <= n1n9_dist <= n1n9_range[1]:
            sequence = res1.base_type + res2.base_type

            # Filter to canonical WC sequences only
            if canonical_only and sequence not in CANONICAL_WC_SEQUENCES:
                continue

            filtered.append(CandidatePair(
                res_id1=res_id1,
                res_id2=res_id2,
                sequence=sequence,
                c1_distance=c1_dist,
                n1n9_distance=n1n9_dist,
            ))

    return filtered


def convert_to_dna_format(res_id: str) -> str:
    """Convert single-letter base to DNA format (A -> DA, etc.)."""
    parts = res_id.split("-")
    if len(parts) >= 2:
        base = parts[1]
        # Convert single-letter to DNA prefix
        if base in ("A", "G", "C", "T"):
            parts[1] = "D" + base
            return "-".join(parts)
    return res_id


def convert_from_dna_format(res_id: str) -> str:
    """Convert DNA format to single-letter (DA -> A, etc.)."""
    parts = res_id.split("-")
    if len(parts) >= 2:
        base = parts[1]
        if base in ("DA", "DG", "DC", "DT"):
            parts[1] = base[1]
            return "-".join(parts)
    return res_id


def load_slot_hbonds(hbond_path: Path) -> Dict[Tuple[str, str], List[Dict]]:
    """Load slot H-bonds from JSON file.

    Stores entries under multiple key formats to handle DNA/RNA naming variations
    and modified bases.
    """
    if not hbond_path.exists():
        return {}

    with open(hbond_path) as f:
        data = json.load(f)

    result = {}
    for entry in data:
        # Handle both naming conventions
        res_id1 = entry.get("res_id_i", entry.get("res_id1", ""))
        res_id2 = entry.get("res_id_j", entry.get("res_id2", ""))
        hbonds = entry.get("hbonds", [])

        # Generate all possible key formats
        ids1 = [res_id1, convert_from_dna_format(res_id1), convert_to_dna_format(res_id1),
                normalize_res_id_base(res_id1)]
        ids2 = [res_id2, convert_from_dna_format(res_id2), convert_to_dna_format(res_id2),
                normalize_res_id_base(res_id2)]

        # Remove duplicates
        ids1 = list(set(ids1))
        ids2 = list(set(ids2))

        # Store under all combinations
        for id1 in ids1:
            for id2 in ids2:
                key = (id1, id2)
                rev_key = (id2, id1)
                if key not in result:
                    result[key] = hbonds
                if rev_key not in result:
                    result[rev_key] = hbonds

    return result


def count_wc_hbonds(
    sequence: str,
    hbonds: List[Dict],
) -> Tuple[int, int]:
    """Count expected and found Watson-Crick H-bonds.

    Returns (expected_count, found_count).
    """
    normalized = normalize_bp_type(sequence)

    # Handle DNA equivalents (T = U for H-bond patterns)
    dna_to_rna = {"AT": "AU", "TA": "UA", "GT": "GU", "TG": "UG"}
    if normalized in dna_to_rna:
        normalized = dna_to_rna[normalized]

    pattern = WC_HBOND_PATTERNS.get(normalized, [])

    if not pattern:
        return 0, 0

    expected = len(pattern)
    found = 0

    # Extract atoms from H-bonds
    hbond_atoms = set()
    for hb in hbonds:
        if hb.get("context") == "base_base":
            donor = hb.get("donor_atom", "")
            acceptor = hb.get("acceptor_atom", "")
            hbond_atoms.add((donor, acceptor))
            hbond_atoms.add((acceptor, donor))  # Check both directions

    for _, donor_atom, _, acceptor_atom in pattern:
        if (donor_atom, acceptor_atom) in hbond_atoms:
            found += 1

    return expected, found


def select_best_pairs(validated: List[CandidatePair]) -> List[CandidatePair]:
    """Select best pair for each residue using greedy mutual selection.

    Each residue can only be in one pair. When a residue appears in multiple
    validated pairs, keep only the one with the highest score.
    """
    # Sort by score descending (best first)
    sorted_pairs = sorted(validated, key=lambda p: p.score, reverse=True)

    # Track which residues are already paired
    paired_residues = set()
    selected = []

    for pair in sorted_pairs:
        # Normalize res_ids for comparison
        res1 = normalize_res_id_base(pair.res_id1)
        res2 = normalize_res_id_base(pair.res_id2)

        # Skip if either residue is already paired
        if res1 in paired_residues or res2 in paired_residues:
            continue

        # Accept this pair
        selected.append(pair)
        paired_residues.add(res1)
        paired_residues.add(res2)

    return selected


def validate_candidates(
    residues: Dict[str, Residue],
    candidates: List[CandidatePair],
    slot_hbonds: Dict[Tuple[str, str], List[Dict]],
    template_aligner: Optional[TemplateAligner] = None,
    greedy_selection: bool = True,
) -> List[CandidatePair]:
    """Validate candidate pairs using H-bonds and template RMSD.

    Args:
        greedy_selection: If True, each residue can only be in one pair (best score wins)
    """
    validated = []

    for cand in candidates:
        res1 = residues[cand.res_id1]
        res2 = residues[cand.res_id2]

        # Get H-bonds for this pair
        key = (cand.res_id1, cand.res_id2)
        hbonds = slot_hbonds.get(key, [])

        # Count WC H-bonds
        expected, found = count_wc_hbonds(cand.sequence, hbonds)
        cand.hbond_count = found

        # Compute template RMSD if aligner available
        if template_aligner:
            result = template_aligner.classify_pair(res1, res2, lw_classes=["cWW"])
            if result.best_lw == "cWW":
                cand.rmsd_cww = result.best_rmsd

        # Compute interbase angle
        angle = compute_interbase_angle(res1, res2)

        # Scoring: prioritize H-bonds and geometry
        # H-bond score: 0.5 for each found
        hbond_score = min(found / max(expected, 1), 1.0)

        # RMSD score: 1.0 if <0.5Å, 0.0 if >1.5Å
        if cand.rmsd_cww is not None:
            if cand.rmsd_cww <= 0.5:
                rmsd_score = 1.0
            elif cand.rmsd_cww >= 1.5:
                rmsd_score = 0.0
            else:
                rmsd_score = 1.0 - (cand.rmsd_cww - 0.5) / 1.0
        else:
            rmsd_score = 0.5  # Unknown

        # Angle score: 1.0 if <15°, 0.0 if >30°
        if angle <= 15:
            angle_score = 1.0
        elif angle >= 30:
            angle_score = 0.0
        else:
            angle_score = 1.0 - (angle - 15) / 15

        # Combined score
        cand.score = 0.4 * hbond_score + 0.3 * rmsd_score + 0.3 * angle_score

        # Valid if:
        # 1. Score >= 0.6 AND at least 1 H-bond, OR
        # 2. RMSD < 0.5 AND at least 2 H-bonds (clear geometry match with H-bonds)
        # Requiring H-bonds prevents false positives from stacked non-WC pairs
        # Note: Lowering to 0.5 was tested but added too many false positives
        cand.is_valid = (
            (cand.score >= 0.6 and found >= 1) or
            (cand.rmsd_cww is not None and cand.rmsd_cww < 0.5 and found >= 2)
        )

        if cand.is_valid:
            validated.append(cand)

    # Apply greedy selection so each residue is in at most one pair
    if greedy_selection:
        validated = select_best_pairs(validated)

    return validated


def load_dssr_cww_pairs(dssr_path: Path) -> Set[Tuple[str, str]]:
    """Load cWW pairs from DSSR JSON as set of (res_id1, res_id2) tuples."""
    if not dssr_path.exists():
        return set()

    with open(dssr_path) as f:
        data = json.load(f)

    pairs = set()
    for p in data.get("pairs", []):
        if p.get("LW") != "cWW":
            continue

        nt1 = p.get("nt1", "")
        nt2 = p.get("nt2", "")

        # Convert DSSR format to our res_id format
        res_id1 = normalize_dssr_nt(nt1)
        res_id2 = normalize_dssr_nt(nt2)

        # Store as sorted tuple for order-independent comparison
        pairs.add(tuple(sorted([res_id1, res_id2])))

    return pairs


def normalize_dssr_nt(nt: str) -> str:
    """Convert DSSR nt format to res_id format.

    Keeps original base names (DA, DG, 2MG, etc.) without normalizing.
    Handles DSSR's special format for modified bases: AF2/103 -> AF2-103
    """
    if "." not in nt:
        return nt

    parts = nt.split(".")
    chain = parts[0]
    rest = parts[1]

    # Handle ^X alternate conformer suffix
    if "^" in rest:
        rest, alt = rest.split("^", 1)
    else:
        alt = ""

    # DSSR uses / before residue number for modified bases (e.g., AF2/103)
    # Split on / to separate base name from residue number
    if "/" in rest:
        base, num = rest.rsplit("/", 1)
        if alt:
            num += alt
        return f"{chain}-{base}-{num}"

    # Find where number starts (handle negative numbers too)
    # Look for last stretch of digits, possibly preceded by minus sign
    i = len(rest) - 1
    # Skip trailing insertion code (single letter at end after digit)
    if i >= 0 and rest[i].isalpha() and i > 0 and rest[i-1].isdigit():
        i -= 1
    # Find start of number
    while i >= 0 and rest[i].isdigit():
        i -= 1
    # Include minus sign if present
    if i >= 0 and rest[i] == '-':
        i -= 1
    i += 1

    if i <= 0:
        base = rest
        num = "0"
    else:
        base = rest[:i]
        num = rest[i:]

    # Keep base as-is (no normalization)

    # Append alt conformer to number
    if alt:
        num += alt

    return f"{chain}-{base}-{num}"


def find_cww_pairs_single(
    pdb_id: str,
    pdb_dir: Path,
    dssr_dir: Path,
    hbond_dir: Path,
    template_dir: Path,
) -> FinderResult:
    """Find cWW pairs for a single PDB."""
    pdb_path = pdb_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"

    result = FinderResult(pdb_id=pdb_id, total_residues=0,
                          initial_candidates=0, n1n9_filtered=0, final_pairs=0)

    if not pdb_path.exists():
        return result

    # Parse PDB
    residues = parse_pdb(pdb_path)
    result.total_residues = len(residues)

    # Find initial candidates
    initial = find_initial_candidates(residues, MAX_C1_DISTANCE)
    result.initial_candidates = len(initial)

    # Filter by N1N9 distance
    candidates = filter_by_n1n9(residues, initial, N1N9_RANGE_CWW)
    result.n1n9_filtered = len(candidates)

    # Load H-bonds
    hbond_path = hbond_dir / f"{pdb_id}.json"
    slot_hbonds = load_slot_hbonds(hbond_path)

    # Create template aligner
    aligner = TemplateAligner(idealized_dir=template_dir)

    # Validate candidates
    validated = validate_candidates(residues, candidates, slot_hbonds, aligner)
    result.final_pairs = len(validated)
    result.pairs = validated

    # Compare with DSSR
    dssr_path = dssr_dir / f"{pdb_id}.json"
    dssr_pairs_raw = load_dssr_cww_pairs(dssr_path)
    result.dssr_pairs = len(dssr_pairs_raw)

    # Normalize both sides to parent bases for comparison
    # This handles modified bases (5MC->C, PSU->U) and DNA (DA->A)
    dssr_pairs_norm = {
        tuple(sorted([normalize_res_id_base(r1), normalize_res_id_base(r2)]))
        for r1, r2 in dssr_pairs_raw
    }

    our_pairs_norm = {
        tuple(sorted([normalize_res_id_base(p.res_id1), normalize_res_id_base(p.res_id2)]))
        for p in validated
    }

    result.true_positives = len(our_pairs_norm & dssr_pairs_norm)
    result.false_positives = len(our_pairs_norm - dssr_pairs_norm)
    result.false_negatives = len(dssr_pairs_norm - our_pairs_norm)

    return result


def find_cww_pairs_batch(
    pdb_ids: List[str],
    pdb_dir: Path,
    dssr_dir: Path,
    hbond_dir: Path,
    template_dir: Path,
    max_workers: Optional[int] = None,
) -> List[FinderResult]:
    """Find cWW pairs for multiple PDBs in parallel."""
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()

    results = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                find_cww_pairs_single, pdb_id, pdb_dir, dssr_dir, hbond_dir, template_dir
            ): pdb_id
            for pdb_id in pdb_ids
        }

        for future in as_completed(futures):
            pdb_id = futures[future]
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                print(f"Error processing {pdb_id}: {e}")

    return results


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Find cWW base pairs using distance filters and validation"
    )

    # Default paths relative to project root
    project_root = Path(__file__).parent.parent.parent

    parser.add_argument(
        "--pdb-dir", type=Path,
        default=project_root / "data" / "pdb",
        help="Directory with PDB files"
    )
    parser.add_argument(
        "--dssr-dir", type=Path,
        default=project_root / "data" / "json_dssr",
        help="Directory with DSSR JSON files"
    )
    parser.add_argument(
        "--hbond-dir", type=Path,
        default=project_root / "data" / "json" / "slot_hbonds",
        help="Directory with slot H-bond JSON files"
    )
    parser.add_argument(
        "--template-dir", type=Path,
        default=project_root / "basepair-idealized",
        help="Directory with idealized templates"
    )
    parser.add_argument(
        "--pdb", type=str,
        help="Single PDB ID to analyze"
    )
    parser.add_argument(
        "--test-set", type=str,
        choices=["10", "50", "100", "fast"],
        default="100",
        help="Test set to use (default: 100)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Verbose output"
    )
    parser.add_argument(
        "--output", "-o", type=Path,
        help="Output JSON file for results"
    )

    args = parser.parse_args()

    # Single PDB mode
    if args.pdb:
        result = find_cww_pairs_single(
            args.pdb, args.pdb_dir, args.dssr_dir, args.hbond_dir, args.template_dir
        )

        print(f"\n{result.pdb_id}: {result.total_residues} residues")
        print(f"  Initial candidates (C1' < {MAX_C1_DISTANCE}Å): {result.initial_candidates}")
        print(f"  After N1N9 filter ({N1N9_RANGE_CWW[0]}-{N1N9_RANGE_CWW[1]}Å): {result.n1n9_filtered}")
        print(f"  Final validated pairs: {result.final_pairs}")
        print(f"\n  DSSR comparison:")
        print(f"    DSSR cWW pairs: {result.dssr_pairs}")
        print(f"    True positives: {result.true_positives}")
        print(f"    False positives: {result.false_positives}")
        print(f"    False negatives: {result.false_negatives}")
        print(f"    Precision: {result.precision:.2%}")
        print(f"    Recall: {result.recall:.2%}")

        if args.verbose and result.pairs:
            print(f"\n  Found pairs:")
            for p in result.pairs[:20]:
                rmsd_str = f"{p.rmsd_cww:.2f}" if p.rmsd_cww else "N/A"
                print(f"    {p.res_id1} - {p.res_id2} ({p.sequence}) "
                      f"N1N9={p.n1n9_distance:.2f}Å RMSD={rmsd_str} "
                      f"H-bonds={p.hbond_count} score={p.score:.2f}")

        return

    # Batch mode
    test_set_file = project_root / "data" / f"test_set_{args.test_set}_list.txt"
    if not test_set_file.exists():
        # Try alternate location
        test_set_file = project_root / "data" / "test_sets" / f"{args.test_set}.txt"

    if test_set_file.exists():
        with open(test_set_file) as f:
            pdb_ids = [line.strip() for line in f if line.strip()]
    else:
        # Fall back to listing PDB files
        pdb_ids = [p.stem.upper() for p in args.pdb_dir.glob("*.pdb")][:100]

    print(f"Finding cWW pairs for {len(pdb_ids)} PDBs...")

    results = find_cww_pairs_batch(
        pdb_ids, args.pdb_dir, args.dssr_dir, args.hbond_dir, args.template_dir
    )

    # Aggregate statistics
    total_dssr = sum(r.dssr_pairs for r in results)
    total_tp = sum(r.true_positives for r in results)
    total_fp = sum(r.false_positives for r in results)
    total_fn = sum(r.false_negatives for r in results)
    total_found = sum(r.final_pairs for r in results)

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / total_dssr if total_dssr > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    print(f"\n{'='*60}")
    print(f"AGGREGATE RESULTS ({len(results)} PDBs)")
    print(f"{'='*60}")
    print(f"DSSR cWW pairs: {total_dssr}")
    print(f"Our pairs found: {total_found}")
    print(f"True positives: {total_tp}")
    print(f"False positives: {total_fp}")
    print(f"False negatives: {total_fn}")
    print(f"\nPrecision: {precision:.2%}")
    print(f"Recall: {recall:.2%}")
    print(f"F1 Score: {f1:.2%}")

    # Save results if output specified
    if args.output:
        output_data = {
            "summary": {
                "pdb_count": len(results),
                "dssr_pairs": total_dssr,
                "found_pairs": total_found,
                "true_positives": total_tp,
                "false_positives": total_fp,
                "false_negatives": total_fn,
                "precision": precision,
                "recall": recall,
                "f1": f1,
            },
            "per_pdb": [
                {
                    "pdb_id": r.pdb_id,
                    "dssr_pairs": r.dssr_pairs,
                    "found_pairs": r.final_pairs,
                    "true_positives": r.true_positives,
                    "false_positives": r.false_positives,
                    "false_negatives": r.false_negatives,
                }
                for r in results
            ]
        }

        with open(args.output, "w") as f:
            json.dump(output_data, f, indent=2)
        print(f"\nResults saved to {args.output}")


if __name__ == "__main__":
    main()

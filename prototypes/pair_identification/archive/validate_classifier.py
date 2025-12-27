#!/usr/bin/env python3
"""Validate pair classifier against DSSR output.

Tests classification on well-matched pair types and compares results
with DSSR's LW classifications.
"""

import json
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import numpy as np

from prototypes.pair_identification.template_overlay import (
    parse_pdb_residues, parse_template_pdb, align_template_to_target,
    find_template, get_base_letter, load_dssr_pairs, Residue
)


@dataclass
class ValidationResult:
    """Result of validating a single pair."""
    pdb_id: str
    res1_id: str
    res2_id: str
    sequence: str
    dssr_lw: str
    predicted_lw: str
    predicted_rmsd: float
    correct: bool
    confidence: float


def get_template_n1n9_distance(template_path: Path) -> Optional[float]:
    """Get N1N9 distance from a template PDB."""
    res_atoms: Dict[int, Dict[str, np.ndarray]] = {}

    with open(template_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                name = line[12:16].strip()
                res = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if res not in res_atoms:
                    res_atoms[res] = {}
                res_atoms[res][name] = np.array([x, y, z])

    res_list = sorted(res_atoms.keys())
    if len(res_list) < 2:
        return None

    res1, res2 = res_list[0], res_list[1]

    n1 = res_atoms[res1].get("N1")
    if n1 is None:
        n1 = res_atoms[res1].get("N9")
    n2 = res_atoms[res2].get("N1")
    if n2 is None:
        n2 = res_atoms[res2].get("N9")

    if n1 is None or n2 is None:
        return None

    return float(np.linalg.norm(n2 - n1))


def get_pair_n1n9_distance(res1: Residue, res2: Residue) -> Optional[float]:
    """Get N1N9 distance from a pair of residues."""
    n1 = res1.atoms.get("N1")
    if n1 is None:
        n1 = res1.atoms.get("N9")
    n2 = res2.atoms.get("N1")
    if n2 is None:
        n2 = res2.atoms.get("N9")

    if n1 is None or n2 is None:
        return None

    return float(np.linalg.norm(n2.coords - n1.coords))


def classify_pair_by_n1n9(
    target_res1: Residue,
    target_res2: Residue,
    sequence: str,
    idealized_dir: Path,
    exemplar_dir: Path,
    lw_classes: List[str],
) -> Tuple[str, float, float]:
    """Classify a pair by matching N1N9 distance to templates.

    Returns:
        Tuple of (best_lw_class, n1n9_diff, confidence)
    """
    # Get measured N1N9 distance
    measured_n1n9 = get_pair_n1n9_distance(target_res1, target_res2)
    if measured_n1n9 is None:
        return None, float('inf'), 0.0

    best_lw = None
    best_diff = float('inf')
    second_best_diff = float('inf')

    for lw in lw_classes:
        # Try forward sequence
        template_path = find_template(sequence, lw, idealized_dir, exemplar_dir)

        # Try reversed sequence if forward not found
        if not template_path:
            rev_seq = sequence[1] + sequence[0]
            template_path = find_template(rev_seq, lw, idealized_dir, exemplar_dir)

        if not template_path:
            continue

        template_n1n9 = get_template_n1n9_distance(template_path)
        if template_n1n9 is None:
            continue

        diff = abs(measured_n1n9 - template_n1n9)

        if diff < best_diff:
            second_best_diff = best_diff
            best_diff = diff
            best_lw = lw
        elif diff < second_best_diff:
            second_best_diff = diff

    # Confidence based on gap between best and second best
    if best_lw and second_best_diff < float('inf'):
        gap = second_best_diff - best_diff
        confidence = min(1.0, gap / 1.0)  # Full confidence if gap > 1 Å
    else:
        confidence = 0.5 if best_lw else 0.0

    return best_lw, best_diff, confidence


def classify_pair_by_rmsd(
    target_res1: Residue,
    target_res2: Residue,
    sequence: str,
    idealized_dir: Path,
    exemplar_dir: Path,
    lw_classes: List[str],
) -> Tuple[str, float, float]:
    """Classify a pair by finding the LW class with lowest RMSD.

    Returns:
        Tuple of (best_lw_class, best_rmsd, confidence)
    """
    best_lw = None
    best_rmsd = float('inf')
    second_best_rmsd = float('inf')

    for lw in lw_classes:
        # Try forward sequence
        template_path = find_template(sequence, lw, idealized_dir, exemplar_dir)
        res1, res2 = target_res1, target_res2

        # Try reversed sequence if forward not found
        if not template_path:
            rev_seq = sequence[1] + sequence[0]
            template_path = find_template(rev_seq, lw, idealized_dir, exemplar_dir)
            if template_path:
                res1, res2 = target_res2, target_res1

        if not template_path:
            continue

        try:
            template_res1, template_res2 = parse_template_pdb(template_path)
            _, _, rmsd = align_template_to_target(
                template_res1, template_res2, res1, res2
            )

            if rmsd < best_rmsd:
                second_best_rmsd = best_rmsd
                best_rmsd = rmsd
                best_lw = lw
            elif rmsd < second_best_rmsd:
                second_best_rmsd = rmsd

        except Exception:
            continue

    # Confidence based on gap between best and second best
    if best_lw and second_best_rmsd < float('inf'):
        gap = second_best_rmsd - best_rmsd
        confidence = min(1.0, gap / 0.5)  # Full confidence if gap > 0.5 Å
    else:
        confidence = 0.5 if best_lw else 0.0

    return best_lw, best_rmsd, confidence


def validate_on_dssr(
    dssr_dir: Path,
    pdb_dir: Path,
    idealized_dir: Path,
    exemplar_dir: Path,
    target_lw_classes: List[str],
    max_pdbs: int = 100,
    max_pairs_per_type: int = 50,
    verbose: bool = False,
) -> Dict[str, List[ValidationResult]]:
    """Validate classifier against DSSR on multiple PDBs.

    Args:
        dssr_dir: Directory with DSSR JSON files
        pdb_dir: Directory with PDB files
        idealized_dir: Directory with idealized templates
        exemplar_dir: Directory with exemplar templates
        target_lw_classes: LW classes to test
        max_pdbs: Maximum PDBs to process
        max_pairs_per_type: Maximum pairs per LW type
        verbose: Print progress

    Returns:
        Dict mapping LW class to list of validation results
    """
    results: Dict[str, List[ValidationResult]] = defaultdict(list)
    pdb_cache: Dict[str, Dict[str, Residue]] = {}
    pairs_per_type: Dict[str, int] = defaultdict(int)

    dssr_files = sorted(dssr_dir.glob("*.json"))[:max_pdbs * 10]  # Over-sample to find enough

    pdbs_processed = 0

    for dssr_path in dssr_files:
        if pdbs_processed >= max_pdbs:
            break

        # Check if we have enough pairs for all types
        all_done = all(
            pairs_per_type[lw] >= max_pairs_per_type
            for lw in target_lw_classes
        )
        if all_done:
            break

        pdb_id = dssr_path.stem
        pairs = load_dssr_pairs(dssr_path)

        found_any = False

        for pair in pairs:
            dssr_lw = pair.get("LW", "")
            if dssr_lw not in target_lw_classes:
                continue

            if pairs_per_type[dssr_lw] >= max_pairs_per_type:
                continue

            bp = pair.get("bp", "")
            if not bp or "-" not in bp:
                continue

            parts = bp.split("-")
            if len(parts) != 2:
                continue

            seq = parts[0][-1].upper() + parts[1][-1].upper()
            res1_id = pair["nt1"]
            res2_id = pair["nt2"]

            # Load PDB if needed
            if pdb_id not in pdb_cache:
                pdb_path = pdb_dir / f"{pdb_id}.pdb"
                if not pdb_path.exists():
                    pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"
                if pdb_path.exists():
                    pdb_cache[pdb_id] = parse_pdb_residues(pdb_path)
                else:
                    pdb_cache[pdb_id] = None

            if pdb_cache[pdb_id] is None:
                continue

            residues = pdb_cache[pdb_id]
            if res1_id not in residues or res2_id not in residues:
                continue

            target_res1 = residues[res1_id]
            target_res2 = residues[res2_id]

            # Classify by N1N9 distance (more discriminative than RMSD)
            predicted_lw, score, confidence = classify_pair_by_n1n9(
                target_res1, target_res2, seq,
                idealized_dir, exemplar_dir,
                target_lw_classes
            )

            if predicted_lw is None:
                continue

            correct = (predicted_lw == dssr_lw)

            result = ValidationResult(
                pdb_id=pdb_id,
                res1_id=res1_id,
                res2_id=res2_id,
                sequence=seq,
                dssr_lw=dssr_lw,
                predicted_lw=predicted_lw,
                predicted_rmsd=score,  # Actually N1N9 diff now
                correct=correct,
                confidence=confidence,
            )

            results[dssr_lw].append(result)
            pairs_per_type[dssr_lw] += 1
            found_any = True

            if verbose:
                status = "✓" if correct else "✗"
                print(f"  {status} {pdb_id} {res1_id}-{res2_id} ({seq}): "
                      f"DSSR={dssr_lw}, pred={predicted_lw}, N1N9_diff={score:.3f}")

        if found_any:
            pdbs_processed += 1
            if verbose:
                print(f"Processed {pdb_id} ({pdbs_processed}/{max_pdbs})")

    return results


def print_validation_summary(results: Dict[str, List[ValidationResult]]) -> None:
    """Print summary of validation results."""
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)

    total_correct = 0
    total_pairs = 0

    print(f"\n{'LW Class':<10} {'Total':>6} {'Correct':>8} {'Accuracy':>10} {'Avg RMSD':>10}")
    print("-" * 50)

    for lw_class in sorted(results.keys()):
        lw_results = results[lw_class]
        total = len(lw_results)
        correct = sum(1 for r in lw_results if r.correct)
        accuracy = 100 * correct / total if total > 0 else 0
        avg_rmsd = np.mean([r.predicted_rmsd for r in lw_results]) if lw_results else 0

        print(f"{lw_class:<10} {total:>6} {correct:>8} {accuracy:>9.1f}% {avg_rmsd:>9.3f} Å")

        total_correct += correct
        total_pairs += total

    print("-" * 50)
    overall_accuracy = 100 * total_correct / total_pairs if total_pairs > 0 else 0
    print(f"{'OVERALL':<10} {total_pairs:>6} {total_correct:>8} {overall_accuracy:>9.1f}%")

    # Confusion analysis
    print("\n" + "=" * 70)
    print("CONFUSION ANALYSIS (misclassifications)")
    print("=" * 70)

    confusions: Dict[Tuple[str, str], int] = defaultdict(int)
    for lw_class, lw_results in results.items():
        for r in lw_results:
            if not r.correct:
                confusions[(r.dssr_lw, r.predicted_lw)] += 1

    if confusions:
        print(f"\n{'DSSR':<10} {'Predicted':<10} {'Count':>6}")
        print("-" * 30)
        for (dssr, pred), count in sorted(confusions.items(), key=lambda x: -x[1]):
            print(f"{dssr:<10} {pred:<10} {count:>6}")
    else:
        print("\nNo misclassifications!")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Validate pair classifier against DSSR")
    parser.add_argument(
        "--dssr-dir", type=Path, default=Path("data/json_dssr"),
        help="Directory with DSSR JSON files"
    )
    parser.add_argument(
        "--pdb-dir", type=Path, default=Path("data/pdb"),
        help="Directory with PDB files"
    )
    parser.add_argument(
        "--idealized-dir", type=Path, default=Path("basepair-idealized"),
        help="Directory with idealized templates"
    )
    parser.add_argument(
        "--exemplar-dir", type=Path, default=Path("basepair-exemplars"),
        help="Directory with exemplar templates"
    )
    parser.add_argument(
        "--max-pdbs", type=int, default=100,
        help="Maximum PDBs to process"
    )
    parser.add_argument(
        "--max-pairs", type=int, default=50,
        help="Maximum pairs per LW type"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Verbose output"
    )
    parser.add_argument(
        "--lw-classes", type=str, nargs="+",
        default=["cWW", "tSH", "tHW", "cSW", "tWH"],
        help="LW classes to test"
    )

    args = parser.parse_args()

    print(f"Validating classifier on LW classes: {args.lw_classes}")
    print(f"Max PDBs: {args.max_pdbs}, Max pairs per type: {args.max_pairs}")

    results = validate_on_dssr(
        dssr_dir=args.dssr_dir,
        pdb_dir=args.pdb_dir,
        idealized_dir=args.idealized_dir,
        exemplar_dir=args.exemplar_dir,
        target_lw_classes=args.lw_classes,
        max_pdbs=args.max_pdbs,
        max_pairs_per_type=args.max_pairs,
        verbose=args.verbose,
    )

    print_validation_summary(results)


if __name__ == "__main__":
    main()

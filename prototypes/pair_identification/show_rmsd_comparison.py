#!/usr/bin/env python3
"""Show detailed RMSD comparison for misclassified pairs."""

import json
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import sys

sys.path.insert(0, str(Path(__file__).parent))
from template_aligner import TemplateAligner, parse_pdb_residues, Residue
from lw_classifier import LWClassifier, normalize_dssr_nt
from hbond_scorer import load_modern_hbonds


def load_template_coords(template_path: Path) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Load template PDB and return atom coords for both residues."""
    res1_atoms = {}
    res2_atoms = {}

    with open(template_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            res_seq = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            if res_seq == 1:
                res1_atoms[atom_name] = np.array([x, y, z])
            else:
                res2_atoms[atom_name] = np.array([x, y, z])

    return res1_atoms, res2_atoms


def align_and_get_details(
    target_res1: Residue,
    target_res2: Residue,
    template_path: Path,
    ring_atoms: List[str]
) -> Dict:
    """Align target to template and return detailed results."""
    template_res1, template_res2 = load_template_coords(template_path)

    # Collect matching atoms
    template_points = []
    target_points = []
    atom_names = []

    for atom_name in ring_atoms:
        if atom_name in template_res1 and atom_name in target_res1.atoms:
            template_points.append(template_res1[atom_name])
            target_points.append(target_res1.atoms[atom_name].coords)
            atom_names.append(f"res1.{atom_name}")
        if atom_name in template_res2 and atom_name in target_res2.atoms:
            template_points.append(template_res2[atom_name])
            target_points.append(target_res2.atoms[atom_name].coords)
            atom_names.append(f"res2.{atom_name}")

    if len(template_points) < 4:
        return {"error": "Too few atoms"}

    template_points = np.array(template_points)
    target_points = np.array(target_points)

    # Kabsch alignment
    centroid_T = np.mean(template_points, axis=0)
    centroid_P = np.mean(target_points, axis=0)

    T_centered = template_points - centroid_T
    P_centered = target_points - centroid_P

    H = T_centered.T @ P_centered
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # Apply transformation to template
    aligned_template = np.array([R @ (p - centroid_T) + centroid_P for p in template_points])

    # Per-atom deviations
    deviations = np.linalg.norm(aligned_template - target_points, axis=1)
    rmsd = float(np.sqrt(np.mean(deviations**2)))

    return {
        "rmsd": rmsd,
        "num_atoms": len(atom_names),
        "atom_names": atom_names,
        "deviations": deviations.tolist(),
        "max_deviation": float(np.max(deviations)),
        "template_aligned": aligned_template.tolist(),
        "target_coords": target_points.tolist(),
    }


def analyze_miss(
    pdb_id: str,
    res_id1: str,
    res_id2: str,
    pdb_dir: Path,
    hbond_dir: Path,
    dssr_dir: Path,
    idealized_dir: Path,
    exemplar_dir: Path,
):
    """Analyze a single misclassified pair in detail."""

    pdb_path = pdb_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"

    residues = parse_pdb_residues(pdb_path)
    hbonds = load_modern_hbonds(hbond_dir / f"{pdb_id}.json")

    if res_id1 not in residues or res_id2 not in residues:
        print(f"Residue not found: {res_id1} or {res_id2}")
        return

    res1 = residues[res_id1]
    res2 = residues[res_id2]
    seq = res1.base_type + res2.base_type

    pair_key = tuple(sorted([res_id1, res_id2]))
    pair_hbonds = hbonds.get(pair_key, [])

    print(f"\n{'='*70}")
    print(f"Pair: {res_id1} - {res_id2} ({seq})")
    print(f"{'='*70}")

    print(f"\nObserved H-bonds ({len(pair_hbonds)}):")
    for hb in pair_hbonds:
        print(f"  {hb.donor_atom} -> {hb.acceptor_atom} ({hb.distance:.2f}Å)")

    # Load DSSR data for this pair
    dssr_path = dssr_dir / f"{pdb_id}.json"
    if dssr_path.exists():
        with open(dssr_path) as f:
            dssr_data = json.load(f)
        for p in dssr_data.get("pairs", []):
            r1 = normalize_dssr_nt(p.get("nt1", ""))
            r2 = normalize_dssr_nt(p.get("nt2", ""))
            if tuple(sorted([r1, r2])) == pair_key:
                print(f"\nDSSR classification: {p.get('LW', 'unknown')}")
                print(f"DSSR H-bonds: {p.get('hbonds_desc', 'N/A')}")
                print(f"DSSR N1N9_dist: {p.get('N1N9_dist', 0):.3f}Å")
                print(f"DSSR interBase_angle: {p.get('interBase_angle', 0):.1f}°")
                break

    # Compare alignments to different templates
    aligner = TemplateAligner(idealized_dir, exemplar_dir)
    ring_atoms = ["C2", "C4", "C5", "C6", "N1", "N3", "N7", "C8", "N9"]

    print(f"\n{'='*70}")
    print(f"Template Alignment Comparison")
    print(f"{'='*70}")

    templates_to_try = [
        ("cWW", seq),
        ("tWW", seq),
        ("cWS", seq),
        ("tHS", seq),
    ]

    for lw_class, template_seq in templates_to_try:
        # Try to find template
        template_path = idealized_dir / lw_class / f"{template_seq}.pdb"
        if not template_path.exists():
            template_path = idealized_dir / lw_class / f"{template_seq[0]}_{template_seq[1].lower()}.pdb"
        if not template_path.exists():
            # Try reversed
            rev_seq = template_seq[1] + template_seq[0]
            template_path = idealized_dir / lw_class / f"{rev_seq}.pdb"
            if not template_path.exists():
                template_path = idealized_dir / lw_class / f"{rev_seq[0]}_{rev_seq[1].lower()}.pdb"

        if not template_path.exists():
            print(f"\n{lw_class} ({template_seq}): No template found")
            continue

        # Get alignment details
        if template_seq == seq:
            details = align_and_get_details(res1, res2, template_path, ring_atoms)
        else:
            details = align_and_get_details(res2, res1, template_path, ring_atoms)

        if "error" in details:
            print(f"\n{lw_class} ({template_seq}): {details['error']}")
            continue

        print(f"\n{lw_class} template ({template_path.name}):")
        print(f"  RMSD: {details['rmsd']:.4f}Å ({details['num_atoms']} atoms)")
        print(f"  Max deviation: {details['max_deviation']:.4f}Å")
        print(f"  Per-atom deviations:")
        for name, dev in zip(details['atom_names'], details['deviations']):
            marker = " ***" if dev > 0.5 else ""
            print(f"    {name:12s}: {dev:.4f}Å{marker}")


def find_rmsd_misses(
    pdb_ids: List[str],
    pdb_dir: Path,
    hbond_dir: Path,
    dssr_dir: Path,
    idealized_dir: Path,
    exemplar_dir: Path,
    max_examples: int = 10,
) -> List[Dict]:
    """Find examples where RMSD prefers non-cWW."""

    classifier = LWClassifier(idealized_dir, exemplar_dir)
    examples = []

    for pdb_id in pdb_ids:
        if len(examples) >= max_examples:
            break

        pdb_path = pdb_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"
        if not pdb_path.exists():
            continue

        dssr_path = dssr_dir / f"{pdb_id}.json"
        if not dssr_path.exists():
            continue

        hbond_path = hbond_dir / f"{pdb_id}.json"
        if not hbond_path.exists():
            continue

        try:
            residues = parse_pdb_residues(pdb_path)
            hbonds = load_modern_hbonds(hbond_path)

            with open(dssr_path) as f:
                dssr_data = json.load(f)

            for pair in dssr_data.get("pairs", []):
                if len(examples) >= max_examples:
                    break

                if pair.get("LW") != "cWW":
                    continue

                bp = pair.get("bp", "")
                if not bp or bp[0] not in "GCAU" or bp[-1] not in "GCAU":
                    continue

                seq = bp[0] + bp[-1]
                if seq not in ["GC", "CG", "AU", "UA"]:
                    continue

                res_id1 = normalize_dssr_nt(pair["nt1"])
                res_id2 = normalize_dssr_nt(pair["nt2"])

                if res_id1 not in residues or res_id2 not in residues:
                    continue

                res1 = residues[res_id1]
                res2 = residues[res_id2]
                pair_key = tuple(sorted([res_id1, res_id2]))
                pair_hbonds = hbonds.get(pair_key, [])

                result = classifier.classify(res1, res2, pair_hbonds)

                # Check if it's an RMSD-prefers-other miss
                if result.best_lw != "cWW" and result.best_hbond_score >= 0.5:
                    examples.append({
                        "pdb_id": pdb_id,
                        "res_id1": res_id1,
                        "res_id2": res_id2,
                        "seq": seq,
                        "predicted": result.best_lw,
                        "pred_rmsd": result.best_rmsd,
                        "pred_hbond": result.best_hbond_score,
                    })
        except Exception as e:
            print(f"Error processing {pdb_id}: {e}")

    return examples


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Show RMSD comparison for misses")
    parser.add_argument("--pdb-dir", type=Path, default=Path("data/pdb"))
    parser.add_argument("--hbond-dir", type=Path, default=Path("data/json/all_hbond_list"))
    parser.add_argument("--dssr-dir", type=Path, default=Path("data/json_dssr"))
    parser.add_argument("--idealized-dir", type=Path, default=Path("basepair-idealized"))
    parser.add_argument("--exemplar-dir", type=Path, default=Path("basepair-exemplars"))
    parser.add_argument("--pdb", type=str, default=None, help="Specific PDB to analyze")
    parser.add_argument("--res1", type=str, default=None, help="First residue ID")
    parser.add_argument("--res2", type=str, default=None, help="Second residue ID")
    parser.add_argument("--find-examples", type=int, default=5, help="Find N examples automatically")
    parser.add_argument("--pdb-list", type=Path, default=None)

    args = parser.parse_args()

    if args.pdb and args.res1 and args.res2:
        # Analyze specific pair
        analyze_miss(
            args.pdb, args.res1, args.res2,
            args.pdb_dir, args.hbond_dir, args.dssr_dir,
            args.idealized_dir, args.exemplar_dir
        )
    else:
        # Find examples automatically
        if args.pdb_list:
            with open(args.pdb_list) as f:
                pdb_list_data = json.load(f)
            if isinstance(pdb_list_data, list):
                pdb_ids = pdb_list_data
            else:
                pdb_ids = list(next(iter(pdb_list_data.values())))
        else:
            pdb_ids = sorted([f.stem for f in args.dssr_dir.glob("*.json")])

        print(f"Finding {args.find_examples} RMSD-prefers-other examples...")
        examples = find_rmsd_misses(
            pdb_ids[:500],  # Limit search
            args.pdb_dir, args.hbond_dir, args.dssr_dir,
            args.idealized_dir, args.exemplar_dir,
            max_examples=args.find_examples
        )

        print(f"\nFound {len(examples)} examples")

        for ex in examples:
            analyze_miss(
                ex["pdb_id"], ex["res_id1"], ex["res_id2"],
                args.pdb_dir, args.hbond_dir, args.dssr_dir,
                args.idealized_dir, args.exemplar_dir
            )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Visualize specific outliers from cww_miss_annotator results."""

import argparse
import json
import sys
from pathlib import Path
from typing import List, Dict, Optional

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

from visualization import PairVisualizer
from cww_miss_annotator.bp_score import compute_bp_score, score_to_grade


def load_outliers(results_dir: Path, reason_filter: Optional[str] = None, pdb_filter: Optional[str] = None, include_dssr_questionable: bool = False) -> List[Dict]:
    """Load outliers from analysis results.

    Args:
        results_dir: Directory containing individual PDB reports
        reason_filter: Optional filter for specific reason code
        pdb_filter: Optional filter for specific PDB ID
        include_dssr_questionable: If False (default), exclude dssr_questionable pairs

    Returns:
        List of outlier dicts with pdb_id, res_id1, res_id2, sequence, reasons, bp_score
    """
    outliers = []

    for json_file in sorted(results_dir.glob("*.json")):
        # Apply PDB filter early
        if pdb_filter and json_file.stem.upper() != pdb_filter.upper():
            continue
        if json_file.name == "aggregate.json":
            continue

        try:
            with open(json_file) as f:
                report = json.load(f)
        except:
            continue

        pdb_id = report.get("pdb_id", json_file.stem)

        for fn in report.get("false_negatives", []):
            reasons = fn.get("reasons", [])

            # Skip dssr_questionable pairs by default (they're DSSR errors, not our misses)
            if not include_dssr_questionable and "dssr_questionable" in reasons:
                continue

            # Apply filter if specified
            if reason_filter and reason_filter not in reasons:
                continue

            # Extract data for BP score
            sequence = fn["sequence"]
            rmsd_cww = fn.get("geometric_diagnostics", {}).get("rmsd_cww", 1.0)
            found_hbonds = fn.get("hbond_diagnostics", {}).get("found_hbonds", [])

            # Compute BP score
            bp_score, bp_components = compute_bp_score(sequence, rmsd_cww, found_hbonds)
            bp_grade = score_to_grade(bp_score)

            outliers.append({
                "pdb_id": pdb_id,
                "res_id1": fn["res_id1"],
                "res_id2": fn["res_id2"],
                "sequence": sequence,
                "reasons": reasons,
                "our_prediction": fn.get("our_prediction", "unknown"),
                "rmsd_cww": rmsd_cww,
                "rmsd_best": fn.get("geometric_diagnostics", {}).get("rmsd_best", 0),
                "best_lw": fn.get("geometric_diagnostics", {}).get("best_lw", "unknown"),
                "bp_score": bp_score,
                "bp_grade": bp_grade,
                "bp_components": bp_components,
            })

    # Sort by BP score (best to worst = highest to lowest)
    outliers.sort(key=lambda x: x["bp_score"], reverse=True)

    return outliers


def list_outliers(outliers: List[Dict], limit: int = 50):
    """Print list of outliers sorted by BP score (best to worst)."""
    print(f"\n{'#':<4} {'Score':<8} {'PDB':<6} {'Pair':<22} {'Seq':<4} {'RMSD':<6} {'Reasons'}")
    print("-" * 110)

    for i, o in enumerate(outliers[:limit]):
        pair = f"{o['res_id1']} - {o['res_id2']}"
        reasons = ", ".join(o['reasons'][:3])
        if len(o['reasons']) > 3:
            reasons += f" (+{len(o['reasons'])-3})"
        rmsd = f"{o['rmsd_cww']:.2f}" if o['rmsd_cww'] else "N/A"
        score_str = f"{o['bp_score']:.2f}({o['bp_grade']})"

        print(f"{i:<4} {score_str:<8} {o['pdb_id']:<6} {pair:<22} {o['sequence']:<4} {rmsd:<6} {reasons}")

    if len(outliers) > limit:
        print(f"\n... and {len(outliers) - limit} more. Use --limit to show more.")

    # Print score explanation
    print("\n" + "=" * 60)
    print("BP SCORE EXPLANATION")
    print("=" * 60)
    print("Score = 0.30×RMSD + 0.40×Coverage + 0.30×Quality")
    print()
    print("Components:")
    print("  RMSD:     1.0 if ≤0.3Å, 0.0 if ≥1.0Å, linear between")
    print("  Coverage: (found H-bonds) / (expected H-bonds)")
    print("            GC/CG expect 3, AU/UA expect 2")
    print("  Quality:  H-bond quality with geometry-adjusted leniency")
    print("            Good geometry (RMSD<0.5Å): accept up to 4.2Å H-bonds")
    print("            Poor geometry (RMSD>0.8Å): strict 3.2Å threshold")
    print()
    print("Grades: A≥0.9, B≥0.8, C≥0.7, D≥0.6, F<0.6")
    print()
    print("NOTE: Scores shown here are from saved H-bonds only.")
    print("      The annotator also uses extended H-bond search (up to 5Å)")
    print("      which may find additional H-bonds for pairs with good geometry.")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="List and visualize outliers from cWW analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # List all outliers
  %(prog)s --results /tmp/cww_test_fixed

  # Filter by reason
  %(prog)s --results /tmp/cww_test_fixed --reason geometric_outlier

  # Visualize specific outlier by index
  %(prog)s --results /tmp/cww_test_fixed --render 5

  # Visualize by PDB and residues
  %(prog)s --render-pair 1EHZ A-G-1 A-C-72
        """,
    )

    parser.add_argument(
        "--results", "-r",
        type=Path,
        default=Path(__file__).parent / "analysis_results" / "cww_analysis_v6",
        help="Results directory from cww_annotate.py",
    )
    parser.add_argument(
        "--reason",
        type=str,
        choices=[
            "no_hbonds", "missing_hbonds", "long_hbonds", "short_hbonds",
            "wrong_hbonds", "extra_hbonds", "poor_planarity",
            "geometric_outlier", "rmsd_prefers_other", "non_canonical", "dssr_questionable"
        ],
        help="Filter by specific reason",
    )
    parser.add_argument(
        "--pdb",
        type=str,
        help="Filter by specific PDB ID (e.g., 1EHZ)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=50,
        help="Max outliers to list (default: 50)",
    )
    parser.add_argument(
        "--render",
        type=int,
        help="Render outlier by index number from list",
    )
    parser.add_argument(
        "--render-pair",
        nargs=3,
        metavar=("PDB", "RES1", "RES2"),
        help="Render specific pair: PDB_ID RES_ID1 RES_ID2",
    )
    # Default paths relative to project root
    project_root = Path(__file__).parent.parent.parent  # find_pair_2/

    parser.add_argument(
        "--pdb-dir",
        type=Path,
        default=project_root / "data" / "pdb",
        help="PDB directory",
    )
    parser.add_argument(
        "--template-dir",
        type=Path,
        default=project_root / "basepair-idealized",
        help="Template directory",
    )
    parser.add_argument(
        "--hbond-dir",
        type=Path,
        default=project_root / "data" / "json" / "slot_hbonds",
        help="H-bond JSON directory",
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=Path("viz_output"),
        help="Output directory for visualizations",
    )
    parser.add_argument(
        "--launch",
        action="store_true",
        help="Launch PyMOL after generating",
    )
    parser.add_argument(
        "--include-questionable",
        action="store_true",
        help="Include dssr_questionable pairs (excluded by default as they're DSSR errors)",
    )

    args = parser.parse_args()

    # Handle direct pair rendering
    if args.render_pair:
        pdb_id, res1, res2 = args.render_pair
        viz = PairVisualizer(
            pdb_dir=args.pdb_dir,
            template_dir=args.template_dir,
            hbond_dir=args.hbond_dir,
            output_dir=args.output_dir,
        )

        print(f"Rendering {pdb_id} {res1} - {res2}...")
        result = viz.visualize_pair(pdb_id, res1, res2)

        if result:
            print(f"\nGenerated: {result.pymol_script}")
            print(f"Best match: {result.best_lw} (RMSD: {result.best_rmsd:.3f} Å)")
            for lw, rmsd in sorted(result.rmsd_values.items(), key=lambda x: x[1]):
                marker = " <--" if lw == result.best_lw else ""
                print(f"  {lw}: {rmsd:.3f} Å{marker}")

            if args.launch:
                import subprocess
                subprocess.run(["pymol", str(result.pymol_script)])
        else:
            print("Failed to generate visualization")
            return 1
        return 0

    # Load outliers
    if not args.results.exists():
        print(f"Results directory not found: {args.results}")
        return 1

    outliers = load_outliers(args.results, args.reason, args.pdb, args.include_questionable)

    if not outliers:
        print("No outliers found")
        return 0

    filters = []
    if args.pdb:
        filters.append(f"PDB '{args.pdb}'")
    if args.reason:
        filters.append(f"reason '{args.reason}'")
    filter_str = " with " + " and ".join(filters) if filters else ""
    print(f"\nFound {len(outliers)} outliers{filter_str}")

    # Render specific outlier
    if args.render is not None:
        if args.render < 0 or args.render >= len(outliers):
            print(f"Invalid index {args.render}. Valid range: 0-{len(outliers)-1}")
            return 1

        o = outliers[args.render]
        viz = PairVisualizer(
            pdb_dir=args.pdb_dir,
            template_dir=args.template_dir,
            hbond_dir=args.hbond_dir,
            output_dir=args.output_dir,
        )

        print(f"\nRendering #{args.render}: {o['pdb_id']} {o['res_id1']} - {o['res_id2']} ({o['sequence']})")
        print(f"Reasons: {', '.join(o['reasons'])}")

        result = viz.visualize_pair(o['pdb_id'], o['res_id1'], o['res_id2'])

        if result:
            print(f"\nGenerated: {result.pymol_script}")
            print(f"Best match: {result.best_lw} (RMSD: {result.best_rmsd:.3f} Å)")
            for lw, rmsd in sorted(result.rmsd_values.items(), key=lambda x: x[1]):
                marker = " <--" if lw == result.best_lw else ""
                print(f"  {lw}: {rmsd:.3f} Å{marker}")

            if args.launch:
                import subprocess
                subprocess.run(["pymol", str(result.pymol_script)])
        else:
            print("Failed to generate visualization")
            return 1
    else:
        # Just list outliers
        list_outliers(outliers, args.limit)
        print(f"\nTo render: python visualize_outlier.py --results {args.results} --render <index>")

    return 0


if __name__ == "__main__":
    sys.exit(main())

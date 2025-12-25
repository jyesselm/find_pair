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


def load_outliers(results_dir: Path, reason_filter: Optional[str] = None, pdb_filter: Optional[str] = None) -> List[Dict]:
    """Load outliers from analysis results.

    Args:
        results_dir: Directory containing individual PDB reports
        reason_filter: Optional filter for specific reason code
        pdb_filter: Optional filter for specific PDB ID

    Returns:
        List of outlier dicts with pdb_id, res_id1, res_id2, sequence, reasons
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

            # Apply filter if specified
            if reason_filter and reason_filter not in reasons:
                continue

            outliers.append({
                "pdb_id": pdb_id,
                "res_id1": fn["res_id1"],
                "res_id2": fn["res_id2"],
                "sequence": fn["sequence"],
                "reasons": reasons,
                "our_prediction": fn.get("our_prediction", "unknown"),
                "rmsd_cww": fn.get("geometric_diagnostics", {}).get("rmsd_cww", 0),
                "rmsd_best": fn.get("geometric_diagnostics", {}).get("rmsd_best", 0),
                "best_lw": fn.get("geometric_diagnostics", {}).get("best_lw", "unknown"),
            })

    return outliers


def list_outliers(outliers: List[Dict], limit: int = 50):
    """Print list of outliers."""
    print(f"\n{'#':<4} {'PDB':<6} {'Pair':<20} {'Seq':<4} {'Reasons':<40} {'Best LW':<8} {'RMSD'}")
    print("-" * 100)

    for i, o in enumerate(outliers[:limit]):
        pair = f"{o['res_id1']} - {o['res_id2']}"
        reasons = ", ".join(o['reasons'][:3])
        if len(o['reasons']) > 3:
            reasons += f" (+{len(o['reasons'])-3})"
        rmsd = f"{o['rmsd_cww']:.2f}" if o['rmsd_cww'] else "N/A"

        print(f"{i:<4} {o['pdb_id']:<6} {pair:<20} {o['sequence']:<4} {reasons:<40} {o['best_lw']:<8} {rmsd}")

    if len(outliers) > limit:
        print(f"\n... and {len(outliers) - limit} more. Use --limit to show more.")


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
        default=Path(__file__).parent / "analysis_results" / "cww_analysis_100",
        help="Results directory from cww_annotate.py",
    )
    parser.add_argument(
        "--reason",
        type=str,
        choices=[
            "geometric_outlier", "rmsd_prefers_other", "distance_issues",
            "wrong_atoms", "extra_hbonds", "missing_hbonds", "no_hbonds", "non_canonical"
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

    outliers = load_outliers(args.results, args.reason, args.pdb)

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

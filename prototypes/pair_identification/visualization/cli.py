#!/usr/bin/env python3
"""Command-line interface for pair visualization."""

import argparse
import sys
from pathlib import Path

from .pair_visualizer import PairVisualizer


def get_issues_table(pdb_id: str, pdb_dir: Path, hbond_dir: Path, dssr_dir: Path,
                     idealized_dir: Path, exemplar_dir: Path):
    """Get table of cWW classification issues for a PDB.

    Returns list of dicts with res_id1, res_id2, sequence, rmsd, best_lw, reasons.
    """
    try:
        from ..cww_miss_annotator.annotator import MissAnnotator

        annotator = MissAnnotator(
            pdb_dir=pdb_dir,
            hbond_dir=hbond_dir,
            dssr_dir=dssr_dir,
            idealized_dir=idealized_dir,
            exemplar_dir=exemplar_dir,
        )

        report = annotator.annotate_pdb(pdb_id)
        if not report:
            return []

        issues = []
        for ann in report.fn_annotations:
            rmsd = ann.geometric_diagnostics.rmsd_cww if ann.geometric_diagnostics else 0.0
            best_lw = ann.geometric_diagnostics.best_lw if ann.geometric_diagnostics else "?"
            bp_score = ann.bp_score if hasattr(ann, 'bp_score') else 0.0
            issues.append({
                'res_id1': ann.res_id1,
                'res_id2': ann.res_id2,
                'sequence': ann.sequence,
                'rmsd': rmsd,
                'bp_score': bp_score,
                'best_lw': best_lw,
                'reasons': ann.reasons,
            })
        return issues
    except Exception as e:
        print(f"Warning: Could not load issues: {e}", file=sys.stderr)
        return []


def print_issues_table(pdb_id: str, issues: list):
    """Print formatted table of issues."""
    print(f"\n{'='*110}")
    print(f"cWW Classification Issues for {pdb_id} ({len(issues)} false negatives)")
    print(f"{'='*110}")
    print(f"{'#':>3} {'Res1':<14} {'Res2':<14} {'Seq':4} {'RMSD':>6} {'Score':>6} {'Best':8} {'Reasons'}")
    print(f"{'-'*110}")
    for i, issue in enumerate(issues, 1):
        reasons_str = ", ".join(issue['reasons'][:3])
        if len(issue['reasons']) > 3:
            reasons_str += "..."
        print(f"{i:3} {issue['res_id1']:<14} {issue['res_id2']:<14} {issue['sequence']:4} "
              f"{issue['rmsd']:6.3f} {issue['bp_score']:6.2f} {issue['best_lw']:8} {reasons_str}")
    print()


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Generate PyMOL visualization for base pairs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Show table of issues and pick one to visualize
  %(prog)s --pdb 1UN6 --table

  # Visualize a specific pair
  %(prog)s --pdb 1EHZ --res1 A-G-1 --res2 A-C-72

  # Visualize with PyMOL auto-launch
  %(prog)s --pdb 1EHZ --res1 A-G-1 --res2 A-C-72 --launch-pymol
        """,
    )

    parser.add_argument("--pdb", required=True, help="PDB ID (e.g., 1EHZ)")
    parser.add_argument(
        "--res1", help="First residue ID (e.g., A-G-1). Optional with --table"
    )
    parser.add_argument(
        "--res2", help="Second residue ID (e.g., A-C-72). Optional with --table"
    )
    parser.add_argument(
        "--table", "-t", action="store_true",
        help="Show table of cWW issues and select which to visualize"
    )
    parser.add_argument(
        "--all", "-a", action="store_true",
        help="With --table, visualize ALL issues (generates multiple files)"
    )

    parser.add_argument(
        "--lw",
        default="cWW,tWW,cWS,tHS",
        help="Comma-separated LW classes to show (default: cWW,tWW,cWS,tHS)",
    )

    parser.add_argument(
        "--pdb-dir",
        type=Path,
        default=Path("data/pdb"),
        help="Directory with PDB files (default: data/pdb)",
    )

    parser.add_argument(
        "--template-dir",
        type=Path,
        default=Path("basepair-idealized"),
        help="Directory with idealized templates (default: basepair-idealized)",
    )

    parser.add_argument(
        "--hbond-dir",
        type=Path,
        default=Path("data/json/slot_hbonds"),
        help="Directory with H-bond JSON (default: data/json/slot_hbonds)",
    )

    parser.add_argument(
        "--dssr-dir",
        type=Path,
        default=Path("data/json_dssr"),
        help="Directory with DSSR JSON (default: data/json_dssr)",
    )

    parser.add_argument(
        "--idealized-dir",
        type=Path,
        default=Path("basepair-idealized"),
        help="Directory with idealized templates (default: basepair-idealized)",
    )

    parser.add_argument(
        "--exemplar-dir",
        type=Path,
        default=Path("basepair-exemplars"),
        help="Directory with exemplar templates (default: basepair-exemplars)",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("viz_output"),
        help="Output directory (default: viz_output)",
    )

    parser.add_argument(
        "--no-hbonds",
        action="store_true",
        help="Don't load or display H-bonds",
    )

    parser.add_argument(
        "--launch-pymol",
        action="store_true",
        help="Launch PyMOL with the generated script",
    )

    args = parser.parse_args()

    # Handle --table mode
    if args.table:
        issues = get_issues_table(
            args.pdb, args.pdb_dir, args.hbond_dir, args.dssr_dir,
            args.idealized_dir, args.exemplar_dir
        )

        if not issues:
            print(f"No cWW classification issues found for {args.pdb}")
            return 0

        print_issues_table(args.pdb, issues)

        if args.all:
            # Visualize all issues
            pairs_to_viz = [(issue['res_id1'], issue['res_id2']) for issue in issues]
        elif args.res1 and args.res2:
            # User specified a pair
            pairs_to_viz = [(args.res1, args.res2)]
        else:
            # Interactive selection
            try:
                selection = input("Enter number to visualize (or 'q' to quit, 'a' for all): ").strip()
                if selection.lower() == 'q':
                    return 0
                elif selection.lower() == 'a':
                    pairs_to_viz = [(issue['res_id1'], issue['res_id2']) for issue in issues]
                else:
                    idx = int(selection) - 1
                    if 0 <= idx < len(issues):
                        pairs_to_viz = [(issues[idx]['res_id1'], issues[idx]['res_id2'])]
                    else:
                        print(f"Invalid selection: {selection}")
                        return 1
            except (ValueError, EOFError):
                print("No selection made")
                return 0
    else:
        # Direct visualization mode - require res1 and res2
        if not args.res1 or not args.res2:
            print("ERROR: --res1 and --res2 required (or use --table to see issues)", file=sys.stderr)
            return 1
        pairs_to_viz = [(args.res1, args.res2)]

    # Parse LW classes
    lw_classes = [lw.strip() for lw in args.lw.split(",")]

    # Initialize visualizer
    visualizer = PairVisualizer(
        pdb_dir=args.pdb_dir,
        template_dir=args.template_dir,
        hbond_dir=args.hbond_dir,
        output_dir=args.output_dir,
    )

    # Generate visualizations
    results = []
    for res1, res2 in pairs_to_viz:
        print(f"\nGenerating visualization for {args.pdb} {res1}-{res2}...")

        result = visualizer.visualize_pair(
            pdb_id=args.pdb,
            res_id1=res1,
            res_id2=res2,
            lw_classes=lw_classes,
            include_hbonds=not args.no_hbonds,
        )

        if result is None:
            print(f"ERROR: Failed to generate visualization for {res1}-{res2}", file=sys.stderr)
            continue

        results.append(result)

        # Print results
        print(f"\n{'='*70}")
        print(f"Visualization Complete: {args.pdb} {res1} - {res2}")
        print(f"{'='*70}")
        print(f"Sequence: {result.sequence}")
        print(f"\nGenerated files:")
        print(f"  Target PDB:   {result.target_pdb}")
        print(f"  PyMOL script: {result.pymol_script}")

        if result.template_pdbs:
            print(f"\nTemplates aligned:")
            for lw_class, path in result.template_pdbs.items():
                print(f"  {lw_class}: {path}")

        if result.rmsd_values:
            print(f"\nRMSD Summary:")
            for lw_class, rmsd in sorted(result.rmsd_values.items(), key=lambda x: x[1]):
                marker = " <-- BEST MATCH" if lw_class == result.best_lw else ""
                print(f"  {lw_class:6s}: {rmsd:6.4f} Å{marker}")

        print(f"\nTo visualize:")
        print(f"  pymol {result.pymol_script}")

    # Launch PyMOL if requested (with first result)
    if args.launch_pymol and results:
        import subprocess

        print(f"\nLaunching PyMOL...")
        try:
            subprocess.run(["pymol", str(results[0].pymol_script)])
        except FileNotFoundError:
            print("ERROR: PyMOL not found in PATH", file=sys.stderr)
            return 1

    return 0 if results else 1


if __name__ == "__main__":
    sys.exit(main())

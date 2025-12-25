#!/usr/bin/env python3
"""Command-line interface for pair visualization."""

import argparse
import sys
from pathlib import Path

from .pair_visualizer import PairVisualizer


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Generate PyMOL visualization for base pairs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Visualize a GC pair from 1EHZ
  %(prog)s --pdb 1EHZ --res1 A-G-1 --res2 A-C-72

  # Custom template classes
  %(prog)s --pdb 1EHZ --res1 A-G-1 --res2 A-C-72 --lw cWW,tWW,cWS

  # Custom directories
  %(prog)s --pdb 1EHZ --res1 A-G-1 --res2 A-C-72 \\
           --pdb-dir /path/to/pdbs \\
           --template-dir /path/to/templates
        """,
    )

    parser.add_argument("--pdb", required=True, help="PDB ID (e.g., 1EHZ)")
    parser.add_argument(
        "--res1", required=True, help="First residue ID (e.g., A-G-1)"
    )
    parser.add_argument(
        "--res2", required=True, help="Second residue ID (e.g., A-C-72)"
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

    # Parse LW classes
    lw_classes = [lw.strip() for lw in args.lw.split(",")]

    # Initialize visualizer
    visualizer = PairVisualizer(
        pdb_dir=args.pdb_dir,
        template_dir=args.template_dir,
        hbond_dir=args.hbond_dir,
        output_dir=args.output_dir,
    )

    # Generate visualization
    print(f"Generating visualization for {args.pdb} {args.res1}-{args.res2}...")

    result = visualizer.visualize_pair(
        pdb_id=args.pdb,
        res_id1=args.res1,
        res_id2=args.res2,
        lw_classes=lw_classes,
        include_hbonds=not args.no_hbonds,
    )

    if result is None:
        print("ERROR: Failed to generate visualization", file=sys.stderr)
        return 1

    # Print results
    print("\n" + "=" * 70)
    print(f"Visualization Complete: {args.pdb} {args.res1} - {args.res2}")
    print("=" * 70)
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

    # Launch PyMOL if requested
    if args.launch_pymol:
        import subprocess

        print(f"\nLaunching PyMOL...")
        try:
            subprocess.run(["pymol", str(result.pymol_script)])
        except FileNotFoundError:
            print("ERROR: PyMOL not found in PATH", file=sys.stderr)
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

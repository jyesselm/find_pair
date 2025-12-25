#!/usr/bin/env python3
"""Example usage of PairVisualizer."""

from pathlib import Path
from pair_visualizer import PairVisualizer


def main():
    """Demonstrate PairVisualizer usage."""
    # Initialize visualizer
    visualizer = PairVisualizer(
        pdb_dir=Path("data/pdb"),
        template_dir=Path("basepair-idealized"),
        hbond_dir=Path("data/json/slot_hbonds"),
        output_dir=Path("viz_output"),
    )

    # Example 1: Standard Watson-Crick GC pair
    print("=" * 60)
    print("Example 1: GC Watson-Crick pair from 1EHZ")
    print("=" * 60)

    result = visualizer.visualize_pair(
        pdb_id="1EHZ",
        res_id1="A-G-1",
        res_id2="A-C-72",
        lw_classes=["cWW", "tWW", "cWS", "tHS"],
        include_hbonds=True,
    )

    if result:
        print(f"\nVisualization complete!")
        print(f"  Sequence: {result.sequence}")
        print(f"  Target PDB: {result.target_pdb}")
        print(f"  PyMOL script: {result.pymol_script}")
        print(f"\nRMSD values:")
        for lw_class, rmsd in sorted(result.rmsd_values.items(), key=lambda x: x[1]):
            marker = " <-- BEST" if lw_class == result.best_lw else ""
            print(f"  {lw_class}: {rmsd:.4f} Å{marker}")
        print(f"\nTo view: pymol {result.pymol_script}")

    # Example 2: GU wobble pair
    print("\n" + "=" * 60)
    print("Example 2: GU wobble pair")
    print("=" * 60)

    result2 = visualizer.visualize_pair(
        pdb_id="1EHZ",
        res_id1="A-G-5",
        res_id2="A-U-68",
        lw_classes=["cWW", "tWW"],
        include_hbonds=True,
    )

    if result2:
        print(f"\nVisualization complete!")
        print(f"  Best match: {result2.best_lw} (RMSD={result2.best_rmsd:.4f} Å)")
        print(f"  PyMOL script: {result2.pymol_script}")


if __name__ == "__main__":
    main()

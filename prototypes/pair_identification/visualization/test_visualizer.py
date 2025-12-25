#!/usr/bin/env python3
"""Quick test of PairVisualizer module."""

import sys
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from visualization import PairVisualizer


def test_basic_imports():
    """Test that all imports work."""
    from visualization import (
        PyMOLScriptGenerator,
        HBondViz,
        TemplateViz,
        VisualizationResult,
    )

    print("✓ All imports successful")


def test_visualizer_init():
    """Test PairVisualizer initialization."""
    visualizer = PairVisualizer(
        pdb_dir=Path("data/pdb"),
        template_dir=Path("basepair-idealized"),
        hbond_dir=Path("data/json/slot_hbonds"),
        output_dir=Path("viz_test_output"),
    )

    assert visualizer.pdb_dir == Path("data/pdb")
    assert visualizer.template_dir == Path("basepair-idealized")
    print("✓ PairVisualizer initialization successful")


def test_template_finding():
    """Test template path resolution."""
    visualizer = PairVisualizer()

    # Test standard sequence
    path = visualizer._find_template("GC", "cWW")
    if path:
        print(f"✓ Found cWW GC template: {path}")
    else:
        print("⚠ cWW GC template not found (may be expected)")

    # Test reversed sequence
    path = visualizer._find_template("CG", "cWW")
    if path:
        print(f"✓ Found cWW CG template: {path}")
    else:
        print("⚠ cWW CG template not found (may be expected)")


def test_full_visualization():
    """Test complete visualization pipeline (if PDB exists)."""
    visualizer = PairVisualizer(output_dir=Path("viz_test_output"))

    pdb_path = Path("data/pdb/1EHZ.pdb")
    if not pdb_path.exists():
        print("⚠ Skipping full test: 1EHZ.pdb not found")
        return

    print("\nAttempting full visualization of 1EHZ A-G-1 / A-C-72...")
    result = visualizer.visualize_pair(
        pdb_id="1EHZ",
        res_id1="A-G-1",
        res_id2="A-C-72",
        lw_classes=["cWW"],
        include_hbonds=False,  # Skip H-bonds for quick test
    )

    if result:
        print("✓ Visualization successful!")
        print(f"  Sequence: {result.sequence}")
        print(f"  Target PDB: {result.target_pdb}")
        print(f"  PyMOL script: {result.pymol_script}")
        if result.rmsd_values:
            print(f"  RMSD values: {result.rmsd_values}")
    else:
        print("⚠ Visualization returned None")


def main():
    """Run all tests."""
    print("=" * 60)
    print("Testing PairVisualizer Module")
    print("=" * 60)

    try:
        test_basic_imports()
        test_visualizer_init()
        test_template_finding()
        test_full_visualization()

        print("\n" + "=" * 60)
        print("All tests completed!")
        print("=" * 60)

    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback

        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

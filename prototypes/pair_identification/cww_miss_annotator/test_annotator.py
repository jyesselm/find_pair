#!/usr/bin/env python3
"""Quick test of the annotator module."""

import sys
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from cww_miss_annotator.annotator import MissAnnotator


def test_annotator():
    """Test the annotator on a sample PDB."""
    # Get absolute paths
    base_dir = Path(__file__).parent.parent.parent.parent

    annotator = MissAnnotator(
        pdb_dir=base_dir / "data" / "pdb",
        hbond_dir=base_dir / "data" / "json" / "slot_hbonds",
        dssr_dir=base_dir / "data" / "json_dssr",
        idealized_dir=base_dir / "basepair-idealized",
        exemplar_dir=base_dir / "basepair-exemplars",
    )

    # Test on a small PDB
    pdb_id = "100D"
    print(f"Testing annotator on {pdb_id}...")

    report = annotator.annotate_pdb(pdb_id)

    if report is None:
        print(f"ERROR: No report generated for {pdb_id}")
        return

    print(f"\nResults for {pdb_id}:")
    print(f"  Total canonical cWW: {report.total_canonical_cww}")
    print(f"  True positives: {report.true_positives}")
    print(f"  False negatives: {report.false_negatives}")
    print(f"  False positives: {report.false_positives}")

    if report.false_negatives > 0:
        print(f"\nFalse negative examples:")
        for i, ann in enumerate(report.fn_annotations[:3]):
            print(f"\n  {i+1}. {ann.res_id1} - {ann.res_id2} ({ann.sequence})")
            print(f"     Saenger: {ann.saenger}")
            print(f"     Our prediction: {ann.our_prediction}")
            print(f"     Reasons: {', '.join(ann.reasons)}")
            print(f"     H-bonds found: {len(ann.hbond_diagnostics.found_hbonds)}")
            print(
                f"     H-bonds missing: {len(ann.hbond_diagnostics.missing_hbonds)}"
            )
            print(f"     RMSD cWW: {ann.geometric_diagnostics.rmsd_cww:.3f}")
            print(f"     RMSD best: {ann.geometric_diagnostics.rmsd_best:.3f}")

    print("\nTest complete!")


if __name__ == "__main__":
    test_annotator()

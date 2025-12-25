#!/usr/bin/env python3
"""CLI interface for cWW miss annotation tool.

This tool annotates Watson-Crick classification differences between our system
and DSSR, providing detailed diagnostics for each miss (false negatives and
false positives).
"""

import argparse
import json
import sys
from multiprocessing import Pool
from pathlib import Path
from typing import List, Optional, Tuple

from .annotator import MissAnnotator
from .diagnostics import PDBReport
from .report import AggregateReporter, save_aggregate_report, save_pdb_report


def process_pdb(args: Tuple[str, dict]) -> Tuple[str, Optional[PDBReport]]:
    """Process a single PDB (for parallel execution).

    Args:
        args: Tuple of (pdb_id, config) where config contains directory paths.

    Returns:
        Tuple of (pdb_id, PDBReport or None if failed).
    """
    pdb_id, config = args

    annotator = MissAnnotator(
        pdb_dir=Path(config["pdb_dir"]),
        hbond_dir=Path(config["hbond_dir"]),
        dssr_dir=Path(config["dssr_dir"]),
        idealized_dir=Path(config["idealized_dir"]),
        exemplar_dir=Path(config["exemplar_dir"]),
    )

    try:
        report = annotator.annotate_pdb(pdb_id)
        return pdb_id, report
    except Exception as e:
        print(f"Error processing {pdb_id}: {e}", file=sys.stderr)
        return pdb_id, None


def load_pdb_list(args: argparse.Namespace) -> List[str]:
    """Load list of PDB IDs to process.

    Args:
        args: Command-line arguments.

    Returns:
        List of PDB IDs to process.
    """
    if args.pdb:
        return [args.pdb]

    if args.pdb_list:
        with open(args.pdb_list) as f:
            data = json.load(f)
        if isinstance(data, list):
            pdb_ids = data
        else:
            # Assume it's a dict with PDB IDs in values
            pdb_ids = list(next(iter(data.values())))
    else:
        # Default: find all PDBs with both DSSR and slot_hbonds
        dssr_ids = set(f.stem for f in args.dssr_dir.glob("*.json"))
        hbond_ids = set(f.stem for f in args.hbond_dir.glob("*.json"))
        pdb_ids = sorted(dssr_ids & hbond_ids)

    if args.max_pdbs:
        pdb_ids = pdb_ids[: args.max_pdbs]

    return pdb_ids


def print_summary(reports: List[PDBReport]):
    """Print summary statistics to console.

    Args:
        reports: List of PDBReport objects.
    """
    if not reports:
        print("No reports to summarize")
        return

    reporter = AggregateReporter()
    aggregate = reporter.aggregate(reports)

    print(f"\n{'='*60}")
    print(f"Aggregate Results ({aggregate.total_pdbs} PDBs)")
    print(f"{'='*60}")
    print(f"Total canonical cWW:     {aggregate.total_canonical_cww:,}")
    print(f"True positives:          {aggregate.total_true_positives:,}")
    print(f"False negatives:         {aggregate.total_false_negatives:,}")
    print(f"False positives:         {aggregate.total_false_positives:,}")
    print(f"Accuracy:                {aggregate.accuracy:.2%}")
    print(f"Sensitivity:             {aggregate.sensitivity:.2%}")
    print(f"Precision:               {aggregate.precision:.2%}")

    if aggregate.reason_counts:
        print(f"\nTop miss reasons:")
        sorted_reasons = sorted(
            aggregate.reason_counts.items(), key=lambda x: -x[1]
        )
        for reason, count in sorted_reasons[:10]:
            pct = 100 * count / aggregate.total_false_negatives
            print(f"  {reason:30s}: {count:5d} ({pct:5.1f}%)")

    if aggregate.sequence_stats:
        print(f"\nPer-sequence breakdown:")
        for seq in ["GC", "CG", "AU", "UA"]:
            if seq in aggregate.sequence_stats:
                stats = aggregate.sequence_stats[seq]
                total = stats["total"]
                misses = stats["misses"]
                acc = (total - misses) / total if total > 0 else 0.0
                print(f"  {seq}: {total-misses}/{total} correct ({acc:.1%})")


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Annotate cWW classification differences between our system and DSSR",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze single PDB
  %(prog)s --pdb 1EHZ

  # Analyze PDBs from list
  %(prog)s --pdb-list data/fast_pdbs.json --max-pdbs 100

  # Use all available cores
  %(prog)s --pdb-list data/fast_pdbs.json --workers 10

  # Verbose output
  %(prog)s --pdb 1EHZ -v
        """,
    )

    # Input options
    parser.add_argument("--pdb", type=str, help="Single PDB ID to analyze")
    parser.add_argument(
        "--pdb-list", type=Path, help="JSON file with list of PDB IDs"
    )

    # Directory options
    parser.add_argument(
        "--pdb-dir",
        type=Path,
        default=Path("data/pdb"),
        help="Directory containing PDB files (default: data/pdb)",
    )
    parser.add_argument(
        "--hbond-dir",
        type=Path,
        default=Path("data/json/slot_hbonds"),
        help="Directory containing slot H-bond JSON files (default: data/json/slot_hbonds)",
    )
    parser.add_argument(
        "--dssr-dir",
        type=Path,
        default=Path("data/json_dssr"),
        help="Directory containing DSSR JSON files (default: data/json_dssr)",
    )
    parser.add_argument(
        "--idealized-dir",
        type=Path,
        default=Path("basepair-idealized"),
        help="Directory containing idealized basepair templates (default: basepair-idealized)",
    )
    parser.add_argument(
        "--exemplar-dir",
        type=Path,
        default=Path("basepair-exemplars"),
        help="Directory containing exemplar basepair templates (default: basepair-exemplars)",
    )

    # Output options
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("results/cww_analysis"),
        help="Output directory for reports (default: results/cww_analysis)",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Print detailed output"
    )

    # Processing options
    parser.add_argument(
        "--workers",
        "-w",
        type=int,
        default=10,
        help="Number of parallel workers (default: 10)",
    )
    parser.add_argument(
        "--max-pdbs",
        type=int,
        default=None,
        help="Maximum number of PDBs to process",
    )

    args = parser.parse_args()

    # Load PDB list
    pdb_ids = load_pdb_list(args)

    if not pdb_ids:
        print("No PDBs to process", file=sys.stderr)
        return 1

    print(f"Processing {len(pdb_ids)} PDBs with {args.workers} workers...")

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Build config for workers
    config = {
        "pdb_dir": str(args.pdb_dir),
        "hbond_dir": str(args.hbond_dir),
        "dssr_dir": str(args.dssr_dir),
        "idealized_dir": str(args.idealized_dir),
        "exemplar_dir": str(args.exemplar_dir),
    }

    # Process PDBs
    work_items = [(pdb_id, config) for pdb_id in pdb_ids]

    if args.workers == 1:
        results = [process_pdb(item) for item in work_items]
    else:
        with Pool(args.workers) as pool:
            results = pool.map(process_pdb, work_items)

    # Collect successful reports
    reports = []
    for pdb_id, report in results:
        if report:
            reports.append(report)

            # Save individual report
            report_path = args.output / f"{pdb_id}.json"
            save_pdb_report(report, report_path)

            if args.verbose:
                tp = report.true_positives
                total = report.total_canonical_cww
                fn = report.false_negatives
                fp = report.false_positives
                acc = tp / total if total > 0 else 0.0
                print(
                    f"{pdb_id}: {tp}/{total} TP ({acc:.1%}), "
                    f"{fn} FN, {fp} FP"
                )

    # Generate and save aggregate report
    if reports:
        reporter = AggregateReporter()
        aggregate = reporter.aggregate(reports)

        save_aggregate_report(aggregate, args.output / "aggregate.json")

        # Print summary
        print_summary(reports)

        print(f"\nResults saved to {args.output}/")
        print(f"  - Individual reports: {len(reports)} files")
        print(f"  - Aggregate report: aggregate.json")
    else:
        print("No successful reports generated", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

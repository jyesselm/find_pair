#!/usr/bin/env python3
"""
Aggregate results from all cluster comparison batches.

This script:
1. Loads results from all batch result files
2. Combines statistics
3. Generates a comprehensive report
4. Identifies PDBs with differences
"""

import argparse
import json
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Any

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))


def load_batch_results(results_dir: Path) -> List[Dict[str, Any]]:
    """Load all batch result files."""
    batch_results = []
    
    for result_file in sorted(results_dir.glob("batch_*_results.json")):
        try:
            with open(result_file) as f:
                batch_results.append(json.load(f))
        except Exception as e:
            print(f"Warning: Failed to load {result_file}: {e}", file=sys.stderr)
    
    return batch_results


def aggregate_statistics(batch_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Aggregate statistics from all batches."""
    aggregated = {
        "total_pdbs": 0,
        "matches": 0,
        "differences": 0,
        "errors": 0,
        "total_duration": 0.0,
        "by_status": defaultdict(int),
        "pdb_status": {},
        "batches": {
            "total": len(batch_results),
            "successful": 0,
            "failed": 0,
        }
    }
    
    for batch in batch_results:
        stats = batch.get("stats", {})
        
        aggregated["total_pdbs"] += stats.get("total", 0)
        aggregated["matches"] += stats.get("matches", 0)
        aggregated["differences"] += stats.get("differences", 0)
        aggregated["errors"] += stats.get("errors", 0)
        aggregated["total_duration"] += batch.get("duration_seconds", 0.0)
        
        aggregated["batches"]["successful"] += 1
        
        # Track individual PDB statuses
        individual = batch.get("individual_results", {})
        for pdb_id, pdb_info in individual.items():
            status = pdb_info.get("status", "unknown")
            aggregated["by_status"][status] += 1
            aggregated["pdb_status"][pdb_id] = {
                "status": status,
                "has_differences": pdb_info.get("has_differences", False),
                "batch_id": batch.get("batch_id"),
            }
    
    return aggregated


def generate_summary_report(aggregated: Dict[str, Any], output_file: Path):
    """Generate a summary report."""
    lines = []
    lines.append("=" * 80)
    lines.append("CLUSTER COMPARISON SUMMARY")
    lines.append("=" * 80)
    lines.append("")
    
    # Overall statistics
    lines.append("OVERALL STATISTICS")
    lines.append("-" * 80)
    lines.append(f"Total PDBs processed: {aggregated['total_pdbs']}")
    lines.append(f"Matches: {aggregated['matches']}")
    lines.append(f"Differences: {aggregated['differences']}")
    lines.append(f"Errors: {aggregated['errors']}")
    lines.append(f"Total duration: {aggregated['total_duration']:.1f} seconds ({aggregated['total_duration']/3600:.2f} hours)")
    lines.append("")
    
    # Status breakdown
    lines.append("STATUS BREAKDOWN")
    lines.append("-" * 80)
    for status, count in sorted(aggregated['by_status'].items()):
        lines.append(f"  {status}: {count}")
    lines.append("")
    
    # Batch statistics
    lines.append("BATCH STATISTICS")
    lines.append("-" * 80)
    batches = aggregated['batches']
    lines.append(f"Total batches: {batches['total']}")
    lines.append(f"Successful: {batches['successful']}")
    lines.append(f"Failed: {batches['failed']}")
    lines.append("")
    
    # PDBs with differences
    pdbs_with_diffs = [
        pdb_id
        for pdb_id, info in aggregated['pdb_status'].items()
        if info.get('has_differences', False)
    ]
    
    if pdbs_with_diffs:
        lines.append("PDBS WITH DIFFERENCES")
        lines.append("-" * 80)
        for pdb_id in sorted(pdbs_with_diffs):
            info = aggregated['pdb_status'][pdb_id]
            lines.append(f"  {pdb_id} (batch {info['batch_id']}, status: {info['status']})")
        lines.append("")
    else:
        lines.append("PDBS WITH DIFFERENCES")
        lines.append("-" * 80)
        lines.append("  None - all comparisons matched!")
        lines.append("")
    
    # PDBs with errors
    pdbs_with_errors = [
        pdb_id
        for pdb_id, info in aggregated['pdb_status'].items()
        if info.get('status') == 'error'
    ]
    
    if pdbs_with_errors:
        lines.append("PDBS WITH ERRORS")
        lines.append("-" * 80)
        for pdb_id in sorted(pdbs_with_errors):
            info = aggregated['pdb_status'][pdb_id]
            lines.append(f"  {pdb_id} (batch {info['batch_id']})")
        lines.append("")
    
    with open(output_file, 'w') as f:
        f.write('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate results from cluster comparison jobs"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory containing batch results"
    )
    parser.add_argument(
        "--summary-file",
        type=Path,
        help="Output file for summary report (default: output-dir/summary_report.txt)"
    )
    parser.add_argument(
        "--json-summary",
        type=Path,
        help="Output file for JSON summary (default: output-dir/aggregated_results.json)"
    )
    
    args = parser.parse_args()
    
    results_dir = args.output_dir / "batch_results"
    
    if not results_dir.exists():
        print(f"Error: Results directory not found: {results_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Load batch results
    print(f"Loading batch results from {results_dir}...")
    batch_results = load_batch_results(results_dir)
    
    if not batch_results:
        print("Error: No batch results found!", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loaded {len(batch_results)} batch result files")
    
    # Aggregate statistics
    print("Aggregating statistics...")
    aggregated = aggregate_statistics(batch_results)
    
    # Generate summary report
    summary_file = args.summary_file or (args.output_dir / "summary_report.txt")
    print(f"Generating summary report: {summary_file}")
    generate_summary_report(aggregated, summary_file)
    
    # Save JSON summary
    json_summary = args.json_summary or (args.output_dir / "aggregated_results.json")
    print(f"Saving JSON summary: {json_summary}")
    with open(json_summary, 'w') as f:
        json.dump(aggregated, f, indent=2)
    
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total PDBs: {aggregated['total_pdbs']}")
    print(f"Matches: {aggregated['matches']}")
    print(f"Differences: {aggregated['differences']}")
    print(f"Errors: {aggregated['errors']}")
    print(f"Total time: {aggregated['total_duration']/3600:.2f} hours")
    print(f"\n📁 Summary report: {summary_file}")
    print(f"📁 JSON summary: {json_summary}")


if __name__ == "__main__":
    main()


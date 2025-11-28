#!/usr/bin/env python3
"""
Run comparison for a single batch of PDBs.

This script is called by cluster job scripts to process one batch.
It reads PDB IDs from a JSON file and runs the comparison.
"""

import argparse
import json
import sys
from pathlib import Path
from datetime import datetime

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from scripts.compare_json import (
    run_comparison,
    create_comparator,
    analyze_results,
    generate_report,
)


def main():
    parser = argparse.ArgumentParser(
        description="Run comparison for a single batch of PDBs"
    )
    parser.add_argument(
        "--batch-file",
        type=Path,
        required=True,
        help="JSON file containing list of PDB IDs"
    )
    parser.add_argument(
        "--batch-id",
        type=int,
        required=True,
        help="Batch ID number"
    )
    parser.add_argument(
        "--comparison-type",
        type=str,
        choices=["compare", "atoms", "frames", "steps"],
        default="compare",
        help="Type of comparison to run"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for results"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads to use"
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Disable result caching"
    )
    
    args = parser.parse_args()
    
    # Load PDB IDs
    with open(args.batch_file) as f:
        pdb_ids = json.load(f)
    
    print(f"Batch {args.batch_id:04d}: Processing {len(pdb_ids)} PDB files")
    print(f"  PDBs: {', '.join(pdb_ids[:5])}{'...' if len(pdb_ids) > 5 else ''}")
    
    # Create comparator based on comparison type
    config_map = {
        "compare": {"compare_atoms": True, "compare_frames": True, "compare_steps": True},
        "atoms": {"compare_atoms": True, "compare_frames": False, "compare_steps": False},
        "frames": {"compare_atoms": False, "compare_frames": True, "compare_steps": False},
        "steps": {"compare_atoms": False, "compare_frames": True, "compare_steps": True},
    }
    
    config = config_map.get(args.comparison_type, config_map["compare"])
    comparator = create_comparator(None, project_root)
    comparator.compare_atoms = config.get("compare_atoms", False)
    comparator.compare_frames = config.get("compare_frames", False)
    comparator.compare_steps = config.get("compare_steps", False)
    comparator.enable_cache = not args.no_cache
    
    # Run comparison
    start_time = datetime.now()
    try:
        results = run_comparison(
            pdb_ids=pdb_ids,
            project_root=project_root,
            use_legacy_mode=False,
            comparator=comparator,
            max_workers=args.threads,
            regenerate=False,
        )
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        # Analyze results
        stats = analyze_results(results)
        
        # Generate report
        report = generate_report(
            stats=stats,
            results=results,
            verbose=False,
            diff_only=True,
        )
        
        # Save results
        result_file = args.output_dir / "batch_results" / f"batch_{args.batch_id:04d}_results.json"
        report_file = args.output_dir / "batch_results" / f"batch_{args.batch_id:04d}_report.txt"
        
        # Save summary JSON
        summary = {
            "batch_id": args.batch_id,
            "pdb_ids": pdb_ids,
            "start_time": start_time.isoformat(),
            "end_time": end_time.isoformat(),
            "duration_seconds": duration,
            "stats": {
                "total": stats.get("total", 0),
                "matches": stats.get("matches", 0),
                "differences": stats.get("diff_count", 0),
                "errors": stats.get("errors", 0),
            },
            "individual_results": {
                pdb_id: {
                    "status": result.status,
                    "has_differences": result.has_differences(),
                }
                for pdb_id, result in results.items()
            }
        }
        
        with open(result_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        with open(report_file, 'w') as f:
            f.write(report)
        
        print(f"✅ Batch {args.batch_id:04d} completed in {duration:.1f}s")
        print(f"   Results: {stats.get('matches', 0)} matches, {stats.get('diff_count', 0)} differences")
        print(f"   Saved to: {result_file}")
        
    except Exception as e:
        error_file = args.output_dir / "batch_results" / f"batch_{args.batch_id:04d}_error.txt"
        with open(error_file, 'w') as f:
            f.write(f"Error processing batch {args.batch_id:04d}\n")
            f.write(f"Time: {datetime.now().isoformat()}\n")
            f.write(f"Error: {str(e)}\n")
            import traceback
            f.write(f"\nTraceback:\n{traceback.format_exc()}\n")
        
        print(f"❌ Batch {args.batch_id:04d} failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()


#!/usr/bin/env python3
"""
Cluster comparison runner - Splits PDBs into batches and runs comparisons.

This script:
1. Loads PDB IDs from test sets or all available PDBs
2. Splits them into batches for parallel cluster jobs
3. Creates job scripts for each batch
4. Submits jobs to the cluster
5. Collects results when jobs complete

Usage:
    python3 scripts/cluster/run_cluster_comparison.py \
        --test-set 100 \
        --jobs-per-node 4 \
        --comparison-type compare \
        --output-dir cluster_results
"""

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional
import subprocess

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from scripts.compare_json import (
    load_test_set,
    find_available_pdbs,
    get_test_set_path,
)


def load_pdb_list(project_root: Path, test_set: Optional[int], legacy_mode: bool) -> List[str]:
    """Load list of PDB IDs to compare."""
    if test_set:
        test_pdb_ids = load_test_set(project_root, test_set)
        if test_pdb_ids is None:
            print(f"Error: Test set of size {test_set} not found.", file=sys.stderr)
            sys.exit(1)
        return test_pdb_ids
    else:
        pdb_ids = find_available_pdbs(project_root, legacy_mode)
        if not pdb_ids:
            print("Error: No PDB files found!", file=sys.stderr)
            sys.exit(1)
        return pdb_ids


def split_into_batches(items: List[str], batch_size: int) -> List[List[str]]:
    """Split list into batches of specified size."""
    batches = []
    for i in range(0, len(items), batch_size):
        batches.append(items[i:i + batch_size])
    return batches


def create_job_script(
    batch_id: int,
    pdb_ids: List[str],
    project_root: Path,
    comparison_type: str,
    output_dir: Path,
    jobs_per_node: int,
    partition: str,
    time_limit: str,
    memory: str,
) -> Path:
    """Create a SLURM job script for a batch of PDBs."""
    job_script = output_dir / "job_scripts" / f"batch_{batch_id:04d}.sh"
    job_script.parent.mkdir(parents=True, exist_ok=True)
    
    # Write PDB list to JSON file
    pdb_list_file = output_dir / "batch_lists" / f"batch_{batch_id:04d}.json"
    pdb_list_file.parent.mkdir(parents=True, exist_ok=True)
    with open(pdb_list_file, 'w') as f:
        json.dump(pdb_ids, f, indent=2)
    
    # Create job script
    script_content = f"""#!/bin/bash
#SBATCH --job-name=comp_batch_{batch_id:04d}
#SBATCH --output={output_dir}/job_logs/batch_{batch_id:04d}.out
#SBATCH --error={output_dir}/job_logs/batch_{batch_id:04d}.err
#SBATCH --partition={partition}
#SBATCH --time={time_limit}
#SBATCH --mem={memory}
#SBATCH --cpus-per-task={jobs_per_node}
#SBATCH --nodes=1
#SBATCH --ntasks=1

# Job information
echo "Job started at: $(date)"
echo "Batch ID: {batch_id:04d}"
echo "PDB IDs in batch: {len(pdb_ids)}"
echo "Comparison type: {comparison_type}"
echo "Working directory: $(pwd)"

# Load environment (modify as needed for your cluster)
# module load python/3.9
# module load gcc/11.2.0

# Change to project directory
cd "{project_root.absolute()}"

# Activate virtual environment if used (uncomment if needed)
# source venv/bin/activate

# Run comparison for this batch
python3 scripts/cluster/run_single_batch.py \\
    --batch-file "{pdb_list_file.absolute()}" \\
    --batch-id {batch_id:04d} \\
    --comparison-type {comparison_type} \\
    --output-dir "{output_dir.absolute()}" \\
    --threads {jobs_per_node}

echo "Job finished at: $(date)"
"""
    
    with open(job_script, 'w') as f:
        f.write(script_content)
    
    job_script.chmod(0o755)
    return job_script


def submit_jobs(job_scripts: List[Path], dry_run: bool = False) -> List[str]:
    """Submit all job scripts to the cluster."""
    job_ids = []
    
    for job_script in job_scripts:
        if dry_run:
            print(f"[DRY RUN] Would submit: {job_script}")
            job_ids.append("DRY_RUN")
        else:
            try:
                result = subprocess.run(
                    ["sbatch", str(job_script)],
                    capture_output=True,
                    text=True,
                    check=True
                )
                # Parse job ID from output: "Submitted batch job 12345"
                job_id = result.stdout.strip().split()[-1]
                job_ids.append(job_id)
                print(f"Submitted batch {job_script.stem}: Job ID {job_id}")
            except subprocess.CalledProcessError as e:
                print(f"Error submitting {job_script}: {e}", file=sys.stderr)
                print(f"  stderr: {e.stderr}", file=sys.stderr)
                job_ids.append(None)
            except FileNotFoundError:
                print("Error: 'sbatch' command not found. Are you on a cluster with SLURM?", file=sys.stderr)
                sys.exit(1)
    
    return job_ids


def main():
    parser = argparse.ArgumentParser(
        description="Run JSON comparison jobs on a cluster",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare test set of 100 PDBs, 10 PDBs per job
  python3 scripts/cluster/run_cluster_comparison.py \\
      --test-set 100 \\
      --batch-size 10 \\
      --comparison-type compare

  # Compare all available PDBs, 5 PDBs per job
  python3 scripts/cluster/run_cluster_comparison.py \\
      --batch-size 5 \\
      --comparison-type steps

  # Dry run to see what would be submitted
  python3 scripts/cluster/run_cluster_comparison.py \\
      --test-set 50 \\
      --batch-size 5 \\
      --dry-run
"""
    )
    
    parser.add_argument(
        "--test-set",
        type=int,
        choices=[10, 50, 100, 500, 1000],
        help="Use a test set of the specified size (10, 50, 100, 500, or 1000)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10,
        help="Number of PDBs per batch (default: 10)"
    )
    parser.add_argument(
        "--comparison-type",
        type=str,
        choices=["compare", "atoms", "frames", "steps"],
        default="compare",
        help="Type of comparison to run (default: compare)"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("cluster_results"),
        help="Output directory for results (default: cluster_results)"
    )
    parser.add_argument(
        "--jobs-per-node",
        type=int,
        default=4,
        help="Number of parallel jobs per node (default: 4)"
    )
    parser.add_argument(
        "--partition",
        type=str,
        default="normal",
        help="SLURM partition name (default: normal)"
    )
    parser.add_argument(
        "--time-limit",
        type=str,
        default="2:00:00",
        help="Time limit per job (default: 2:00:00)"
    )
    parser.add_argument(
        "--memory",
        type=str,
        default="8G",
        help="Memory per job (default: 8G)"
    )
    parser.add_argument(
        "--legacy-mode",
        action="store_true",
        help="Use legacy mode JSON files"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Don't submit jobs, just create scripts"
    )
    
    args = parser.parse_args()
    
    # Setup directories
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "job_scripts").mkdir(exist_ok=True)
    (args.output_dir / "job_logs").mkdir(exist_ok=True)
    (args.output_dir / "batch_lists").mkdir(exist_ok=True)
    (args.output_dir / "batch_results").mkdir(exist_ok=True)
    
    # Load PDB list
    print(f"Loading PDB list...")
    pdb_ids = load_pdb_list(project_root, args.test_set, args.legacy_mode)
    print(f"Found {len(pdb_ids)} PDB files to compare")
    
    # Split into batches
    batches = split_into_batches(pdb_ids, args.batch_size)
    print(f"Split into {len(batches)} batches of ~{args.batch_size} PDBs each")
    
    # Create job scripts
    print(f"Creating job scripts...")
    job_scripts = []
    for batch_id, batch_pdb_ids in enumerate(batches, start=1):
        job_script = create_job_script(
            batch_id=batch_id,
            pdb_ids=batch_pdb_ids,
            project_root=project_root,
            comparison_type=args.comparison_type,
            output_dir=args.output_dir,
            jobs_per_node=args.jobs_per_node,
            partition=args.partition,
            time_limit=args.time_limit,
            memory=args.memory,
        )
        job_scripts.append(job_script)
        print(f"  Created batch_{batch_id:04d}.sh ({len(batch_pdb_ids)} PDBs)")
    
    # Save job metadata
    metadata = {
        "total_pdbs": len(pdb_ids),
        "batch_size": args.batch_size,
        "num_batches": len(batches),
        "comparison_type": args.comparison_type,
        "test_set": args.test_set,
        "legacy_mode": args.legacy_mode,
        "job_ids": [],
    }
    
    # Submit jobs
    if args.dry_run:
        print(f"\n[DRY RUN] Would submit {len(job_scripts)} jobs")
        print("Use --no-dry-run to actually submit jobs")
    else:
        print(f"\nSubmitting {len(job_scripts)} jobs...")
        job_ids = submit_jobs(job_scripts, dry_run=False)
        metadata["job_ids"] = job_ids
        
        # Save metadata
        metadata_file = args.output_dir / "job_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        print(f"\n✅ Submitted {len([j for j in job_ids if j])} jobs successfully")
        print(f"📁 Results will be in: {args.output_dir}")
        print(f"📋 Metadata saved to: {metadata_file}")
        print(f"\nMonitor jobs with: squeue -u $USER")
        print(f"After jobs complete, run:")
        print(f"  python3 scripts/cluster/aggregate_results.py --output-dir {args.output_dir}")


if __name__ == "__main__":
    main()


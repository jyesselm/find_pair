# Quick Start: Running Comparisons on SLURM Cluster

## Step-by-Step Guide

### 0. Download PDB Files (if needed)

Before running comparisons, ensure PDB files are available:

```bash
# Download PDBs from test set
python3 scripts/download_pdbs.py --test-set 100

# Or download specific PDBs
python3 scripts/download_pdbs.py 1H4S 2BNA 3DNA

# Skip existing files (resume interrupted downloads)
python3 scripts/download_pdbs.py --test-set 1000 --skip-existing
```

### 1. Generate Test Sets (if needed)

```bash
python3 scripts/compare_json.py generate-test-sets
```

This creates test sets of 10, 50, 100, 500, and 1000 PDBs in `resources/test_sets/`.

### 2. Submit Jobs to Cluster

```bash
python3 scripts/cluster/run_cluster_comparison.py \
    --test-set 100 \
    --batch-size 10 \
    --comparison-type compare \
    --output-dir cluster_results
```

**Options explained:**
- `--test-set 100` - Use the 100-PDB test set
- `--batch-size 10` - 10 PDBs per batch (adjust based on PDB size/complexity)
- `--comparison-type compare` - Full comparison (atoms, frames, steps)
- `--output-dir cluster_results` - Where results will be saved

### 3. Monitor Jobs

```bash
# Check job status
squeue -u $USER

# Monitor progress
bash scripts/cluster/monitor_jobs.sh cluster_results
```

### 4. Wait for Jobs to Complete

Jobs will run in parallel on cluster nodes. Check periodically:

```bash
squeue -u $USER  # Should show no jobs when done
```

### 5. Aggregate Results

Once all jobs complete, combine the results:

```bash
python3 scripts/cluster/aggregate_results.py \
    --output-dir cluster_results
```

This creates:
- `cluster_results/summary_report.txt` - Human-readable summary
- `cluster_results/aggregated_results.json` - JSON summary with statistics

## Common Examples

### Small test (10 PDBs, 5 per batch)
```bash
python3 scripts/cluster/run_cluster_comparison.py \
    --test-set 10 \
    --batch-size 5 \
    --comparison-type compare
```

### Large comparison (1000 PDBs, 25 per batch)
```bash
python3 scripts/cluster/run_cluster_comparison.py \
    --test-set 1000 \
    --batch-size 25 \
    --comparison-type compare \
    --time-limit 4:00:00 \
    --memory 16G
```

### Steps-only comparison
```bash
python3 scripts/cluster/run_cluster_comparison.py \
    --test-set 500 \
    --batch-size 20 \
    --comparison-type steps \
    --output-dir steps_results
```

### Test first (dry run)
```bash
python3 scripts/cluster/run_cluster_comparison.py \
    --test-set 100 \
    --batch-size 10 \
    --dry-run
```

## Adjusting for Your Cluster

You may need to modify cluster settings in the script or job templates:

1. **Partition name**: Change `--partition normal` to your cluster's partition
2. **Module loading**: Uncomment/modify module load lines in generated job scripts
3. **Virtual environment**: Uncomment venv activation if using one

## Troubleshooting

### Jobs not submitting
- Check if `sbatch` is available: `which sbatch`
- Check partition name: `sinfo`
- Use `--dry-run` to test first

### Jobs failing
- Check logs: `cat cluster_results/job_logs/batch_*.err`
- Verify JSON files exist: `ls data/json/`
- Check if test set exists: `cat resources/test_sets/test_set_100.json`

### Out of memory
- Increase memory: `--memory 32G`

### Timeouts
- Increase time limit: `--time-limit 8:00:00`
- Reduce batch size: `--batch-size 5`

## Canceling Jobs

```bash
# Cancel all your jobs
scancel -u $USER

# Cancel specific job
scancel JOBID  # Get JOBID from squeue output
```

## Next Steps

See [README.md](README.md) for detailed documentation.


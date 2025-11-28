# Cluster Comparison Scripts

Scripts for running large-scale JSON comparisons on cluster computing systems.

## Overview

These scripts allow you to:
1. Split PDB comparisons into batches
2. Submit jobs to a cluster (SLURM or PBS/Torque)
3. Run comparisons in parallel across multiple nodes
4. Aggregate results from all batches

## Quick Start

```bash
# 1. Submit jobs to cluster
python3 scripts/cluster/run_cluster_comparison.py \
    --test-set 100 \
    --batch-size 10 \
    --comparison-type compare \
    --output-dir cluster_results

# 2. Wait for jobs to complete (monitor with squeue)
squeue -u $USER

# 3. Aggregate results
python3 scripts/cluster/aggregate_results.py \
    --output-dir cluster_results
```

## Scripts

### 1. `run_cluster_comparison.py`
**Main script** that splits PDBs into batches and submits cluster jobs.

**Usage:**
```bash
python3 scripts/cluster/run_cluster_comparison.py [OPTIONS]
```

**Options:**
- `--test-set SIZE` - Use test set (10, 50, 100, 500, 1000)
- `--batch-size N` - PDBs per batch (default: 10)
- `--comparison-type TYPE` - compare/atoms/frames/steps (default: compare)
- `--output-dir DIR` - Output directory (default: cluster_results)
- `--jobs-per-node N` - Parallel jobs per node (default: 4)
- `--partition NAME` - SLURM partition (default: normal)
- `--time-limit TIME` - Time limit per job (default: 2:00:00)
- `--memory MEM` - Memory per job (default: 8G)
- `--legacy-mode` - Use legacy mode JSON files
- `--dry-run` - Create scripts without submitting

**Example:**
```bash
python3 scripts/cluster/run_cluster_comparison.py \
    --test-set 500 \
    --batch-size 20 \
    --comparison-type steps \
    --partition gpu \
    --time-limit 4:00:00 \
    --memory 16G
```

### 2. `run_single_batch.py`
**Worker script** that runs comparisons for one batch. Called automatically by cluster jobs.

### 3. `aggregate_results.py`
**Aggregation script** that combines results from all batches.

**Usage:**
```bash
python3 scripts/cluster/aggregate_results.py \
    --output-dir cluster_results
```

**Output:**
- `summary_report.txt` - Human-readable summary
- `aggregated_results.json` - JSON summary with all statistics

### 4. `monitor_jobs.sh`
Helper script to monitor job status and check for failures.

**Usage:**
```bash
bash scripts/cluster/monitor_jobs.sh cluster_results
```

## Cluster System: SLURM

These scripts are designed for **SLURM** cluster systems. All job submission uses `sbatch` commands.

**SLURM commands:**
- `sbatch` - Submit jobs
- `squeue` - Check job status
- `scancel` - Cancel jobs
- `sacct` - View job accounting


## Workflow

1. **Prepare**: Ensure you have test sets generated
   ```bash
   python3 scripts/compare_json.py generate-test-sets
   ```

2. **Submit**: Run the main script to submit jobs
   ```bash
   python3 scripts/cluster/run_cluster_comparison.py \
       --test-set 100 \
       --batch-size 10
   ```

3. **Monitor**: Check job status
   ```bash
   squeue -u $USER
   bash scripts/cluster/monitor_jobs.sh cluster_results
   ```

4. **Wait**: Jobs run in parallel on cluster nodes

5. **Aggregate**: Combine results when jobs complete
   ```bash
   python3 scripts/cluster/aggregate_results.py \
       --output-dir cluster_results
   ```

## Output Structure

```
cluster_results/
├── job_metadata.json          # Job submission metadata
├── summary_report.txt          # Aggregated summary (after aggregation)
├── aggregated_results.json     # JSON summary (after aggregation)
├── job_scripts/               # SLURM job scripts
│   ├── batch_0001.sh
│   ├── batch_0002.sh
│   └── ...
├── job_logs/                  # Job output logs
│   ├── batch_0001.out
│   ├── batch_0001.err
│   └── ...
├── batch_lists/               # PDB lists per batch
│   ├── batch_0001.json
│   ├── batch_0002.json
│   └── ...
└── batch_results/             # Individual batch results
    ├── batch_0001_results.json
    ├── batch_0001_report.txt
    └── ...
```

## Configuration

### Adjusting Batch Size

Larger batches = fewer jobs but longer runtime per job:
```bash
--batch-size 50  # 50 PDBs per batch
```

Smaller batches = more jobs but shorter runtime per job:
```bash
--batch-size 5   # 5 PDBs per batch
```

### Memory Requirements

Adjust based on PDB file sizes:
```bash
--memory 16G     # For larger PDBs
--memory 8G      # For smaller PDBs (default)
```

### Time Limits

Adjust based on batch size and PDB complexity:
```bash
--time-limit 4:00:00  # 4 hours
--time-limit 1:00:00  # 1 hour (default: 2:00:00)
```

## Examples

### Compare 1000 PDBs
```bash
python3 scripts/cluster/run_cluster_comparison.py \
    --test-set 1000 \
    --batch-size 25 \
    --comparison-type compare \
    --partition normal \
    --time-limit 4:00:00 \
    --memory 16G
```

### Quick test with 10 PDBs
```bash
python3 scripts/cluster/run_cluster_comparison.py \
    --test-set 10 \
    --batch-size 5 \
    --comparison-type frames \
    --dry-run  # Test first
```

### Steps only comparison
```bash
python3 scripts/cluster/run_cluster_comparison.py \
    --test-set 500 \
    --batch-size 20 \
    --comparison-type steps \
    --output-dir steps_results
```

## Troubleshooting

### Jobs failing
1. Check job logs: `cat cluster_results/job_logs/batch_*.err`
2. Check if JSON files exist: `ls data/json/`
3. Check if test set exists: `cat resources/test_sets/test_set_100.json`

### Out of memory
Increase memory: `--memory 32G`

### Timeouts
Increase time limit: `--time-limit 8:00:00`

### Job submission errors
- Check if SLURM is available: `which sbatch`
- Check partition name: `sinfo`
- Use `--dry-run` to test script generation

## Canceling Jobs

```bash
# Cancel all jobs
scancel -u $USER

# Cancel specific batch
scancel JOBID  # Get JOBID from squeue
```

## Performance Tips

1. **Batch size**: Start with 10-20 PDBs per batch
2. **Parallel jobs**: Match `--jobs-per-node` to available CPUs
3. **Memory**: Larger PDBs may need more memory
4. **Partition**: Use GPU partition if available for faster processing

## Related Documentation

- [TESTING_GUIDE.md](../../docs/TESTING_GUIDE.md) - Comparison workflows
- [COMPARISON_WORKFLOW.md](../../docs/COMPARISON_WORKFLOW.md) - Frame-first approach


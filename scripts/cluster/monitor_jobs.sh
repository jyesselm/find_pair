#!/bin/bash
# Monitor cluster comparison jobs

OUTPUT_DIR="${1:-cluster_results}"

if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Error: Output directory not found: $OUTPUT_DIR" >&2
    exit 1
fi

echo "=== Cluster Comparison Job Monitor ==="
echo "Output directory: $OUTPUT_DIR"
echo ""

# Check if metadata exists
if [ -f "$OUTPUT_DIR/job_metadata.json" ]; then
    echo "Job metadata found"
    python3 -c "
import json
import sys
with open('$OUTPUT_DIR/job_metadata.json') as f:
    meta = json.load(f)
    print(f\"  Total PDBs: {meta.get('total_pdbs', 0)}\")
    print(f\"  Batches: {meta.get('num_batches', 0)}\")
    print(f\"  Comparison type: {meta.get('comparison_type', 'unknown')}\")
"
    echo ""
fi

# Count result files
RESULTS_DIR="$OUTPUT_DIR/batch_results"
if [ -d "$RESULTS_DIR" ]; then
    TOTAL=$(find "$RESULTS_DIR" -name "batch_*_results.json" | wc -l | tr -d ' ')
    ERRORS=$(find "$RESULTS_DIR" -name "batch_*_error.txt" | wc -l | tr -d ' ')
    echo "Batch results:"
    echo "  Completed: $TOTAL"
    echo "  Errors: $ERRORS"
    echo ""
fi

# Check SLURM jobs
if command -v squeue &> /dev/null; then
    echo "Active SLURM jobs:"
    USER=$(whoami)
    JOBS=$(squeue -u "$USER" -h | wc -l | tr -d ' ')
    if [ "$JOBS" -gt 0 ]; then
        squeue -u "$USER" --format="%.18i %.9P %.20j %.8u %.8T %.10M %.9l %.6D %R"
    else
        echo "  No active jobs"
    fi
    echo ""
fi

# Check for recent errors
if [ -d "$OUTPUT_DIR/job_logs" ]; then
    RECENT_ERRORS=$(find "$OUTPUT_DIR/job_logs" -name "*.err" -mmin -60 -size +0 | head -5)
    if [ -n "$RECENT_ERRORS" ]; then
        echo "Recent errors (last 60 minutes):"
        while IFS= read -r err_file; do
            echo "  $(basename $err_file):"
            tail -3 "$err_file" | sed 's/^/    /'
        done <<< "$RECENT_ERRORS"
        echo ""
    fi
fi

# Progress estimate
if [ -f "$OUTPUT_DIR/job_metadata.json" ] && [ -d "$RESULTS_DIR" ]; then
    python3 -c "
import json
from pathlib import Path
meta_file = Path('$OUTPUT_DIR/job_metadata.json')
results_dir = Path('$RESULTS_DIR')
if meta_file.exists() and results_dir.exists():
    with open(meta_file) as f:
        meta = json.load(f)
    total_batches = meta.get('num_batches', 0)
    completed = len(list(results_dir.glob('batch_*_results.json')))
    if total_batches > 0:
        progress = (completed / total_batches) * 100
        print(f\"Progress: {completed}/{total_batches} batches ({progress:.1f}%)\")
"
fi


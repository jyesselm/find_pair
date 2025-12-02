#!/bin/bash
#
# Run full index validation on all PDBs
#
# This script will:
# - Automatically resume from where it left off (skips already passed PDBs)
# - Test all untested, failed, skipped, or timed-out PDBs
# - Use 8 parallel threads
# - Process in batches of 100
# - Clean up generated files for successful validations
# - Stop on first failure for investigation
#
# Usage:
#   ./RUN_FULL_VALIDATION.sh
#
# To revalidate ALL PDBs (including those that passed):
#   ./RUN_FULL_VALIDATION.sh --revalidate
#

set -e  # Exit on error

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Check if build exists
if [ ! -f "build/generate_modern_json" ]; then
    echo "Error: build/generate_modern_json not found"
    echo "Please build the project first:"
    echo "  cd build && cmake .. && make"
    exit 1
fi

# Run validation
echo "Starting index validation..."
echo "This will test all PDBs that haven't passed validation yet."
echo ""

python3 scripts/validate_all_indices.py \
    --threads 8 \
    --batch-size 100 \
    "$@"

echo ""
echo "Validation complete! See data/index_validation_status.csv for results."


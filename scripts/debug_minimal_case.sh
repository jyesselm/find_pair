#!/bin/bash
# Comprehensive debugging workflow for minimal test cases
# Automates the entire process: extract, generate, compare

set -e

PDB_ID="${1:-1AQ4}"
MINIMAL_DIR="data/pdb/minimal"
LEGACY_JSON_DIR="data/json_legacy"
MODERN_JSON_DIR="data/json"
NUM_PAIRS="${2:-2}"

if [ -z "$PDB_ID" ]; then
    echo "Usage: $0 <PDB_ID> [num_pairs]"
    echo "Example: $0 1AQ4 2"
    exit 1
fi

PDB_FILE="data/pdb/${PDB_ID}.pdb"
LEGACY_SELECTION_JSON="${LEGACY_JSON_DIR}/find_bestpair_selection/${PDB_ID}.json"

echo "=========================================="
echo "Step-by-Step Debugging Workflow"
echo "=========================================="
echo "PDB ID: $PDB_ID"
echo "Minimal pairs: $NUM_PAIRS"
echo ""

# Step 1: Check if PDB exists
if [ ! -f "$PDB_FILE" ]; then
    echo "ERROR: PDB file not found: $PDB_FILE"
    exit 1
fi

if [ ! -f "$LEGACY_SELECTION_JSON" ]; then
    echo "ERROR: Legacy JSON not found: $LEGACY_SELECTION_JSON"
    exit 1
fi

# Step 2: Extract minimal test case
echo "Step 1: Extracting minimal test case..."
python3 scripts/extract_minimal_pairs.py \
    --pdb "$PDB_FILE" \
    --legacy-json "$LEGACY_SELECTION_JSON" \
    --output-dir "$MINIMAL_DIR" \
    --num-pairs "$NUM_PAIRS" \
    --max-fragments 1

MINIMAL_PDB="${MINIMAL_DIR}/${PDB_ID}_minimal_pairs_1_${NUM_PAIRS}.pdb"
if [ ! -f "$MINIMAL_PDB" ]; then
    echo "ERROR: Minimal PDB not created: $MINIMAL_PDB"
    exit 1
fi

MINIMAL_PDB_ID=$(basename "$MINIMAL_PDB" .pdb)
echo "Created: $MINIMAL_PDB_ID"
echo ""

# Step 3: Generate legacy JSON
echo "Step 2: Generating legacy JSON..."
rm -f "${LEGACY_JSON_DIR}/best_partner_candidates/${MINIMAL_PDB_ID}.json"
rm -f "${LEGACY_JSON_DIR}/mutual_best_decisions/${MINIMAL_PDB_ID}.json"
rm -f "${LEGACY_JSON_DIR}/iteration_states/${MINIMAL_PDB_ID}.json"

./org/build/bin/find_pair_original "$MINIMAL_PDB" > /dev/null 2>&1

echo "Legacy JSON generated"
echo ""

# Step 4: Generate modern JSON
echo "Step 3: Generating modern JSON..."
MODERN_JSON_DIR_FULL="${MODERN_JSON_DIR}/${MINIMAL_PDB_ID}.json"
rm -rf "$MODERN_JSON_DIR_FULL"

./build/generate_modern_json "$MINIMAL_PDB" "$MODERN_JSON_DIR_FULL" --fix-indices > /dev/null 2>&1

echo "Modern JSON generated"
echo ""

# Step 5: Compare results
echo "Step 4: Comparing results..."
echo ""

# Find residues to check
LEGACY_CANDIDATES="${LEGACY_JSON_DIR}/best_partner_candidates/${MINIMAL_PDB_ID}.json"
if [ -f "$LEGACY_CANDIDATES" ]; then
    # Extract first residue that has candidates
    FIRST_RES=$(python3 -c "
import json
with open('$LEGACY_CANDIDATES') as f:
    data = json.load(f)
    if data and len(data) > 0:
        print(data[0]['res_i'])
    else:
        print('1')
")
else
    FIRST_RES=1
fi

echo "--- Best Partner Comparison (residue $FIRST_RES) ---"
python3 scripts/compare_best_partner.py "$MINIMAL_PDB_ID" "$FIRST_RES" \
    --legacy-dir "$LEGACY_JSON_DIR" \
    --modern-dir "$MODERN_JSON_DIR_FULL" \
    --verbose 2>&1 | head -30

echo ""
echo "--- Iteration States Comparison ---"
python3 scripts/compare_iteration.py "$MINIMAL_PDB_ID" \
    --legacy-dir "$LEGACY_JSON_DIR" \
    --modern-dir "$MODERN_JSON_DIR_FULL" \
    --verbose 2>&1

echo ""
echo "--- Mutual Best Decisions Comparison ---"
python3 scripts/compare_mutual_best.py "$MINIMAL_PDB_ID" \
    --legacy-dir "$LEGACY_JSON_DIR" \
    --modern-dir "$MODERN_JSON_DIR_FULL" \
    --verbose 2>&1

echo ""
echo "=========================================="
echo "Debugging complete!"
echo "=========================================="
echo ""
echo "Minimal PDB: $MINIMAL_PDB"
echo "Legacy JSON: $LEGACY_JSON_DIR"
echo "Modern JSON: $MODERN_JSON_DIR_FULL"
echo ""
echo "To investigate further:"
echo "  python3 scripts/compare_best_partner.py $MINIMAL_PDB_ID <res_i> --verbose"
echo "  python3 scripts/compare_iteration.py $MINIMAL_PDB_ID --verbose"
echo "  python3 scripts/compare_mutual_best.py $MINIMAL_PDB_ID --verbose"


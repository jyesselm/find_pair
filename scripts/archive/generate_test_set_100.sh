#!/bin/bash
# Generate modern JSON for all PDBs in test_set_100

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

# Load PDB list from test_set_100.json
PDBS=$(python3 -c "
import json
with open('data/test_sets/test_set_100.json') as f:
    data = json.load(f)
    pdbs = data.get('pdb_ids', data.get('pdbs', []))
    print(' '.join(pdbs))
")

PDB_DIR="data/pdb"
OUTPUT_DIR="data/json"

TOTAL=$(echo "$PDBS" | wc -w | tr -d ' ')
COUNT=0
SUCCESS=0
FAILED=0
SKIPPED=0

echo "Generating modern JSON for test_set_100 ($TOTAL PDBs)..."
echo ""

for PDB in $PDBS; do
    COUNT=$((COUNT + 1))
    PDB_FILE="$PDB_DIR/${PDB}.pdb"
    JSON_FILE="$OUTPUT_DIR/${PDB}.json"
    
    if [ ! -f "$PDB_FILE" ]; then
        SKIPPED=$((SKIPPED + 1))
        if [ $((COUNT % 20)) -eq 0 ]; then
            echo "[$COUNT/$TOTAL] $PDB: SKIP (PDB file not found)"
        fi
        continue
    fi
    
    if ./build/generate_modern_json "$PDB_FILE" "$JSON_FILE" > /dev/null 2>&1; then
        SUCCESS=$((SUCCESS + 1))
        if [ $((COUNT % 10)) -eq 0 ]; then
            echo "[$COUNT/$TOTAL] $PDB: OK (Progress: $SUCCESS success, $FAILED failed, $SKIPPED skipped)"
        fi
    else
        FAILED=$((FAILED + 1))
        echo "[$COUNT/$TOTAL] $PDB: FAILED"
    fi
done

echo ""
echo "Summary:"
echo "  Total: $TOTAL"
echo "  Success: $SUCCESS"
echo "  Failed: $FAILED"
echo "  Skipped: $SKIPPED"


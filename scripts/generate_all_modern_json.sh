#!/bin/bash
# Generate modern JSON for all PDBs that have legacy base_pair files

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

# Find all PDBs with legacy base_pair files
LEGACY_DIR="data/json_legacy"
PDB_DIR="data/pdb"
OUTPUT_DIR="data/json"

# Get list of PDBs
PDBS=$(find "$LEGACY_DIR" -name "*_base_pair.json" | sed 's|.*/||' | sed 's|_base_pair.json||' | sort)

TOTAL=$(echo "$PDBS" | wc -l | tr -d ' ')
COUNT=0
SUCCESS=0
FAILED=0

echo "Generating modern JSON for $TOTAL PDBs..."
echo ""

for PDB in $PDBS; do
    COUNT=$((COUNT + 1))
    PDB_FILE="$PDB_DIR/${PDB}.pdb"
    JSON_FILE="$OUTPUT_DIR/${PDB}.json"
    
    if [ ! -f "$PDB_FILE" ]; then
        echo "[$COUNT/$TOTAL] $PDB: SKIP (PDB file not found)"
        continue
    fi
    
    if ./build/generate_modern_json "$PDB_FILE" "$JSON_FILE" > /dev/null 2>&1; then
        SUCCESS=$((SUCCESS + 1))
        if [ $((COUNT % 10)) -eq 0 ]; then
            echo "[$COUNT/$TOTAL] $PDB: OK (Progress: $SUCCESS success, $FAILED failed)"
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


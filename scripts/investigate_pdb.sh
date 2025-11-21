#!/bin/bash
# Helper script to investigate a single PDB file
# Usage: ./scripts/investigate_pdb.sh <pdb_name>

if [ $# -ne 1 ]; then
    echo "Usage: $0 <pdb_name>"
    echo "Example: $0 1AQ3"
    exit 1
fi

PDB_NAME=$1
PDB_FILE="data/pdb/${PDB_NAME}.pdb"
GEN_JSON="data/json/${PDB_NAME}.json"
LEG_JSON="data/json_legacy/${PDB_NAME}.json"

echo "=" | tr -d '\n' | head -c 70; echo
echo "Investigating: $PDB_NAME"
echo "=" | tr -d '\n' | head -c 70; echo
echo ""

# Check if files exist
if [ ! -f "$PDB_FILE" ]; then
    echo "ERROR: PDB file not found: $PDB_FILE"
    exit 1
fi

if [ ! -f "$LEG_JSON" ]; then
    echo "ERROR: Legacy JSON not found: $LEG_JSON"
    exit 1
fi

# Generate JSON if needed
if [ ! -f "$GEN_JSON" ]; then
    echo "Generating JSON..."
    ./build/tests/integration/test_single_pdb "$PDB_NAME" > /dev/null 2>&1
fi

# Run comparison
echo "Running detailed comparison..."
./build/tests/integration/test_single_pdb "$PDB_NAME"

echo ""
echo "=" | tr -d '\n' | head -c 70; echo
echo "PDB File Analysis"
echo "=" | tr -d '\n' | head -c 70; echo
echo ""

# Check for MODEL/ENDMDL
if grep -q "^MODEL" "$PDB_FILE"; then
    echo "⚠️  WARNING: PDB file contains MODEL/ENDMDL records"
    echo "   Model count: $(grep -c "^MODEL" "$PDB_FILE")"
fi

# Check for alternate conformations
ALT_COUNT=$(grep -E "^ATOM|^HETATM" "$PDB_FILE" | awk '{print substr($0, 17, 1)}' | grep -v " " | sort -u | wc -l | tr -d ' ')
if [ "$ALT_COUNT" -gt 1 ]; then
    echo "⚠️  WARNING: Multiple alternate conformations found"
    echo "   Unique alt_loc values: $(grep -E "^ATOM|^HETATM" "$PDB_FILE" | awk '{print substr($0, 17, 1)}' | grep -v " " | sort -u | tr '\n' ',' | sed 's/,$//')"
fi

# Check occupancy distribution
echo ""
echo "Occupancy distribution:"
grep -E "^ATOM|^HETATM" "$PDB_FILE" | awk '{print substr($0, 55, 6)}' | sort | uniq -c | sort -rn | head -5

# Check chain IDs
echo ""
echo "Chain IDs found:"
grep -E "^ATOM|^HETATM" "$PDB_FILE" | awk '{print substr($0, 22, 1)}' | sort -u | tr '\n' ' '
echo ""


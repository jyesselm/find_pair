#!/bin/bash
# Simple comparison script for legacy vs modern

PDB=$1
if [ -z "$PDB" ]; then
    echo "Usage: $0 <pdb_id>"
    exit 1
fi

PDB_FILE="data/pdb/${PDB}.pdb"

echo "=========================================="
echo "Comparing $PDB"
echo "=========================================="
echo ""

# Step 1: Run find_pair
echo "1. Running find_pair..."
echo "   Legacy:"
org/build/bin/find_pair_original "$PDB_FILE" "${PDB}_legacy.inp" 2>&1 | tail -2
echo "   Modern:"
./build/find_pair_app "$PDB_FILE" "${PDB}_modern.inp" 2>&1 | tail -2

echo ""
echo "2. Base pair comparison:"
echo "   Legacy base pairs:"
grep -E "^[[:space:]]*[0-9]+[[:space:]]+[0-9]+[[:space:]]+[0-9]+" "${PDB}_legacy.inp" | awk '{printf "  %d-%d\n", $1, $2}' | head -10
echo "   Modern base pairs:"
grep -E "^[[:space:]]*[0-9]+[[:space:]]+[0-9]+[[:space:]]+[0-9]+" "${PDB}_modern.inp" | awk '{printf "  %d-%d\n", $2, $3}' | head -10

echo ""
echo "3. Running analyze..."
cd /tmp
cp "$OLDPWD/${PDB}_legacy.inp" .
cp "$OLDPWD/${PDB}_modern.inp" .

echo "   Legacy analyze:"
cd "$OLDPWD/org" && ../org/build/bin/analyze_original "../${PDB}_legacy.inp" 2>&1 | tail -2
cd "$OLDPWD"
if [ -f "bp_step.par" ]; then
    echo "   Step params found: $(wc -l < bp_step.par) lines"
    head -5 bp_step.par
fi

echo "   Modern analyze:"
./build/analyze_app "${PDB}_modern.inp" 2>&1 | grep -A 8 "Step Parameters" | head -10


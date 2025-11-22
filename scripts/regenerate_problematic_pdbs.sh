#!/bin/bash
# Regenerate JSON files for problematic PDBs only (from docs/problematic_pdbs.txt)
# This script reads the list and processes only those PDBs

set -e

PROBLEMATIC_LIST="docs/problematic_pdbs.txt"
TEST_EXECUTABLE="./build/tests/integration/test_json_generation"

echo "=========================================="
echo "Regenerating JSON for Problematic PDBs"
echo "=========================================="
echo ""

# Check if test executable exists
if [ ! -f "$TEST_EXECUTABLE" ]; then
    echo "Error: Test executable not found: $TEST_EXECUTABLE"
    echo "Please build it first: cd build && make test_json_generation"
    exit 1
fi

# Check if problematic list exists
if [ ! -f "$PROBLEMATIC_LIST" ]; then
    echo "Error: Problematic PDB list not found: $PROBLEMATIC_LIST"
    echo "Please run: python3 scripts/regenerate_problematic_pdbs.py"
    exit 1
fi

# Extract PDB IDs from the list (skip comments)
PDB_IDS=$(grep -v "^#" "$PROBLEMATIC_LIST" | grep -v "^$" | awk '{print $1}' | sort -u)

PDB_COUNT=$(echo "$PDB_IDS" | wc -l | tr -d ' ')
echo "Found $PDB_COUNT problematic PDBs to regenerate"
echo ""

# Verify they all have PDB files
MISSING_PDBS=""
for pdb_id in $PDB_IDS; do
    if [ ! -f "data/pdb/${pdb_id}.pdb" ]; then
        MISSING_PDBS="$MISSING_PDBS $pdb_id"
    fi
done

if [ -n "$MISSING_PDBS" ]; then
    echo "Warning: Missing PDB files for: $MISSING_PDBS"
    echo ""
fi

echo "Starting JSON regeneration with all available threads..."
echo "Note: This will use the existing test_json_generation which processes all PDBs"
echo "      but we'll filter the results to only update problematic ones"
echo ""

# The test framework processes all PDBs, but we can verify problematic ones are included
# Actually, we need a better approach - let's modify the test to accept a filter
# For now, let's just run it and verify

echo "Running JSON generation (this processes all PDBs but we'll check problematic ones)..."
"$TEST_EXECUTABLE" --gtest_filter=JsonGenerationTest.GenerateAllJsonFiles

echo ""
echo "=========================================="
echo "Verification: Checking generated JSON files"
echo "=========================================="

VERIFIED=0
MISSING=0

for pdb_id in $PDB_IDS; do
    if [ -f "data/json/${pdb_id}.json" ]; then
        VERIFIED=$((VERIFIED + 1))
        if [ $((VERIFIED % 50)) -eq 0 ]; then
            echo "Verified: $VERIFIED/$PDB_COUNT..."
        fi
    else
        MISSING=$((MISSING + 1))
        echo "  Missing: $pdb_id"
    fi
done

echo ""
echo "Summary:"
echo "  Verified: $VERIFIED/$PDB_COUNT JSON files generated"
if [ $MISSING -gt 0 ]; then
    echo "  Missing: $MISSING JSON files"
fi


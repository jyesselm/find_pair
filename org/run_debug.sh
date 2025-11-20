#!/bin/bash
# Script to run original x3dna with debug output

PDB_FILE="$1"
if [ -z "$PDB_FILE" ]; then
    echo "Usage: $0 <pdb_file>"
    exit 1
fi

BUILD_DIR="build"
if [ ! -d "$BUILD_DIR" ]; then
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    cmake ..
    make
    cd ..
fi

echo "Running find_pair_original with debug output..."
echo "PDB file: $PDB_FILE"
echo ""

# Run find_pair and capture both stdout and stderr
"$BUILD_DIR/bin/find_pair_original" "$PDB_FILE" > "${PDB_FILE%.pdb}.inp" 2> "${PDB_FILE%.pdb}_find_pair_debug.txt"

echo "find_pair output saved to: ${PDB_FILE%.pdb}.inp"
echo "find_pair debug output saved to: ${PDB_FILE%.pdb}_find_pair_debug.txt"

if [ -f "${PDB_FILE%.pdb}.inp" ]; then
    echo ""
    echo "Running analyze_original with debug output..."
    "$BUILD_DIR/bin/analyze_original" "${PDB_FILE%.pdb}.inp" > "${PDB_FILE%.pdb}_analyze_output.txt" 2> "${PDB_FILE%.pdb}_analyze_debug.txt"
    
    echo "analyze output saved to: ${PDB_FILE%.pdb}_analyze_output.txt"
    echo "analyze debug output saved to: ${PDB_FILE%.pdb}_analyze_debug.txt"
fi

echo ""
echo "Debug files created:"
ls -lh "${PDB_FILE%.pdb}"*debug*.txt 2>/dev/null || true


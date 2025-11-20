#!/bin/bash
# Quick test script to run the original code and capture debug output

PDB_FILE="${1:-legacy_archive/fixtures/pdb/1EHZ.pdb}"

if [ ! -f "$PDB_FILE" ]; then
    echo "Error: PDB file not found: $PDB_FILE"
    echo "Usage: $0 [pdb_file]"
    exit 1
fi

echo "Testing original x3dna with debug output..."
echo "PDB file: $PDB_FILE"
echo ""

# Check if built
if [ ! -f "build/bin/find_pair_original" ]; then
    echo "Building first..."
    mkdir -p build && cd build
    cmake .. && make
    cd ..
fi

# Run find_pair
echo "Running find_pair_original..."
./build/bin/find_pair_original "$PDB_FILE" > test_output.inp 2> test_debug.txt

echo ""
echo "Results:"
echo "  Standard output: test_output.inp ($(wc -l < test_output.inp) lines)"
echo "  Debug output: test_debug.txt ($(wc -l < test_debug.txt) lines)"
echo ""

# Show key debug info
echo "Key debug information:"
echo "----------------------"
grep -A 5 "\[DEBUG calculate_more_bppars\]" test_debug.txt | head -20
echo ""
grep "\[DEBUG\] Result:" test_debug.txt | head -5


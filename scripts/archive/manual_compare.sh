#!/bin/bash
# Manual comparison - run from project root

PDB=$1
[ -z "$PDB" ] && PDB="6V9Q"

echo "=== Comparing $PDB ==="
echo ""

# Legacy find_pair
echo "Legacy find_pair:"
org/build/bin/find_pair_original "data/pdb/${PDB}.pdb" "${PDB}_legacy.inp" 2>&1 | grep -E "(base pairs|Time used)" 
echo ""

# Modern find_pair  
echo "Modern find_pair:"
./build/find_pair_app "data/pdb/${PDB}.pdb" "${PDB}_modern.inp" 2>&1 | grep -E "(base pairs|Output file)"
echo ""

# Show base pairs
echo "Base Pairs:"
echo "  Legacy: $(grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9]+[[:space:]]+[0-9]+' ${PDB}_legacy.inp | wc -l) pairs"
echo "  Modern: $(grep -E '^[[:space:]]*[0-9]+' ${PDB}_modern.inp | tail -n +5 | wc -l) pairs"
echo ""

# Run analyze if .inp files exist
if [ -f "${PDB}_legacy.inp" ] && [ -f "${PDB}_modern.inp" ]; then
    echo "Running analyze..."
    echo ""
    
    # Legacy analyze (from project root)
    echo "Legacy analyze:"
    cd org && ../org/build/bin/analyze_original "../${PDB}_legacy.inp" 2>&1 | tail -3 && cd ..
    [ -f "bp_step.par" ] && echo "  Created bp_step.par ($(wc -l < bp_step.par) lines)" || echo "  No bp_step.par found"
    
    echo ""
    echo "Modern analyze:"
    ./build/analyze_app "${PDB}_modern.inp" 2>&1 | grep -A 3 "Calculated"
fi


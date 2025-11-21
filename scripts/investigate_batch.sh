#!/bin/bash
# Investigate multiple PDB files from a list
# Usage: ./scripts/investigate_batch.sh <file_with_pdb_names>

if [ $# -ne 1 ]; then
    echo "Usage: $0 <file_with_pdb_names>"
    echo "Example: $0 docs/representative_pdbs.txt"
    exit 1
fi

INPUT_FILE=$1

if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

echo "Investigating PDB files from: $INPUT_FILE"
echo "=" | tr -d '\n' | head -c 70; echo
echo ""

COUNT=0
TOTAL=$(wc -l < "$INPUT_FILE")

while IFS= read -r pdb_name || [ -n "$pdb_name" ]; do
    COUNT=$((COUNT + 1))
    echo "[$COUNT/$TOTAL] Investigating: $pdb_name"
    echo ""
    
    ./scripts/investigate_pdb.sh "$pdb_name" 2>&1 | tail -20
    
    echo ""
    echo "---"
    echo ""
    
    # Pause every 5 files to review
    if [ $((COUNT % 5)) -eq 0 ]; then
        echo "Paused after $COUNT files. Press Enter to continue..."
        read
    fi
done < "$INPUT_FILE"

echo "Investigation complete!"


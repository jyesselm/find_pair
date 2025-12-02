#!/bin/bash
# Regenerate modern JSON for mismatched PDBs with --fix-indices

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

# PDBs with confirmed index mismatches
PDBS=("1TTT" "9CF3" "1TN1" "1TN2" "3F2T" "5V0O")

echo "=========================================="
echo "Regenerating Modern JSON with --fix-indices"
echo "=========================================="
echo ""

for pdb in "${PDBS[@]}"; do
    echo "Processing $pdb..."
    
    pdb_file="data/pdb/${pdb}.pdb"
    output_file="/tmp/${pdb}.inp"
    legacy_json="data/json_legacy/base_frame_calc/${pdb}.json"
    
    # Check if PDB file exists
    if [ ! -f "$pdb_file" ]; then
        echo "  ⚠️  PDB file not found: $pdb_file"
        continue
    fi
    
    # Check if legacy JSON exists
    if [ ! -f "$legacy_json" ]; then
        echo "  ⚠️  Legacy JSON not found: $legacy_json"
        echo "  Attempting auto-detect..."
    fi
    
    # Regenerate with --fix-indices
    echo "  Running: ./build/find_pair_app --fix-indices $pdb_file $output_file"
    if ./build/find_pair_app --fix-indices "$pdb_file" "$output_file" 2>&1 | tee "/tmp/${pdb}_regenerate.log"; then
        echo "  ✅ $pdb regenerated successfully"
        
        # Verify indices match
        echo "  Verifying residue indices..."
        python3 -c "
import json
import sys

pdb_id = '$pdb'

try:
    with open(f'data/json_legacy/base_frame_calc/{pdb_id}.json') as f:
        legacy = json.load(f)
    
    with open(f'data/json/base_frame_calc/{pdb_id}.json') as f:
        modern = json.load(f)
    
    mismatches = []
    for i in [1, 10, 25, 50, 100]:
        legacy_res = next((r for r in legacy if r.get('residue_idx') == i), None)
        modern_res = next((r for r in modern if r.get('residue_idx') == i or r.get('legacy_residue_idx') == i), None)
        
        if legacy_res and modern_res:
            legacy_name = legacy_res.get('residue_name', '?')
            modern_name = modern_res.get('residue_name', '?')
            legacy_seq = legacy_res.get('residue_seq', '?')
            modern_seq = modern_res.get('residue_seq', '?')
            
            if legacy_name != modern_name or legacy_seq != modern_seq:
                mismatches.append(f'Residue {i}: Legacy={legacy_name} seq{legacy_seq}, Modern={modern_name} seq{modern_seq}')
    
    if mismatches:
        print('  ❌ Index mismatches found:')
        for m in mismatches:
            print(f'    {m}')
        sys.exit(1)
    else:
        print('  ✅ Residue indices match')
        sys.exit(0)
except Exception as e:
    print(f'  ⚠️  Verification error: {e}')
    sys.exit(1)
" || echo "  ⚠️  Verification failed (check manually)"
    else
        echo "  ❌ $pdb regeneration failed (check /tmp/${pdb}_regenerate.log)"
    fi
    
    echo ""
done

echo "=========================================="
echo "Regeneration Complete"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Verify residue indices match for all PDBs"
echo "2. Re-run comparisons: python3 scripts/compare_json.py compare <PDB_ID> --verbose"
echo "3. Check if pair differences are resolved"


# Action Plan: Regenerate Modern JSON with --fix-indices

**Date**: 2025-01-XX  
**Status**: ⚠️ **URGENT** - Multiple PDBs have residue index mismatches

---

## Summary

**Critical Finding**: Modern JSON files for mismatched PDBs were **NOT generated with `--fix-indices`**, causing residue indices to not match legacy. This means we're comparing **different pairs**, not real validation differences!

---

## Affected PDBs

### Confirmed Index Mismatches

1. **1TTT** ⚠️ **CRITICAL**
   - Residue 162: Legacy = 2MG, Modern = C (DIFFERENT!)
   - Residue 177: Legacy = C, Modern = M2G (DIFFERENT!)
   - **Status**: Needs regeneration

2. **9CF3** ⚠️
   - Residue 10: Legacy = G seq13, Modern = A seq12 (MISMATCH)
   - Residue 50: Legacy = G seq53, Modern = A seq52 (MISMATCH)
   - **Status**: Needs regeneration

3. **1TN1** ⚠️
   - Residue 1: Legacy = G seq1, Modern = C seq2 (MISMATCH)
   - Residue 10: Legacy = 2MG seq10, Modern = C seq11 (MISMATCH)
   - Residue 50: Legacy = U seq50, Modern = G seq51 (MISMATCH)
   - **Status**: Needs regeneration

4. **1TN2** ⚠️
   - Same mismatches as 1TN1
   - **Status**: Needs regeneration

5. **3F2T** ⚠️
   - Residue 10: Legacy = G seq10, Modern = G seq9 (MISMATCH - seq off by 1)
   - Residue 50: Legacy = A seq50, Modern = A seq49 (MISMATCH - seq off by 1)
   - **Status**: Needs regeneration

6. **5V0O** ✅
   - Sample residues match (1, 10, 25)
   - **Status**: May be OK, but should verify

---

## Regeneration Commands

### For Each Affected PDB

```bash
# 1TTT
./build/find_pair_app --fix-indices data/pdb/1TTT.pdb /tmp/1TTT.inp

# 9CF3
./build/find_pair_app --fix-indices data/pdb/9CF3.pdb /tmp/9CF3.inp

# 1TN1
./build/find_pair_app --fix-indices data/pdb/1TN1.pdb /tmp/1TN1.inp

# 1TN2
./build/find_pair_app --fix-indices data/pdb/1TN2.pdb /tmp/1TN2.inp

# 3F2T
./build/find_pair_app --fix-indices data/pdb/3F2T.pdb /tmp/3F2T.inp

# 5V0O (verify)
./build/find_pair_app --fix-indices data/pdb/5V0O.pdb /tmp/5V0O.inp
```

### Batch Regeneration Script

```bash
#!/bin/bash
# regenerate_mismatched_pdbs.sh

PDBS=("1TTT" "9CF3" "1TN1" "1TN2" "3F2T" "5V0O")

for pdb in "${PDBS[@]}"; do
    echo "Regenerating $pdb with --fix-indices..."
    ./build/find_pair_app --fix-indices data/pdb/${pdb}.pdb /tmp/${pdb}.inp
    if [ $? -eq 0 ]; then
        echo "✅ $pdb regenerated successfully"
    else
        echo "❌ $pdb regeneration failed"
    fi
done
```

---

## Verification Steps

### Step 1: Verify Residue Indices Match

After regeneration, verify indices match:

```bash
# For each PDB
python3 -c "
import json

pdb_id = '1TTT'  # Change for each PDB

# Check legacy
with open(f'data/json_legacy/base_frame_calc/{pdb_id}.json') as f:
    legacy = json.load(f)

# Check modern
with open(f'data/json/base_frame_calc/{pdb_id}.json') as f:
    modern = json.load(f)

# Compare sample residues
for i in [1, 10, 25, 50, 100]:
    legacy_res = next((r for r in legacy if r.get('residue_idx') == i), None)
    modern_res = next((r for r in modern if r.get('residue_idx') == i or r.get('legacy_residue_idx') == i), None)
    
    if legacy_res and modern_res:
        legacy_name = legacy_res.get('residue_name', '?')
        modern_name = modern_res.get('residue_name', '?')
        legacy_seq = legacy_res.get('residue_seq', '?')
        modern_seq = modern_res.get('residue_seq', '?')
        
        if legacy_name == modern_name and legacy_seq == modern_seq:
            print(f'✅ Residue {i}: Match - {legacy_name} seq{legacy_seq}')
        else:
            print(f'❌ Residue {i}: MISMATCH - Legacy: {legacy_name} seq{legacy_seq}, Modern: {modern_name} seq{modern_seq}')
"
```

### Step 2: Re-compare Pairs

```bash
# Compare find_bestpair_selection
python3 scripts/compare_json.py compare 1TTT --record-type find_bestpair_selection --verbose

# Full comparison
python3 scripts/compare_json.py compare 1TTT --verbose
```

### Step 3: Check if Differences Resolved

After regeneration, check if differences are resolved:

```bash
# For each PDB
python3 -c "
import json

pdb_id = '1TTT'

# Load pairs
with open(f'data/json_legacy/find_bestpair_selection/{pdb_id}.json') as f:
    legacy_data = json.load(f)
legacy_pairs = set()
for r in legacy_data:
    if 'pairs' in r:
        for p in r['pairs']:
            if len(p) >= 2:
                legacy_pairs.add((min(p[0], p[1]), max(p[0], p[1])))

with open(f'data/json/find_bestpair_selection/{pdb_id}.json') as f:
    modern_data = json.load(f)
modern_pairs = set()
for r in modern_data:
    if 'pairs' in r:
        for p in r['pairs']:
            if len(p) >= 2:
                modern_pairs.add((min(p[0], p[1]), max(p[0], p[1])))

missing = legacy_pairs - modern_pairs
extra = modern_pairs - legacy_pairs

print(f'{pdb_id}:')
print(f'  Legacy: {len(legacy_pairs)} pairs')
print(f'  Modern: {len(modern_pairs)} pairs')
print(f'  Missing in modern: {len(missing)} - {sorted(missing)}')
print(f'  Extra in modern: {len(extra)} - {sorted(extra)}')
"
```

---

## Expected Results

### Before Regeneration
- ❌ Residue indices don't match
- ❌ Comparing different pairs
- ❌ False validation differences

### After Regeneration
- ✅ Residue indices match exactly
- ✅ Comparing same pairs
- ✅ Real validation differences (if any)

---

## Priority Order

1. **1TTT** - Most critical (clear index mismatch)
2. **9CF3** - Has pair differences
3. **1TN1, 1TN2** - Same pattern, likely same issue
4. **3F2T** - Sequence number offset
5. **5V0O** - Verify (may be OK)

---

## Notes

- The `--fix-indices` option auto-detects legacy JSON from `data/json_legacy/base_frame_calc/<PDB_ID>.json`
- If auto-detection fails, specify explicitly: `--fix-indices=data/json_legacy/base_frame_calc/<PDB_ID>.json`
- After regeneration, all modern JSON files will use legacy indices
- Re-run all comparisons after regeneration

---

## Related Documentation

- [CRITICAL_FINDING_RESIDUE_INDEX_MISMATCH.md](CRITICAL_FINDING_RESIDUE_INDEX_MISMATCH.md) - Critical finding details
- [LEGACY_INDICES_GUIDE.md](LEGACY_INDICES_GUIDE.md) - Complete guide on using legacy indices
- [FIX_INDICES_OPTION.md](FIX_INDICES_OPTION.md) - --fix-indices option documentation

---

*Regenerate all affected PDBs with --fix-indices before continuing investigation.*


# Legacy Indices Guide - Ensuring Correct Comparisons

**Date**: 2025-01-XX  
**Status**: Critical for accurate comparisons

---

## ⚠️ CRITICAL: Always Use Legacy Indices

**Why**: Modern code may assign different residue indices than legacy due to parsing differences. To ensure we're comparing the same residues, we must use legacy indices.

---

## Quick Reference

### Step 1: Generate Modern JSON with Legacy Indices

**ALWAYS do this first** when preparing for comparison:

```bash
# Auto-detect legacy JSON
./build/find_pair_app --fix-indices data/pdb/<PDB_ID>.pdb output.inp

# Or specify legacy JSON explicitly
./build/find_pair_app --fix-indices=data/json_legacy/base_frame_calc/<PDB_ID>.json data/pdb/<PDB_ID>.pdb output.inp
```

**What this does**:
- Matches residues by PDB properties: `(residue_name, chain_id, residue_seq, insertion)`
- Assigns legacy indices from legacy JSON
- Ensures modern JSON uses same indices as legacy

### Step 2: Use Legacy Indices in All Tools

All investigation tools expect **legacy indices** (the indices from legacy JSON files):

```bash
# Quality score comparison (use legacy indices)
./build/compare_quality_scores <PDB_ID> <legacy_residue1> <legacy_residue2>

# H-bond detection (use legacy indices)
./build/detect_hbonds_standalone <pdb_file> <legacy_residue1> <legacy_residue2>
./org/build/bin/test_hbond_detection <pdb_file> <legacy_residue1> <legacy_residue2>

# JSON comparison (assumes modern JSON was generated with --fix-indices)
python3 scripts/compare_json.py compare <PDB_ID> --verbose
```

---

## How to Find Legacy Indices

### From Legacy JSON Files

Legacy indices are in legacy JSON files. For example, in `data/json_legacy/find_bestpair_selection/<PDB_ID>.json`:

```json
{
  "type": "find_bestpair_selection",
  "pairs": [
    [162, 177],  // These are legacy indices
    [16, 59],
    ...
  ]
}
```

### From Legacy Input Files

Legacy input files (`.inp`) also contain legacy indices:

```
# Base pairs (legacy indices)
162 177
16 59
...
```

### From Legacy Base Pair JSON

In `data/json_legacy/base_pair/<PDB_ID>.json`:

```json
{
  "type": "base_pair",
  "base_i": 162,  // Legacy index
  "base_j": 177,  // Legacy index
  ...
}
```

---

## Common Workflow

### Investigating a Specific Pair

**Example**: Investigating pair (162, 177) in 1TTT

```bash
# Step 1: Generate modern JSON with legacy indices
./build/find_pair_app --fix-indices data/pdb/1TTT.pdb /tmp/1TTT.inp

# Step 2: Compare quality scores (using legacy indices 162, 177)
./build/compare_quality_scores 1TTT 162 177

# Step 3: Compare H-bond detection (using legacy indices)
./build/detect_hbonds_standalone data/pdb/1TTT.pdb 162 177
./org/build/bin/test_hbond_detection data/pdb/1TTT.pdb 162 177

# Step 4: Compare JSON files
python3 scripts/compare_json.py compare 1TTT --verbose
```

### Comparing All PDBs in Test Set

```bash
# Step 1: Generate modern JSON for all PDBs with --fix-indices
for pdb in 1Q96 1VBY 3AVY 3G8T 3KNC 4AL5 5UJ2 6CAQ 6LTU 8J1J; do
    ./build/find_pair_app --fix-indices data/pdb/${pdb}.pdb /tmp/${pdb}.inp
done

# Step 2: Compare all
python3 scripts/compare_json.py compare --test-set 10
```

---

## Verification

### Check if Modern JSON Has Correct Indices

Compare residue indices between modern and legacy:

```bash
# Compare residue indices
python3 scripts/compare_json.py compare <PDB_ID> --record-type residue_indices
```

**Expected**: Should show perfect match if `--fix-indices` was used.

### Check if Pairs Match

```bash
# Compare find_bestpair_selection
python3 scripts/compare_json.py compare <PDB_ID> --record-type find_bestpair_selection
```

**Expected**: Pairs should match by indices if `--fix-indices` was used.

---

## Common Mistakes

### ❌ Wrong: Using Modern Indices

```bash
# DON'T do this - modern indices may differ from legacy
./build/compare_quality_scores 1TTT 150 165  # Wrong indices
```

### ✅ Correct: Using Legacy Indices

```bash
# DO this - use indices from legacy JSON files
./build/compare_quality_scores 1TTT 162 177  # Correct legacy indices
```

### ❌ Wrong: Generating Modern JSON Without --fix-indices

```bash
# DON'T do this - indices won't match legacy
./build/find_pair_app data/pdb/1TTT.pdb output.inp
python3 scripts/compare_json.py compare 1TTT  # Will show false differences
```

### ✅ Correct: Generating Modern JSON With --fix-indices

```bash
# DO this - indices will match legacy
./build/find_pair_app --fix-indices data/pdb/1TTT.pdb output.inp
python3 scripts/compare_json.py compare 1TTT  # Will show real differences
```

---

## Troubleshooting

### Issue: Pairs Don't Match Even After --fix-indices

**Possible causes**:
1. Legacy JSON file not found or incorrect path
2. Residue matching failed (PDB properties don't match)
3. Legacy JSON file is corrupted or incomplete

**Solution**:
```bash
# Check if legacy JSON exists
ls -la data/json_legacy/base_frame_calc/<PDB_ID>.json

# Check fix-indices output for warnings
./build/find_pair_app --fix-indices data/pdb/<PDB_ID>.pdb output.inp 2>&1 | grep -i "fix\|index\|match"
```

### Issue: Quality Scores Don't Match

**Possible causes**:
1. Using wrong indices (modern instead of legacy)
2. Modern JSON not generated with --fix-indices
3. Comparing different pairs

**Solution**:
1. Verify indices are from legacy JSON files
2. Regenerate modern JSON with --fix-indices
3. Double-check pair indices match between legacy and modern

---

## Related Documentation

- [FIX_INDICES_OPTION.md](FIX_INDICES_OPTION.md) - Detailed --fix-indices option documentation
- [ALL_DIFFERENCES_SUMMARY.md](ALL_DIFFERENCES_SUMMARY.md) - All differences summary
- [INVESTIGATION_FINDINGS.md](INVESTIGATION_FINDINGS.md) - Investigation findings

---

*Always use legacy indices when comparing modern and legacy code to ensure accurate comparisons.*


# Investigation Findings - Remaining Differences

**Date**: 2025-01-XX  
**Status**: Ongoing investigation of 6 mismatched PDBs

---

## Summary

After reviewing all documentation and code, here are the key findings about remaining differences between modern and legacy find_pair:

---

## 1. H-Bond Type vs Conflict Flag in adjust_pairQuality

### Potential Issue Identified

**Legacy `adjust_pairQuality`** (`org/src/cmn_fncs.c:4559-4578`):
```c
for (k = 1; k <= num_hb; k++) {
    if (num_list[k][0])  // Check conflict FLAG (not type!)
        continue;
    dval = num_list[k][3] / MFACTOR;
    if (dval_in_range(dval, 2.5, 3.5))
        num_good_hb++;
}
```

**Modern `adjust_pair_quality`** (`src/x3dna/algorithms/base_pair_finder.cpp:582-605`):
```cpp
for (const auto& hbond : hbonds) {
    // Skip non-standard hydrogen bonds (type != '-')
    if (hbond.type != '-') {  // Check TYPE (not conflict flag!)
        continue;
    }
    if (hbond.distance >= 2.5 && hbond.distance <= 3.5) {
        num_good_hb++;
    }
}
```

### Key Difference

- **Legacy**: Checks `num_list[k][0]` (conflict flag) - skips H-bonds with conflict flag set
- **Modern**: Checks `hbond.type != '-'` - skips H-bonds that are not type '-'

### Are These Equivalent?

**Not necessarily!** The conflict flag and H-bond type are determined differently:

1. **Conflict Flag** (`num_list[k][0]`):
   - Set during `hb_numlist()` based on linkage types
   - Marks H-bonds that share atoms with other H-bonds (linkage_type != 18)
   - Can be set even if H-bond is valid

2. **H-Bond Type** (`hbond.type`):
   - Determined by `donor_acceptor()` function
   - Only assigned to H-bonds with NEGATIVE distance (conflicts marked by negative)
   - H-bonds with positive distance remain type=' ' (space) and are NOT processed

**Key Insight**: In legacy `hb_numlist`, the conflict flag `num_list[k][0]` is set based on linkage types (from `hb_atompair`), not just whether distance is negative. This means:
- An H-bond can have conflict flag set (`num_list[k][0] != 0`) but still have positive distance
- Legacy `adjust_pairQuality` skips such H-bonds
- Modern `adjust_pair_quality` would count them if they have type='-'

**However**: After `validate_hbonds`, only H-bonds with negative distance get their type assigned. H-bonds with positive distance remain type=' '. So in practice:
- Modern receives H-bonds after validation
- Only H-bonds with type != ' ' are in the list (those that had negative distance)
- So checking `type != '-'` might be equivalent to checking conflict flag, but need to verify

### Impact

If an H-bond:
- Has conflict flag set (`num_list[k][0] != 0`) but type='-'
- Legacy: Skipped (conflict flag check)
- Modern: Counted (type='-' check)

This could cause modern to count more "good" H-bonds than legacy, leading to different `adjust_pairQuality` values.

### Investigation Needed

1. Check if H-bonds with conflict flag set can still have type='-'
2. Compare conflict flag assignment between legacy and modern
3. Verify if this causes quality score differences for mismatched pairs

---

## 2. H-Bond Validation Logic

### Legacy Behavior

**`validate_hbonds`** (`org/src/cmn_fncs.c:3989-4019`):
- Only processes H-bonds with **NEGATIVE distance** (conflicts)
- If `hb_dist[k] > 0.0`, it continues (skips)
- Only conflicted H-bonds get their type assigned

### Modern Behavior

**`validate_hbonds`** (`src/x3dna/algorithms/hydrogen_bond_finder.cpp:240-297`):
- Also only processes H-bonds with **NEGATIVE distance**
- If `hbond.distance > 0.0`, it continues (skips)
- Matches legacy behavior

### Status

✅ **Matches** - Both only process negative distances (conflicts)

---

## 3. adjust_pairQuality Logic

### Formula

Both implementations use the same formula:
```cpp
if (num_good_hb >= 2)
    return -3.0;
else
    return -num_good_hb;
```

### Good H-Bond Criteria

**Legacy**:
- Conflict flag NOT set (`num_list[k][0] == 0`)
- Distance in [2.5, 3.5] Angstroms

**Modern**:
- Type is '-' (`hbond.type == '-'`)
- Distance in [2.5, 3.5] Angstroms

### Potential Issue

The criteria differ:
- Legacy: Conflict flag check
- Modern: Type check

These may not be equivalent (see Section 1).

---

## 4. Specific Pairs to Investigate

### 1TTT

**Missing in Modern**: (162, 177)
- Legacy: Validated and selected
- Modern: Invalid (`is_valid=0`)
- Need to check: Why does modern invalidate this pair?

**Extra in Modern**: (16, 59), (177, 197)
- (16, 59): Modern validates, legacy doesn't
- (177, 197): Both validate, but legacy selects (162, 177) instead
- Need to check: Quality score differences

### 9CF3

**Missing in Modern**: (25, 27)
- Legacy: Validated and selected
- Modern: Invalid (`is_valid=0`)
- Need to check: Why does modern invalidate this pair?

**Extra in Modern**: (27, 29)
- Modern: Invalid (`is_valid=0`) but STILL SELECTED ⚠️
- Legacy: Not validated
- **Suspicious**: Pairs should not be selected if invalid!

### Other PDBs

- **1TN1, 1TN2**: Extra pair (45, 77)
- **3F2T**: Extra pair (96, 110)
- **5V0O**: Extra pair (9, 34)

---

## 5. Validation Thresholds

### Need to Verify

All validation parameters should match exactly:
- `min_dorg`, `max_dorg`
- `min_dv`, `max_dv`
- `min_dNN`, `max_dNN`
- `min_plane_angle`, `max_plane_angle`
- `min_base_hb`
- `overlap_threshold`

### Code Locations

- **Legacy**: `org/src/cmn_fncs.c` - `check_pair()`, validation thresholds
- **Modern**: `include/x3dna/algorithms/base_pair_validator.hpp` - Validation thresholds

### Action Needed

Compare all threshold values to ensure exact match.

---

## 6. Quality Score Calculation Flow

### Legacy Flow

1. Calculate base score: `rtn_val[5] = dorg + 2.0*d_v + plane_angle/20.0`
2. Apply `adjust_pairQuality()`: `rtn_val[5] += adjust_pairQuality(...)`
3. Apply `bp_type_id == 2` adjustment: `rtn_val[5] -= 2.0` (if applicable)
4. Use `rtn_val[5]` for pair selection

### Modern Flow

1. Calculate base score: `quality_score = dorg + 2.0*d_v + plane_angle/20.0`
2. Apply `adjust_pair_quality()`: `adjusted_quality_score = quality_score + adjust_pair_quality(...)`
3. Apply `bp_type_id == 2` adjustment: `adjusted_quality_score -= 2.0` (if applicable)
4. Use `adjusted_quality_score` for pair selection

### Status

✅ **Flow matches** - Same order of operations

---

## 7. Next Steps

### Priority 1: Fix H-Bond Conflict Flag Check

**Action**: Update `adjust_pair_quality()` to check conflict flag instead of type

**Current Code**:
```cpp
if (hbond.type != '-') {
    continue;
}
```

**Proposed Fix**:
```cpp
// Check conflict flag (matches legacy num_list[k][0] check)
// Need to add conflict flag to HydrogenBond structure
if (hbond.has_conflict) {  // or similar field
    continue;
}
```

**Challenge**: Need to ensure conflict flag is properly set during H-bond detection

### Priority 2: Investigate Invalid Pair Selection

**Issue**: Pair (27, 29) in 9CF3 is marked invalid but still selected

**Action**: 
1. Check why validation records show `is_valid=0` but pair is selected
2. Verify Phase 1 validation results match what's written to JSON
3. Ensure pairs cannot be selected if they fail validation

### Priority 3: Compare Validation Thresholds

**Action**: 
1. Extract all validation thresholds from legacy code
2. Compare with modern thresholds
3. Fix any mismatches

### Priority 4: Debug Specific Pairs

**Action**: For each mismatched pair:
1. **Ensure modern JSON was generated with `--fix-indices`**:
   ```bash
   ./build/find_pair_app --fix-indices data/pdb/<PDB_ID>.pdb /tmp/output.inp
   ```
2. **Use legacy indices when looking up pairs** (indices from legacy JSON files)
3. Compare H-bond detection results using legacy indices
4. Compare quality scores using legacy indices
5. Compare validation results using legacy indices
6. Identify root cause

**Example for 1TTT pair (162, 177)**:
```bash
# Step 1: Generate modern JSON with legacy indices
./build/find_pair_app --fix-indices data/pdb/1TTT.pdb /tmp/1TTT.inp

# Step 2: Compare quality scores (using legacy indices 162, 177)
./build/compare_quality_scores 1TTT 162 177

# Step 3: Compare H-bond detection (using legacy indices)
./build/detect_hbonds_standalone data/pdb/1TTT.pdb 162 177
./org/build/bin/test_hbond_detection data/pdb/1TTT.pdb 162 177
```

---

## 8. Tools Available

### ⚠️ CRITICAL: Use Legacy Indices for All Comparisons

**IMPORTANT**: All tools assume indices are in **legacy format**. When generating modern JSON for comparison, you MUST use the `--fix-indices` option:

```bash
# Generate modern JSON with legacy indices
./build/find_pair_app --fix-indices data/pdb/<PDB_ID>.pdb output.inp

# Or specify legacy JSON explicitly
./build/find_pair_app --fix-indices=data/json_legacy/base_frame_calc/<PDB_ID>.json data/pdb/<PDB_ID>.pdb output.inp
```

**Why**: Modern code may assign different residue indices than legacy due to parsing differences. The `--fix-indices` option matches residues by PDB properties `(residue_name, chain_id, residue_seq, insertion)` and assigns legacy indices, ensuring correct comparison.

### Quality Score Comparison
```bash
# Use LEGACY indices (from legacy JSON files)
./build/compare_quality_scores <PDB_ID> <legacy_residue1> <legacy_residue2>

# Example: Compare pair (162, 177) in 1TTT using legacy indices
./build/compare_quality_scores 1TTT 162 177
```

**Note**: The tool looks up pairs by indices in JSON files. These indices must be legacy indices for correct comparison.

### H-Bond Detection
```bash
# Modern H-bond detection (use legacy indices)
./build/detect_hbonds_standalone <pdb_file> <legacy_residue1> <legacy_residue2>

# Legacy H-bond detection (use legacy indices)
./org/build/bin/test_hbond_detection <pdb_file> <legacy_residue1> <legacy_residue2>
```

**Note**: Residue indices passed to these tools should be legacy indices to match the residues correctly.

### JSON Comparison
```bash
# Compare JSON files (assumes modern JSON was generated with --fix-indices)
python3 scripts/compare_json.py compare <PDB_ID> --verbose
```

**Note**: This tool compares indices directly from JSON files. Modern JSON must be generated with `--fix-indices` to ensure indices match legacy.

### Pair Investigation
```bash
# Investigate pair differences (uses indices from JSON files)
python3 scripts/investigate_pair_differences.py
```

**Note**: This script reads pairs from JSON files. Modern JSON must be generated with `--fix-indices` to ensure correct pair matching.

### Generating Modern JSON for Comparison

**ALWAYS use `--fix-indices` when generating modern JSON for comparison**:

```bash
# Step 1: Generate modern JSON with legacy indices
./build/find_pair_app --fix-indices data/pdb/<PDB_ID>.pdb /tmp/output.inp

# Step 2: Compare (indices will match)
python3 scripts/compare_json.py compare <PDB_ID> --verbose
```

**Without `--fix-indices`**: Indices may not match, leading to incorrect comparisons and false differences.

---

## 9. Code Locations

### Legacy
- `org/src/cmn_fncs.c:4559-4578` - `adjust_pairQuality()`
- `org/src/cmn_fncs.c:4580-4652` - `check_pair()` (validation)
- `org/src/cmn_fncs.c:3932-4019` - H-bond conflict resolution and validation

### Modern
- `src/x3dna/algorithms/base_pair_finder.cpp:582-605` - `adjust_pair_quality()`
- `src/x3dna/algorithms/base_pair_validator.cpp` - Validation logic
- `src/x3dna/algorithms/hydrogen_bond_finder.cpp` - H-bond detection and validation

---

## 10. Summary

### Key Findings

1. **H-Bond Conflict Flag vs Type**: Potential mismatch in `adjust_pairQuality` logic
2. **Invalid Pair Selection**: Suspicious behavior in 9CF3
3. **Validation Differences**: Some pairs validated differently between legacy and modern
4. **Quality Score Differences**: When both validate, different "best partners" selected

### Root Causes

1. **H-Bond Counting**: Conflict flag check vs type check may not be equivalent
2. **Validation Logic**: May have subtle differences in threshold enforcement
3. **Quality Score Calculation**: H-bond adjustments may differ due to counting differences

### Next Actions

1. Fix H-bond conflict flag check in `adjust_pair_quality()`
2. Investigate invalid pair selection bug
3. Compare and verify all validation thresholds
4. Debug specific mismatched pairs

---

## 11. Legacy Indices - Critical for Accurate Comparisons

### ⚠️ CRITICAL FINDING: Residue Index Mismatch in 1TTT

**DISCOVERED**: Modern JSON for 1TTT was **NOT generated with `--fix-indices`**, causing residue indices to not match legacy!

**Evidence**:
- Legacy residue 162: **2MG** (modified guanine), chain F, seq 10
- Modern residue 162: **C** (cytosine), chain F, seq 11
- **These are completely different residues!**

**Impact**: All differences in 1TTT are likely due to comparing different pairs, not real validation differences.

**Solution**: Regenerate modern JSON with `--fix-indices`:
```bash
./build/find_pair_app --fix-indices data/pdb/1TTT.pdb /tmp/1TTT.inp
```

See [CRITICAL_FINDING_RESIDUE_INDEX_MISMATCH.md](CRITICAL_FINDING_RESIDUE_INDEX_MISMATCH.md) for complete details.

### ⚠️ IMPORTANT

**All comparisons must use legacy indices** to ensure we're comparing the same residues. See [LEGACY_INDICES_GUIDE.md](LEGACY_INDICES_GUIDE.md) for complete details.

### Quick Checklist

- [ ] Generate modern JSON with `--fix-indices` before comparing
- [ ] Use legacy indices (from legacy JSON files) in all tools
- [ ] Verify indices match between legacy and modern JSON files
- [ ] Double-check pair indices when investigating specific pairs

### Example Workflow

```bash
# Step 1: Generate modern JSON with legacy indices
./build/find_pair_app --fix-indices data/pdb/1TTT.pdb /tmp/1TTT.inp

# Step 2: Use legacy indices in investigation tools
./build/compare_quality_scores 1TTT 162 177  # 162, 177 are legacy indices
```

---

*Last Updated: 2025-01-XX*


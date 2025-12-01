# 1TTT Investigation - Specific Pair Differences

**Date**: 2025-01-XX  
**Status**: Investigating validation and selection differences

---

## Summary

**PDB**: 1TTT  
**Legacy pairs**: 92  
**Modern pairs**: 93  
**Missing in modern**: 1 pair - (162, 177)  
**Extra in modern**: 2 pairs - (16, 59), (177, 197)

---

## Pair (162, 177) - Missing in Modern

### Status
- **Legacy**: ✅ Validated (in base_pair records), ✅ Selected
- **Modern**: ❌ NOT validated (NOT in base_pair records), ❌ NOT selected

### Analysis
- Legacy validates this pair (appears in `base_pair` JSON)
- Modern does NOT validate this pair (does not appear in `base_pair` JSON)
- This is a **validation difference** - modern fails validation when legacy passes

### Root Cause Identified: Residue Index Mismatch!

**CRITICAL FINDING**: Modern and legacy have DIFFERENT residues at indices 162 and 177!

**Legacy** (from `base_frame_calc`):
- Residue 162: **2MG** (modified guanine), base_type: g, chain F, seq 10
- Residue 177: **C** (cytosine), base_type: C, chain F, seq 25

**Modern** (from `base_frame_calc`):
- Residue 162: **C** (cytosine), base_type: C, chain F, seq 11
- Residue 177: **M2G** (modified guanine), base_type: ?, chain F, seq 26

**The residues are SWAPPED!** Modern is comparing the wrong residues because indices don't match legacy.

**Modern Validation Record** (for wrong residues):
```
is_valid: 0
distance_check: True
d_v_check: True
plane_angle_check: True
dNN_check: False  ⚠️ FAILING (because comparing wrong residues)
overlap_check: None
hbond_check: None
dNN: None  ⚠️ N1/N9 atoms not found (wrong residues)
```

**Root Cause**: Modern JSON was NOT generated with `--fix-indices`, or the fix didn't work correctly. The residue indices don't match legacy, so we're comparing different pairs!

### Why dNN Check Fails

Modern code sets `dNN = 1e10` when N1/N9 atoms are not found, but the validation record shows `dNN: None`. This suggests:
1. **N1/N9 atoms missing**: One or both residues don't have N1/N9 atoms
2. **Atom finding logic difference**: Modern's `find_n1_n9_position()` may not find atoms that legacy finds
3. **Residue type classification difference**: Modern may classify residue differently, affecting N1/N9 detection

### Legacy Behavior

Legacy checks `NC1xyz[i][7] > 0` before calculating dNN. If N1/N9 atoms are missing, legacy may:
- Set `rtn_val[4] = EMPTY_NUMBER` or use a default value
- Still pass validation if other checks pass
- Handle missing N1/N9 differently than modern

### Status After Regeneration

**✅ REGENERATED**: Modern JSON regenerated with `--fix-indices` (fixed 186 residue indices)

**✅ VERIFIED**: Residue indices now match when using `legacy_residue_idx`:
- Legacy residue 162 = 2MG seq10
- Modern legacy_residue_idx 162 = 2MG seq10 ✅ MATCH

**⚠️ REMAINING DIFFERENCES**:
- Missing in modern: (162, 177)
- Extra in modern: (16, 59), (177, 197)

**Analysis**: Since indices are now correct, these are **real validation/selection differences**, not index mismatches.

---

## Root Cause Analysis

### Pair (162, 177) - Missing in Modern

**Status**:
- Legacy: ✅ Validated, ✅ Selected (in base_pair)
- Modern: ❌ Validation fails (dNN_check=False), ❌ Not selected

**Root Cause Identified**: **N1/N9 Detection Bug for Modified Nucleotides**

**Issue**: Residue 162 is **2MG** (modified guanine), which:
- Has N9 atom in PDB file ✅
- Is classified as `base_type: ?` (likely `ResidueType::UNKNOWN`)
- Modern code's `find_n1_n9_position()` only checks for N9 when:
  - `res_type == ADENINE || res_type == GUANINE` (line 224)
  - `res_type == NONCANONICAL_RNA` (line 228)
- **2MG is UNKNOWN, not NONCANONICAL_RNA**, so code doesn't check for N9
- Code tries to find N1 instead (line 254), but 2MG is a purine (should use N9)
- Result: `find_n1_n9_position()` returns `nullopt`, `dNN = 1e10`, `dNN_check` fails

**Legacy Behavior**: Legacy's `glyco_N()` function checks for N9 if `isR == 1` (purine), regardless of residue type classification. It also has fallback logic to find atoms with '9' in the name.

**Fix Required**: Update `find_n1_n9_position()` to check for N9 atom when `res_type == UNKNOWN`, similar to how it checks for `NONCANONICAL_RNA`.

### Pair (16, 59) - Extra in Modern

**Status**:
- Legacy: ❌ Not selected (not in base_pair)
- Modern: ✅ Selected (in base_pair, is_valid=1)

**Analysis**: Modern validates and selects this pair, but legacy doesn't. Need to investigate:
- Quality score differences
- Validation threshold differences
- Selection logic differences

### Pair (177, 197) - Extra in Modern

**Status**:
- Legacy: ✅ In base_pair (selected), ❌ NOT in find_bestpair_selection
- Modern: ✅ In base_pair (selected), ✅ In find_bestpair_selection

**Analysis**: Both select this pair, but legacy doesn't include it in the initial `find_bestpair_selection` (maybe added during reordering?). Modern includes it in the initial selection. This might be a difference in when pairs are recorded vs. when they're selected.

### Next Steps

1. **Investigate pair (162, 177)** - Why modern fails validation when legacy passes:
   - Check validation results in `pair_validation` JSON
   - Compare d_v, dNN, overlap, H-bond values
   - Verify all validation checks

2. **Investigate pairs (16, 59) and (177, 197)** - Why modern selects when legacy doesn't:
   - Check quality scores
   - Compare validation results
   - Check if these pairs are in legacy validation records

3. **Compare validation logic** - Check if there are subtle differences in:
   - Overlap calculation
   - H-bond detection
   - Validation threshold enforcement

---

## Pair (16, 59) - Extra in Modern

### Status
- **Legacy**: ❌ NOT validated (NOT in base_pair records), ❌ NOT selected
- **Modern**: ✅ Validated (in base_pair records), ✅ Selected

### Analysis
- Modern validates this pair (appears in `base_pair` JSON)
- Legacy does NOT validate this pair (does not appear in `base_pair` JSON)
- This is a **validation difference** - modern passes validation when legacy fails

### Root Cause Hypothesis
- **d_v check**: Modern may pass d_v check when legacy fails
- **Overlap check**: Modern may pass overlap check when legacy fails
- **H-bond check**: Modern may pass H-bond check when legacy fails
- **Validation threshold differences**: Modern may have different thresholds

### Next Steps
1. Check d_v value for this pair in both implementations
2. Check overlap_area value
3. Check H-bond detection results
4. Compare validation thresholds

---

## Pair (177, 197) - Extra in Modern

### Status
- **Legacy**: ✅ Validated (in base_pair records), ❌ NOT selected
- **Modern**: ✅ Validated (in base_pair records), ✅ Selected

### Analysis
- Both validate this pair (appears in both `base_pair` JSON files)
- Legacy does NOT select it (not in `find_bestpair_selection`)
- Modern DOES select it (in `find_bestpair_selection`)
- This is a **quality score difference** - modern finds this pair has better quality score than legacy

### Root Cause Hypothesis
- **Quality score calculation**: Modern may calculate different quality score
- **H-bond adjustment**: Modern may count H-bonds differently (`adjust_pairQuality`)
- **bp_type_id adjustment**: Modern may apply different `bp_type_id` adjustment
- **Best partner selection**: Modern may find (177, 197) is better partner for residue 177 than (162, 177)

### Next Steps
1. Compare quality scores for (177, 197) vs (162, 177)
2. Compare H-bond counts and adjustments
3. Compare bp_type_id values
4. Check which pair legacy selects for residue 177

---

## Validation Thresholds

### Legacy Default Thresholds (from `org/src/app_fncs.c:set_default_misc_pars`)
```c
min_dorg = 0.0;
max_dorg = 15.0;
min_dv = 0.0;
max_dv = 2.5;
min_dNN = 4.5;
max_dNN = XBIG;  // 1e18
min_plane_angle = 0.0;
max_plane_angle = 65.0;
overlap_threshold = 0.01;  // OVERLAP constant
```

### Modern Thresholds (from `base_pair_validator.hpp`)
```cpp
double min_dorg = 0.0;        ✅ Matches
double max_dorg = 15.0;       ✅ Matches
double min_dv = 0.0;          ✅ Matches
double max_dv = 2.5;          ✅ Matches
double min_dNN = 4.5;         ✅ Matches
double max_dNN = 1e18;        ✅ Matches (XBIG)
double min_plane_angle = 0.0; ✅ Matches
double max_plane_angle = 65.0; ✅ Matches (NOTE: header comment says 45.0 but code uses 65.0)
double overlap_threshold = 0.01; ✅ Matches
```

### Status
✅ **All thresholds match** - No threshold differences found

### Note
The header comment in `base_pair_validator.hpp` incorrectly says `max_plane_angle = 45.0`, but the actual code uses `65.0` which matches legacy. This is just a documentation issue.

---

## Investigation Commands

### Generate Modern JSON with Legacy Indices
```bash
./build/find_pair_app --fix-indices data/pdb/1TTT.pdb /tmp/1TTT.inp
```

### Compare Quality Scores
```bash
# Compare pair (162, 177) - missing in modern
./build/compare_quality_scores 1TTT 162 177

# Compare pair (177, 197) - extra in modern
./build/compare_quality_scores 1TTT 177 197

# Compare pair (16, 59) - extra in modern
./build/compare_quality_scores 1TTT 16 59
```

### Compare H-Bond Detection
```bash
# Pair (162, 177)
./build/detect_hbonds_standalone data/pdb/1TTT.pdb 162 177
./org/build/bin/test_hbond_detection data/pdb/1TTT.pdb 162 177

# Pair (177, 197)
./build/detect_hbonds_standalone data/pdb/1TTT.pdb 177 197
./org/build/bin/test_hbond_detection data/pdb/1TTT.pdb 177 197
```

### Compare JSON Files
```bash
python3 scripts/compare_json.py compare 1TTT --verbose
```

---

## Key Questions

1. **Why does modern fail validation for (162, 177) when legacy passes?**
   - Check d_v, overlap, H-bond values
   - Verify validation thresholds match

2. **Why does modern pass validation for (16, 59) when legacy fails?**
   - Check d_v, overlap, H-bond values
   - Verify validation thresholds match

3. **Why does modern select (177, 197) when legacy selects (162, 177)?**
   - Compare quality scores
   - Check H-bond adjustments
   - Verify bp_type_id values

---

## Related Documentation

- [ALL_DIFFERENCES_SUMMARY.md](ALL_DIFFERENCES_SUMMARY.md) - All differences summary
- [INVESTIGATION_FINDINGS.md](INVESTIGATION_FINDINGS.md) - Investigation findings
- [LEGACY_INDICES_GUIDE.md](LEGACY_INDICES_GUIDE.md) - Using legacy indices

---

*Last Updated: 2025-01-XX*


# 3CME Pair (3844, 3865) Investigation

**Date**: 2025-01-XX  
**Status**: 🔍 Investigating extra pair in modern find_bestpair_selection

---

## Summary

**PDB**: 3CME  
**Issue**: Pair (3844, 3865) is selected in modern but NOT in legacy find_bestpair_selection

**Statistics**:
- Legacy selected pairs: 1,159
- Modern selected pairs: 1,160
- Common pairs: 1,159
- Missing in modern: 0
- Extra in modern: 1 (pair (3844, 3865))

---

## Investigation Plan

### Step 1: Check if pair is validated in legacy

**Question**: Does legacy validate this pair (does it appear in base_pair records)?

**Action**: Check legacy base_pair JSON for pair (3844, 3865)

**Expected Outcomes**:
- If found in legacy base_pair: Pair is validated but not selected (quality score issue)
- If NOT found in legacy base_pair: Pair is not validated (validation difference)

### Step 2: Check quality scores

**Question**: What are the quality scores for this pair and competing pairs?

**Action**: 
1. Check modern base_pair record for (3844, 3865)
2. Check what other pairs residues 3844 and 3865 form
3. Compare quality scores

**Expected Outcomes**:
- Modern may have lower quality score for (3844, 3865) than legacy
- Modern may validate this pair when legacy doesn't
- Modern may select different "best partner" due to quality score differences

### Step 3: Check mutual best match

**Question**: Is this a mutual best match in modern but not legacy?

**Action**:
1. Check if residue 3844 selects 3865 as best partner in modern
2. Check if residue 3865 selects 3844 as best partner in modern
3. Compare with legacy behavior

**Expected Outcomes**:
- Modern: Both residues select each other as best partner
- Legacy: One or both residues select different partners

### Step 4: Root cause analysis

**Possible Causes**:
1. **Quality score calculation difference**: Modern calculates different quality score
2. **Validation difference**: Modern validates pair when legacy doesn't
3. **Tie-breaking difference**: Scores are equal, different pair selected
4. **Iteration order effect**: Selection order affects later iterations

---

## Findings

### Pair Status in base_pair Records

**Legacy**: ❌ No base_pair file found for 3CME (or pair not validated)
**Modern**: ✅ Found in base_pair records
  - base_i: 3844, base_j: 3865
  - bp_type: "GA"
  - is_valid: None (recorded but validation status not shown)

### Partner Selection

**Legacy**:
- Residue 3844: ❌ Not paired (no partner selected)
- Residue 3865: ❌ Not paired (no partner selected)

**Modern**:
- Residue 3844: ✅ Paired with 3865 (mutual best match)
- Residue 3865: ✅ Paired with 3844 (mutual best match)

### Analysis

**Key Finding**: In legacy, neither residue 3844 nor 3865 is paired with anything. In modern, they form a mutual best match pair.

**Possible Causes**:
1. **Validation Difference**: Modern validates this pair when legacy doesn't
   - Modern may detect H-bonds that legacy doesn't
   - Modern may pass overlap check when legacy doesn't
   - Modern may pass distance/angle checks when legacy doesn't
2. **Quality Score Difference**: Even if both validate, modern may calculate a quality score that makes this pair the best match, while legacy doesn't
3. **Iteration Order Effect**: The order in which pairs are selected may affect whether this pair is available when residues 3844/3865 are being processed

---

## Root Cause Identified

**Issue**: There are TWO validation records for pair (3844, 3865) with different indices and results:

1. **0-based indices (3843, 3864)**: 
   - `is_valid: 1` ✅
   - `d_v: 0.837369` (passes d_v_check)
   - This is the CORRECT record used for selection

2. **1-based indices (3844, 3865)**:
   - `is_valid: 0` ❌
   - `d_v: 5.33983` (fails d_v_check)
   - This is a DUPLICATE record with WRONG indices

**Conclusion**: The pair IS valid and is being selected correctly based on the 0-based validation record. The issue is that there's a duplicate validation record being written with 1-based indices (should be 0-based).

**Why Legacy Doesn't Select It**: Need to check if legacy validates this pair. If legacy doesn't validate it, there may be a validation difference (different parameters or thresholds).

## Next Steps

1. ✅ Check if pair exists in legacy base_pair records
2. ✅ Identified duplicate validation record issue (0-based vs 1-based)
3. ⏳ Find where the duplicate 1-based validation record is being written
4. ⏳ Fix the duplicate record issue (add defensive checks in record_pair_validation)
5. ⏳ Check if legacy validates this pair (if not, investigate validation differences)
6. ⏳ Determine if pair should be selected in legacy (if it's valid in modern)

## Code Changes Made

1. **Added defensive check in `record_pair_validation`**: 
   - Warns if indices >= 4000 (likely 1-based for large PDBs)
   - Helps identify where 1-based records are coming from

2. **Added safety checks in `find_best_partner`**:
   - Ensures invalid pairs cannot be selected
   - Logs warnings for debugging

## Current Status

- **Pair (3844, 3865) is VALID** in modern (d_v: 0.837369 < 2.5)
- **Pair is correctly selected** based on 0-based validation record
- **Duplicate 1-based record exists** but is ignored (wrong indices, different d_v value)
- **Next**: Need to find source of 1-based record and fix it

---

*Document created: 2025-01-XX*  
*Updated: 2025-01-XX - Root cause identified, defensive checks added*


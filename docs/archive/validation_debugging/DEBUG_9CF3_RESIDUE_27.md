# Debug Report: 9CF3 Residue 27 Pair Selection Issue

**Date**: 2025-12-01  
**Status**: 🔴 **CRITICAL BUG FOUND** - Off-by-One Error in Residue Index Mapping

---

## Problem

**PDB**: 9CF3  
**Issue**: Modern code selects different best partner for residue 27 than legacy

### Legacy Behavior
- Residue 27 pairs with residue 25 → Pair (25, 27) selected
- Pair (25, 27) is validated and selected as mutual best match

### Modern Behavior  
- Residue 27 pairs with residue 29 → Pair (27, 29) selected
- Pair (27, 29) is selected instead of (25, 27)

### Impact
- **Missing in Modern**: (25, 27)
- **Extra in Modern**: (27, 29)
- Both legacy and modern find 40 base pairs total, but one pair is different

---

## Root Cause Hypothesis

According to documentation:
1. **Quality Score Calculation Differences**: Different quality scores for (25, 27) vs (27, 29)
2. **Validation Differences**: Possible differences in H-bond counting or validation thresholds
3. **Best Partner Selection**: Residue 27's best partner is different between legacy and modern

---

## Investigation Steps

### Step 1: Compare Best Partner Candidates

Use step-by-step debugging infrastructure:

```bash
# Compare best partner candidates for residue 27
python3 scripts/compare_best_partner.py 9CF3 27 --verbose
```

**Expected**: Should show:
- Which candidates legacy considers for residue 27
- Which candidates modern considers for residue 27
- Quality scores for each candidate
- Which candidate is selected as "best partner"

### Step 2: Compare Mutual Best Decisions

```bash
# Check mutual best decisions
python3 scripts/compare_mutual_best.py 9CF3 --verbose | grep -E "27|25|29"
```

**Expected**: Should show:
- Whether (25, 27) is mutual best in legacy
- Whether (27, 29) is mutual best in modern
- Quality scores for both pairs

### Step 3: Extract Minimal Test Case

If needed, extract minimal PDB around residues 25, 27, 29:

```bash
# Extract residues involved
python3 scripts/extract_residues.py data/pdb/9CF3.pdb 25 27 29 --output data/pdb/minimal/9CF3_res_25_27_29.pdb

# Then debug the minimal case
./scripts/debug_minimal_case.sh 9CF3_res_25_27_29
```

---

## Root Cause Identified

### Key Finding

**Modern code has BOTH candidates with correct scores:**
- Candidate 25: score = 9.894551 (same as legacy), eligible = 1
- Candidate 29: score = 7.532966, eligible = 1

**But modern selects candidate 29 (score 7.53) instead of candidate 25 (score 9.89)!**

**Legacy behavior:**
- Candidate 25: score = 9.894551, selected as best
- Candidate 29: score = 1e18 (invalid/ineligible), not considered

### Issue

Modern considers pair (27, 29) as **valid** with a quality score of 7.53, while legacy considers it **invalid** (score 1e18 = sentinel for ineligible pairs).

**This suggests a validation difference:**
- Modern validates (27, 29) → quality_score = 7.53 → selected as best
- Legacy does NOT validate (27, 29) → never considered

## Root Cause: Validation Threshold Mismatch

### Critical Finding

**Modern validation REJECTS pair (25, 27), but legacy ACCEPTS it:**

| Pair | Modern | Legacy |
|------|--------|--------|
| (25, 27) | is_valid=0, quality_score=27.33 | is_valid=1, quality_score=9.89 |
| (27, 29) | is_valid=0, quality_score=25.71 | is_valid=0, quality_score=7.53 |

### The Problem

1. **Modern rejects (25, 27)** → quality_score becomes 27.33 (high penalty)
2. **Legacy accepts (25, 27)** → quality_score is 9.89 (valid pair)
3. Modern then selects (27, 29) with score 7.53 because (25, 27) is marked invalid
4. Legacy selects (25, 27) with score 9.89 because it's valid

### Why the Scores Differ

The quality_score in `pair_validation` records is **after validation**. When a pair fails validation:
- Modern: Sets quality_score to a high value (27.33) indicating rejection
- Legacy: Still shows the calculated quality_score (9.89) in some contexts, but validation may use different thresholds

**But in `best_partner_candidates`, modern shows candidate 25 with score 9.894551 (matching legacy),** which suggests the **base quality score calculation is correct**, but the **validation decision differs**.

## Critical Discovery: Geometry Values Are Completely Different

### The Real Problem

Modern and legacy calculate **completely different geometry values** for pair (25, 27):

| Parameter | Modern | Legacy | Threshold | Status |
|-----------|--------|--------|-----------|--------|
| dorg | 18.49 | 4.87 | 0-15.0 | Modern exceeds ❌ |
| d_v | 2.73 | 1.33 | 0-2.5 | Modern exceeds ❌ |
| plane_angle | 67.63 | 47.28 | 0-65.0 | Modern exceeds ❌ |
| dNN | 12.77 | 6.70 | 4.5-1e18 | Both pass ✅ |

### Direction Vectors Also Different

- **Modern**: dir_x=-0.830, dir_y=-0.232, dir_z=0.381
- **Legacy**: dir_x=0.035, dir_y=0.313, dir_z=0.678

This suggests **different reference frames** are being used, not just validation threshold differences.

### Possible Causes

1. **Frame calculation differences**: Reference frames for residues 25 and/or 27 differ between legacy and modern
2. **Residue index mapping issue**: Wrong residues being compared
3. **Residue order swap**: (25, 27) vs (27, 25) - frame order matters for calculation
4. **Frame orientation differences**: Same frames but different orientations

### Next Steps

1. ✅ Generate modern JSON for 9CF3 (done)
2. ✅ Confirm pair difference exists (done)
3. ✅ Compare best partner candidates (done)
4. ✅ Compare pair validation (done)
5. ✅ Compare calculated geometry values (done - **MAJOR DIFFERENCE FOUND**)
6. ✅ **Compare frame origins**: Extracting frame origins from base_frame_calc JSON (in progress)
7. 🔄 **Verify residue mapping**: Verify residues 25 and 27 are the same physical residues
8. 🔄 **Check frame calculation**: Compare how frames are calculated for these residues

## CRITICAL DISCOVERY: Frame Origins Are Correct!

### Finding

**Modern frame origins (from `frame_calc` JSON) yield the CORRECT dorg:**
- Calculated dorg from modern origins: **4.874563 Å** ✅
- Legacy validation dorg: **4.874563 Å** ✅
- **MATCH PERFECTLY!**

**BUT modern validation reports:**
- Modern validation dorg: **18.490534 Å** ❌

### Root Cause Identified

The **frame origins themselves are correct**, but during validation:
- **Wrong frames are being retrieved** from residues
- **OR wrong residues are being passed** to `validate()`
- **OR frames aren't being set correctly** on residue objects during validation

### Validation Order

The validation is called as `validator_.validate(*res1, *res2)` where:
- `res1` = residue with `legacy_idx1` (25)
- `res2` = residue with `legacy_idx2` (27)
- dorg calculation: `frame1.origin() - frame2.origin()` = `frame_25.origin() - frame_27.origin()`

### Next Steps

1. ✅ Frame origins verified as correct (done)
2. 🔄 **Check if frames are correctly set on residue objects during validation**
3. 🔄 **Check if correct residue objects are being retrieved from `residue_by_legacy_idx`**
4. 🔄 **Check if `residue.reference_frame()` returns the correct frame during validation**
5. 🔄 **Add debug logging to see which frames are actually used during validation**

This suggests a **frame retrieval/setting bug** rather than a frame calculation bug!

---

## 🚨 CRITICAL DISCOVERY: Residue Mapping Issue (NOT 0-index vs 1-index)

**Date**: 2025-12-01

### The Bug

Validation for pair (25, 27) is using frames from residues (26, 28)!

**User Question**: Is this a 0-index vs 1-index issue?  
**Answer**: No, both use 1-based legacy indices. The issue is **iteration order**!

**Evidence**:
- Frame origins for pair (25, 27): dorg = 4.874563 Å ✅ (correct)
- Frame origins for pair (26, 28): dorg = 18.490534 Å
- Validation dorg for pair (25, 27): dorg = 18.490534 Å ❌ (wrong!)

**Conclusion**: When validating pair (25, 27), the code retrieves residues with legacy indices 25 and 27, but those residue objects have frames from residues 26 and 28.

### Root Cause

**NOT a 0-index vs 1-index issue** - both use 1-based legacy indices!

The problem is **iteration order mismatch**:
- When building `residue_by_legacy_idx`, we iterate: `for (const auto& chain : structure.chains()) for (const auto& residue : chain.residues())`
- This gives residues in **chain/residue order**, not **legacy order** (PDB file order)
- Legacy indices are assigned in **PDB file order** during parsing
- When we use `legacy_residue_idx` from atoms as keys, we might be mapping to wrong residue objects

**Solution**: Use `structure.get_residue_by_legacy_idx()` or iterate using `structure.residues_in_legacy_order()` to ensure correct mapping!

### Next Steps

1. **Fix residue mapping**: Use `structure.get_residue_by_legacy_idx(legacy_idx)` to get residues in legacy order
2. **OR**: Iterate using `structure.residues_in_legacy_order()` when building the mapping
3. **Verify**: Test that residue pointers match legacy index expectations

---

## Related Documentation

- `docs/ALL_DIFFERENCES_SUMMARY.md` - Original issue documentation
- `docs/DEBUGGING_WORKFLOW.md` - Step-by-step debugging guide
- `docs/MATCHING_PLAN.md` - Overall matching strategy


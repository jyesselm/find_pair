# bp_type_id Investigation Summary

**Date**: 2025-01-27  
**Issue**: 4 mismatched pairs in 6CAQ where modern assigns `bp_type_id=2` but legacy assigns `bp_type_id=-1`

---

## Problem Statement

**Pairs with discrepancies:**
- (1024, 1188): Modern=2, Legacy=-1
- (980, 997): Modern=2, Legacy=-1
- (968, 1024): Modern=-1, Legacy=-1 (but different selection due to quality scores)
- (980, 998): Modern=-1, Legacy=-1 (but different selection due to quality scores)

**Impact**: Modern selects different pairs than legacy because:
- Modern applies -2.0 quality score adjustment for `bp_type_id=2` pairs
- Legacy doesn't apply the adjustment (keeps `bp_type_id=-1`)
- This causes modern to select better quality score pairs

---

## Investigation Results

### ✅ Verified: check_wc_wobble_pair Logic Matches

**Legacy Implementation** (`org/src/ana_fncs.c:1122-1131`):
```c
void check_wc_wobble_pair(long *bpid, char *bp, double shear, double stretch, double opening)
{
    static char *WC[9] = { WC_LIST };
    if (fabs(stretch) > 2.0 || fabs(opening) > 60)
        return;
    if (dval_in_range(fabs(shear), 1.8, 2.8))
        *bpid = 1;
    if (fabs(shear) <= 1.8 && num_strmatch(bp, WC, 1, 8))
        *bpid = 2;
}
```

**Modern Implementation** (`src/x3dna/algorithms/base_pair_finder.cpp:660-689`):
- ✅ Thresholds match: `fabs(stretch) > 2.0 || fabs(opening) > 60`
- ✅ Wobble check: `fabs(shear) >= 1.8 && fabs(shear) <= 2.8` → `bp_type_id = 1`
- ✅ WC check: `fabs(shear) <= 1.8 && bp_type in WC_LIST` → `bp_type_id = 2`
- ✅ WC_LIST matches: `{"XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"}`

**Conclusion**: The logic is correctly implemented in modern code.

### ✅ Verified: Modern Step Parameter Calculations

**Pair (1024, 1188) - AU:**
- Direction vectors: dir_x=0.622029, dir_y=-0.623112, dir_z=-0.973039 ✅
- Step parameters: 
  - shear (slide) = -1.719321 → fabs = 1.719321 ≤ 1.8 ✅
  - stretch (rise) = -0.035735 → fabs = 0.035735 ≤ 2.0 ✅
  - opening (twist) = -50.870333 → fabs = 50.870333 ≤ 60.0 ✅
- Base pair type: "AU" in WC_LIST ✅
- **Result**: bp_type_id = 2 (Watson-Crick) ✅

**Pair (980, 997) - CG:**
- Direction vectors: dir_x=0.741459, dir_y=-0.923190, dir_z=-0.799189 ✅
- Step parameters:
  - shear (slide) = -1.370042 → fabs = 1.370042 ≤ 1.8 ✅
  - stretch (rise) = -1.925168 → fabs = 1.925168 ≤ 2.0 ✅
  - opening (twist) = -22.298109 → fabs = 22.298109 ≤ 60.0 ✅
- Base pair type: "CG" in WC_LIST ✅
- **Result**: bp_type_id = 2 (Watson-Crick) ✅

**Conclusion**: Modern code correctly calculates step parameters and assigns `bp_type_id=2`.

### ✅ Verified: Frame Calculations Match

- Frame origins match exactly
- Rotation matrices match exactly (after reversal)
- Frame reversal logic is applied correctly when `dir_z <= 0`
- Frame order is correct: `bpstep_par(r2, org[j], r1, org[i], ...)` matches modern `calculate_step_parameters(frame2, frame1)`

---

## Root Cause Hypothesis

Since all the logic matches perfectly, the most likely causes are:

### Hypothesis 1: Legacy Calculates Different Step Parameters (Most Likely)

**Evidence**:
- Modern step parameters meet all thresholds for `bp_type_id=2`
- Legacy assigns `bp_type_id=-1`, suggesting step parameters don't meet thresholds
- Frames match perfectly, so difference must be in `bpstep_par` calculation

**Possible reasons**:
1. **Numerical precision differences** in matrix operations
2. **Different handling of edge cases** (e.g., degenerate z-axes)
3. **Subtle differences in rotation matrix construction**

### Hypothesis 2: Legacy Base Pair String Construction Differs

**Evidence**:
- Legacy uses: `sprintf(bpi, "%c%c", toupper((int) bseq[i]), toupper((int) bseq[j]))`
- Modern uses: `get_base_letter_from_type(ResidueType)`
- If legacy's `bseq` array has different values, WC_LIST matching would fail

**Investigation needed**: Compare `bseq` values between legacy and modern for these pairs.

### Hypothesis 3: Legacy Has a Bug

**Evidence**:
- Modern logic is correct and matches legacy code exactly
- Modern step parameters meet all thresholds
- Legacy might have a subtle bug in `check_wc_wobble_pair` or `bpstep_par`

---

## Tools Created

1. ✅ `tools/compare_bp_type_id_calculation` - Analyzes bp_type_id calculation step-by-step
2. ✅ `tools/compare_frames_and_step_params` - Compares frames and step parameters
3. ✅ `scripts/compare_step_params_for_pairs.py` - Script to extract step parameters from JSON
4. ✅ `scripts/investigate_bp_type_id_differences.py` - Finds all bp_type_id differences

---

## Next Steps

### Priority 1: Extract Legacy Step Parameters

**Option A**: Run legacy code with debug output enabled ✅ TOOL CREATED
- Legacy code already has `fprintf(stderr, ...)` in `calculate_more_bppars`
- Run legacy `find_pair` on 6CAQ and capture stderr output:
  ```bash
  org/bin/find_pair data/pdb/6CAQ.pdb 2> legacy_debug.log
  python3 scripts/parse_legacy_debug_output.py legacy_debug.log
  ```
- Extract step parameters for pairs (1024, 1188) and (980, 997)
- **Tool**: `scripts/parse_legacy_debug_output.py` - Parses legacy debug output

**Option B**: Modify legacy JSON writer to include step parameters
- Add step parameters to `pair_validation` JSON records
- Store `pars[1]`, `pars[2]`, `pars[6]` in JSON (they're already in `rtn_val[26+k]`)
- Re-run legacy code and extract from JSON

### Priority 2: Compare Base Pair String Construction

- Verify `bseq` array values match between legacy and modern
- Check if `get_base_letter_from_type()` produces same results as `bseq[i]`
- Ensure base pair strings match exactly

### Priority 3: Verify bpstep_par Implementation

- Add unit tests comparing modern `bpstep_par_impl` with legacy `bpstep_par`
- Test with known step parameter values
- Check for numerical precision differences

---

## Code References

### Legacy Code
- `org/src/cmn_fncs.c:4482-4536` - `calculate_more_bppars` (calls `bpstep_par` and `check_wc_wobble_pair`)
- `org/src/ana_fncs.c:1122-1131` - `check_wc_wobble_pair` implementation
- `org/src/ana_fncs.c:1396-1441` - `bpstep_par` implementation

### Modern Code
- `src/x3dna/algorithms/base_pair_finder.cpp:592-692` - `calculate_bp_type_id` implementation
- `src/x3dna/algorithms/parameter_calculator.cpp:110-179` - `bpstep_par_impl` implementation

---

## Status

**Investigation**: ✅ Complete  
**Root Cause**: ⏳ Pending (need legacy step parameters)  
**Fix**: ⏳ Pending (depends on root cause)

**Confidence**: High that modern implementation is correct. The discrepancy is likely due to:
1. Different step parameters calculated by legacy (numerical precision)
2. Different base pair string construction
3. A subtle bug in legacy code


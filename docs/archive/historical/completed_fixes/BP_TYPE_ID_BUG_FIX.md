# bp_type_id Bug Fix - Legacy Parameter Order Issue

## Summary

**Status**: ✅ **FIXED** - 100% match achieved for 6CAQ

**Date**: 2025-01-27

**Issue**: 4 mismatched pairs in 6CAQ due to `bp_type_id` calculation differences.

## Root Cause

Legacy code has a bug where it passes **wrong parameters** to `check_wc_wobble_pair`:

### Legacy Code (Buggy)
```c
// In org/src/cmn_fncs.c:4533
check_wc_wobble_pair(bpid, bpi, pars[1], pars[2], pars[6]);
```

Where:
- `pars[1]` = **Shift** (x-displacement)
- `pars[2]` = **Slide** (y-displacement, should be "shear")
- `pars[6]` = **Twist** (rotation about z-axis, should be "opening")

### Function Signature
```c
void check_wc_wobble_pair(long *bpid, char *bp, double shear, double stretch, double opening)
```

**The Bug**: Legacy passes:
- `pars[1]` (Shift) as `shear` ❌ **WRONG** - should be Slide
- `pars[2]` (Slide) as `stretch` ❌ **WRONG** - should be Rise
- `pars[6]` (Twist) as `opening` ✅ **CORRECT**

### Impact

For pairs like (1024, 1188):
- **Legacy**: Uses `fabs(Shift) = 5.397193` as shear → `5.397193 > 1.8` → check fails → `bp_type_id = -1`
- **Modern (before fix)**: Used `fabs(Slide) = 1.719321` as shear → `1.719321 <= 1.8` → check passes → `bp_type_id = 2`
- **Result**: Different `bp_type_id` values → different quality score adjustments → different pair selections

## The Fix

Updated modern code to match legacy's buggy parameter order:

```cpp
// In src/x3dna/algorithms/base_pair_finder.cpp:654-656
// CRITICAL: Legacy has a bug - it passes wrong parameters to check_wc_wobble_pair!
// To match legacy output exactly, we must replicate this bug:
double shear = params.shift;      // BUG: Should be params.slide
double stretch = params.slide;    // BUG: Should be params.rise
double opening = params.twist;    // Correct
```

## Verification

### Before Fix
- **6CAQ**: 619/623 pairs (99.4% match)
- **4 mismatched pairs**: (1024, 1188), (980, 997), (1024, 1190), (980, 998)

### After Fix
- **6CAQ**: ✅ **100% match** - All pairs in `find_bestpair_selection` now match perfectly!

### Test Results
```bash
$ python3 scripts/analyze_mismatched_pairs.py 6CAQ
✅ Perfect match! No mismatched pairs for 6CAQ
```

## Files Modified

1. **`src/x3dna/algorithms/base_pair_finder.cpp`**
   - Updated `calculate_bp_type_id()` to use buggy parameter order
   - Added detailed comments explaining the legacy bug

2. **`docs/LEGACY_STEP_PARAMETER_ANALYSIS.md`**
   - Documented the bug discovery and analysis

3. **`docs/BP_TYPE_ID_BUG_FIX.md`** (this file)
   - Complete documentation of the fix

## Key Learnings

1. **Legacy code has bugs**: Not all legacy behavior is intentional - some are actual bugs that must be replicated for 100% compatibility.

2. **Parameter order matters**: The order of parameters passed to functions is critical, even when the parameter names suggest a different mapping.

3. **Step parameters**: The 6 base pair step parameters are:
   - `pars[1]` = Shift (x-displacement)
   - `pars[2]` = Slide (y-displacement, "shear")
   - `pars[3]` = Rise (z-displacement, "stretch")
   - `pars[4]` = Tilt
   - `pars[5]` = Roll
   - `pars[6]` = Twist ("opening")

4. **Quality score impact**: `bp_type_id=2` (Watson-Crick) triggers a -2.0 quality score adjustment, which can change pair selection.

## Next Steps

- ✅ Fix verified for 6CAQ
- ⏳ Test other PDBs to ensure no regressions
- ⏳ Update overall match rate statistics


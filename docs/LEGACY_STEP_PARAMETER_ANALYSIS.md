# Legacy Step Parameter Analysis

## Summary

Extracted step parameters from legacy debug output for the 4 mismatched pairs in 6CAQ.

## Legacy Step Parameters (from debug output)

### Pair (1024, 1188)
- **Shift**: -5.397193
- **Slide (shear)**: -1.719321
- **Rise (stretch)**: -0.035735
- **Tilt**: -9.237754
- **Roll**: 9.616516
- **Twist (opening)**: -50.870333

**Threshold Checks:**
- `fabs(slide) = 1.719321 <= 1.8` ✅
- `fabs(rise) = 0.035735 <= 2.0` ✅
- `fabs(twist) = 50.870333 <= 60.0` ✅

**Base Pair Type**: AU (in WC_LIST)

**Expected bp_type_id**: 2 (Watson-Crick)

### Pair (980, 997)
- **Shift**: 3.641891
- **Slide (shear)**: -1.370042
- **Rise (stretch)**: -1.925168
- **Tilt**: -8.053285
- **Roll**: -36.058950
- **Twist (opening)**: -22.298109

**Threshold Checks:**
- `fabs(slide) = 1.370042 <= 1.8` ✅
- `fabs(rise) = 1.925168 <= 2.0` ✅
- `fabs(twist) = 22.298109 <= 60.0` ✅

**Base Pair Type**: (need to verify)

**Expected bp_type_id**: 2 (Watson-Crick)

## Critical Discovery

### Legacy Code Bug

Legacy code calls `check_wc_wobble_pair` with:
```c
check_wc_wobble_pair(bpid, bpi, pars[1], pars[2], pars[6]);
```

Where:
- `pars[1]` = **Shift** (not shear!)
- `pars[2]` = **Slide** (shear)
- `pars[6]` = **Twist** (opening)

But `check_wc_wobble_pair` expects:
```c
void check_wc_wobble_pair(long *bpid, char *bp, double shear, double stretch, double opening)
```

**This means legacy is passing:**
- `pars[1]` (Shift) as `shear` ❌ **WRONG PARAMETER!**
- `pars[2]` (Slide) as `stretch` ❌ **WRONG PARAMETER!**
- `pars[6]` (Twist) as `opening` ✅ Correct

**The correct call should be:**
```c
check_wc_wobble_pair(bpid, bpi, pars[2], pars[3], pars[6]);
//                                    ^      ^      ^
//                                 Slide   Rise   Twist
//                                 shear  stretch opening
```

## Impact

This bug means legacy is using **Shift** instead of **Slide** for the shear check, and **Slide** instead of **Rise** for the stretch check. This could explain why legacy assigns `bp_type_id=-1` while modern correctly assigns `bp_type_id=2`.

## Resolution

**FIXED**: Updated modern code to match legacy's buggy parameter order:

```cpp
// In src/x3dna/algorithms/base_pair_finder.cpp
double shear = params.shift;      // BUG: Should be params.slide
double stretch = params.slide;    // BUG: Should be params.rise
double opening = params.twist;    // Correct
```

This ensures modern code produces the same `bp_type_id=-1` for pairs that legacy assigns `-1`, matching legacy output exactly.

## Verification

After this fix:
- Modern code now uses Shift (instead of Slide) for shear threshold check
- Modern code now uses Slide (instead of Rise) for stretch threshold check
- This matches legacy's buggy behavior exactly
- Expected result: 4 mismatched pairs in 6CAQ should now match


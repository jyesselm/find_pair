# Reference Frame Differences Investigation

**Date**: 2025-01-XX  
**Status**: 🔍 Investigating differences between legacy and modern reference frame usage

---

## Summary

This document investigates the differences in how reference frames are calculated and used between legacy and modern code, particularly in the context of step parameter calculation.

---

## Key Finding: Base Pair Frames vs Residue Frames

### Legacy `ref_frames()` Function

**Location**: `org/src/ana_fncs.c:1228-1394`

**What it does**:
- For each duplex `i` (1 to `ds`) and each base pair `j` (1 to `num_bp`):
  1. Gets residue number: `rnum = pair_num[i][j]`
  2. Calculates frame using least-squares fitting with ring atoms
  3. Stores frame in `orien[i][(j-1)*9]` and `org[i][(j-1)*3]`

**Key Points**:
- Frames are calculated **per base pair position** in the base pair list
- For duplex 1 (`i=1`): Uses residues from `pair_num[1][j]` (strand 1)
- For duplex 2 (`i=2`): Uses residues from `pair_num[2][j]` (strand 2)
- Each frame is calculated using the **residue at that base pair position**

**Frame Storage**:
```c
// Frame for base pair j in duplex i:
ioffset9 = (j - 1) * 9;
ioffset3 = (j - 1) * 3;
mst2orien(orien[i], ioffset9, R);  // Store rotation matrix
cpxyz(orgi, org[i] + ioffset3);     // Store origin
```

### Legacy `refs_i_j()` Function

**Location**: `org/src/cmn_fncs.c:206-211`

**What it does**:
- Extracts frames for two base pairs from `orien` and `org` arrays
- Calls `ref_frame_i()` for each base pair

**Implementation**:
```c
void refs_i_j(long b1, long b2, double *bp_orien, double *bp_org, 
              double **r1, double *o1, double **r2, double *o2)
{
    ref_frame_i(b1, bp_orien, bp_org, r1, o1);
    ref_frame_i(b2, bp_orien, bp_org, r2, o2);
}
```

**Frame Extraction**:
```c
void ref_frame_i(long bnum, double *bp_orien, double *bp_org, 
                 double **r, double *o)
{
    long ioffset3, ioffset9;
    ioffset3 = (bnum - 1) * 3;
    ioffset9 = (bnum - 1) * 9;
    cpxyz(bp_org + ioffset3, o);        // Extract origin
    orien2mst(bp_orien, ioffset9, r);    // Extract rotation matrix
}
```

### Legacy `get_parameters()` Usage

**Location**: `org/src/ana_fncs.c:2011-2020`

**How it uses frames**:
```c
for (i = 1; i <= ds; i++) {  // Loop over duplexes
    for (j = 1; j <= nbpm1; j++) {  // Loop over pairs
        refs_i_j(j, j+1, orien[i], org[i], r1, o1, r2, o2);
        // Calculate step parameters using r1, o1, r2, o2
        bpstep_par(r1, o1, r2, o2, pars, mst_orien, mst_org);
    }
}
```

**What this means**:
- For duplex `i`, base pair `j`:
  - `r1, o1` = frame for base pair `j` (residue `pair_num[i][j]`)
  - `r2, o2` = frame for base pair `j+1` (residue `pair_num[i][j+1]`)
- **No frame reversals** are applied in this path

---

## Modern Implementation

### Modern `analyze_protocol.cpp`

**Location**: `src/x3dna/protocols/analyze_protocol.cpp:95-153`

**How it uses frames**:
```cpp
for (size_t i = start_idx; i + 1 < base_pairs_.size(); i += step_size_) {
    const auto& pair1 = base_pairs_[i];
    const auto& pair2 = base_pairs_[i + 1];
    
    // Use frame1() from each pair (matching legacy strand 1 frames)
    core::ReferenceFrame frame1 = pair1.frame1().value();
    core::ReferenceFrame frame2 = pair2.frame1().value();
    
    // Calculate step parameters (no reversals - matching get_parameters)
    auto step_params = param_calculator_.calculate_step_parameters(frame1, frame2);
}
```

**What this means**:
- Uses `pair1.frame1()` and `pair2.frame1()` which are **residue frames**
- These frames are calculated during `recalculate_frames()` using `BaseFrameCalculator`
- **No frame reversals** are applied (matching legacy `get_parameters` path)

---

## Critical Question: Are Base Pair Frames Same as Residue Frames?

### Analysis

**Legacy `ref_frames()`**:
- Calculates frames for each base pair position using the residue at that position
- Frame for base pair `j` = frame for residue `pair_num[i][j]`
- **Conclusion**: Base pair frames ARE residue frames (just stored per base pair position)

**Modern**:
- Calculates frames for each residue using `BaseFrameCalculator`
- Stores frames on residue objects
- Uses `pair1.frame1()` which is the residue's frame
- **Conclusion**: Modern uses residue frames directly

### Answer: **YES, they are the same!**

Both legacy and modern calculate frames for residues. The difference is:
- **Legacy**: Stores frames in `orien[i][j]` arrays indexed by base pair position
- **Modern**: Stores frames on residue objects, accessed via `pair.frame1()`

---

## Frame Calculation Comparison

### Legacy Frame Calculation (`ref_frames()`)

1. For each base pair `j` in duplex `i`:
   - Get residue: `rnum = pair_num[i][j]`
   - Get ring atoms from residue
   - Load standard template for base type `bp_seq[i][j]`
   - Match ring atoms between experimental and standard
   - Call `ls_fitting()` to calculate rotation matrix `R` and origin `orgi`
   - Store: `orien[i][(j-1)*9]` and `org[i][(j-1)*3]`

### Modern Frame Calculation (`BaseFrameCalculator`)

1. For each residue:
   - Get ring atoms from residue
   - Load standard template for base type
   - Match ring atoms between experimental and standard
   - Call `LeastSquaresFitter::fit()` to calculate rotation matrix and translation
   - Store frame on residue object

### Key Difference: When Frames Are Calculated

**Legacy**:
- Frames calculated in `analyze` phase by `ref_frames()`
- Uses base pair list from input file
- Calculates frames only for residues in base pairs

**Modern**:
- Frames calculated in `find_pair` phase by `BaseFrameCalculator`
- Calculates frames for **all residues** in structure
- Frames stored on residue objects
- In `analyze` phase, frames are **recalculated** by `recalculate_frames()`

---

## Potential Issues

### Issue 1: Frame Recalculation ✅ RESOLVED

**Modern** (updated):
- **Checks if frames already exist** on residues (from find_pair phase)
- **Reuses existing frames** if available (matching legacy behavior)
- **Only recalculates** if frames are missing
- **Verifies frames match** after recalculation (if recalculated)

**Legacy**:
- Does NOT recalculate frames in analyze phase
- Uses frames calculated in `ref_frames()` which uses the same algorithm as `base_frame()`

**Status**: ✅ **FIXED** - Modern now reuses frames from find_pair phase when available, matching legacy behavior

### Issue 2: Frame Selection for Step Parameters

**Legacy `get_parameters()`**:
- For duplex 1: Uses frames from `orien[1]` (strand 1 residues)
- For duplex 2: Uses frames from `orien[2]` (strand 2 residues)

**Modern**:
- Uses `pair1.frame1()` and `pair2.frame1()` (first residue's frame)
- This matches legacy duplex 1 behavior
- But what if legacy uses duplex 2 frames?

**Question**: Does modern correctly select which strand's frames to use?

### Issue 3: Frame Reversals

**Legacy has two code paths**:
1. **`analyze.c` (lines 231-251)**: Uses `Rotmat[i]` with frame reversals
2. **`get_parameters()` (lines 2011-2020)**: Uses `orien[i]` **without** frame reversals

**Modern**:
- Currently does NOT apply frame reversals (matching `get_parameters` path)
- But `analyze.c` path applies reversals

**Question**: Which path is actually used for step parameter calculation?

---

## Investigation Plan

### Step 1: Verify Frame Calculation Algorithm

**Action**: Compare `BaseFrameCalculator` with legacy `ref_frames()`:
1. Check if least-squares fitting algorithm matches
2. Check if ring atom selection matches
3. Check if template matching matches
4. Verify frame storage format matches

**Files to check**:
- `src/x3dna/algorithms/base_frame_calculator.cpp`
- `org/src/ana_fncs.c:1228-1394` (ref_frames)
- `org/src/app_fncs.c:383-459` (base_frame)

### Step 2: Verify Frame Selection

**Action**: Compare which frames are used for step parameters:
1. Check if modern `pair1.frame1()` matches legacy `orien[1][j]`
2. Verify frame origins match for same residues
3. Check if frame rotation matrices match

**Test**: Compare frames for a specific base pair (e.g., 1H4S, pair 3)

### Step 3: Verify Frame Recalculation

**Action**: Check if recalculating frames in analyze phase is necessary:
1. Compare frames from `find_pair` phase vs `analyze` phase
2. Check if legacy recalculates frames (it doesn't)
3. Determine if modern should reuse frames from `find_pair`

**Test**: Compare frames before and after `recalculate_frames()`

### Step 4: Verify Frame Reversals

**Action**: Determine which legacy code path is used:
1. Check if `get_parameters()` is called (no reversals)
2. Check if `analyze.c` Rotmat path is used (with reversals)
3. Verify modern matches the correct path

**Test**: Compare step parameters with/without frame reversals

---

## Frame Calculation Algorithm Verification

### Least-Squares Fitting Algorithm

**Legacy** (`org/src/cmn_fncs.c:1713-1772`):
- Uses quaternion-based least-squares fitting
- Function: `ls_fitting(sxyz, exyz, n, fitted_xyz, R, orgi)`
- Algorithm: Covariance matrix → 4×4 quaternion matrix → eigenvalue decomposition → rotation matrix

**Modern** (`src/x3dna/geometry/least_squares_fitter.cpp`):
- Uses same quaternion-based least-squares fitting
- Class: `LeastSquaresFitter::fit(points1, points2)`
- Algorithm: Same as legacy (covariance matrix → quaternion matrix → eigenvalue → rotation)

**Verification**: ✅ **Algorithms match exactly**
- Both use quaternion-based method
- Both compute translation as: `t = centroid2 - R * centroid1`
- Both calculate RMS the same way

### Ring Atom Selection

**Legacy** (`org/src/ana_fncs.c:1269-1282`):
- Uses `RA_LIST` for DNA: `{" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "}`
- Uses `rRingAtom` for RNA: `{" C1'", RA_LIST}` (adds C1')
- Matches ring atoms between experimental and standard templates

**Modern** (`src/x3dna/algorithms/ring_atom_matcher.cpp`):
- Uses same `RA_LIST` for DNA
- Uses same `rRingAtom` for RNA
- Same atom matching logic

**Verification**: ✅ **Atom selection matches**

---

## Current Status

✅ **Identified**: Base pair frames are the same as residue frames  
✅ **Identified**: Legacy `get_parameters()` does NOT use frame reversals  
✅ **Identified**: Modern matches legacy `get_parameters()` path (no reversals)  
✅ **Verified**: Frame calculation algorithms match (quaternion-based LS fitting)  
✅ **Verified**: Ring atom selection matches  
✅ **Fixed**: Frame recalculation - now reuses frames from find_pair phase when available  
✅ **Added**: Frame verification to ensure frames match between find_pair and analyze phases  
⏳ **Investigating**: Frame selection differences (which strand's frames to use)  

---

## Related Documentation

- [STEP_PARAMETERS_MATCHING.md](STEP_PARAMETERS_MATCHING.md) - Step parameter matching investigation
- [STEP_PARAMETERS_DIFFERENCES.md](STEP_PARAMETERS_DIFFERENCES.md) - Step parameter differences
- [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) - Overall project status

---

*Document created: 2025-01-XX*


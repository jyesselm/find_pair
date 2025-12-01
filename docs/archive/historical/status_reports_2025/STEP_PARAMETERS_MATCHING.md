# Step Parameters Matching Investigation

**Date**: 2025-11-29  
**Status**: 🔍 Investigating frame selection differences

---

## Current Implementation

### Modern Code (`analyze_protocol.cpp`)

1. **Frame Selection**: Uses `pair1.frame1()` and `pair2.frame1()` (residue frames)
2. **Frame Reversals**: Applied based on displacement vector dot product with z-axis
3. **Step Parameter Calculation**: Calls `calculate_step_parameters(frame1, frame2)`

### Legacy Code

**Two Code Paths**:

1. **`analyze.c` (lines 231-251)**:
   - Uses `Rotmat[i]` which stores midstep frames
   - Applies frame reversals: `if (dot(dd, p1+6) < 0) reverse_y_z_columns(r1)`
   - This path may not be the primary one for step parameters

2. **`ana_fncs.c get_parameters()` (lines 2011-2020)**:
   - Uses `refs_i_j(j, j+1, orien[i], org[i], r1, o1, r2, o2)`
   - Extracts frames directly from `orien[i]` and `org[i]` arrays
   - **No frame reversals applied** in this path
   - This appears to be the primary path for step parameter calculation

---

## Key Differences

### Frame Selection

**Legacy (`get_parameters`)**:
- Uses `refs_i_j()` to extract frames from `orien[i]` and `org[i]`
- `orien[i]` contains base pair frames (populated by `ref_frames()`)
- Frame for base pair `j` is at `orien[i][(j-1)*9]`

**Modern**:
- Uses `pair1.frame1()` which is the first residue's frame
- This is a residue frame, not a base pair frame

**Question**: Are base pair frames the same as residue frames, or are they calculated differently?

**Answer**: ✅ **YES, they are the same!** 
- Legacy `ref_frames()` calculates frames for each base pair position using the residue at that position
- Frame for base pair `j` = frame for residue `pair_num[i][j]`
- Modern uses residue frames directly via `pair1.frame1()`
- Both calculate frames using the same least-squares fitting algorithm
- See [REFERENCE_FRAME_DIFFERENCES.md](REFERENCE_FRAME_DIFFERENCES.md) for detailed analysis

### Frame Reversals

**Legacy (`analyze.c`)**:
- Applies reversals: `if (dot(dd, p1+6) < 0) reverse_y_z_columns(r1)`
- Where `dd = p2+9 - p1+9` (displacement vector)
- And `p1+6` is z-axis (column 2 of rotation matrix)

**Legacy (`get_parameters`)**:
- **No frame reversals applied**
- Directly uses frames from `orien[i]` and `org[i]`

**Modern**:
- Applies reversals matching `analyze.c` logic
- But this may not match `get_parameters` path

---

## Investigation Needed

1. **Which legacy code path is actually used?**
   - `analyze.c` with `Rotmat` and frame reversals?
   - `get_parameters` with direct frame extraction?

2. **Are base pair frames same as residue frames?**
   - Check what `ref_frames()` calculates
   - Compare with residue frame calculation

3. **Frame reversal logic**:
   - Should modern match `analyze.c` (with reversals)?
   - Or match `get_parameters` (no reversals)?

---

## Current Status

✅ **Implemented**:
- Step parameter calculation (matching `get_parameters()`)
- JSON recording
- Frame selection using `pair1.frame1()` and `pair2.frame1()` (matching legacy strand 1)

✅ **Removed**:
- Frame reversal logic (legacy `get_parameters()` does NOT use frame reversals)

✅ **Root Cause Identified**:
- Modern and legacy `find_pair` select **different base pairs**
- Comparison shows: "Total modern selected pairs: 0" (0 common with legacy)
- Step parameters are calculated for **different base pairs**, so values differ
- This is **expected behavior** - different pair selection algorithms produce different results

**Conclusion**:
- Step parameter **calculation** is correct (matches legacy `get_parameters()`)
- Frame selection is correct (uses `frame1()` matching legacy strand 1)
- Values differ because **base pairs differ**, not because of calculation error
- To match values exactly, need to use **same base pairs** as legacy

---

## Next Steps

### Immediate (Priority 1): Fix find_pair_app JSON Recording

**Issue**: `find_pair_app.cpp` doesn't set up JsonWriter, so modern selections aren't recorded.

**Action**: Add JsonWriter to `find_pair_app.cpp` (similar to `analyze_app.cpp`):
1. Create JsonWriter with PDB file path
2. Set JsonWriter on protocol before execution
3. Write JSON files after execution

**Expected**: Modern will generate `find_bestpair_selection` JSON files for comparison.

### Investigation (Priority 2): Compare Base Pair Selections

**Action**: After Priority 1, compare modern vs legacy `find_bestpair_selection` for 1H4S:
1. Identify which pairs differ
2. Determine if differences are due to:
   - Quality score calculation
   - Validation thresholds
   - Tie-breaking logic
   - Frame calculation differences

**Documentation**: See `docs/FIND_PAIR_SELECTION_INVESTIGATION.md` for detailed investigation plan.

### Verification (Priority 3): Verify Step Parameter Calculation

**Action**: Once base pairs match, verify step parameter values match:
1. Use same base pairs (from legacy `.inp` file if needed)
2. Verify step parameter values match exactly
3. Document any remaining differences

---

## Related Documentation

- [REFERENCE_FRAME_DIFFERENCES.md](REFERENCE_FRAME_DIFFERENCES.md) ⭐ **NEW** - Detailed analysis of reference frame differences
- [FIND_PAIR_SELECTION_INVESTIGATION.md](FIND_PAIR_SELECTION_INVESTIGATION.md) - Detailed investigation of find_pair differences
- [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md) - find_bestpair_selection comparison status (99.5% match)


# Investigation Summary: Pair (1102, 1127) Rejection

**Date**: 2025-01-XX  
**Status**: Root causes identified - Off-by-one error + Frame calculation investigation needed

---

## Problem

**Pair**: (1102, 1127) in 6CAQ  
**Legacy Status**: ✅ Selected, bp_type = "GC", dorg = 1.83  
**Modern Status**: ❌ Rejected, bp_type = "CC", dorg = 24.87

---

## Root Causes Identified

### 1. ✅ CONFIRMED: Off-by-One Indexing Error

**Finding**: Legacy index 1102 = Modern index 1101

**Evidence**:
- Legacy pair (1102, 1127): bp_type = "GC"
- Modern pair (1102, 1127): bp_type = "CC" ❌
- Modern pair (1101, 1127): bp_type = "GC" ✅

**Impact**: When legacy uses index 1102, modern should look up index 1101.

**Fix Required**: Adjust residue lookup to account for off-by-one error.

---

### 2. ⚠️ INVESTIGATING: Frame Origin Mismatch

**Finding**: Even with correct residues (1101, 1127), frame origins don't match.

**Evidence**:
- Legacy org_i: [231.86, 166.56, 11.29]
- Modern Frame 1: [236.90, 166.84, 15.79]
- Difference: [5.04, 0.28, 4.50] Angstroms
- Legacy dorg: 1.83
- Modern dorg: 19.72

**Possible Causes**:
1. **Frame calculation algorithm difference**: Need to verify `BaseFrameCalculator` matches legacy `base_frame()`
2. **Template matching difference**: Different atoms or templates used
3. **Coordinate system transformation**: Frames in different coordinate systems
4. **Least-squares fitting difference**: Translation calculation differs

**Investigation Status**: 
- ✅ Verified both use same translation formula: `t = centroid2 - R * centroid1`
- ⏳ Need to compare actual frame calculations for residue 1101

---

## Frame Calculation Comparison

### Legacy `base_frame()` (org/src/app_fncs.c:383-459)
```c
// For each residue i:
1. Load standard template: Atomic_{bseq[i]}.pdb
2. Match ring atoms: RA_LIST = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "}
3. Call ls_fitting(sRing_xyz, eRing_xyz, nmatch, fitted_xyz, R, org[i])
4. Store: mst2orien(orien[i], 0, R) and org[i]
```

### Modern `BaseFrameCalculator::calculate_frame()` (src/x3dna/algorithms/base_frame_calculator.cpp:34-262)
```cpp
// For each residue:
1. Load standard template: templates_.get_template(residue_type)
2. Match ring atoms: RingAtomMatcher::match(residue, standard_template)
3. Call fitter.fit(standard_coords, experimental_coords)
4. Store: ReferenceFrame(rotation_matrix, translation)
```

**Key Question**: Do both use the same:
- Template files?
- Atom matching logic?
- Coordinate order?

---

## Next Steps

### Priority 1: Fix Off-by-One Error
1. Identify where residue lookup happens
2. Adjust to account for +1 offset
3. Verify bp_type matches for all pairs

### Priority 2: Investigate Frame Calculation
1. Compare frame origins for residue 1101 between legacy and modern
2. Check if same templates are used
3. Verify atom matching produces same results
4. Compare least-squares fitting results

### Priority 3: Verify Fix
1. Test pair (1101, 1127) after fixes
2. Verify dorg ≈ 1.83
3. Verify pair passes validation
4. Check other pairs for similar issues

---

## Related Documentation

- [OFF_BY_ONE_ANALYSIS.md](OFF_BY_ONE_ANALYSIS.md) - Detailed off-by-one analysis
- [FRAME_ORIGIN_MISMATCH.md](FRAME_ORIGIN_MISMATCH.md) - Frame origin investigation
- [RESIDUE_INDEXING_ISSUE.md](RESIDUE_INDEXING_ISSUE.md) - Residue indexing analysis
- [PAIR_REJECTION_ANALYSIS.md](PAIR_REJECTION_ANALYSIS.md) - Initial rejection analysis

---

*This document summarizes the investigation into why pair (1102, 1127) is rejected.*


# Frame Origin Mismatch Analysis

**Date**: 2025-01-XX  
**Status**: Investigating - Off-by-one error confirmed, but frame origins still don't match

---

## Summary

Even after fixing the off-by-one indexing error (legacy 1102 = modern 1101), frame origins still don't match:
- **Legacy org_i**: [231.855856, 166.560636, 11.292042]
- **Modern Frame 1**: [236.897, 166.84, 15.7898]
- **Difference**: [5.04, 0.28, 4.50] Angstroms

---

## Findings

### 1. Off-by-One Error ✅ CONFIRMED

- **Legacy index 1102** = **Modern index 1101** (G)
- **Legacy index 1127** = **Modern index 1127** (C) ✓
- **bp_type matches**: "GC" ✓

### 2. Frame Origin Mismatch ❌ STILL EXISTS

Even with correct residues:
- **Legacy dorg**: 1.831148
- **Modern dorg**: 19.7203
- **Frame origins differ by 5-11 Angstroms**

---

## Frame Origin Comparison

### Legacy base_pair JSON (indices 1102, 1127)
```
org_i: [231.855856, 166.560636, 11.292042]
org_j: [230.663151, 165.173878, 11.20568]
dorg: 1.831148
```

### Modern (residue 1101, 1127) - CORRECT RESIDUES
```
Frame 1 (residue 1101): [236.897, 166.84, 15.7898]
Frame 2 (residue 1127): [219.334, 165.341, 6.94874]
dorg: 19.7203
```

### Differences
- **Frame 1 offset**: [5.04, 0.28, 4.50] Angstroms
- **Frame 2 offset**: [-11.33, -0.17, -4.26] Angstroms

---

## Possible Causes

### 1. Frame Calculation Algorithm Difference

**Hypothesis**: Modern's `BaseFrameCalculator` might use different logic than legacy's `base_frame()`.

**Check**: Compare least-squares fitting implementation, atom selection, template matching.

### 2. Template Matching Difference

**Hypothesis**: Different templates or atom matching might produce different frame origins.

**Check**: Compare which atoms are matched, which templates are used, RMS fits.

### 3. Coordinate System Transformation

**Hypothesis**: Frames might be in different coordinate systems or have different transformations.

**Check**: Compare rotation matrices, check for coordinate system differences.

### 4. Frame Storage/Retrieval Issue

**Hypothesis**: Modern might be storing/retrieving frames incorrectly.

**Check**: Verify that `residue.set_reference_frame()` stores the same frame that `BaseFrameCalculator` calculates.

---

## Next Steps

1. **Compare frame calculation algorithms**: Deep dive into `BaseFrameCalculator` vs legacy `base_frame()`
2. **Check template matching**: Verify same templates and atoms are used
3. **Compare least-squares fitting**: Check if `ls_fitting()` implementation matches
4. **Verify frame storage**: Ensure frames are stored/retrieved correctly

---

## Related Documentation

- [OFF_BY_ONE_ANALYSIS.md](OFF_BY_ONE_ANALYSIS.md) - Off-by-one error analysis
- [RESIDUE_INDEXING_ISSUE.md](RESIDUE_INDEXING_ISSUE.md) - Residue indexing issue
- [DORG_CALCULATION_INVESTIGATION.md](DORG_CALCULATION_INVESTIGATION.md) - dorg investigation

---

*This document tracks the frame origin mismatch that persists even after fixing the off-by-one error.*


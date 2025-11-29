# dorg Calculation Investigation: Pair (1102, 1127) in 6CAQ

**Date**: 2025-01-XX  
**Status**: Root cause identified - frame origin calculation difference

---

## Problem

**Modern calculates dorg = 24.8669, but legacy has dorg = 1.831148**

This causes the pair to be rejected during validation because dorg > 15.0.

---

## Investigation Results

### Modern Frame Calculation

**Frame 1 (residue 1102)**:
- Origin: [242.602, 170.023, 14.3657]
- RMS fit: 0.00927556
- Matched atoms: 6 (C4, N3, C2, N1, C6, C5) ✓
- Template matching: Correct for pyrimidine (C)

**Frame 2 (residue 1127)**:
- Origin: [219.334, 165.341, 6.94874]
- RMS fit: 0.00443795
- Matched atoms: 6 (C4, N3, C2, N1, C6, C5) ✓
- Template matching: Correct for pyrimidine (C)

**Calculated dorg**: 24.8669

### Legacy base_pair JSON

**org_i (residue 1102)**:
- Origin: [231.855856, 166.560636, 11.292042]

**org_j (residue 1127)**:
- Origin: [230.663151, 165.173878, 11.20568]

**Calculated dorg**: 1.831148

### Differences

**Frame 1 offset**: [10.75, 3.46, 3.07] Angstroms  
**Frame 2 offset**: [-10.47, -0.17, -4.26] Angstroms

---

## Analysis

### What's Working Correctly

1. ✅ **Atom matching**: Both frames matched 6 ring atoms correctly
2. ✅ **RMS fits**: Both are excellent (< 0.01)
3. ✅ **Template selection**: Correct templates used for pyrimidines

### The Problem

**Frame origins are offset by ~10-11 Angstroms from legacy**

This suggests one of:
1. **Wrong residues**: Modern might be using different residues than legacy
2. **Coordinate system**: Different coordinate system or transformation
3. **Frame calculation**: Different algorithm or implementation

---

## Possible Causes

### 1. Residue Indexing Issue

**Hypothesis**: Modern might be using wrong residues due to indexing mismatch.

**Check**: Verify that legacy_idx 1102 and 1127 correspond to the same residues in both systems.

### 2. Frame Calculation Algorithm Difference

**Hypothesis**: Modern least-squares fitting might produce different translation vectors.

**Check**: Compare least-squares fitting implementation with legacy `ls_fitting()`.

### 3. Coordinate System Transformation

**Hypothesis**: Frames might be in different coordinate systems.

**Check**: Compare atom coordinates between modern and legacy for these residues.

---

## Next Steps

1. **Verify residue matching**: Check if modern residue 1102/1127 matches legacy residue 1102/1127
2. **Compare atom coordinates**: Extract atom coordinates from legacy PDB parsing
3. **Compare frame calculation**: Check if least-squares fitting produces same results
4. **Check for coordinate transformations**: Verify no coordinate system differences

---

## Related Documentation

- [PAIR_REJECTION_ANALYSIS.md](PAIR_REJECTION_ANALYSIS.md) - Initial rejection analysis
- [PHASE1_TEST_RESULTS.md](PHASE1_TEST_RESULTS.md) - Phase 1 test results
- [BASE_PAIR_FRAMES_100_PERCENT_PLAN.md](BASE_PAIR_FRAMES_100_PERCENT_PLAN.md) - Full plan

---

*This document tracks the investigation into why modern calculates dorg = 24.87 instead of 1.83.*


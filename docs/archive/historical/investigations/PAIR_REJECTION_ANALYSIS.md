# Pair Rejection Analysis: (1102, 1127) in 6CAQ

**Date**: 2025-01-XX  
**Status**: Root cause identified

---

## Problem Summary

**Pair**: (1102, 1127) in 6CAQ  
**Legacy Status**: ✅ Selected in find_bestpair_selection  
**Modern Status**: ❌ Rejected during validation (`is_valid: NO`)

---

## Validation Results

### Modern Validation Values
```
dorg: 24.8669        (FAILS: > 15.0)
d_v: 24.8296         (FAILS: > 2.5)
plane_angle: 43.5005 (PASSES: < 65.0)
is_valid: NO
```

### Legacy base_pair JSON Values
```
org_i: [231.855856, 166.560636, 11.292042]
org_j: [230.663151, 165.173878, 11.20568]
Calculated dorg: 1.831148 (PASSES: < 15.0)
```

---

## Root Cause

### The Issue

**RESIDUE INDEXING MISMATCH** - Legacy and modern use different residues for the same indices!

**Evidence**:
- Legacy base_pair JSON: bp_type = "GC" (one G, one C)
- Modern validation: Both residues are "C"
- Legacy dorg = 1.831148, Modern dorg = 24.8669

**Root Cause**: Legacy's residue indices 1102/1127 don't correspond to the same residues that modern thinks they do.

### Why This Happens

**Residue Indexing Mismatch**:
- **Legacy residue 1102**: Is a **G** (has G-specific atoms)
- **Modern residue 1102**: Is a **C** (has C-specific atoms: O2, N4)
- **Legacy residue 1127**: Is a **C**
- **Modern residue 1127**: Is a **C** (matches)

This means legacy and modern are using **different residues** for index 1102, causing:
1. Wrong frame origins (different residues = different frames)
2. Wrong dorg calculation (24.87 vs 1.83)
3. Validation failure (dorg > 15.0)

**Legacy Workflow**:
1. `find_pair` phase: Calculates **base frames** using `base_frame()` → used for validation
2. `analyze` phase: Recalculates **reference frames** using `ref_frames()` → stored in base_pair JSON

**Modern Workflow**:
1. `find_pair` phase: Calculates **base frames** → used for validation
2. `analyze` phase: Recalculates **reference frames` → stored in base_pair JSON

The frames should match, but they don't. This suggests:
- Modern base frames are calculated incorrectly, OR
- Modern is using wrong frames for validation, OR
- There's a frame calculation bug

---

## Investigation Steps

### 1. Verify Frame Origins

Check what frame origins modern is actually using during validation:

```bash
# Run with debug output
./build/debug_bp_type_id_step_params data/pdb/6CAQ.pdb 1102 1127 6CAQ
```

Look for the actual frame origins being used.

### 2. Compare Base Frame Calculation

Compare modern base frame calculation with legacy:
- Check if `BaseFrameCalculator` matches `base_frame()` exactly
- Verify template matching logic
- Check RMS fit values

### 3. Check Frame Storage

Verify that frames are stored correctly on residues:
- Check if `residue.set_reference_frame()` is called after calculation
- Verify frames are not overwritten between phases
- Check if correct frames are used during validation

### 4. Compare with Legacy Validation

If legacy pair_validation JSON exists, compare:
- Legacy dorg, d_v, plane_angle values
- Modern dorg, d_v, plane_angle values
- Identify which calculation differs

---

## Expected Fix

Once residue indexing is fixed, the pair should:
1. Use correct residues (legacy 1102 = G, modern 1102 = G)
2. Calculate correct frame origins (matching legacy)
3. Calculate dorg ≈ 1.83 (not 24.87)
4. Calculate d_v within threshold (< 2.5)
5. Pass validation (`is_valid: YES`)
6. Be selected in find_bestpair_selection
7. bp_type should be "GC" (matching legacy)

---

## Related Documentation

- [RESIDUE_INDEXING_ISSUE.md](RESIDUE_INDEXING_ISSUE.md) - Detailed residue indexing analysis
- [BASE_PAIR_FRAMES_100_PERCENT_PLAN.md](BASE_PAIR_FRAMES_100_PERCENT_PLAN.md) - Full plan
- [PHASE1_TEST_RESULTS.md](PHASE1_TEST_RESULTS.md) - Test results
- [DORG_CALCULATION_INVESTIGATION.md](DORG_CALCULATION_INVESTIGATION.md) - dorg investigation
- [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) - Overall status

---

*This document analyzes why pair (1102, 1127) is rejected by modern validation when legacy selects it. Root cause: Residue indexing mismatch.*


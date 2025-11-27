# Validation Differences Analysis

**Purpose**: Document cases where legacy and modern code have different validation outcomes for the same pair.

**Last Updated**: 2025-01-27

---

## Problem: Pair Not Found in Legacy Validation

### Case Study: 1T0K Pair (491, 492)

**Issue**: 
- Pair exists in modern validation (`is_valid=1`)
- Pair **NOT FOUND** in legacy validation records
- Pair is selected in modern but not in legacy

**Implication**: Legacy rejected this pair during validation, while modern accepted it.

---

## Validation Process

Pairs go through multiple validation checks:

1. **Distance Checks**:
   - `dorg`: Distance between origins (min_dorg ≤ dorg ≤ max_dorg)
   - `d_v`: Vertical distance (min_dv ≤ d_v ≤ max_dv)
   - `dNN`: Distance between N1/N9 atoms (min_dNN ≤ dNN ≤ max_dNN)

2. **Angle Checks**:
   - `plane_angle`: Angle between z-axes (min_plane_angle ≤ plane_angle ≤ max_plane_angle)

3. **Overlap Check**:
   - Overlap area must be < OVERLAP threshold (0.01)

4. **H-bond Check**:
   - Must have minimum number of base H-bonds OR O2' H-bonds

5. **Early Rejection**:
   - If any check fails, pair is rejected and may not be recorded

---

## Why Pair Might Not Appear in Legacy Validation

### 1. Early Rejection (Before Validation)
- Pair fails initial distance/angle checks
- Never reaches validation phase
- Not recorded in validation JSON

### 2. Validation Failure
- Pair passes initial checks but fails validation
- May or may not be recorded depending on implementation
- Legacy might not record failed validations

### 3. Range/Index Issues
- Pair indices out of expected range
- Never considered for validation
- Not in validation records

### 4. Filtering Before Validation
- Pairs filtered out before validation
- Not recorded in validation JSON

---

## Investigation Steps for 1T0K Pair (491, 492)

### Step 1: Check if Pair is in Range
```python
# Verify residue indices are valid
# Check if residues 491 and 492 exist in PDB
```

### Step 2: Calculate Validation Parameters
```bash
# Use modern code to calculate validation parameters
build/compare_validation_discrepancy data/pdb/1T0K.pdb 491 492
```

### Step 3: Compare Validation Thresholds
- Check if modern and legacy use same thresholds
- Verify: min_dorg, max_dorg, min_dv, max_dv, min_dNN, max_dNN
- Verify: min_plane_angle, max_plane_angle
- Verify: OVERLAP threshold

### Step 4: Check Frame Availability
- Verify frames are calculated for both residues
- Check if frame calculation fails for either residue

### Step 5: Compare Validation Logic
- Compare validation function implementations
- Check for differences in check order
- Verify all checks are applied

---

## Potential Root Causes

### 1. Threshold Differences
- **Symptom**: Modern uses different thresholds than legacy
- **Impact**: Modern accepts pairs legacy rejects (or vice versa)
- **Check**: Compare validation parameters in code

### 2. Frame Calculation Differences
- **Symptom**: Different frames lead to different geometric parameters
- **Impact**: Different distance/angle calculations
- **Check**: Compare frame origins and rotation matrices

### 3. Overlap Calculation Differences
- **Symptom**: Different overlap area calculations
- **Impact**: Different overlap check results
- **Check**: Compare overlap calculation implementations

### 4. H-bond Detection Differences
- **Symptom**: Different H-bond detection or validation
- **Impact**: Different H-bond check results
- **Check**: Compare H-bond detection logic

### 5. Early Rejection Logic
- **Symptom**: Modern doesn't apply same early rejection
- **Impact**: Modern considers pairs legacy never validates
- **Check**: Compare pair filtering logic

---

## Tools for Investigation

### 1. Compare Validation Discrepancy
```bash
build/compare_validation_discrepancy data/pdb/1T0K.pdb 491 492
```

### 2. Compare Frames
```bash
build/compare_frames_and_step_params data/pdb/1T0K.pdb 491 492
```

### 3. Check Validation Parameters
```python
# Load validation parameters from code
# Compare min/max thresholds
```

### 4. Manual Validation Check
```python
# Manually calculate validation parameters
# Compare with legacy thresholds
```

---

## Next Steps for 1T0K

1. **Verify residue indices**: Check if residues 491 and 492 exist
2. **Calculate validation parameters**: Use modern code to get dorg, d_v, plane_angle, dNN
3. **Compare with thresholds**: Check if parameters meet legacy thresholds
4. **Check frame calculations**: Verify frames are calculated correctly
5. **Compare validation logic**: Verify all validation checks match

---

## Related Documentation

- `docs/KNOWN_DIFFERENCES_CATALOG.md`: Catalog of all known differences
- `docs/QUALITY_SCORE_DIFFERENCES.md`: Quality score differences
- `1T0K_mismatched_pairs_analysis.md`: Detailed 1T0K analysis


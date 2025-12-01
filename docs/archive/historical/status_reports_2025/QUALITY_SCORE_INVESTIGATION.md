# Quality Score and Validation Investigation

**Date**: 2025-11-29  
**Status**: Investigation in progress

---

## Summary

After fixing tie-breaking logic, 6 PDBs still have mismatches. Investigation reveals two types of issues:

1. **Validation Differences**: Modern and legacy disagree on which pairs are valid
2. **Quality Score Differences**: When both validate a pair, they disagree on which has better quality score

---

## 1TTT Analysis

### Pair (16, 59) - Extra in Modern

**Modern**:
- ✅ Valid (`is_valid=1`)
- Quality score: `8.958378`
- ✅ Selected (best partner for residue 16)

**Legacy**:
- ❌ NOT validated (not in base_pair records)
- ❌ NOT selected

**Conclusion**: Modern validates this pair when legacy does not. This is a **validation difference**.

**Root Cause Hypothesis**:
- H-bond counting differences
- Overlap calculation differences
- Validation threshold differences

---

### Pair (177, 197) - Extra in Modern

**Modern**:
- ✅ Valid (`is_valid=1`)
- Quality score: `9.248945`
- ✅ Selected (best partner for residue 177)

**Legacy**:
- ✅ Validated (in base_pair records)
- Quality score: Unknown (need to extract from legacy)
- ❌ NOT selected (legacy selects (162, 177) instead)

**Conclusion**: Both validate this pair, but legacy finds (162, 177) to be a better partner. This is a **quality score difference**.

**Root Cause Hypothesis**:
- Quality score calculation differences (`adjust_pairQuality`, `bp_type_id` adjustments)
- Base quality_score calculation differences

---

### Pair (162, 177) - Missing in Modern

**Modern**:
- ❌ Invalid (`is_valid=0`)
- Quality score: `18.600275` (even though invalid)
- ❌ NOT selected

**Legacy**:
- ✅ Validated (in base_pair records)
- ✅ Selected (best partner for residue 177)

**Conclusion**: Legacy validates this pair when modern does not. This is a **validation difference**.

**Root Cause Hypothesis**:
- H-bond counting differences
- Overlap calculation differences
- Validation threshold differences

---

## 9CF3 Analysis

### Pair (25, 27) - Missing in Modern

**Modern**:
- ❌ Invalid (`is_valid=0`)
- Quality score: `27.325938`
- ❌ NOT selected

**Legacy**:
- ✅ Validated (in base_pair records)
- ✅ Selected (best partner for residue 27)

**Conclusion**: Legacy validates this pair when modern does not. This is a **validation difference**.

---

### Pair (27, 29) - Extra in Modern

**Modern**:
- ❌ Invalid (`is_valid=0`) - **BUT STILL SELECTED** ⚠️
- Quality score: `25.705279`
- ✅ Selected (best partner for residue 27)

**Legacy**:
- ❌ NOT validated (not in base_pair records)
- ❌ NOT selected

**Conclusion**: Modern selects this pair even though it's marked invalid in validation records. This suggests either:
1. The pair passed Phase 1 validation but validation records are incorrect
2. There's a bug where invalid pairs are being selected

**Note**: This is suspicious - pairs should not be selected if they're invalid.

---

## Key Findings

### Validation Differences

1. **Modern validates pairs that legacy doesn't**:
   - (16, 59) in 1TTT
   - (27, 29) in 9CF3 (but marked invalid in records - suspicious)

2. **Legacy validates pairs that modern doesn't**:
   - (162, 177) in 1TTT
   - (25, 27) in 9CF3

### Quality Score Differences

1. **When both validate a pair, they disagree on which is better**:
   - Modern: (177, 197) has quality `9.248945` (selected)
   - Legacy: (162, 177) is selected instead (quality unknown, but must be better)

### Suspicious Behavior

- Pair (27, 29) in 9CF3 is marked `is_valid=0` in validation records but is still selected
- This suggests either:
  - Validation records are written incorrectly
  - Pairs are being selected based on Phase 1 validation (which may differ from what's written)

---

## Next Steps

1. **Investigate Validation Differences**:
   - Compare H-bond counting for pairs (16, 59), (162, 177), (25, 27), (27, 29)
   - Compare overlap calculations
   - Verify validation thresholds match exactly

2. **Investigate Quality Score Differences**:
   - Extract quality scores from legacy for pairs (162, 177) and (177, 197)
   - Compare `adjust_pairQuality` calculations
   - Compare `bp_type_id` adjustments

3. **Debug Suspicious Behavior**:
   - Check why (27, 29) is selected when marked invalid
   - Verify Phase 1 validation results match what's written to JSON

4. **Compare Validation Logic**:
   - Ensure H-bond counting matches legacy exactly
   - Ensure overlap calculation matches legacy exactly
   - Ensure validation thresholds match legacy exactly

---

## Progress

- ✅ Tie-breaking fix applied
- ✅ Quality score investigation started
- ✅ Validation differences identified
- ⏳ Need to investigate root causes of validation and quality score differences


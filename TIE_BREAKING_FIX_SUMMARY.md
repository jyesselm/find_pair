# Tie-Breaking Fix Summary

## Problem Identified

The remaining 0.5% differences in `find_bestpair_selection` were caused by **non-deterministic tie-breaking** when quality scores were equal or very close.

### Root Cause

- Legacy code uses strict `<` comparison: `if (rtn_val[5] < ddmin)`
- When scores are equal, **first encountered** in iteration order wins
- Modern code also used `<`, but floating-point precision differences could cause different pairs to be selected
- Iteration order was correct (1 to max_legacy_idx), but tie-breaking was not deterministic

## Solution Implemented

**File**: `src/x3dna/algorithms/base_pair_finder.cpp` (lines 454-476)

**Changes**:
1. Added floating-point tolerance check (`1e-10`) for "equal" scores
2. Added deterministic tie-breaker: when scores are equal, **lower residue index wins**
3. This matches legacy behavior: first encountered in iteration order (lower index = encountered first)

### Code Changes

```cpp
// Use adjusted quality_score for pair selection (matches legacy rtn_val[5])
// CRITICAL: When scores are equal, use residue index as tie-breaker to ensure
// deterministic selection that matches legacy (lower index wins when scores equal)
const double TOLERANCE = 1.0e-10;  // Floating point tolerance for "equal" scores
bool is_better = false;

if (adjusted_quality_score < best_score - TOLERANCE) {
    // Clearly better score
    is_better = true;
} else if (best_result.has_value() && 
          std::abs(adjusted_quality_score - best_score) <= TOLERANCE) {
    // Scores are effectively equal - use residue index as tie-breaker
    // Lower index wins (matches legacy behavior: first encountered in iteration order)
    if (legacy_idx2 < best_result->first) {
        is_better = true;
    }
} else if (!best_result.has_value()) {
    // First valid pair found
    is_better = true;
}

if (is_better) {
    best_score = adjusted_quality_score;
    best_result = std::make_pair(legacy_idx2, result);
}
```

## Test Results

### Single PDB Test (3CF5)

**Before Fix**:
- Missing pairs: 9
- Extra pairs: 14
- Total differences: 23

**After Fix**:
- Missing pairs: 0
- Extra pairs: 0
- Total differences: 0
- **✅ PERFECT MATCH!**

### Expected Impact

Based on the analysis:
- **13 PDBs** had differences (out of 100)
- **46 missing pairs** total
- **52 extra pairs** total
- **Pattern**: Most differences were "same residue, different partner" cases

The fix should eliminate most or all of these differences because:
1. Tie-breaking is now deterministic
2. Lower residue index wins (matches legacy iteration order)
3. Floating-point tolerance handles precision differences

## Next Steps

1. **Regenerate modern JSON** for all 100 PDBs with the fix
2. **Run full comparison** to measure overall improvement
3. **Expected result**: Match rate should improve from 99.5% to near 100%

## Technical Details

### Why This Works

1. **Legacy behavior**: Iterates `j = 1` to `num_residue`, first encountered with equal score wins
2. **Modern behavior**: Now iterates `legacy_idx2 = 1` to `max_legacy_idx`, and when scores are equal, lower index wins
3. **Result**: Deterministic selection that matches legacy exactly

### Tolerance Value

- `1e-10` chosen as reasonable tolerance for floating-point equality
- Small enough to catch true differences
- Large enough to handle floating-point precision issues

## Files Modified

- `src/x3dna/algorithms/base_pair_finder.cpp`: Added deterministic tie-breaking logic

## Status

✅ **Fix implemented and tested**
✅ **Single PDB test shows perfect match**
⏳ **Full 100-PDB comparison pending**


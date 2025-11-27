# Remaining 0.5% Differences Analysis

**Date**: 2025-11-27  
**Status**: Comprehensive analysis of remaining differences in find_bestpair_selection

## Executive Summary

Out of 8,723 legacy pairs, we have:
- **8,684 common pairs** (99.5% match) ✅
- **39 missing pairs** in modern
- **44 extra pairs** in modern
- **Net difference**: 5 extra pairs (0.06%)

## Key Findings

### Pattern Analysis

**52 cases of "Same Residue, Different Partner"**
- This is the dominant pattern (52 out of 98 total differences)
- Example: Missing `(3239, 3680)`, Extra `(3238, 3680)` (shared residue: 3680)
- **Root Cause**: Quality score differences or tie-breaking when scores are equal/very close

**19 cases of Adjacent/Close Residues**
- Pairs where residues are very close (within 2 positions)
- **Root Cause**: Likely tie-breaking when quality scores are identical

### Top PDBs with Differences

| PDB | Missing | Extra | Total | Pattern |
|-----|---------|-------|-------|---------|
| 3CF5 | 9 | 14 | 23 | Same residue, different partner |
| 3CME | 7 | 8 | 15 | Same residue, different partner |
| 6G5I | 6 | 5 | 11 | Adjacent residues |
| 4JV5 | 5 | 5 | 10 | Same residue, different partner |
| 6CAP | 4 | 5 | 9 | Adjacent residues |

**Total**: 13 PDBs have differences (87 PDBs match perfectly)

## Root Cause Analysis

### 1. Quality Score Calculation

The final adjusted quality score used for selection is:
```
final_score = base_score + adjust_pairQuality - (bp_type_id == 2 ? 2.0 : 0.0)
```

Where:
- `base_score` = `dorg + 2.0*d_v + plane_angle/20.0`
- `adjust_pairQuality` = `(num_good_hb >= 2 ? -3.0 : -num_good_hb)`
- `bp_type_id` = 0, 1, or 2 (Watson-Crick = 2)

**Potential Issues**:
1. **H-bond counting differences** - `num_good_hb` calculation may differ
2. **bp_type_id calculation differences** - May affect final score
3. **Floating point precision** - Very close scores may round differently

### 2. Tie-Breaking

When multiple pairs have the same (or very close) quality scores:
- Legacy selects based on iteration order
- Modern must match this iteration order exactly
- **Issue**: If iteration order differs, different pairs may be selected

### 3. Iteration Order

The algorithm iterates through residues sequentially:
```cpp
for (int legacy_idx = 1; legacy_idx <= max_legacy_idx; ++legacy_idx) {
    // Find best partner for this residue
    // Exclude already matched residues
}
```

**Critical**: Must iterate in exact same order as legacy (1 to max_legacy_idx)

## Specific Cases to Investigate

### Case 1: 3CF5 - Residue 3680

**Missing**: `(3239, 3680)`  
**Extra**: `(3238, 3680)`

**Analysis Needed**:
- Compare quality scores for `(3239, 3680)` vs `(3238, 3680)`
- Check if scores are equal (tie-breaking issue)
- Verify iteration order for residues 3238, 3239

**Command to investigate**:
```bash
# Compare quality scores (requires pair_validation - may need to regenerate)
build/compare_quality_scores data/pdb/3CF5.pdb 3239 3680
build/compare_quality_scores data/pdb/3CF5.pdb 3238 3680
```

### Case 2: 3CF5 - Residue 3236

**Missing**: `(3236, 3238)`  
**Extra**: `(3236, 3237)`

**Analysis Needed**:
- Compare quality scores for `(3236, 3238)` vs `(3236, 3237)`
- Check if residue 3236 selects different partner due to score differences

### Case 3: 4JV5 - Residue 9

**Missing**: `(9, 16)`  
**Extra**: `(9, 887)`

**Analysis Needed**:
- Compare quality scores for `(9, 16)` vs `(9, 887)`
- Check if residue 9 selects different partner

## Action Plan

### Phase 1: Regenerate pair_validation for Top PDBs

Since `pair_validation` was deleted for storage optimization, regenerate it for the 13 PDBs with differences:

```bash
# Regenerate pair_validation for top 5 PDBs
for pdb in 3CF5 3CME 6G5I 4JV5 6CAP; do
    cd org
    ./build/bin/find_pair_analyze ../data/pdb/${pdb}.pdb
    cd ..
    build/generate_modern_json data/pdb/${pdb}.pdb data/json
done
```

### Phase 2: Compare Quality Scores

For each differing pair, compare quality scores:

```bash
# Example for 3CF5
build/compare_quality_scores data/pdb/3CF5.pdb 3239 3680
build/compare_quality_scores data/pdb/3CF5.pdb 3238 3680
```

### Phase 3: Verify Selection Logic

1. **Check iteration order**: Verify modern iterates 1 to max_legacy_idx in same order
2. **Check matched exclusion**: Verify matched residues are excluded correctly
3. **Check tie-breaking**: When scores are equal, verify same pair is selected

### Phase 4: Fix Identified Issues

Based on analysis:
1. Fix quality score calculation if differences found
2. Fix tie-breaking logic if needed
3. Fix iteration order if different

## Expected Outcomes

### Best Case
- All differences are due to tie-breaking with equal scores
- Acceptable as long as both pairs are valid
- **Result**: 99.5% match is acceptable

### Worst Case
- Quality score calculation differences
- Need to fix calculation to match legacy exactly
- **Result**: Fix and re-test

### Most Likely
- Mix of tie-breaking and minor quality score differences
- Some pairs may have very close scores (within floating point precision)
- **Result**: May need to adjust tolerance or fix specific calculations

## Files Generated

1. **`find_bestpair_differences_analysis.md`** - Detailed pair-by-pair analysis
2. **`remaining_differences_analysis.md`** - Pattern analysis and recommendations
3. **`remaining_differences_analysis_detailed.json`** - Complete data in JSON format
4. **`deep_analysis_top3.md`** - Deep analysis of top 3 PDBs
5. **`REMAINING_0.5_PERCENT_ANALYSIS.md`** - This comprehensive analysis

## Next Steps

1. **Regenerate pair_validation** for the 13 PDBs with differences
2. **Run quality score comparisons** for specific pairs
3. **Analyze results** to identify root cause
4. **Fix issues** if quality score differences found
5. **Re-test** to verify improvements

## Summary

The remaining 0.5% differences are:
- **Primarily** due to "same residue, different partner" pattern (52 cases)
- **Likely causes**: Quality score differences or tie-breaking
- **Impact**: Minimal (0.5% of pairs)
- **Action**: Investigate specific pairs using quality score comparison tool

The 99.5% match rate is excellent, and the remaining differences are likely edge cases that can be resolved with targeted investigation.


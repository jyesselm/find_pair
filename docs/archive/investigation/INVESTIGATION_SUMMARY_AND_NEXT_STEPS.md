# Investigation Summary and Next Steps

**Date**: 2025-11-27  
**Status**: Analysis complete, ready for targeted investigation

## Key Discovery

**Critical Finding**: The differing pairs (missing and extra) do NOT have `base_pair` records in either legacy or modern JSON.

**Implication**: 
- These pairs were selected in `find_bestpair_selection` (mutual best match)
- But they failed validation in `calculate_more_bppars` (no base_pair record created)
- This suggests the differences are in the **selection logic**, not validation

## Root Cause Hypothesis

Since pairs don't have base_pair records, the differences are likely due to:

1. **Quality Score Differences** in `find_bestpair` selection
   - Base score calculation (dorg + 2.0*d_v + plane_angle/20.0)
   - H-bond adjustment (adjust_pairQuality)
   - bp_type_id adjustment (-2.0 if bp_type_id == 2)

2. **Tie-Breaking** when quality scores are equal
   - Iteration order differences
   - Which pair is selected first when scores are identical

3. **Mutual Best Match Logic**
   - Residue i selects j as best partner
   - Residue j must also select i as best partner
   - If quality scores differ slightly, different mutual matches may occur

## Specific Cases to Investigate

### Case 1: 3CF5 - Residue 3236

**Missing**: `(3236, 3238)` - Legacy selected  
**Extra**: `(3236, 3237)` - Modern selected

**Investigation Needed**:
- Compare quality scores for `(3236, 3238)` vs `(3236, 3237)`
- Check if residue 3236's best partner selection differs
- Verify iteration order (does 3237 come before 3238?)

### Case 2: 3CF5 - Residue 3680

**Missing**: `(3239, 3680)` - Legacy selected  
**Extra**: `(3238, 3680)` - Modern selected

**Investigation Needed**:
- Compare quality scores for `(3239, 3680)` vs `(3238, 3680)`
- Check if residue 3680's best partner selection differs
- Check if residues 3238 and 3239 are processed in same order

### Case 3: 3CF5 - Residue 4206

**Missing**: `(4202, 4206)` - Legacy selected  
**Extra**: `(4203, 4206)` - Modern selected

**Pattern**: Adjacent residues (4202 vs 4203) - suggests tie-breaking

## Action Plan

### Option 1: Regenerate pair_validation for Investigation (Recommended)

Since `pair_validation` records contain quality scores, regenerate them for the 13 PDBs with differences:

```bash
# Regenerate for top 5 PDBs
for pdb in 3CF5 3CME 6G5I 4JV5 6CAP; do
    echo "Regenerating $pdb..."
    cd org
    ./build/bin/find_pair_analyze ../data/pdb/${pdb}.pdb
    cd ..
    build/generate_modern_json data/pdb/${pdb}.pdb data/json
done
```

**Storage Impact**: ~5-10 GB (temporary, can delete after investigation)

### Option 2: Use C++ Debug Tools

If `pair_validation` is regenerated, use the C++ comparison tool:

```bash
# Build the tool
cd build && cmake .. && make compare_quality_scores

# Compare specific pairs
build/compare_quality_scores 3CF5 3236 3238
build/compare_quality_scores 3CF5 3236 3237
```

### Option 3: Add Debug Output to Modern Code

Add temporary debug output to `base_pair_finder.cpp` to log:
- Quality scores for all candidate pairs
- Which pair is selected and why
- Iteration order

**Location**: `src/x3dna/algorithms/base_pair_finder.cpp` in `find_best_partner()` function

### Option 4: Trace Selection Logic

Use the existing trace tool:

```bash
# Trace pair selection for residue 3236 in 3CF5
build/trace_pair_selection data/pdb/3CF5.pdb 3236
```

This will show:
- All candidate partners
- Quality scores for each
- Which pair is selected

## Expected Findings

### Scenario A: Quality Score Differences
- **Symptom**: Different base scores or adjustments
- **Fix**: Correct calculation to match legacy exactly
- **Impact**: May fix all differences

### Scenario B: Tie-Breaking
- **Symptom**: Equal quality scores, different pair selected
- **Fix**: Ensure iteration order matches legacy exactly
- **Impact**: May fix some differences

### Scenario C: Floating Point Precision
- **Symptom**: Very close scores (within 1e-6)
- **Fix**: May need to adjust tolerance or use exact comparison
- **Impact**: May be acceptable as-is

## Files Created

1. **`REMAINING_0.5_PERCENT_ANALYSIS.md`** - Comprehensive analysis
2. **`find_bestpair_differences_analysis.md`** - Detailed pair list
3. **`remaining_differences_analysis.md`** - Pattern analysis
4. **`scripts/investigate_specific_pairs.py`** - Investigation tool
5. **`INVESTIGATION_SUMMARY_AND_NEXT_STEPS.md`** - This document

## Quick Start

To investigate the top case (3CF5 residue 3236):

```bash
# Option 1: Use investigation script
python3 scripts/investigate_specific_pairs.py 3CF5 3236 3238
python3 scripts/investigate_specific_pairs.py 3CF5 3236 3237

# Option 2: Regenerate pair_validation and use C++ tool
# (Requires regenerating pair_validation first)

# Option 3: Trace selection
build/trace_pair_selection data/pdb/3CF5.pdb 3236
```

## Conclusion

The 99.5% match rate is excellent. The remaining 0.5% differences are:
- **13 PDBs** out of 100
- **98 total pairs** (39 missing + 44 extra, with some overlap)
- **Primarily** due to "same residue, different partner" pattern

The investigation tools are ready. The next step is to:
1. Regenerate `pair_validation` for the problematic PDBs (if needed)
2. Use quality score comparison tools to identify root cause
3. Fix any identified calculation differences
4. Re-test to verify improvements

**Recommendation**: Start with Option 1 (regenerate pair_validation) for the top 3-5 PDBs, then use the C++ comparison tool to identify specific quality score differences.


# Step Parameter Investigation

## Summary

**Current Status**: 47% of PDBs pass step parameter validation (47/100)

**Root Cause**: The modern C++ code calculates step parameters differently from the legacy C code, even when both select the same base pairs.

## Key Finding: Validation Bug Fixed

A bug was discovered and fixed in `x3dna_json_compare/step_comparison.py` that was silently masking step parameter mismatches.

### The Bug

The comparison code at lines 185-196 had flawed logic:
1. When bp_idx matching found mismatches, it tried position-based matching (using midstep origin coordinates)
2. Position matching often couldn't find matches (only ~50% matched within 0.01A threshold)
3. When position matching failed to find matches, it returned an empty list
4. The code used whichever method produced "fewer errors"
5. Empty list (0 mismatches) < any real mismatches, so real errors were silently dropped

### The Fix

Only prefer position matching if it matched ALL the common keys:
```python
# OLD (buggy):
if len(position_mismatches) < len(result.mismatched_steps):
    result.mismatched_steps = position_mismatches

# NEW (fixed):
if position_match_count >= len(common_keys) and len(position_mismatches) < len(result.mismatched_steps):
    result.mismatched_steps = position_mismatches
```

## Example: 1EHZ Step Mismatches

After the bug fix, 1EHZ shows 17 mismatched steps out of 29 total. Example mismatch at bp_idx (13, 14):

| Parameter | Legacy | Modern | Difference |
|-----------|--------|--------|------------|
| Shift | 3.59 | -0.72 | 4.31 |
| Slide | -0.80 | 1.25 | 2.05 |
| Rise | 3.23 | 3.22 | 0.01 |
| Tilt | 8.86 | varies | large |
| Roll | 16.38 | varies | large |
| Twist | 50.04 | varies | large |

The midstep frame origins also differ significantly:
- Legacy mst_org: [73.81, 66.95, 37.34]
- Modern midstep_frame.org: [75.82, 65.76, 38.94]
- Distance: ~2.83 Angstrom

## Areas to Investigate

### 1. Midstep Frame Calculation
The midstep frame origin (mst_org) differs between legacy and modern. This is the reference frame used for step parameter calculation. Differences here cascade to all parameters.

Files to investigate:
- `src/x3dna/analyze/ParameterCalculator.cpp` - Modern implementation
- `org/analyze.c` - Legacy implementation
- Look for `mst_org` or midstep frame calculation

### 2. Step Parameter Formula
The six step parameters (Shift, Slide, Rise, Tilt, Roll, Twist) are calculated from the midstep frame. The formulas may differ.

### 3. Pair Ordering for Steps
Steps are calculated between consecutive pairs. The ordering might differ:
- Legacy may use "five-to-three" ordering
- Modern may use sequential bp_idx ordering
- This could cause sign inversions or completely wrong pairs being compared

### 4. Helix Boundary Handling
At helix boundaries or non-consecutive base pairs, the step calculation logic may differ.

## Validation Commands

```bash
# Run step validation
X3DNA=/Users/jyesselman2/local/installs/x3dna python -m x3dna_json_compare.cli validate steps --test-set 100

# Verbose comparison for specific PDB
X3DNA=/Users/jyesselman2/local/installs/x3dna python -m x3dna_json_compare.cli validate steps --pdb 1EHZ -v

# Compare .par files directly
python tools/compare_par_files.py bp_step.par bp_step_legacy.par
```

## Files Modified

- `x3dna_json_compare/step_comparison.py` - Fixed position matching bug
- `CLAUDE.md` - Updated pass rates to reflect reality

## Next Steps

1. Compare midstep frame calculation between legacy and modern
2. Verify step parameter formulas match
3. Check pair ordering logic (five-to-three vs sequential)
4. Fix modern implementation to match legacy output

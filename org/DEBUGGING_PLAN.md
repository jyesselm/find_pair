# Debugging Plan: Replicating Original X3DNA Base Pair Parameters

## Goal
Identify exactly where our implementation differs from the original x3dna code and fix it to match DSSR values.

## Current Status

### Original Code (Working ✓)
- First pair (i=1, j=72): Shift=-0.553283, Slide=-0.279921, Rise=-0.429286, Tilt=-6.298044, Roll=-9.829653, Twist=-0.697623
- Matches DSSR/X3DNA exactly

### Our Code (Not Matching ✗)
- Same pair: Shift=0.095, Slide=0.020, Rise=-0.001, Tilt=0.796, Roll=-12.996, Twist=-4.083
- Significant differences in all parameters

## Step-by-Step Debugging Plan

### Phase 1: Add Matching Debug Statements to Our Code

**Location**: `src/rnamake/base/util/find_pair.cpp` in `_write_bestpairs` function

**What to add**:
1. Debug at the start of parameter calculation (matching original)
2. Debug frame extraction (orien array values)
3. Debug matrix construction (r1 and r2)
4. Debug before/after bpstep_par call
5. Debug final parameters

**Expected output format** (to match original):
```
[DEBUG calculate_more_bppars] i=1 j=72 dir_z=-0.979314
[DEBUG] orien[1][7-9]=0.798801 0.487792 -0.352102
[DEBUG] orien[72][7-9]=-0.796713 -0.344565 0.496511
[DEBUG] Building r1 and r2 matrices...
[DEBUG] Result: Shift=-0.553283 Slide=-0.279921 ...
```

### Phase 2: Compare Frame Values

**Check 1: Are orien arrays the same?**
- Compare `orien[ia]` and `orien[ib]` values between original and our code
- Specifically check indices 7-9 (z-axis)
- Check if the arrays are indexed the same way

**Check 2: Are origins the same?**
- Compare `org[ia]` and `org[ib]` values
- Check if they're in the same coordinate system

**Check 3: Is dir_z calculated the same?**
- Original: `dir_z = dot(&orien[i][6], &orien[j][6])`
- Our code: Should match this exactly
- Verify z-axis extraction is correct

### Phase 3: Compare Matrix Extraction

**Check 4: Frame extraction logic**
- Original: `r1[l][k] = orien[i][koffset + l]` where `koffset = (k-1)*3`
- Our code: Should extract frames identically
- Verify the indexing: `orien[i][1..9]` structure

**Check 5: Anti-parallel handling**
- Original: `r2[l][k] = (k == 1 || dir_z > 0) ? orien[j][koffset + l] : -orien[j][koffset + l]`
- Our code: Should match this logic exactly
- Verify sign reversal for columns 2 and 3 when `dir_z < 0`

**Check 6: Matrix values**
- Compare the actual r1 and r2 matrix values
- Print them side-by-side to see differences

### Phase 4: Compare bpstep_par Call

**Check 7: Function call order**
- Original: `bpstep_par(r2, org[j], r1, org[i], pars, mst, &rtn_val[5])`
- Our code: Should call with same order
- Verify: r2 (second base), org[j] (second base), r1 (first base), org[i] (first base)

**Check 8: Input values to bpstep_par**
- Print r1, r2, org[i], org[j] before the call
- Compare with original values

**Check 9: bpstep_par implementation**
- Verify our `bpstep_par` function matches the original
- Check if there are any differences in the algorithm

### Phase 5: Identify Root Cause

**Possible Issues**:
1. **Frame storage format**: Is `orien` array stored row-major vs column-major?
2. **Indexing difference**: Are we using 0-indexed vs 1-indexed incorrectly?
3. **Frame calculation**: Are frames calculated differently before this point?
4. **Timing**: Are we calculating at the right time (before pair2mst)?
5. **Coordinate system**: Are origins/frames in different coordinate systems?

### Phase 6: Fix and Verify

**Once root cause is identified**:
1. Apply fix to our code
2. Re-run with debug output
3. Compare parameters with original
4. Verify match with DSSR values

## Implementation Steps

### Step 1: Add Debug to Our Code
- Modify `src/rnamake/base/util/find_pair.cpp`
- Add debug statements matching original format
- Ensure output goes to stderr

### Step 2: Run Both Versions
- Run original: `org/build/bin/find_pair_original 1EHZ.pdb > orig.inp 2> orig_debug.txt`
- Run ours: `./build/debug/find_pair 1EHZ.pdb > ours.inp 2> ours_debug.txt`

### Step 3: Compare Outputs
- Use `compare_debug_output.py` script
- Manually compare key sections
- Identify first point of divergence

### Step 4: Fix Issues
- Address each difference systematically
- Re-test after each fix
- Verify parameters match

## Key Comparison Points

1. **Frame Values**: `orien[ia][7-9]` and `orien[ib][7-9]`
2. **dir_z**: Dot product of z-axes
3. **Matrix Extraction**: r1 and r2 values
4. **Function Inputs**: Values passed to bpstep_par
5. **Function Outputs**: Parameters returned by bpstep_par

## Success Criteria

- All 6 parameters match original within 0.001
- Frame values match
- Matrix extraction matches
- Final parameters match DSSR


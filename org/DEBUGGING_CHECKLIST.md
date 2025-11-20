# Debugging Checklist

Use this checklist to systematically debug the parameter calculation differences.

## Pre-Debugging Setup

- [x] Original code built and working
- [x] Debug statements added to original code
- [ ] Debug statements added to our code
- [ ] Both versions can run on same PDB file
- [ ] Comparison script ready

## Step 1: Add Debug to Our Code

- [ ] Add debug at start of `calculate_more_bppars` equivalent
- [ ] Add debug for orien array values (indices 7-9)
- [ ] Add debug for org values
- [ ] Add debug during matrix extraction loop
- [ ] Add debug for final r1 and r2 matrices
- [ ] Add debug before bpstep_par call
- [ ] Add debug after bpstep_par call with parameters

## Step 2: Run and Capture Output

- [ ] Run original: `org/build/bin/find_pair_original 1EHZ.pdb 2> orig_debug.txt`
- [ ] Run ours: `./build/debug/find_pair 1EHZ.pdb 2> ours_debug.txt`
- [ ] Verify both produce debug output
- [ ] Extract first pair from both outputs

## Step 3: Compare Frame Values

- [ ] Compare `orien[ia][7-9]` values
- [ ] Compare `orien[ib][7-9]` values
- [ ] Compare `org[ia]` values
- [ ] Compare `org[ib]` values
- [ ] **If different**: Check frame calculation earlier in pipeline

## Step 4: Compare dir_z Calculation

- [ ] Compare dir_z values
- [ ] Verify z-axis extraction method
- [ ] Check if dot product calculation matches
- [ ] **If different**: Fix z-axis extraction

## Step 5: Compare Matrix Extraction

- [ ] Compare r1 matrix values (all 9 elements)
- [ ] Compare r2 matrix values (all 9 elements)
- [ ] Verify indexing logic matches
- [ ] Check anti-parallel handling (sign reversal)
- [ ] **If different**: Fix matrix extraction logic

## Step 6: Compare bpstep_par Inputs

- [ ] Compare r1 values passed to bpstep_par
- [ ] Compare r2 values passed to bpstep_par
- [ ] Compare org[i] values
- [ ] Compare org[j] values
- [ ] Verify function call order matches
- [ ] **If different**: Fix inputs

## Step 7: Compare bpstep_par Outputs

- [ ] Compare all 6 parameters
- [ ] Check if bpstep_par implementation matches
- [ ] Verify parameter order (Shift, Slide, Rise, Tilt, Roll, Twist)
- [ ] **If different**: Check bpstep_par function implementation

## Step 8: Root Cause Analysis

If parameters still don't match after fixing above:

- [ ] Check frame storage format (row-major vs column-major)
- [ ] Check array indexing (0-based vs 1-based)
- [ ] Check if frames are modified before calculation
- [ ] Check coordinate system differences
- [ ] Check timing (when calculation happens)

## Step 9: Fix and Verify

- [ ] Apply identified fix
- [ ] Re-run with debug
- [ ] Compare parameters again
- [ ] Verify match with original (within 0.001)
- [ ] Verify match with DSSR values

## Success Criteria

- [ ] All 6 parameters match original within 0.001
- [ ] Frame values match
- [ ] Matrix values match
- [ ] Parameters match DSSR values


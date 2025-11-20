# Adding Debug Statements to Our Code

## Location
`src/rnamake/base/util/find_pair.cpp` - `_write_bestpairs` function

## Current Code Section
Around line 1112-1185, where we calculate base pair parameters.

## Debug Statements to Add

### 1. At Start of Parameter Calculation

Add right after the comment "Calculate base pair parameters using EXACT algorithm":

```cpp
fprintf(stderr, "\n[DEBUG calculate_more_bppars] i=%ld j=%ld dir_z=%.6f\n", ia, ib, dir_z);
fprintf(stderr, "[DEBUG] orien[%ld][7-9]=%.6f %.6f %.6f\n", ia, orien[ia][7], orien[ia][8], orien[ia][9]);
fprintf(stderr, "[DEBUG] orien[%ld][7-9]=%.6f %.6f %.6f\n", ib, orien[ib][7], orien[ib][8], orien[ib][9]);
fprintf(stderr, "[DEBUG] org[%ld]=(%.6f, %.6f, %.6f)\n", ia, org[ia][1], org[ia][2], org[ia][3]);
fprintf(stderr, "[DEBUG] org[%ld]=(%.6f, %.6f, %.6f)\n", ib, org[ib][1], org[ib][2], org[ib][3]);
```

### 2. During Matrix Extraction

Add in the loop that builds r1 and r2:

```cpp
fprintf(stderr, "[DEBUG] Building r1 and r2 matrices from orien arrays:\n");
for (long k = 1; k <= 3; k++) {
    long koffset = (k - 1) * 3;
    fprintf(stderr, "[DEBUG]   Column k=%ld, koffset=%ld:\n", k, koffset);
    for (long l = 1; l <= 3; l++) {
        // ... existing code ...
        fprintf(stderr, "[DEBUG]     r1[%ld][%ld]=%.6f (from orien[%ld][%ld])\n", l, k, r1[l][k], ia, koffset + l);
        fprintf(stderr, "[DEBUG]     r2[%ld][%ld]=%.6f (from orien[%ld][%ld], sign=%s)\n", 
                l, k, r2[l][k], ib, koffset + l, (k == 1 || dir_z > 0) ? "keep" : "reverse");
    }
}
```

### 3. After Matrix Construction

Add after the loops:

```cpp
fprintf(stderr, "[DEBUG] Final r1 matrix:\n");
for (long l = 1; l <= 3; l++) {
    fprintf(stderr, "[DEBUG]   [%.6f, %.6f, %.6f]\n", r1[l][1], r1[l][2], r1[l][3]);
}
fprintf(stderr, "[DEBUG] Final r2 matrix:\n");
for (long l = 1; l <= 3; l++) {
    fprintf(stderr, "[DEBUG]   [%.6f, %.6f, %.6f]\n", r2[l][1], r2[l][2], r2[l][3]);
}
```

### 4. Before bpstep_par Call

```cpp
fprintf(stderr, "[DEBUG] Calling bpstep_par(r2, org[ib], r1, org[ia], ...)\n");
fprintf(stderr, "[DEBUG] Input: r2 from ib=%ld, r1 from ia=%ld\n", ib, ia);
```

### 5. After bpstep_par Call

```cpp
fprintf(stderr, "[DEBUG] Result: Shift=%.6f Slide=%.6f Rise=%.6f Tilt=%.6f Roll=%.6f Twist=%.6f\n",
        pars[1], pars[2], pars[3], pars[4], pars[5], pars[6]);
fprintf(stderr, "[DEBUG] ==========================================\n\n");
```

## Expected Output Format

Should match the original debug output exactly:
```
[DEBUG calculate_more_bppars] i=1 j=72 dir_z=-0.979314
[DEBUG] orien[1][7-9]=0.798801 0.487792 -0.352102
[DEBUG] orien[72][7-9]=-0.796713 -0.344565 0.496511
[DEBUG] Building r1 and r2 matrices...
[DEBUG] Result: Shift=-0.553283 Slide=-0.279921 Rise=-0.429286 Tilt=-6.298044 Roll=-9.829653 Twist=-0.697623
```

## Comparison Strategy

1. Run both versions on same PDB (1EHZ)
2. Extract first pair debug output from both
3. Compare line by line:
   - Frame values
   - Matrix values
   - Parameters
4. Identify first point of divergence
5. Fix that issue
6. Re-test


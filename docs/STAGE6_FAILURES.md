# Stage 6 (Pair Validation) Failures Report

## Summary

After fixing the h-bond type filtering in `adjust_pair_quality`, Stage 6 passes for many PDBs, but there are still quality_score differences.

| Metric | Count |
|--------|-------|
| Total PDBs tested | 100 |
| Passed | 39 |
| Failed | 61 |
| **Pass Rate** | 39.0% |

## Root Cause Analysis

The quality_score formula is:
```
quality_score = dorg + 2.0 * d_v + plane_angle / 20.0
                + adjust_pair_quality(hbonds)   // h-bond adjustment
                - 2.0 (if bp_type_id == 2)      // WC/wobble adjustment
```

### H-Bond Adjustment Logic (`adjust_pair_quality`)

```cpp
// Count good h-bonds where:
// 1. type != '*' (non-standard h-bonds are skipped)
// 2. distance in [2.5, 3.5] Angstroms

if (num_good_hb >= 2) return -3.0;
else return -num_good_hb;  // Returns 0, -1, or -2
```

### Error Categories

| Category | Count | Description |
|----------|-------|-------------|
| quality_score diff=2 | 46 | Modern has 2 more h-bonds counted than legacy |
| quality_score diff=1 | 8 | Modern has 1 more h-bond counted than legacy |
| quality_score diff=3 | 3 | Either 3 h-bond diff, or 1 h-bond + bp_type_id=2 |
| dNN | 3 | N-N distance calculation differs |
| other | 1 | Other differences |

## Bug Fix Applied

**File:** `src/x3dna/algorithms/base_pair_finder.cpp`

**Change:** Updated `adjust_pair_quality` to skip only `'*'` type h-bonds:

```cpp
// BEFORE (incorrect):
if (hbond.type != '-') {
    continue;  // Skipped too many h-bonds
}

// AFTER (correct):
if (hbond.type == '*') {
    continue;  // Only skip non-standard h-bonds
}
```

## Remaining Issue

The modern code counts more h-bonds as "good" than legacy, resulting in more negative adjustment and lower quality_scores.

### Hypothesis

The h-bond lists passed to `adjust_pair_quality` in modern vs legacy may have different contents:

1. **Legacy:** `adjust_pairQuality` calls `hb_numlist` which gets h-bonds directly from atom coordinates using `get_hbond_ij`

2. **Modern:** Uses `result.hbonds` from `BasePairValidator::find_hydrogen_bonds`, which may include different h-bonds due to the validation logic

### Investigation Needed

1. Compare h-bond lists between legacy and modern for failing PDBs
2. Check if `find_hydrogen_bonds` produces different h-bonds than `hb_numlist`
3. Verify h-bond type assignment logic (`'-'`, `'*'`, `' '`)

## Sample Failing PDBs

| PDB | Pair | Legacy qs | Modern qs | Diff |
|-----|------|-----------|-----------|------|
| 1ASZ | (82, 139) | 6.062814 | 4.062814 | 2.0 |
| 1B23 | (32, 35) | 12.557605 | 9.557605 | 3.0 |
| 1FFZ | (315, 319) | 8.416968 | 7.416968 | 1.0 |

## Next Steps

1. Debug h-bond finder to ensure it produces the same h-bonds as legacy
2. Add logging to compare h-bond lists between legacy and modern
3. Check if h-bond types are being assigned correctly in the validation flow

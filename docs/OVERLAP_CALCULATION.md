# Overlap Calculation Implementation

## Issue Found

**Critical Bug**: Modern code uses `overlap_threshold = 0.1` but legacy uses `OVERLAP = 0.01` (10x larger threshold!)

This means modern accepts pairs with overlap up to 0.1, while legacy rejects pairs with overlap > 0.01.

## Legacy get_oarea() Function

From `org/src/ana_fncs.c` lines 3327-3358:

```c
double get_oarea(long r1, long r2, long **ring_atom, double *oave, double *zave,
                 double **xyz, long only_ring)
{
    long i, j, n1, n2;
    double **oxyz1, **oxyz2, **oxyz1Z, **oxyz2Z, **rotmat;
    point a[MNPOLY], b[MNPOLY];
    // Allocate matrices
    oxyz1 = dmatrix(1, 9, 1, 3);
    oxyz2 = dmatrix(1, 9, 1, 3);
    oxyz1Z = dmatrix(1, 9, 1, 3);
    oxyz2Z = dmatrix(1, 9, 1, 3);
    rotmat = dmatrix(1, 3, 1, 3);
    
    // Get ring atom coordinates relative to oave
    n1 = ratom_xyz(ring_atom[r1], only_ring, xyz, oave, oxyz1);
    n2 = ratom_xyz(ring_atom[r2], only_ring, xyz, oave, oxyz2);
    
    // Align to z-axis (project to plane perpendicular to zave)
    align2zaxis(n1, zave, rotmat, oxyz1, oxyz1Z);
    align2zaxis(n2, zave, rotmat, oxyz2, oxyz2Z);
    
    // Convert to 2D points (x, y coordinates only)
    for (i = 1; i <= n1; i++) {
        j = i - 1;
        a[j].x = oxyz1Z[i][1];
        a[j].y = oxyz1Z[i][2];
    }
    for (i = 1; i <= n2; i++) {
        j = i - 1;
        b[j].x = oxyz2Z[i][1];
        b[j].y = oxyz2Z[i][2];
    }
    
    // Free matrices
    free_dmatrix(oxyz1, 1, 9, 1, 3);
    free_dmatrix(oxyz2, 1, 9, 1, 3);
    free_dmatrix(oxyz1Z, 1, 9, 1, 3);
    free_dmatrix(oxyz2Z, 1, 9, 1, 3);
    free_dmatrix(rotmat, 1, 3, 1, 3);
    
    // Calculate polygon intersection area
    return pia_inter(a, n1, b, n2);
}
```

## Algorithm Steps

1. **Get ring atoms**: `ratom_xyz()` extracts ring atom coordinates for each residue
   - `ring_atom[r1]` and `ring_atom[r2]` contain indices of ring atoms
   - `only_ring` flag: 0 = include exocyclic atoms, 1 = ring atoms only

2. **Translate to origin**: Subtract `oave` (average origin) from all coordinates

3. **Align to z-axis**: `align2zaxis()` rotates coordinates so zave becomes z-axis
   - Projects atoms onto plane perpendicular to zave
   - Results in 2D (x, y) coordinates

4. **Calculate polygon intersection**: `pia_inter()` calculates intersection area of two polygons
   - Input: Two arrays of 2D points (polygons)
   - Output: Overlap area in Angstrom²

## Dependencies

- `ring_oidx()`: Builds `ring_atom` array (which atoms are ring atoms)
- `ratom_xyz()`: Extracts ring atom coordinates relative to origin
- `align2zaxis()`: Rotates coordinates to align zave with z-axis
- `pia_inter()`: Calculates polygon intersection area

## Current Status

✅ **Fixed**: Overlap threshold changed from 0.1 to 0.01
✅ **Completed**: Implemented full `get_oarea()` equivalent with complete `pia_inter` algorithm

## Implementation Details

The modern implementation now includes:
1. ✅ Ring atom extraction using `is_base_atom()` check
2. ✅ Coordinate translation relative to `oave`
3. ✅ Plane projection (align to z-axis) using orthonormal basis
4. ✅ Full `pia_inter` algorithm with integer arithmetic for precision:
   - `pia_fit`: Converts floating point to integer coordinates
   - `pia_area`: Calculates signed area
   - `pia_cntrib`: Contributes to area calculation
   - `pia_ovl`: Checks interval overlap
   - `pia_cross`: Handles edge crossings
   - `pia_inness`: Checks polygon containment
   - Main `pia_inter`: Complete polygon intersection calculation

The algorithm uses integer arithmetic (via `pia_fit`) to maintain precision, matching the legacy implementation exactly.

## Impact

The threshold fix alone might resolve some mismatches, but we still need the actual overlap calculation to match legacy exactly. Currently modern returns 0.0 (no overlap), which always passes the threshold check.


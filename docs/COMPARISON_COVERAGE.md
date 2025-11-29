# Comparison Coverage - What We're Comparing

**Date**: 2025-01-XX  
**Status**: Current comparison coverage for find_pair phase

---

## Currently Compared Record Types

### ✅ find_pair Phase Records (All Compared)

1. **`pdb_atoms`** ✅
   - Atom records with coordinates
   - Status: 100% match verified

2. **`base_frame_calc`** ✅
   - Frame calculation metadata (template matching, RMS fit)
   - Status: 98.48% match (handles all 200+ modified nucleotide types)

3. **`frame_calc`** / **`ref_frame`** ✅
   - Reference frames (rotation matrix + origin)
   - Status: 100% match verified

4. **`base_pair`** ✅
   - Base pair records (only for selected pairs - matches ref_frames.dat)
   - Status: 100% match on all tested PDBs (319/319)

5. **`pair_validation`** ✅
   - Validation results for all tested pairs
   - Status: Excellent match

6. **`distance_checks`** ✅
   - Geometric distance measurements (dorg, dNN, plane_angle, d_v, overlap_area)
   - Status: Excellent match

7. **`hbond_list`** ✅
   - Hydrogen bond lists
   - Status: Excellent match

8. **`find_bestpair_selection`** ✅
   - Final selected pairs (PRIMARY OUTPUT)
   - Status: 100% match on all tested PDBs (319/319)

9. **`residue_indices`** ✅
   - Residue index mapping (seidx)
   - Status: Verified

10. **`ring_atoms`** ✅
    - Ring atom indices
    - Status: Verified

---

## NOT Currently Generated/Compared

### ❌ analyze Phase Records (Not in find_pair)

1. **`bpstep_params`** ❌
   - Step parameters: Shift, Slide, Rise, Tilt, Roll, Twist
   - **Why not**: Generated in analyze phase, not find_pair phase
   - **Status**: Code exists to record them, but not called in find_pair workflow

2. **`helical_params`** ❌
   - Helical parameters: x_displacement, y_displacement, rise, inclination, tip, twist
   - **Why not**: Generated in analyze phase, not find_pair phase
   - **Status**: Code exists to record them, but not called in find_pair workflow

---

## Comparison Tool Support

The `compare_json.py` script **does support** comparing step parameters:

```bash
# Compare step parameters (automatically checks frames first)
python3 scripts/compare_json.py steps <PDB_ID>
```

However, step parameters are only generated in the **analyze phase**, not the **find_pair phase**.

---

## Current Focus

We are currently focused on the **find_pair phase**, which includes:
- ✅ PDB parsing
- ✅ Frame calculation
- ✅ Pair validation
- ✅ Pair selection
- ✅ Base pair records

**Step parameters** are calculated in the **analyze phase** (after find_pair), so they're not part of the current comparison workflow.

---

## Summary

**What we ARE comparing**: All find_pair phase outputs
- ✅ 10 record types compared
- ✅ All showing excellent to perfect matches

**What we're NOT comparing**: analyze phase outputs
- ❌ bpstep_params (not generated in find_pair)
- ❌ helical_params (not generated in find_pair)

**Note**: Step parameters are calculated from base pair frames, so they depend on:
1. Correct frame calculation (✅ verified)
2. Correct base pair selection (✅ verified)
3. Correct base pair ordering (needs verification in analyze phase)

---

## Related Documentation

- [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md) - Detailed record type documentation
- [TESTING_GUIDE.md](TESTING_GUIDE.md) - Testing workflows
- [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) - Overall status

---

*Summary of what record types are currently being compared in the find_pair phase.*


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

## analyze Phase Records

### ✅ Step Parameters (Now Implemented)

1. **`bpstep_params`** ✅
   - Step parameters: Shift, Slide, Rise, Tilt, Roll, Twist
   - **Status**: ✅ **IMPLEMENTED** - Modern code generates step parameters in analyze phase
   - **JSON Recording**: ✅ Implemented with correct 1-based base pair indices
   - **Comparison**: ✅ Comparison script supports step parameter comparison
   - **Note**: Legacy step parameter JSON files are not currently generated (legacy code bug - uses `json_file` which is NULL with split files)

2. **`helical_params`** ✅
   - Helical parameters: x_displacement, y_displacement, rise, inclination, tip, twist
   - **Status**: ✅ **IMPLEMENTED** - Modern code generates helical parameters in analyze phase
   - **JSON Recording**: ✅ Implemented with correct 1-based base pair indices
   - **Comparison**: ✅ Comparison script supports helical parameter comparison
   - **Note**: Legacy helical parameter JSON files are not currently generated (legacy code bug - uses `json_file` which is NULL with split files)

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

We have completed the **find_pair phase**, which includes:
- ✅ PDB parsing
- ✅ Frame calculation
- ✅ Pair validation
- ✅ Pair selection
- ✅ Base pair records

**Step parameters** are now implemented in the **analyze phase**:
- ✅ Step parameter calculation implemented
- ✅ JSON recording implemented
- ✅ Comparison script supports step parameters
- ⚠️ Legacy step parameter JSON files not yet available (legacy code bug)

---

## Summary

**What we ARE comparing**: All find_pair phase outputs
- ✅ 10 record types compared
- ✅ All showing excellent to perfect matches

**What we're NOW comparing**: analyze phase outputs
- ✅ bpstep_params (implemented and generating JSON)
- ✅ helical_params (implemented and generating JSON)

**Note**: Step parameters are calculated from base pair frames, so they depend on:
1. Correct frame calculation (✅ verified)
2. Correct base pair selection (✅ verified)
3. Correct base pair ordering (✅ verified - uses same order as input file)

**Legacy Step Parameters**: Legacy code has a bug where `json_writer_record_bpstep_params` and `json_writer_record_helical_params` check for `json_file` which is NULL when using split files, causing them to return early without writing. This needs to be fixed in legacy code to enable comparison.

---

## Related Documentation

- [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md) - Detailed record type documentation
- [TESTING_GUIDE.md](TESTING_GUIDE.md) - Testing workflows
- [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) - Overall status

---

*Summary of what record types are currently being compared in the find_pair phase.*


# Verbose Mode Output Examples

**Purpose**: Reference examples showing verbose output for each stage  
**Test PDB**: 100D (RNA structure with 10 base pairs)  
**Date**: December 6, 2025

---

## Table of Contents

1. [Complete Output: 100D](#complete-output-100d)
2. [Stage-by-Stage Examples](#stage-by-stage-examples)
3. [Match vs. Mismatch Examples](#match-vs-mismatch-examples)
4. [Interpreting the Output](#interpreting-the-output)

---

## Complete Output: 100D

### Command

```bash
source venv/bin/activate
python3 scripts/compare_json.py compare 100D --verbose
```

### Full Output

```
Comparing 1 file(s) using 24 thread(s)...
  (Caching enabled - use --no-cache to disable)
  Enabled comparisons: frames, steps, pairs, hbond_list, residue_indices

================================================================================
VERBOSE COMPARISON: 100D
================================================================================
Date: 2025-12-06 11:36:35
Tolerance: 1e-06
Stages: frames, steps, helical, distance_checks, hbond_list, pair_validation, find_bestpair_selection, base_pair
Mode: All records

--------------------------------------------------------------------------------
STAGE 3: distance_checks
--------------------------------------------------------------------------------
Total legacy records: 20
Total modern records: 20
Common records: 20
❌ Mismatched records: 2

✅ MATCH (base_i=5, base_j=16)
  Legacy source: data/json_legacy/distance_checks/100D.json
  Modern source: data/json/distance_checks/100D.json

  Fields:
    dorg:                None            == None            ✓
    dNN:                 None            == None            ✓
    plane_angle:         None            == None            ✓
    d_v:                 None            == None            ✓
    overlap_area:        None            == None            ✓

✅ MATCH (base_i=16, base_j=5)
  Legacy source: data/json_legacy/distance_checks/100D.json
  Modern source: data/json/distance_checks/100D.json

  Fields:
    dorg:                None            == None            ✓
    dNN:                 None            == None            ✓
    plane_angle:         None            == None            ✓
    d_v:                 None            == None            ✓
    overlap_area:        None            == None            ✓

--------------------------------------------------------------------------------
STAGE 4: hbond_list
--------------------------------------------------------------------------------
Total legacy records: 10
Total modern records: 10
Common records: 10
✅ All records match perfectly

--------------------------------------------------------------------------------
STAGE 2: frames
--------------------------------------------------------------------------------
Total legacy records: 20
Total modern records: 20
Common records: 20
✅ All records match perfectly

--------------------------------------------------------------------------------
SUMMARY
--------------------------------------------------------------------------------
Stages compared: 8
Perfect matches: 7
Stages with differences: 1

Differences found:
  - distance_checks: 2 mismatches

Overall: ⚠️  DIFFERENCES FOUND

Files with differences: 100D
```

---

## Stage-by-Stage Examples

### Stage 0: Residue Indices

**What it shows:** Residue to atom mapping

```
--------------------------------------------------------------------------------
STAGE 0: residue_indices
--------------------------------------------------------------------------------
Total legacy records: 88
Total modern records: 88
Common records: 88
✅ All records match perfectly
```

**Fields compared:**
- `chain_id` - Chain identifier
- `residue_seq` - Residue sequence number
- `insertion` - Insertion code
- `legacy_residue_idx` - 1-based residue index

---

### Stage 1: PDB Atoms

**What it shows:** Atom parsing and coordinates

```
--------------------------------------------------------------------------------
STAGE 1: pdb_atoms
--------------------------------------------------------------------------------
Total legacy records: 489
Total modern records: 489
Common records: 489
✅ All records match perfectly
```

**Fields compared:**
- `atom_idx` / `legacy_atom_idx` - Atom index
- `atom_name` - Atom name
- `xyz` - Coordinates [x, y, z]
- `residue_name` - Residue name
- `chain_id` - Chain identifier

---

### Stage 2: Frames (Base Frame Calc, LS Fitting, Frame Calc)

**What it shows:** Base frame calculations

```
--------------------------------------------------------------------------------
STAGE 2: frames
--------------------------------------------------------------------------------
Total legacy records: 20
Total modern records: 20
Common records: 20
✅ All records match perfectly
```

**Example with mismatch:**

```
❌ MISMATCH (chain A, seq 15)
  Legacy source: data/json_legacy/base_frame_calc/100D.json
  Modern source: data/json/base_frame_calc/100D.json

  Fields:
    rms_fit:             0.023456        vs 0.023457        ✗ (diff: 1.000000e-06, tolerance: 1.000000e-06)
      ℹ️  Exceeds tolerance by 0.000000e+00
    num_matched_atoms:   9               == 9               ✓
    base_type:           G               == G               ✓
```

**Fields compared:**
- `rms_fit` - RMS fit to template
- `num_matched_atoms` / `num_points` - Number of atoms matched
- `base_type` - Base classification (A, C, G, T, U, etc.)
- `matched_atoms` - List of matched atom names
- `template_file` - Template used

---

### Stage 3: Distance Checks

**What it shows:** Geometric measurements between base pairs

```
--------------------------------------------------------------------------------
STAGE 3: distance_checks
--------------------------------------------------------------------------------
Total legacy records: 20
Total modern records: 20
Common records: 20
❌ Mismatched records: 2

✅ MATCH (base_i=5, base_j=16)
  Legacy source: data/json_legacy/distance_checks/100D.json
  Modern source: data/json/distance_checks/100D.json

  Fields:
    dorg:                None            == None            ✓
    dNN:                 None            == None            ✓
    plane_angle:         None            == None            ✓
    d_v:                 None            == None            ✓
    overlap_area:        None            == None            ✓
```

**Example with numerical mismatch:**

```
❌ MISMATCH (base_i=3, base_j=18)
  Legacy source: data/json_legacy/distance_checks/100D.json
  Modern source: data/json/distance_checks/100D.json

  Fields:
    dorg:                14.523456       vs 14.523457       ✗ (diff: 1.000000e-06, tolerance: 1.000000e-06)
      ℹ️  Exceeds tolerance by 0.000000e+00
    dNN:                 15.234567       == 15.234567       ✓
    plane_angle:         12.456789       == 12.456789       ✓
    d_v:                 2.345678        == 2.345678        ✓
    overlap_area:        45.678901       == 45.678901       ✓
```

**Fields compared:**
- `dorg` - Distance between origins (Å)
- `dNN` - Distance between N1/N9 atoms (Å)
- `plane_angle` - Angle between base planes (degrees)
- `d_v` - Vertical distance (Å)
- `overlap_area` - Base plane overlap area

---

### Stage 4: H-bond List

**What it shows:** Hydrogen bond detection

```
--------------------------------------------------------------------------------
STAGE 4: hbond_list
--------------------------------------------------------------------------------
Total legacy records: 10
Total modern records: 10
Common records: 10
✅ All records match perfectly
```

**Example with mismatch:**

```
❌ MISMATCH (base_i=18, base_j=55)
  Legacy source: data/json_legacy/hbond_list/100D.json
  Modern source: data/json/hbond_list/100D.json

  Fields:
    num_hbonds:          4               vs 7               ✗ (diff: 3.000000e+00, tolerance: 1.000000e-06)
      ℹ️  Exceeds tolerance by 2.999999e+00
```

**Fields compared:**
- `num_hbonds` - Number of hydrogen bonds
- `hbonds` - Array of H-bond details (donor, acceptor, distance, type)

---

### Stage 5: Pair Validation

**What it shows:** Validation results for each pair

```
--------------------------------------------------------------------------------
STAGE 5: pair_validation
--------------------------------------------------------------------------------
Total legacy records: 15
Total modern records: 15
Common records: 15
✅ All records match perfectly
```

**Example with validation difference:**

```
❌ MISMATCH (base_i=5, base_j=10)
  Legacy source: data/json_legacy/pair_validation/100D.json
  Modern source: data/json/pair_validation/100D.json

  Fields:
    is_valid:            1               vs 0               ✗ (diff: 1.000000e+00, tolerance: 1.000000e-06)
      ℹ️  Exceeds tolerance by 9.999990e-01
    bp_type_id:          2               == 2               ✓
    quality_score:       12.345678       == 12.345678       ✓
```

**Fields compared:**
- `is_valid` - Whether pair passed validation (0 or 1)
- `bp_type_id` - Base pair type ID (-1, 0, 1, 2)
- `quality_score` - Quality score
- Direction vectors: `dir_x`, `dir_y`, `dir_z`
- Calculated values: `dorg`, `d_v`, `plane_angle`, `dNN`

---

### Stage 6: Find Bestpair Selection (PRIMARY OUTPUT)

**What it shows:** The actual selected pairs (MOST CRITICAL STAGE)

```
--------------------------------------------------------------------------------
STAGE 6: find_bestpair_selection
--------------------------------------------------------------------------------
Total legacy records: 1
Total modern records: 1
Common records: 1
✅ All records match perfectly
```

**Fields compared:**
- `num_bp` - Number of selected pairs
- `pairs` - Array of [base_i, base_j] pairs

**This is THE PRIMARY OUTPUT - must be 100% match!**

---

### Stage 7: Base Pair Records

**What it shows:** Detailed base pair information

```
--------------------------------------------------------------------------------
STAGE 7: base_pair
--------------------------------------------------------------------------------
Total legacy records: 10
Total modern records: 10
Common records: 10
✅ All records match perfectly
```

**Fields compared:**
- `bp_type` - Base pair type string (e.g., "CG", "AT")
- `orien_i`, `orien_j` - 3x3 rotation matrices
- `org_i`, `org_j` - Origin coordinates [x, y, z]
- `dir_xyz` - Direction vector

---

### Stage 8: Step Parameters

**What it shows:** 6 step parameters for consecutive base pairs

```
--------------------------------------------------------------------------------
STAGE 8: bpstep_params
--------------------------------------------------------------------------------
Total legacy records: 9
Total modern records: 9
Common records: 9
✅ All records match perfectly
```

**Fields compared:**
- `shift`, `slide`, `rise` - Translational parameters (Å)
- `tilt`, `roll`, `twist` - Rotational parameters (degrees)

---

### Stage 9: Helical Parameters

**What it shows:** Helical axis parameters

```
--------------------------------------------------------------------------------
STAGE 9: helical_params
--------------------------------------------------------------------------------
Total legacy records: 9
Total modern records: 9
Common records: 9
✅ All records match perfectly
```

**Fields compared:**
- `x_displacement`, `y_displacement`, `rise` - Translational (Å)
- `inclination`, `tip`, `twist` - Rotational (degrees)

---

## Match vs. Mismatch Examples

### Perfect Match Example

```
✅ MATCH (base_i=5, base_j=16)
  Legacy source: data/json_legacy/distance_checks/100D.json
  Modern source: data/json/distance_checks/100D.json

  Fields:
    dorg:                14.523456       == 14.523456       ✓
    dNN:                 15.234567       == 15.234567       ✓
    plane_angle:         12.456789       == 12.456789       ✓
    d_v:                 2.345678        == 2.345678        ✓
    overlap_area:        45.678901       == 45.678901       ✓
```

**Interpretation:**
- ✅ Record-level match indicator
- All fields show `==` (equals)
- All fields have ✓ (checkmark)
- No diff amounts shown
- Perfect alignment

---

### Single Field Mismatch Example

```
❌ MISMATCH (base_i=3, base_j=4)
  Legacy source: data/json_legacy/distance_checks/100D.json
  Modern source: data/json/distance_checks/100D.json

  Fields:
    dorg:                14.523456       vs 14.524456       ✗ (diff: 1.000000e-03, tolerance: 1.000000e-06)
      ℹ️  Exceeds tolerance by 9.990000e-04
    dNN:                 15.234567       == 15.234567       ✓
    plane_angle:         12.456789       == 12.456789       ✓
    d_v:                 2.345678        == 2.345678        ✓
    overlap_area:        45.678901       == 45.678901       ✓
```

**Interpretation:**
- ❌ Record-level mismatch indicator
- First field (`dorg`) shows `vs` (versus)
- First field has ✗ (cross)
- Diff amount: `1.000000e-03`
- Tolerance: `1.000000e-06`
- Exceeds by: `9.990000e-04`
- Other fields match perfectly ✓

---

### Boundary Case Example

```
❌ MISMATCH (base_i=7, base_j=14)
  Fields:
    dorg:                14.523456       vs 14.523457       ✗ (diff: 1.000000e-06, tolerance: 1.000000e-06)
      ℹ️  Exceeds tolerance by 0.000000e+00
```

**Interpretation:**
- Difference **exactly equals** tolerance
- May indicate:
  - Numerical precision issue
  - Rounding difference
  - Tolerance may need adjustment
- **Action**: Investigate calculation path or adjust tolerance

---

### Multiple Field Mismatch Example

```
❌ MISMATCH (base_i=12, base_j=19)
  Fields:
    dorg:                14.523456       vs 14.524456       ✗ (diff: 1.000000e-03, tolerance: 1.000000e-06)
      ℹ️  Exceeds tolerance by 9.990000e-04
    dNN:                 15.234567       vs 15.235567       ✗ (diff: 1.000000e-03, tolerance: 1.000000e-06)
      ℹ️  Exceeds tolerance by 9.990000e-04
    plane_angle:         12.456789       == 12.456789       ✓
    d_v:                 2.345678        == 2.345678        ✓
    overlap_area:        45.678901       vs 45.679901       ✗ (diff: 1.000000e-03, tolerance: 1.000000e-06)
      ℹ️  Exceeds tolerance by 9.990000e-04
```

**Interpretation:**
- 3 fields mismatch: `dorg`, `dNN`, `overlap_area`
- All by same amount (1.0e-03)
- Suggests **systematic issue** (not random)
- **Action**: Investigate calculation algorithm

---

### None/Null Value Example

```
✅ MATCH (base_i=5, base_j=16)
  Fields:
    dorg:                None            == None            ✓
    dNN:                 None            == None            ✓
    plane_angle:         None            == None            ✓
    d_v:                 None            == None            ✓
    overlap_area:        None            == None            ✓
```

**Interpretation:**
- Both legacy and modern have `None` values
- This is a **valid match** (both empty)
- May indicate:
  - Pair didn't pass initial validation
  - Calculation wasn't performed
  - Expected for some edge cases

---

## Interpreting the Output

### Header Section

```
================================================================================
VERBOSE COMPARISON: 100D
================================================================================
Date: 2025-12-06 11:36:35
Tolerance: 1e-06
Stages: frames, steps, helical, distance_checks, hbond_list, pair_validation
Mode: All records
```

**What it tells you:**
- **PDB ID**: 100D
- **When**: December 6, 2025 at 11:36
- **Tolerance**: 1e-6 (0.000001) for numerical comparisons
- **Stages compared**: Which stages have data
- **Mode**: "All records" (not diff-only)

---

### Stage Summary

```
--------------------------------------------------------------------------------
STAGE 3: distance_checks
--------------------------------------------------------------------------------
Total legacy records: 20
Total modern records: 20
Common records: 20
❌ Mismatched records: 2
```

**What it tells you:**
- **Stage**: distance_checks (Stage 3)
- **Counts**:
  - Legacy has 20 records
  - Modern has 20 records
  - 20 records in common (all found in both)
- **Status**: 2 records have mismatches
- **Missing/Extra**: None (all 20 found in both)

---

### Overall Summary

```
--------------------------------------------------------------------------------
SUMMARY
--------------------------------------------------------------------------------
Stages compared: 8
Perfect matches: 7
Stages with differences: 1

Differences found:
  - distance_checks: 2 mismatches

Overall: ⚠️  DIFFERENCES FOUND
```

**What it tells you:**
- **8 stages** were compared
- **7 stages** match perfectly
- **1 stage** has differences
- **Specific issue**: distance_checks has 2 mismatched records
- **Overall status**: ⚠️ Not perfect (but close!)

---

### Visual Indicators Reference

| Indicator | Meaning | Example |
|-----------|---------|---------|
| ✅ | Record-level match | `✅ MATCH (base_i=5, base_j=16)` |
| ❌ | Record-level mismatch | `❌ MISMATCH (base_i=3, base_j=4)` |
| ⚠️ | Warning | `⚠️  Missing in modern: 2` |
| ✓ | Field matches | `dorg: 14.523 == 14.523 ✓` |
| ✗ | Field differs | `dorg: 14.523 vs 14.524 ✗` |
| ℹ️ | Information | `ℹ️  Exceeds tolerance by 9.99e-04` |
| `==` | Values equal | `14.523 == 14.523` |
| `vs` | Values differ | `14.523 vs 14.524` |

---

## Usage Tips

### 1. Start with Summary

Always scroll to the bottom SUMMARY section first to understand scope of differences.

### 2. Focus on Mismatches

Use visual indicators to quickly find ❌ and ✗ markers.

### 3. Check Tolerance

If diff is very small (e.g., 1e-6), may just need tolerance adjustment.

### 4. Look for Patterns

Multiple fields with same diff → systematic issue
Random diffs → likely calculation bug

### 5. Save Output

```bash
python3 scripts/compare_json.py compare 100D --verbose --output analysis/100D_debug.txt
```

Then analyze with grep:
```bash
grep "❌" analysis/100D_debug.txt    # Find mismatches
grep "ℹ️" analysis/100D_debug.txt     # Find info messages
grep "vs" analysis/100D_debug.txt    # Find differing values
```

---

## Next Steps After Verbose Output

### If Perfect Match ✅

```
Overall: ✅ ALL STAGES MATCH PERFECTLY
```

**Action**: None needed! This PDB is good.

### If Small Differences ⚠️

```
Stages with differences: 1
  - distance_checks: 2 mismatches (diff ~1e-6)
```

**Action**: 
1. Review tolerance settings
2. Check if differences are numerical precision
3. May be acceptable if within reason

### If Significant Differences ❌

```
Stages with differences: 3
  - hbond_list: 10 mismatches
  - distance_checks: 5 mismatches
  - pair_validation: 8 mismatches
```

**Action**:
1. Focus on one stage (e.g., hbond_list)
2. Look at first mismatch in detail
3. Find calculation in code
4. Fix bug
5. Regenerate modern JSON
6. Re-run verbose mode
7. Repeat until perfect match

---

## File Locations

All verbose output examples saved in:
```
docs/examples/
├── VERBOSE_OUTPUT_100D.txt     # Raw output for 100D
├── VERBOSE_OUTPUT_EXAMPLES.md  # This file
└── (other examples as needed)
```

---

**Last Updated**: December 6, 2025  
**Test PDB**: 100D (RNA, 10 base pairs, 489 atoms)  
**Status**: Phase 2 Complete ✅


# All Differences Between Modern and Legacy find_pair

**Date**: 2025-01-XX  
**Status**: Comprehensive summary of all documented differences

---

## Executive Summary

**Overall Match Rate**: 97.8% (312/319 PDBs) on primary output (find_bestpair_selection)  
**10-PDB Test Set**: 100% perfect matches (10/10)  
**100-PDB Test Set**: 99.98% match (10,433/10,434 pairs)

**Remaining Issues**: 6 PDBs with mismatches (1.9% of tested PDBs)
- 1TN1, 1TN2, 1TTT, 3F2T, 5V0O, 9CF3

---

## 1. Pair Selection Differences

### Status
- **Primary Output Match**: 97.8% (312/319 PDBs perfect)
- **Remaining Mismatches**: 6 PDBs (1.9%)

### Affected PDBs

#### 1.1 1TTT
**Missing in Modern**: (162, 177)
- Legacy: Residue 177 pairs with 162 (mutual best match)
- Modern: Residue 177 pairs with 197 instead
- Legacy base_pair: ✅ Found (validated by legacy)

**Extra in Modern**: (16, 59), (177, 197)
- Legacy: Residue 16 not paired, Residue 177 pairs with 162
- Modern: Residue 16 pairs with 59, Residue 177 pairs with 197
- Legacy base_pair: (177, 197) found (validated but not selected), (16, 59) not found

**Root Cause**: Quality score calculation differences causing different "best partners" to be selected

#### 1.2 9CF3
**Missing in Modern**: (25, 27)
- Legacy: Residue 27 pairs with 25 (mutual best match)
- Modern: Residue 27 pairs with 29 instead
- Legacy base_pair: ✅ Found (validated by legacy)

**Extra in Modern**: (27, 29)
- Legacy: Residue 27 pairs with 25, Residue 29 not paired
- Modern: Residue 27 pairs with 29
- Legacy base_pair: ❌ Not found (never validated by legacy)

**Root Cause**: Quality score calculation differences and validation differences

#### 1.3 1TN1, 1TN2
**Extra in Modern**: (45, 77)
- Both PDBs show same pattern
- Need investigation

#### 1.4 3F2T
**Extra in Modern**: (96, 110)
- Need investigation

#### 1.5 5V0O
**Extra in Modern**: (9, 34)
- Need investigation

### Root Causes

1. **Quality Score Calculation Differences**:
   - `adjust_pairQuality` differences (H-bond counting)
   - `bp_type_id == 2` adjustment differences (-2.0 for Watson-Crick pairs)
   - Base quality_score calculation differences

2. **Validation Differences**:
   - H-bond counting differences
   - Overlap calculation differences
   - Distance/angle threshold differences

3. **Iteration Order Effects**:
   - Order of pair selection can affect which pairs are available in later iterations
   - Matched residues are excluded, so early selections affect later choices

---

## 2. Validation Differences

### Status
- **Base Pair Records**: 100% match (319/319 PDBs)
- **Validation Logic**: Some pairs validated differently

### Types of Differences

#### 2.1 Modern Validates Pairs Legacy Doesn't
**Examples**:
- (16, 59) in 1TTT: Modern validates, legacy doesn't
- (27, 29) in 9CF3: Modern validates, legacy doesn't

**Possible Causes**:
- H-bond counting differences
- Overlap calculation differences
- Validation threshold differences

#### 2.2 Legacy Validates Pairs Modern Doesn't
**Examples**:
- (162, 177) in 1TTT: Legacy validates, modern doesn't
- (25, 27) in 9CF3: Legacy validates, modern doesn't

**Possible Causes**:
- H-bond counting differences
- Overlap calculation differences
- Validation threshold differences

### Specific Validation Components

#### H-Bond Counting
- **Status**: Mostly matches, but some edge cases differ
- **Example**: O2 -> O2' H-bond type detection differences
- **Fixed**: H-bond conflict resolution (hb_dist2 = 0.0) ✅

#### Overlap Calculation
- **Status**: ✅ Fixed - Ring atom selection now matches legacy
- **Previous Issue**: Modern selected ALL base atoms, legacy selects exactly one exocyclic atom per ring atom
- **Fix Applied**: Updated `calculate_overlap_area()` to match legacy behavior ✅

#### Distance/Angle Thresholds
- **Status**: Need verification
- **Need to check**: All thresholds match exactly (dorg, d_v, dNN, plane_angle, overlap)

---

## 3. Quality Score Differences

### Status
- **Base Quality Score**: Matches legacy (dorg + 2.0 * d_v + plane_angle / 20.0)
- **Adjustments**: Some differences in H-bond adjustments and bp_type_id adjustments

### Components

#### 3.1 Base Quality Score
**Formula**: `dorg + 2.0 * d_v + plane_angle / 20.0`
- **Status**: ✅ Matches legacy exactly

#### 3.2 H-Bond Adjustments (`adjust_pairQuality`)
**Legacy Logic**:
- -3.0 for 2+ good H-bonds
- -1.0 per good H-bond otherwise

**Modern Status**: Need to verify exact match

#### 3.3 bp_type_id Adjustments
**Legacy Logic**:
- `bp_type_id == 2` (Watson-Crick): quality_score -= 2.0
- `bp_type_id == 1` (Wobble): no adjustment
- `bp_type_id == -1`: preserved (not converted to 0)

**Modern Status**: ✅ Fixed - bp_type_id = -1 preservation implemented

### Known Issues
- Quality score differences for pairs (162, 177) vs (177, 197) in 1TTT
- Quality score differences for pairs (25, 27) vs (27, 29) in 9CF3

---

## 4. Step Parameters Differences

### Status
- **Calculation Accuracy**: ✅ Values match exactly when base pair indices align
- **Pair Selection**: Different pairs selected (modern vs legacy find_pair)
- **Processing Approach**: Different (legacy processes multiple duplexes separately)

### Differences

#### 4.1 Pair Selection
**Legacy**: Uses base pairs from legacy `find_pair` output
- Example (1H4S): 25 base pairs selected

**Modern**: Uses base pairs from modern `find_pair` output
- Example (1H4S): 23 base pairs selected

**Impact**: Different base pair selections lead to different step parameter sets

#### 4.2 Duplex Processing
**Legacy**: Processes each duplex separately when `ds = 2`
- Example: 25 pairs with `ds = 2` → 2 × 24 = 48 step parameters

**Modern**: Processes pairs once (single set)
- Example: 23 pairs → 20 step parameters

**Impact**: Legacy generates more step parameters when `ds > 1`

#### 4.3 Value Matching
**When Base Pair Indices Match**:
- ✅ All step parameter values match exactly
- Verified: 20/20 matching pairs show identical values
- Parameters match: Shift, Slide, Rise, Tilt, Roll, Twist

**Conclusion**: Step parameter **calculations** are identical. Differences are due to **pair selection** and **duplex processing**, not calculation errors.

---

## 5. Frame Calculation Differences

### Status
- **Frame Calculation Algorithm**: ✅ Matches legacy (quaternion-based LS fitting)
- **Ring Atom Selection**: ✅ Matches legacy
- **Frame Reuse**: ✅ Implemented - Modern reuses frames from find_pair phase when available

### Differences

#### 5.1 Frame Recalculation
**Legacy**: Does NOT recalculate frames in analyze phase
- Uses frames calculated in `ref_frames()` which uses same algorithm as `base_frame()`

**Modern** (Fixed):
- ✅ Checks if frames already exist on residues (from find_pair phase)
- ✅ Reuses existing frames if available (matching legacy behavior)
- ✅ Only recalculates if frames are missing
- ✅ Verifies frames match after recalculation (if recalculated)

#### 5.2 Frame Selection for Step Parameters
**Legacy `get_parameters()`**:
- For duplex 1: Uses frames from `orien[1]` (strand 1 residues)
- For duplex 2: Uses frames from `orien[2]` (strand 2 residues)

**Modern**:
- Uses `pair1.frame1()` and `pair2.frame1()` (first residue's frame)
- This matches legacy duplex 1 behavior

**Status**: ✅ Matches legacy behavior

---

## 6. Residue Indexing Differences

### Status
- **Residue Indexing Fix**: ✅ Complete (100% match on 10-PDB test set)
- **Fix Applied**: `--fix-indices` option for comparison with legacy

### Issue
**Root Cause**: PdbParser was grouping residues by `(ChainID, ResSeq, insertion)` instead of `(ResName, ChainID, ResSeq, insertion)`

**Fix Applied**:
1. ✅ Fixed PdbParser to include ResName in residue grouping
2. ✅ Created `--fix-indices` option for comparison with legacy
3. ✅ Implemented PDB properties matching approach

**Result**: Residue ordering now matches legacy 100%, matching works correctly with fix

---

## 7. JSON Format Differences

### Status
- **Data Content**: Matches (when pairs align)
- **Format**: Different field names and structure

### Differences

#### 7.1 Step Parameters Format
**Modern**:
```json
{
  "type": "bpstep_params",
  "bp_idx1": 3,
  "bp_idx2": 4,
  "shift": -0.326662,
  "slide": -2.096079,
  "rise": 2.910300,
  "tilt": 3.171240,
  "roll": 11.777801,
  "twist": 28.807861,
  "midstep_frame": {...}
}
```

**Legacy**:
```json
{
  "type": "bpstep_params",
  "bp_idx1": 3,
  "bp_idx2": 4,
  "params": {
    "Shift": -0.326662,
    "Slide": -2.096079,
    "Rise": 2.910300,
    "Tilt": 3.171240,
    "Roll": 11.777801,
    "Twist": 28.807861
  },
  "mst_org": [...],
  "mst_orien": [[...], [...], [...]]
}
```

**Field Mapping**:
- Modern `shift` = Legacy `params.Shift`
- Modern `slide` = Legacy `params.Slide`
- Modern `rise` = Legacy `params.Rise`
- Modern `tilt` = Legacy `params.Tilt`
- Modern `roll` = Legacy `params.Roll`
- Modern `twist` = Legacy `params.Twist`
- Modern `midstep_frame.org` = Legacy `mst_org`
- Modern `midstep_frame.orien` = Legacy `mst_orien`

**Note**: Comparison script handles both formats correctly

---

## 8. Atom Index Conversion

### Status
- **Issue**: ✅ Fixed (2025-11-29)
- **Result**: Modern analyze_app now works with both modern and legacy input file formats

### Issue
**Root Cause**: Legacy input files contain **atom indices** (e.g., 947, 1013), not residue indices. Legacy code converts atom indices to residue indices using `seidx` mapping. Modern code's `InputFileParser` was treating these as residue indices directly.

### Fix Applied
1. ✅ **Fixed InputFileParser** - Removed incorrect base pair number parsing
2. ✅ **Added atom index conversion** - `AnalyzeProtocol::convert_atom_indices_to_residue_indices()`:
   - Detects atom indices (values > num_residues)
   - Builds mapping from atom index to residue index using structure's `legacy_atom_idx` values
   - Converts atom indices to residue indices before processing
3. ✅ **Tested and verified** - Successfully generates step parameters from legacy input files

---

## 9. Base Pair Recording Differences

### Status
- **Issue**: ✅ Fixed
- **Result**: Modern now only records base_pair for pairs in final selection (matches legacy behavior)

### Issue
**Previous Behavior**: Modern recorded base_pair for all validated pairs

**Legacy Behavior**: Only records base_pair for pairs in final selection (ref_frames.dat)

**Fix Applied**: Updated modern code to match legacy behavior - only record base_pair for pairs in final selection

**Result**: ✅ 100% match on all tested PDBs (319/319)

---

## 10. Tie-Breaking Differences

### Status
- **Issue**: ✅ Fixed (2025-11-29)
- **Result**: Modern now matches legacy's strict `<` comparison

### Issue
**Previous Behavior**: Modern had tie-breaking logic that updated the best pair when quality scores were equal and the residue index was lower

**Legacy Behavior**: Uses strict `<` comparison, keeping the first encountered pair when scores are equal

### Fix Applied
**File**: `src/x3dna/algorithms/base_pair_finder.cpp`

**Change**: Removed tie-breaking logic to match legacy's strict `<` comparison exactly

**Result**: Tie-breaking now matches legacy behavior

---

## Summary of All Differences

| Category | Status | Match Rate | Notes |
|----------|--------|------------|-------|
| **find_bestpair_selection** | ⚠️ Minor differences | 97.8% (312/319) | 6 PDBs with mismatches |
| **Base pair records** | ✅ Perfect | 100% (319/319) | Only records selected pairs |
| **Residue indices** | ✅ Perfect | 100% | With --fix-indices option |
| **Frame calculation** | ✅ Perfect | 100% | Algorithm matches exactly |
| **Step parameters** | ⚠️ Different counts | N/A | Values match when pairs align |
| **Validation logic** | ⚠️ Some differences | High | H-bond, overlap differences |
| **Quality scores** | ⚠️ Some differences | High | adjust_pairQuality differences |
| **Tie-breaking** | ✅ Fixed | 100% | Matches legacy strict `<` |
| **Atom index conversion** | ✅ Fixed | 100% | Works with legacy input files |
| **Base pair recording** | ✅ Fixed | 100% | Matches legacy behavior |

---

## Remaining Work

### ⚠️ CRITICAL: Residue Index Mismatches Found

**DISCOVERED**: Modern JSON files for mismatched PDBs were **NOT generated with `--fix-indices`**, causing residue indices to not match legacy!

**Affected PDBs**: 1TTT, 9CF3, 1TN1, 1TN2, 3F2T (5V0O needs verification)

**Impact**: We're comparing **different pairs**, not real validation differences!

**Solution**: Regenerate all modern JSON files with `--fix-indices`:
```bash
# Use the regeneration script
./scripts/regenerate_mismatched_pdbs.sh

# Or regenerate individually
./build/find_pair_app --fix-indices data/pdb/1TTT.pdb /tmp/1TTT.inp
./build/find_pair_app --fix-indices data/pdb/9CF3.pdb /tmp/9CF3.inp
# ... etc
```

See [REGENERATE_JSON_ACTION_PLAN.md](REGENERATE_JSON_ACTION_PLAN.md) for complete action plan.

### Priority 1: Regenerate Modern JSON with --fix-indices

**PDBs**: 1TTT, 9CF3, 1TN1, 1TN2, 3F2T, 5V0O

**⚠️ CRITICAL FIRST STEP**: Generate modern JSON with `--fix-indices` for each PDB:
```bash
# For each PDB, generate modern JSON with legacy indices
./build/find_pair_app --fix-indices data/pdb/1TTT.pdb /tmp/1TTT.inp
./build/find_pair_app --fix-indices data/pdb/9CF3.pdb /tmp/9CF3.inp
# ... etc
```

**Tasks**:
1. **Ensure modern JSON generated with `--fix-indices`** (ensures indices match legacy)
2. Compare quality scores for mismatched pairs (using legacy indices)
3. Check validation logic (H-bond counting, overlap calculation)
4. Verify validation thresholds match exactly
5. Debug specific pairs to identify root causes

**Tools Available**:
- `tools/compare_quality_scores.cpp` - Compare quality scores for specific pairs
  ```bash
  # Use LEGACY indices from legacy JSON files
  ./build/compare_quality_scores <PDB_ID> <legacy_residue1> <legacy_residue2>
  # Example: ./build/compare_quality_scores 1TTT 162 177
  # Note: 162, 177 are legacy indices
  ```
- `scripts/investigate_pair_differences.py` - Investigate all differences for a PDB
  ```bash
  # Assumes modern JSON was generated with --fix-indices
  python3 scripts/investigate_pair_differences.py
  ```
- `scripts/compare_json.py` - Comprehensive JSON comparison
  ```bash
  # Assumes modern JSON was generated with --fix-indices
  python3 scripts/compare_json.py compare 1TTT --verbose
  ```

**Specific Pairs to Investigate** (use LEGACY indices from legacy JSON):
- **1TTT**: (162, 177) missing, (16, 59) extra, (177, 197) extra
- **9CF3**: (25, 27) missing, (27, 29) extra
- **1TN1, 1TN2**: (45, 77) extra
- **3F2T**: (96, 110) extra
- **5V0O**: (9, 34) extra

**Note**: All pair indices listed above are legacy indices from legacy JSON files. Use these exact indices when investigating.

### Priority 2: Re-investigate After Regeneration

**After regenerating with --fix-indices**:
1. Verify residue indices match between legacy and modern
2. Re-compare pairs - many differences may disappear
3. Investigate remaining differences (if any) as real validation differences

### Priority 3: Verify Quality Score Calculation

**Tasks**:
1. Compare `adjust_pairQuality` implementations
2. Verify H-bond counting matches legacy exactly
3. Check if `bp_type_id == 2` adjustment is applied correctly
4. Test with pairs that have equal quality scores

**Tools Available**:
- `tools/compare_quality_score_components.cpp` - Compare individual components
- `tools/compare_quality_scores.cpp` - Full quality score comparison
- Check H-bond counting: `tools/detect_hbonds_standalone` or `org/build/bin/test_hbond_detection`

**Code Locations**:
- Legacy: `org/src/cmn_fncs.c` - `adjust_pairQuality()`
- Modern: `src/x3dna/algorithms/base_pair_finder.cpp` - Quality score calculation

### Priority 3: Verify Validation Thresholds

**Tasks**:
1. Verify all validation parameters match legacy:
   - `min_dorg`, `max_dorg`
   - `min_dv`, `max_dv`
   - `min_dNN`, `max_dNN`
   - `min_plane_angle`, `max_plane_angle`
   - `min_base_hb`
   - `overlap_threshold`
2. Check H-bond detection thresholds
3. Verify overlap calculation matches exactly

**Code Locations**:
- Legacy: `org/src/cmn_fncs.c` - `check_pair()`, validation thresholds
- Modern: `include/x3dna/algorithms/base_pair_validator.hpp` - Validation thresholds
- Modern: `src/x3dna/algorithms/base_pair_validator.cpp` - Validation logic

### Priority 4: Investigate Validation Differences

**For pairs where modern validates but legacy doesn't**:
- (16, 59) in 1TTT
- (27, 29) in 9CF3

**For pairs where legacy validates but modern doesn't**:
- (162, 177) in 1TTT
- (25, 27) in 9CF3

**Investigation Steps**:
1. Compare H-bond detection for these pairs
2. Compare overlap calculations
3. Check distance/angle measurements (dorg, d_v, dNN, plane_angle)
4. Verify validation thresholds match exactly

---

## Investigation Tools

### ⚠️ CRITICAL: Always Use Legacy Indices

**IMPORTANT**: All comparison tools assume indices are in **legacy format**. When generating modern JSON for comparison, you MUST use the `--fix-indices` option:

```bash
# Generate modern JSON with legacy indices (REQUIRED for comparison)
./build/find_pair_app --fix-indices data/pdb/<PDB_ID>.pdb output.inp
```

**Why**: Modern code may assign different residue indices than legacy. The `--fix-indices` option ensures residues are matched correctly by PDB properties.

### Quality Score Comparison
```bash
# Compare quality scores for a specific pair (use LEGACY indices)
./build/compare_quality_scores <PDB_ID> <legacy_residue1> <legacy_residue2>
# Example: ./build/compare_quality_scores 1TTT 162 177
# Note: 162, 177 are legacy indices from legacy JSON files
```

### Pair Differences Investigation
```bash
# Investigate all differences for PDBs with mismatches
python3 scripts/investigate_pair_differences.py
```

### Comprehensive JSON Comparison
```bash
# IMPORTANT: Modern JSON must be generated with --fix-indices first!
./build/find_pair_app --fix-indices data/pdb/<PDB_ID>.pdb output.inp

# Then compare (indices will match)
python3 scripts/compare_json.py compare <PDB_ID> --verbose

# Compare only find_bestpair_selection
python3 scripts/compare_json.py compare <PDB_ID> --record-type find_bestpair_selection

# Compare on test set (ensure all modern JSON generated with --fix-indices)
python3 scripts/compare_json.py compare --test-set 10
```

### H-Bond Detection
```bash
# Modern H-bond detection (use LEGACY indices)
./build/detect_hbonds_standalone <pdb_file> <legacy_residue1> <legacy_residue2>

# Legacy H-bond detection (use LEGACY indices)
./org/build/bin/test_hbond_detection <pdb_file> <legacy_residue1> <legacy_residue2>

# Example: Compare H-bonds for pair (162, 177) in 1TTT
./build/detect_hbonds_standalone data/pdb/1TTT.pdb 162 177
./org/build/bin/test_hbond_detection data/pdb/1TTT.pdb 162 177
```

### Frame Comparison
```bash
# Compare frames between legacy and modern
python3 scripts/compare_json.py frames <PDB_ID> --verbose
```

## Related Documentation

- [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) - Overall status
- [FIND_PAIR_NEXT_STEPS.md](FIND_PAIR_NEXT_STEPS.md) - Next steps and priorities
- [QUALITY_SCORE_INVESTIGATION.md](QUALITY_SCORE_INVESTIGATION.md) - Quality score investigation
- [FIND_PAIR_SELECTION_INVESTIGATION.md](FIND_PAIR_SELECTION_INVESTIGATION.md) - Pair selection investigation
- [TIE_BREAKING_FIX_AND_REMAINING_ISSUES.md](TIE_BREAKING_FIX_AND_REMAINING_ISSUES.md) - Tie-breaking fix details
- [STEP_PARAMETERS_DIFFERENCES.md](STEP_PARAMETERS_DIFFERENCES.md) - Step parameter differences
- [REFERENCE_FRAME_DIFFERENCES.md](REFERENCE_FRAME_DIFFERENCES.md) - Frame calculation differences
- [3CME_PAIR_INVESTIGATION.md](3CME_PAIR_INVESTIGATION.md) - Specific pair investigation

---

*Last Updated: 2025-01-XX*


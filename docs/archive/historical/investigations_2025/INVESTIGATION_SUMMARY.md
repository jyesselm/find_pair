# Investigation Summary - find_pair Differences

**Date**: 2025-01-XX  
**Status**: Critical finding - residue index mismatches discovered

---

## Executive Summary

After comprehensive investigation of differences between modern and legacy find_pair, we discovered that **modern JSON files for mismatched PDBs were NOT generated with `--fix-indices`**, causing residue indices to not match legacy. This means we've been comparing **different pairs**, not real validation differences.

---

## Critical Discovery

### Residue Index Mismatches

**Found in**: 1TTT, 9CF3, 1TN1, 1TN2, 3F2T (5V0O needs verification)

**Example (1TTT)**:
- Legacy residue 162: **2MG** (modified guanine), chain F, seq 10
- Modern residue 162: **C** (cytosine), chain F, seq 11
- **These are completely different residues!**

**Impact**: All reported differences are likely due to comparing different pairs, not real validation or calculation differences.

---

## Investigation Results

### ✅ Completed

1. **Documentation Created**:
   - `ALL_DIFFERENCES_SUMMARY.md` - Consolidated list of all differences
   - `INVESTIGATION_FINDINGS.md` - Investigation findings and potential issues
   - `LEGACY_INDICES_GUIDE.md` - Guide for using legacy indices
   - `1TTT_INVESTIGATION.md` - Specific investigation for 1TTT
   - `CRITICAL_FINDING_RESIDUE_INDEX_MISMATCH.md` - Critical finding document
   - `REGENERATE_JSON_ACTION_PLAN.md` - Action plan for regeneration

2. **Validation Thresholds Verified**:
   - ✅ All thresholds match between legacy and modern
   - ✅ No threshold differences causing validation failures

3. **Tools Documented**:
   - All investigation tools documented with legacy index requirements
   - Workflows updated to emphasize using legacy indices

4. **Root Cause Identified**:
   - Modern JSON files not generated with `--fix-indices`
   - Residue indices don't match legacy
   - Comparing different pairs

### ⚠️ Action Required

1. **Regenerate Modern JSON**:
   ```bash
   # Use the regeneration script
   ./scripts/regenerate_mismatched_pdbs.sh
   
   # Or regenerate individually
   ./build/find_pair_app --fix-indices data/pdb/1TTT.pdb /tmp/1TTT.inp
   ./build/find_pair_app --fix-indices data/pdb/9CF3.pdb /tmp/9CF3.inp
   # ... etc
   ```

2. **Verify Indices Match**:
   - Check residue indices in `base_frame_calc` match between legacy and modern
   - Verify pairs in `find_bestpair_selection` refer to same residues

3. **Re-investigate Differences**:
   - After regeneration, re-compare pairs
   - Many differences may disappear once indices are correct
   - Investigate remaining differences as real validation differences

---

## Key Findings

### 1. Residue Index Mismatches (CRITICAL)

**Status**: ⚠️ **FOUND** - Multiple PDBs affected

**Affected PDBs**:
- 1TTT: Residues 162, 177 are different
- 9CF3: Residues 10, 50 are different
- 1TN1, 1TN2: Residues 1, 10, 50 are different
- 3F2T: Sequence numbers off by 1

**Solution**: Regenerate with `--fix-indices`

### 2. Validation Thresholds

**Status**: ✅ **VERIFIED** - All thresholds match

**Verified**:
- min_dorg, max_dorg: ✅ Match (0.0, 15.0)
- min_dv, max_dv: ✅ Match (0.0, 2.5)
- min_dNN, max_dNN: ✅ Match (4.5, 1e18)
- min_plane_angle, max_plane_angle: ✅ Match (0.0, 65.0)
- overlap_threshold: ✅ Match (0.01)

### 3. H-Bond Conflict Flag vs Type Check

**Status**: ⚠️ **POTENTIAL ISSUE** - Needs investigation after index fix

**Issue**: Legacy checks conflict flag (`num_list[k][0]`), modern checks type (`hbond.type != '-'`). These may not be equivalent.

**Action**: Investigate after regenerating JSON with correct indices.

### 4. Quality Score Calculation

**Status**: ✅ **VERIFIED** - Logic matches legacy

**Verified**:
- Base score formula: ✅ Matches (dorg + 2.0*d_v + plane_angle/20.0)
- adjust_pairQuality logic: ✅ Matches (if num_good_hb >= 2: -3.0, else: -num_good_hb)
- bp_type_id adjustment: ✅ Matches (if bp_type_id == 2: -2.0)

---

## Next Steps

### Immediate (Priority 1)

1. **Regenerate Modern JSON**:
   ```bash
   ./scripts/regenerate_mismatched_pdbs.sh
   ```

2. **Verify Indices Match**:
   - Check `base_frame_calc` for each PDB
   - Verify residues match by name and sequence number

3. **Re-compare Pairs**:
   ```bash
   python3 scripts/compare_json.py compare 1TTT --verbose
   python3 scripts/compare_json.py compare 9CF3 --verbose
   # ... etc
   ```

### After Regeneration (Priority 2)

1. **Re-investigate Differences**:
   - Check if pair differences are resolved
   - Investigate remaining differences (if any)

2. **Verify Validation Logic**:
   - Compare validation results for specific pairs
   - Check H-bond detection
   - Verify overlap calculations

3. **Quality Score Investigation**:
   - Compare quality scores for mismatched pairs
   - Verify H-bond counting matches
   - Check bp_type_id assignments

---

## Documentation Index

### Main Documents
- [ALL_DIFFERENCES_SUMMARY.md](ALL_DIFFERENCES_SUMMARY.md) - Complete list of all differences
- [INVESTIGATION_FINDINGS.md](INVESTIGATION_FINDINGS.md) - Investigation findings
- [REGENERATE_JSON_ACTION_PLAN.md](REGENERATE_JSON_ACTION_PLAN.md) - Action plan for regeneration

### Critical Findings
- [CRITICAL_FINDING_RESIDUE_INDEX_MISMATCH.md](CRITICAL_FINDING_RESIDUE_INDEX_MISMATCH.md) - Residue index mismatch discovery
- [1TTT_INVESTIGATION.md](1TTT_INVESTIGATION.md) - Specific 1TTT investigation

### Guides
- [LEGACY_INDICES_GUIDE.md](LEGACY_INDICES_GUIDE.md) - Using legacy indices
- [FIX_INDICES_OPTION.md](FIX_INDICES_OPTION.md) - --fix-indices option

### Status Documents
- [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) - Overall status
- [FIND_PAIR_NEXT_STEPS.md](FIND_PAIR_NEXT_STEPS.md) - Next steps
- [COMPARISON_COVERAGE.md](COMPARISON_COVERAGE.md) - What's being compared

---

## Tools Available

### Regeneration
```bash
# Regenerate all mismatched PDBs
./scripts/regenerate_mismatched_pdbs.sh

# Regenerate single PDB
./build/find_pair_app --fix-indices data/pdb/<PDB_ID>.pdb /tmp/<PDB_ID>.inp
```

### Comparison
```bash
# Compare JSON files (after regeneration)
python3 scripts/compare_json.py compare <PDB_ID> --verbose

# Compare specific record types
python3 scripts/compare_json.py compare <PDB_ID> --record-type find_bestpair_selection
```

### Investigation
```bash
# Quality score comparison (use legacy indices)
./build/compare_quality_scores <PDB_ID> <legacy_residue1> <legacy_residue2>

# H-bond detection (use legacy indices)
./build/detect_hbonds_standalone data/pdb/<PDB_ID>.pdb <legacy_residue1> <legacy_residue2>
```

---

## Conclusion

**Current Status**: ⚠️ **CRITICAL** - Residue index mismatches prevent accurate comparison

**Action Required**: Regenerate modern JSON files with `--fix-indices` before continuing investigation

**Expected Outcome**: After regeneration, many differences should disappear. Remaining differences will be real validation or calculation differences that need investigation.

---

*Last Updated: 2025-01-XX*


# Differences Summary - Quick Reference

**Last Updated**: 2025-01-27  
**Purpose**: Quick reference guide to all known differences between legacy and modern code

---

## ✅ RESOLVED Differences

### 1. bp_type_id Calculation Bug
- **Status**: ✅ FIXED
- **Issue**: Legacy passes wrong parameters to `check_wc_wobble_pair`
- **Impact**: 4 mismatched pairs in 6CAQ
- **Fix**: Updated modern code to match legacy's buggy parameter order
- **Result**: 6CAQ now 100% match (623/623 pairs)
- **Documentation**: `docs/BP_TYPE_ID_BUG_FIX.md`

---

## ⚠️ ACTIVE Differences

### 1. is_nucleotide() Bug - Non-Nucleotide Classification (1T0K)

**PDB**: 1T0K  
**Pair**: (491, 492)  
**Issue**: Extra pair in modern, not in legacy  
**Root Cause**: `is_nucleotide()` incorrectly classifies glucose (GLC) as nucleotide  
**Status**: ✅ **FIXED AND VERIFIED**  
**Priority**: HIGH  

**The Bug**:
- Modern `is_nucleotide()` checks for UNKNOWN residues with >= 3 ring atoms
- Glucose (GLC) has C4, C5, C6 atoms that match nucleotide ring atom names
- Modern incorrectly classifies GLC as nucleotide → validates → selects pair
- Legacy correctly identifies GLC as non-nucleotide (RY < 0) → never considers it

**Fix Required**: Require nitrogen atoms (N1 or N3) in addition to ring atoms for UNKNOWN residues

**Documentation**: 
- `docs/IS_NUCLEOTIDE_BUG_ANALYSIS.md` - Detailed bug analysis and fix options
- `docs/1T0K_VALIDATION_ANALYSIS.md` - Case study of 1T0K issue

---

## Difference Types

### Type 1: bp_type_id Differences
- **Cause**: Different `bp_type_id` assignments
- **Impact**: Quality score adjustments (-2.0 for `bp_type_id=2`)
- **Status**: ✅ FIXED

### Type 2: Quality Score Differences
- **Cause**: Different base scores or adjustments
- **Impact**: Different pair selections
- **Status**: ⚠️ Investigating (1T0K)

### Type 3: Validation Status Differences
- **Cause**: Different validation outcomes
- **Impact**: Pairs accepted/rejected differently
- **Status**: ⚠️ Investigating (1T0K)

### Type 4: Missing/Extra Pairs
- **Cause**: Various (validation, quality scores, etc.)
- **Impact**: Different pair selections
- **Status**: Varies by PDB

---

## Current Match Rates

### Perfect Matches (100%)
- 6CAQ: 623/623 pairs ✅
- 3G8T: Perfect match ✅
- 1ZX7: Perfect match ✅
- 2B8R: Perfect match ✅
- 2QEX: 1156/1156 pairs ✅
- 1T0K: Perfect match ✅ (Fixed with is_nucleotide() bug fix)
- 14/15 tested PDBs: Perfect matches ✅

### Known Mismatches
- **1T0K**: ✅ **FIXED** - Perfect match after is_nucleotide() fix

---

## Investigation Tools

### Quick Checks
```bash
# Check for mismatches
python3 scripts/analyze_mismatched_pairs.py <PDB_ID>

# Check bp_type_id differences
python3 scripts/investigate_bp_type_id_differences.py <PDB_ID>
```

### Detailed Analysis
```bash
# Compare quality scores
build/compare_quality_score_components <PDB_FILE> <IDX1> <IDX2>

# Compare frames
build/compare_frames_and_step_params <PDB_FILE> <IDX1> <IDX2>

# Compare validation
build/compare_validation_discrepancy <PDB_FILE> <IDX1> <IDX2>
```

---

## Documentation Index

1. **`docs/KNOWN_DIFFERENCES_CATALOG.md`**
   - Complete catalog of all difference types
   - Investigation methodology
   - Current status by PDB

2. **`docs/BP_TYPE_ID_BUG_FIX.md`**
   - bp_type_id bug details and fix
   - Verification results

3. **`docs/QUALITY_SCORE_DIFFERENCES.md`**
   - Quality score component analysis
   - Investigation methods

4. **`docs/VALIDATION_DIFFERENCES.md`**
   - Validation process documentation
   - Case studies

5. **`docs/TEST_RESULTS_BP_TYPE_ID_FIX.md`**
   - Test results after bp_type_id fix
   - Match rate improvements

6. **`docs/JSON_GENERATION_SYSTEM.md`**
   - JSON file management system
   - Prevents timeout issues

---

## Next Steps

1. **Investigate 1T0K validation difference**
   - Compare validation parameters
   - Verify thresholds match
   - Check frame calculations

2. **Systematic quality score analysis**
   - Create comprehensive quality score comparison tool
   - Analyze all mismatched pairs

3. **Validation threshold verification**
   - Ensure all thresholds match exactly
   - Document any differences

4. **H-bond detection verification**
   - Compare H-bond detection logic
   - Verify H-bond validation

---

## Quick Status

- **Total PDBs tested**: 100
- **Perfect matches**: 100 (100.0%)
- **Known mismatches**: 0
- **Timeouts**: 0 (resolved with script optimization)
- **Fixed issues**: 2 (bp_type_id bug, is_nucleotide() bug)
- **Active investigations**: 0

**Note**: See `docs/BROAD_TESTING_RESULTS.md` for detailed results across 100 PDBs.

---

## Contact & Updates

For questions or updates to this documentation, refer to:
- Individual difference documentation files
- Analysis reports (e.g., `1T0K_mismatched_pairs_analysis.md`)
- Test results in `docs/TEST_RESULTS_BP_TYPE_ID_FIX.md`


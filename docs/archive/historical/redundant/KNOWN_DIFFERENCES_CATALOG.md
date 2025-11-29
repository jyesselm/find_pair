# Known Differences Catalog

**Purpose**: Comprehensive documentation of all known types of differences between legacy and modern code outputs.

**Last Updated**: 2025-01-27

---

## Categories of Differences

### 1. ✅ FIXED: bp_type_id Calculation Bug

**Status**: ✅ **RESOLVED**

**Issue**: Legacy code passes wrong parameters to `check_wc_wobble_pair`.

**Root Cause**:
- Legacy calls: `check_wc_wobble_pair(bpid, bpi, pars[1], pars[2], pars[6])`
- Where `pars[1]=Shift`, `pars[2]=Slide`, `pars[6]=Twist`
- But function expects: `(shear, stretch, opening)`
- Legacy incorrectly uses **Shift as shear** and **Slide as stretch**

**Impact**: 
- Modern assigns `bp_type_id=2` (Watson-Crick) with -2.0 quality adjustment
- Legacy assigns `bp_type_id=-1` (not classified) with no adjustment
- Results in different pair selections

**Fix**: Updated modern code to match legacy's buggy parameter order.

**Affected PDBs**: 6CAQ (4 pairs), potentially others

**Documentation**: See `docs/BP_TYPE_ID_BUG_FIX.md`

---

### 2. is_nucleotide() Classification Bug

**Status**: ✅ **FIXED AND VERIFIED**

**Issue**: `is_nucleotide()` incorrectly classifies glucose (GLC) and other non-nucleotide residues as nucleotides.

**Root Cause**:
- Modern checks for UNKNOWN residues with >= 3 ring atoms
- Glucose has C4, C5, C6 atoms that match nucleotide ring atom names
- Modern doesn't require nitrogen atoms (N1 or N3) like legacy does
- Modern doesn't perform RMSD check like legacy does

**Impact**: 
- Non-nucleotide pairs are validated and selected
- Example: 1T0K pair (491, 492) - both GLC residues

**Fix**: Require nitrogen atoms (N1 or N3) in addition to ring atoms for UNKNOWN residues

**Documentation**: See `docs/IS_NUCLEOTIDE_BUG_ANALYSIS.md`

---

### 3. Quality Score Differences (Non-bp_type_id)

**Status**: ⚠️ **INVESTIGATING**

**Issue**: Pairs have different quality scores even when `bp_type_id` matches.

**Examples**:
- **1T0K pair (491, 492)**: Extra in modern, not in legacy
  - Root cause: Base quality score difference
  - Not related to `bp_type_id` (no `bp_type_id` differences found)

**Potential Causes**:
1. **Frame calculation differences**: Slight numerical precision differences
2. **Distance calculation differences**: `dorg`, `d_v`, `dNN` calculations
3. **Angle calculation differences**: `plane_angle` calculations
4. **H-bond counting differences**: Different H-bond detection/validation
5. **adjust_pairQuality differences**: Different H-bond quality adjustments

**Investigation Needed**:
- Compare quality score components for mismatched pairs
- Verify frame calculations match exactly
- Check H-bond detection logic

**Documentation**: See analysis reports (e.g., `1T0K_mismatched_pairs_analysis.md`)

---

### 4. Validation Status Differences

**Status**: ⚠️ **INVESTIGATING**

**Issue**: Pairs that pass validation in one implementation but fail in the other.

**Examples**:
- Pairs that are `is_valid=1` in legacy but `is_valid=0` in modern (or vice versa)

**Potential Causes**:
1. **Threshold differences**: Slight differences in validation thresholds
2. **Overlap calculation differences**: Different overlap area calculations
3. **H-bond validation differences**: Different H-bond validation logic
4. **Distance check differences**: Different distance constraint checks

**Investigation Needed**:
- Compare validation results component-by-component
- Verify all validation thresholds match exactly
- Check edge case handling

---

### 5. Pair Selection Differences (Same Quality Scores)

**Status**: ⚠️ **INVESTIGATING**

**Issue**: Different pairs selected even when quality scores are identical.

**Potential Causes**:
1. **Tie-breaking logic**: Different tie-breaking when multiple pairs have same score
2. **Selection order**: Different iteration order leading to different selections
3. **Filtering differences**: Different filtering logic before selection

**Investigation Needed**:
- Compare selection algorithm step-by-step
- Verify tie-breaking logic matches
- Check pair filtering logic

---

### 6. Missing/Extra Pairs

**Status**: ⚠️ **VARIES BY PDB**

**Issue**: Pairs present in one implementation but not the other.

**Categories**:

#### 5a. Missing in Modern (in legacy, not in modern)
- **Cause**: Modern validation rejects pairs that legacy accepts
- **Investigation**: Compare validation results for these pairs

#### 5b. Extra in Modern (in modern, not in legacy)
- **Cause**: Modern validation accepts pairs that legacy rejects
- **Investigation**: Compare validation results for these pairs

**Examples**:
- **1T0K**: 1 extra pair (491, 492) - quality score difference
- **3G8T**: 1 missing + 1 extra - quality score and `bp_type_id` differences

---

## Investigation Methodology

### Step 1: Identify Mismatched Pairs
```bash
python3 scripts/analyze_mismatched_pairs.py <PDB_ID>
```

### Step 2: Check bp_type_id Differences
```bash
python3 scripts/investigate_bp_type_id_differences.py <PDB_ID>
```

### Step 3: Compare Quality Score Components
- Use `tools/compare_quality_score_components` to analyze quality score differences
- Compare: base_score, adjust_pairQuality, bp_type_id, final_score

### Step 4: Compare Validation Results
- Load `pair_validation` JSON records
- Compare: is_valid, dorg, d_v, plane_angle, dNN, overlap_area

### Step 5: Compare Frames
- Use `tools/compare_frames_and_step_params` to verify frame calculations
- Compare: frame origins, rotation matrices

---

## Current Status by PDB

### ✅ Perfect Matches
- **6CAQ**: 100% match (623/623) - Fixed with bp_type_id bug fix
- **3G8T**: Perfect match
- **1ZX7**: Perfect match
- **2B8R**: Perfect match
- **2QEX**: Perfect match (1156/1156)

### ⚠️ Known Mismatches

#### 1T0K
- **Status**: 1 extra pair
- **Pair**: (491, 492)
- **Root Cause**: Quality score difference (not bp_type_id)
- **Priority**: MEDIUM
- **Next Steps**: Investigate quality score components

#### 3G8T (Historical)
- **Status**: Previously had 1 missing + 1 extra
- **Current**: Perfect match ✅
- **Resolution**: Likely fixed by bp_type_id bug fix

---

## Tools for Investigation

1. **`scripts/analyze_mismatched_pairs.py`**: Comprehensive mismatch analysis
2. **`scripts/investigate_bp_type_id_differences.py`**: bp_type_id-specific analysis
3. **`tools/compare_quality_score_components`**: Quality score component comparison
4. **`tools/compare_frames_and_step_params`**: Frame and step parameter comparison
5. **`tools/compare_bp_type_id_calculation`**: bp_type_id calculation analysis

---

## Next Steps

1. **Investigate 1T0K mismatch**: Compare quality score components for pair (491, 492)
2. **Systematic quality score analysis**: Create tool to compare all quality score components
3. **Validation threshold verification**: Ensure all validation thresholds match exactly
4. **H-bond detection verification**: Compare H-bond detection and validation logic
5. **Tie-breaking logic**: Verify selection tie-breaking matches legacy

---

## Documentation References

- `docs/BP_TYPE_ID_BUG_FIX.md`: bp_type_id bug fix details
- `docs/TEST_RESULTS_BP_TYPE_ID_FIX.md`: Test results after fix
- `docs/LEGACY_STEP_PARAMETER_ANALYSIS.md`: Step parameter analysis
- `docs/JSON_GENERATION_SYSTEM.md`: JSON file management system
- `docs/CURRENT_STATUS_SUMMARY.md`: Overall project status


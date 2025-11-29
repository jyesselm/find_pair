# 100% Match Plan: Base Pair Frames

**Last Updated**: 2025-01-XX  
**Goal**: Achieve 100% match on base pair frames - same pairs found, same unique pairs with same parameters  
**Current Status**: 
- **find_bestpair_selection**: 99.6% match (1044/1048 common pairs)
- **base_pair records**: 51.9% match (1321/2548 common pairs)
- **10-PDB test set**: 90% perfect matches (9/10 PDBs)

---

## Executive Summary

### Current Match Status

**find_bestpair_selection** (Selected Pairs):
- ✅ **99.6% match** (1044/1048 common pairs)
- ⚠️ **4 missing pairs**: From 3G8T and 6CAQ
- ⚠️ **5 extra pairs**: From 3G8T and 6CAQ
- **10-PDB test set**: 90% perfect matches (9/10 PDBs)

**base_pair records** (All Validated Pairs with Frames):
- ⚠️ **51.9% match** (1321/2548 common pairs)
- ⚠️ **1227 missing pairs** in modern
- ⚠️ **107 extra pairs** in modern
- **Critical**: These records contain the frame information (orien_i, orien_j, org_i, org_j)

### Root Causes Identified

1. **Stage 6 Frame Accessibility** (HIGH PRIORITY)
   - Frames may not be accessible when `calculate_bp_type_id()` is called
   - Pair (1102, 1127) in 6CAQ: Legacy `bp_type_id=2`, Modern `bp_type_id=-1`
   - Blocks correct `bp_type_id = 2` assignment and quality score adjustment

2. **Pairs Not Found in Validation** (HIGH PRIORITY)
   - 6CAQ: 9 pairs not found in validation
   - Pairs: (495, 498), (501, 506), (939, 1378), (1029, 1184), (1380, 1473), (1382, 1470), (1385, 1467), (1489, 1492)
   - These pairs should be validated but aren't reaching validation

3. **Quality Score Differences** (MEDIUM PRIORITY)
   - Affects pair selection when scores are close
   - May cause tie-breaking differences
   - Related to `adjust_pairQuality()` and `bp_type_id` calculation

4. **Base Pair Record Generation** (MEDIUM PRIORITY)
   - Modern may not be recording all validated pairs as `base_pair` records
   - Legacy records `base_pair` for ALL pairs that pass validation in `calculate_more_bppars`
   - Modern should match this behavior exactly

---

## Phase 1: Debug Stage 6 Frame Accessibility (HIGH PRIORITY)

**Goal**: Ensure frames are accessible when `calculate_bp_type_id()` is called

### Issue
- Pair (1102, 1127) in 6CAQ: Legacy `bp_type_id=2`, Modern `bp_type_id=-1`
- Direction vectors match condition (`dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0`)
- Validation passes (requires frames), but frames may not be accessible during `bp_type_id` calculation

### Investigation Steps

1. **Verify Frame Storage** ✅
   - Location: `src/x3dna/algorithms/base_frame_calculator.cpp`
   - Verify: `residue.set_reference_frame()` is called after calculation
   - Check: Frames are stored on residue objects

2. **Debug Frame Access in calculate_bp_type_id()** ✅ COMPLETE
   - Location: `src/x3dna/algorithms/base_pair_finder.cpp` (line 608)
   - Added debug output: Checks `res1->reference_frame().has_value()` and `res2->reference_frame().has_value()`
   - Added debug output for step parameters and bp_type_id assignment
   - Enable with: `cmake -DDEBUG_BP_TYPE_ID=ON ..` and rebuild
   - Status: ✅ Debug output added, ready for testing

3. **Test Step Parameter Calculation** ⏳
   - Location: `src/x3dna/algorithms/parameter_calculator.cpp`
   - Compare: Modern step parameters vs legacy `bpstep_params` JSON records
   - Test pair: (1102, 1127) in 6CAQ (if frames become accessible)

4. **Verify bp_type_id Assignment** ⏳
   - Once frames are accessible, verify classification works
   - Expected: Pair (1102, 1127) should get `bp_type_id = 2` (Watson-Crick)
   - Check: All conditions met (direction vectors, stretch ≤ 2.0, opening ≤ 60°, shear ≤ 1.8, in WC_LIST)

### Files to Modify
- `src/x3dna/algorithms/base_pair_finder.cpp`: Add debug output, verify frame access
- `src/x3dna/algorithms/base_frame_calculator.cpp`: Verify frames are stored on residues
- `src/x3dna/algorithms/parameter_calculator.cpp`: Verify step parameter calculation

### Success Criteria
- ✅ Frames accessible in `calculate_bp_type_id()`
- ✅ Step parameters calculated correctly
- ✅ `bp_type_id = 2` assigned for pair (1102, 1127) in 6CAQ
- ✅ Quality score matches legacy (should be 2.0 lower due to `bp_type_id = 2`)

---

## Phase 2: Fix Pairs Not Found in Validation (HIGH PRIORITY)

**Goal**: Ensure all pairs that should be validated are actually validated

### Issue
- 6CAQ: 9 pairs not found in validation
- Pairs: (495, 498), (501, 506), (939, 1378), (1029, 1184), (1380, 1473), (1382, 1470), (1385, 1467), (1489, 1492)
- These pairs should be validated but aren't reaching validation

### Investigation Steps

1. **Check Residue Recognition** ⏳
   - Verify: Both residues in each pair are recognized as nucleotides
   - Check: Frames are calculated for both residues
   - Location: `src/x3dna/algorithms/base_pair_finder.cpp` (Phase 1 validation loop)

2. **Check Iteration Order** ⏳
   - Verify: All residue pairs are being checked
   - Check: Iteration order matches legacy (sequential from 1 to max_legacy_idx)
   - Location: `src/x3dna/algorithms/base_pair_finder.cpp` (line 100-150)

3. **Check Early Exit Conditions** ⏳
   - Verify: No early exits preventing validation
   - Check: All pairs reach `validator_.validate()` call
   - Location: `src/x3dna/algorithms/base_pair_finder.cpp` (Phase 1)

4. **Compare with Legacy Validation** ⏳
   - Load legacy `pair_validation` records for 6CAQ
   - Check: Which pairs have validation records in legacy
   - Identify: Why modern is missing these pairs

### Files to Modify
- `src/x3dna/algorithms/base_pair_finder.cpp`: Add debug output for pair iteration
- `tools/compare_validation_coverage.cpp`: Create tool to compare validation coverage

### Success Criteria
- ✅ All 9 pairs in 6CAQ are validated
- ✅ Validation records match legacy for these pairs
- ✅ Base pair records created for validated pairs

---

## Phase 3: Ensure Base Pair Record Generation Matches Legacy (MEDIUM PRIORITY)

**Goal**: Ensure modern records `base_pair` for ALL validated pairs, matching legacy behavior

### Issue
- Legacy records `base_pair` for ALL pairs that pass validation in `calculate_more_bppars`
- Modern may not be recording all validated pairs
- Current: 1227 missing base_pair records, 107 extra base_pair records

### Investigation Steps

1. **Verify Recording Logic** ⏳
   - Location: `src/x3dna/algorithms/base_pair_finder.cpp` (line 535-562)
   - Check: `writer->record_base_pair()` is called for ALL `result.is_valid == true` pairs
   - Verify: No conditions skip recording

2. **Compare Recording Conditions** ⏳
   - Legacy: Records in `calculate_more_bppars` for all pairs with `is_valid == 1`
   - Modern: Should match this exactly
   - Check: No additional filters or conditions

3. **Verify Frame Availability** ⏳
   - Check: Frames are available when recording `base_pair`
   - Verify: `pair.frame1()` and `pair.frame2()` are set correctly
   - Location: `src/x3dna/algorithms/base_pair_finder.cpp` (line 542-548)

4. **Compare Base Pair Counts** ⏳
   - Count: Legacy `base_pair` records vs modern `base_pair` records
   - Identify: Which pairs are missing/extra
   - Analyze: Why these pairs differ

### Files to Modify
- `src/x3dna/algorithms/base_pair_finder.cpp`: Verify recording logic
- `src/x3dna/io/json_writer.cpp`: Verify `record_base_pair()` implementation
- `scripts/analyze_base_pair_coverage.py`: Create script to compare base_pair coverage

### Success Criteria
- ✅ All validated pairs have `base_pair` records
- ✅ Base pair count matches legacy
- ✅ Frame information (orien_i, orien_j, org_i, org_j) matches legacy exactly

---

## Phase 4: Fix Quality Score Differences (MEDIUM PRIORITY)

**Goal**: Ensure quality scores match legacy exactly for pair selection

### Issue
- Quality score differences affect pair selection
- May cause tie-breaking differences
- Related to `adjust_pairQuality()` and `bp_type_id` calculation

### Investigation Steps

1. **Compare Quality Score Calculation** ⏳
   - Base score: `dorg + 2.0 * d_v + plane_angle / 20.0`
   - Adjustment: `adjust_pairQuality()` (from H-bonds)
   - Watson-Crick adjustment: `-2.0` if `bp_type_id == 2`
   - Final score: `base_score + adjustment - (2.0 if bp_type_id == 2)`
   - Location: `src/x3dna/algorithms/base_pair_finder.cpp` (line 583-606, 608-650)

2. **Verify adjust_pairQuality()** ⏳
   - Count good H-bonds (distance in [2.5, 3.5])
   - If `num_good_hb >= 2`: return `-3.0`
   - Otherwise: return `-num_good_hb`
   - Location: `src/x3dna/algorithms/base_pair_finder.cpp` (line 583-606)

3. **Compare with Legacy Quality Scores** ⏳
   - Load legacy `pair_validation` records
   - Compare: `rtn_val[5]` (final adjusted quality score)
   - Identify: Differences in calculation

4. **Test on Mismatched Pairs** ⏳
   - Test pairs: (92, 160), (946, 947) in 3G8T
   - Test pairs: (75, 78), (968, 1024) in 6CAQ
   - Compare: Quality scores step-by-step

### Files to Modify
- `src/x3dna/algorithms/base_pair_finder.cpp`: Fix quality score calculation
- `tools/compare_quality_scores.cpp`: Create tool to compare quality scores

### Success Criteria
- ✅ Quality scores match legacy exactly for all pairs
- ✅ Pair selection matches legacy (100% match on find_bestpair_selection)
- ✅ No tie-breaking differences

---

## Phase 5: Verify Frame Parameters Match Legacy (MEDIUM PRIORITY)

**Goal**: Ensure frame parameters (orien_i, orien_j, org_i, org_j) match legacy exactly

### Issue
- Base pair frames must match legacy exactly
- Rotation matrices and origins must be identical
- Current: Frame calculations match, but need to verify for all base_pair records

### Investigation Steps

1. **Compare Frame Values** ⏳
   - Load legacy and modern `base_pair` records
   - Compare: `orien_i`, `orien_j`, `org_i`, `org_j` for each pair
   - Tolerance: < 0.0001 per element

2. **Verify Frame Extraction** ⏳
   - Check: Frames are extracted correctly from residue objects
   - Verify: `pair.frame1()` and `pair.frame2()` match residue frames
   - Location: `src/x3dna/algorithms/base_pair_finder.cpp` (line 542-548)

3. **Compare Direction Vectors** ⏳
   - Check: `dir_xyz` matches legacy (including legacy's 1-based indexing bug)
   - Legacy bug: Stores `[dir_y, dir_z, 0.0]` instead of `[dir_x, dir_y, dir_z]`
   - Modern should replicate this bug to match legacy

4. **Test on All Base Pair Records** ⏳
   - Compare: All common base_pair records
   - Verify: Frame parameters match exactly
   - Identify: Any differences

### Files to Modify
- `src/x3dna/core/base_pair.cpp`: Verify `to_json_legacy()` implementation
- `x3dna_json_compare/base_pair_comparison.py`: Enhance frame comparison

### Success Criteria
- ✅ All frame parameters match legacy exactly
- ✅ Rotation matrices match (< 0.0001 per element)
- ✅ Origins match (< 0.0001)
- ✅ Direction vectors match (including legacy bug)

---

## Phase 6: Final Verification and Testing

**Goal**: Achieve 100% match on base pair frames

### Testing Steps

1. **Run Full Comparison** ⏳
   ```bash
   python3 scripts/compare_json.py compare --test-set 100
   ```

2. **Verify find_bestpair_selection** ⏳
   - Target: 100% match (all pairs match exactly)
   - Current: 99.6% match (1044/1048 common pairs)

3. **Verify base_pair Records** ⏳
   - Target: 100% match (all pairs match exactly)
   - Current: 51.9% match (1321/2548 common pairs)

4. **Verify Frame Parameters** ⏳
   - Target: All frame parameters match exactly
   - Check: orien_i, orien_j, org_i, org_j for all base_pair records

5. **Test on 10-PDB Test Set** ⏳
   - Target: 100% perfect matches (10/10 PDBs)
   - Current: 90% perfect matches (9/10 PDBs)

### Success Criteria
- ✅ 100% match on find_bestpair_selection
- ✅ 100% match on base_pair records
- ✅ All frame parameters match legacy exactly
- ✅ 100% perfect matches on 10-PDB test set

---

## Implementation Checklist

### Phase 1: Stage 6 Frame Accessibility
- [x] Add debug output in `calculate_bp_type_id()` to check frame access ✅
- [x] Add CMake option to enable debug flag (DEBUG_BP_TYPE_ID) ✅
- [ ] Verify frames are stored on residues after calculation
- [ ] Test step parameter calculation for pair (1102, 1127) in 6CAQ
- [ ] Verify `bp_type_id = 2` assignment works correctly
- [ ] Verify quality score adjustment matches legacy

### Phase 2: Pairs Not Found in Validation
- [ ] Check residue recognition for 9 missing pairs in 6CAQ
- [ ] Verify iteration order matches legacy
- [ ] Check for early exit conditions
- [ ] Compare validation coverage with legacy
- [ ] Fix any issues found

### Phase 3: Base Pair Record Generation
- [ ] Verify recording logic matches legacy
- [ ] Check frame availability when recording
- [ ] Compare base pair counts with legacy
- [ ] Fix any recording issues

### Phase 4: Quality Score Differences
- [ ] Compare quality score calculation step-by-step
- [ ] Verify `adjust_pairQuality()` matches legacy
- [ ] Test on mismatched pairs
- [ ] Fix any calculation differences

### Phase 5: Frame Parameters
- [ ] Compare frame values for all base_pair records
- [ ] Verify frame extraction from residues
- [ ] Compare direction vectors (including legacy bug)
- [ ] Fix any frame parameter differences

### Phase 6: Final Verification
- [ ] Run full comparison on 100-set
- [ ] Verify 100% match on find_bestpair_selection
- [ ] Verify 100% match on base_pair records
- [ ] Verify all frame parameters match
- [ ] Test on 10-PDB test set

---

## Tools and Scripts

### Comparison Tools
- `scripts/compare_json.py`: Main comparison tool
- `scripts/analyze_mismatched_pairs.py`: Detailed analysis of mismatched pairs
- `x3dna_json_compare/base_pair_comparison.py`: Base pair comparison

### Debug Tools
- `tools/compare_validation_coverage.cpp`: Compare validation coverage (to be created)
- `tools/compare_quality_scores.cpp`: Compare quality scores (to be created)
- `tools/compare_base_pair_frames.cpp`: Compare base pair frames (to be created)

### Analysis Scripts
- `scripts/analyze_base_pair_coverage.py`: Analyze base_pair coverage (to be created)
- `scripts/analyze_frame_parameters.py`: Analyze frame parameters (to be created)

---

## Success Metrics

**Target**: 100% match on base pair frames

**Current Metrics**:
- find_bestpair_selection: 99.6% match (1044/1048 common pairs)
- base_pair records: 51.9% match (1321/2548 common pairs)
- 10-PDB test set: 90% perfect matches (9/10 PDBs)

**Target Metrics**:
- find_bestpair_selection: 100% match
- base_pair records: 100% match
- Frame parameters: 100% match
- 10-PDB test set: 100% perfect matches (10/10 PDBs)

---

## Related Documentation

- [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) - Current status and completed fixes
- [100_PERCENT_MATCH_PLAN.md](100_PERCENT_MATCH_PLAN.md) - Overall 100% match plan
- [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md) - JSON record types and comparison details
- [TESTING_GUIDE.md](TESTING_GUIDE.md) - Testing and comparison workflows

---

*This plan focuses specifically on achieving 100% match on base pair frames - ensuring we find the same pairs and have the same unique pairs with the same parameters.*


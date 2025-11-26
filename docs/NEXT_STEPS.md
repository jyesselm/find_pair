# Next Steps & Action Plan

**Last Updated**: 2025-11-26  
**Current Focus**: Debug Stage 6 (ParameterCalculator) frame accessibility

---

## Immediate Priority: Debug Stage 6 ⏳ IN PROGRESS

### Problem
Pair (1102, 1127) in 6CAQ shows:
- **Legacy**: `bp_type_id = 2` (Watson-Crick), quality_score = -1.000905
- **Modern**: `bp_type_id = -1` (not classified), quality_score = 2.999095
- **Difference**: 4.0 in quality score (2.0 from base score + 2.0 from `bp_type_id = 2` adjustment)

Direction vectors match condition (`dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0`), so step parameters should be calculated.

### Root Cause Hypothesis
Reference frames may not be accessible when `calculate_bp_type_id()` is called, even though:
- Validation passed (which requires frames)
- Frames are used in `BasePairValidator::validate()`

### Investigation Steps

#### Step 1: Verify Frame Storage
**File**: `src/x3dna/algorithms/base_frame_calculator.cpp`

Check if frames are stored on residue objects:
```cpp
// After calculating frame:
residue.set_reference_frame(frame);
```

**Action**: Verify `set_reference_frame()` is called after frame calculation.

#### Step 2: Debug Frame Access
**File**: `src/x3dna/algorithms/base_pair_finder.cpp`

Add temporary debug output in `calculate_bp_type_id()`:
```cpp
if (result.dir_x > 0.0 && result.dir_y < 0.0 && result.dir_z < 0.0) {
    bool has_frame1 = res1->reference_frame().has_value();
    bool has_frame2 = res2->reference_frame().has_value();
    
    // Debug output for pair (1102, 1127)
    if (/* check if this is the problematic pair */) {
        std::cerr << "[DEBUG] Pair (1102, 1127): "
                  << "frame1=" << has_frame1 << ", frame2=" << has_frame2 << std::endl;
    }
    
    if (!has_frame1 || !has_frame2) {
        return bp_type_id; // Keep -1 if frames not available
    }
    // ... rest of calculation
}
```

**Action**: Add debug output, rebuild, and check logs for pair (1102, 1127).

#### Step 3: Test Step Parameter Calculation
**File**: `src/x3dna/algorithms/parameter_calculator.cpp`

Once frames are accessible, verify step parameters match legacy:
- Compare modern `bpstep_par_impl()` output with legacy `bpstep_params` JSON
- Test pair: (1102, 1127) in 6CAQ
- Expected parameters:
  - Slide (shear): Should be ≤ 1.8 for Watson-Crick
  - Rise (stretch): Should be ≤ 2.0
  - Twist (opening): Should be ≤ 60°

**Action**: Add debug output to log calculated step parameters.

#### Step 4: Verify bp_type_id Assignment
**File**: `src/x3dna/algorithms/base_pair_finder.cpp`

Once step parameters are calculated, verify classification:
1. Check stretch ≤ 2.0 and opening ≤ 60° ✓
2. Check shear ≤ 1.8 ✓
3. Check base pair type in WC_LIST ✓
4. Set `bp_type_id = 2` ✓

**Action**: Verify all conditions are met and `bp_type_id = 2` is assigned.

### Success Criteria
- ✅ Frames accessible in `calculate_bp_type_id()`
- ✅ Step parameters calculated correctly (match legacy within tolerance)
- ✅ `bp_type_id = 2` assigned for pair (1102, 1127)
- ✅ Quality score matches legacy (should be 2.0 lower: -1.000905 vs 2.999095)

---

## Secondary Priority: Investigate Missing Pairs in 6CAQ

### Problem
9 pairs are missing in modern but exist in legacy:
- (495, 498), (501, 506), (939, 1378), (1029, 1184)
- (1380, 1473), (1382, 1470), (1385, 1467), (1489, 1492)

### Investigation Steps

1. **Check if pairs are validated**: Are these pairs being validated at all?
   - Check: `data/json/pair_validation/6CAQ.json`
   - Look for: Any records with these residue indices

2. **Check validation failures**: If validated, why do they fail?
   - Check: `is_valid` field
   - Check: Which validation checks fail (distance, angle, overlap, H-bonds)

3. **Check residue indices**: Verify residue indices match between legacy and modern
   - Legacy uses 1-based indices from PDB parsing
   - Modern may use different indexing

4. **Check frame availability**: Do these residues have reference frames?
   - If no frames, pairs can't be validated
   - Check: Frame calculation logs for these residues

### Files to Check
- `data/json/pair_validation/6CAQ.json`: Modern validation records
- `data/json_legacy/6CAQ.json`: Legacy validation records
- `data/json/base_frame_calc/6CAQ.json`: Frame calculation records

---

## Tertiary Priority: Test Stage 6 on Additional PDBs

### Goal
Verify Stage 6 works across different structures, not just 6CAQ.

### Test Plan

1. **Find PDBs with `bp_type_id = 2` pairs**:
   ```bash
   # Search legacy JSON for pairs with bp_type_id = 2
   python3 scripts/find_bp_type_id_pairs.py --type 2
   ```

2. **Generate modern JSON** for these PDBs:
   ```bash
   ./build/generate_modern_json data/pdb/<PDB_ID>.pdb data/json
   ```

3. **Compare bp_type_id assignments**:
   ```bash
   python3 scripts/compare_json.py \
     data/json_legacy/<PDB_ID>.json \
     data/json/pair_validation/<PDB_ID>.json \
     --field bp_type_id
   ```

4. **Verify quality score adjustments**:
   - Check that pairs with `bp_type_id = 2` have quality_score -= 2.0
   - Compare final quality scores with legacy

### Success Criteria
- ✅ All PDBs with `bp_type_id = 2` pairs show correct classification
- ✅ Quality score adjustments match legacy
- ✅ No regressions in other pairs

---

## Long-Term: Tie-Breaking Logic

### Issue
When quality scores are identical, selection order may differ between legacy and modern.

### Current Status
- 3 extra pairs in 6CAQ have matching quality scores
- These are likely acceptable differences (tie-breaking)

### Investigation (Low Priority)
1. Document current tie-breaking behavior
2. Compare iteration order between legacy and modern
3. Determine if differences are acceptable

### Action
- Document in `docs/TIE_BREAKING.md` (if needed)
- Mark as acceptable if scores match exactly

---

## Testing Checklist

### Stage 6 Debugging
- [ ] Verify frames are stored on residues after calculation
- [ ] Add debug output to `calculate_bp_type_id()`
- [ ] Rebuild and test on 6CAQ
- [ ] Check debug logs for pair (1102, 1127)
- [ ] Verify step parameters are calculated
- [ ] Verify `bp_type_id = 2` is assigned
- [ ] Verify quality score matches legacy

### Missing Pairs Investigation
- [ ] Check if 9 missing pairs are validated
- [ ] Check validation failure reasons
- [ ] Verify residue indices match
- [ ] Check frame availability for these residues

### Broader Testing
- [ ] Find PDBs with `bp_type_id = 2` pairs
- [ ] Test Stage 6 on these PDBs
- [ ] Verify quality score adjustments
- [ ] Run 100-set test to measure overall improvement

---

## Files Modified/To Modify

### Already Modified ✅
- `include/x3dna/algorithms/parameter_calculator.hpp`
- `src/x3dna/algorithms/parameter_calculator.cpp`
- `include/x3dna/algorithms/base_pair_finder.hpp`
- `src/x3dna/algorithms/base_pair_finder.cpp`
- `CMakeLists.txt`
- `docs/100_PERCENT_MATCH_STATUS.md`

### To Modify ⏳
- `src/x3dna/algorithms/base_pair_finder.cpp`: Add debug output
- `src/x3dna/algorithms/base_frame_calculator.cpp`: Verify frame storage (if needed)

### Tools to Create (if needed)
- `scripts/find_bp_type_id_pairs.py`: Find pairs with specific `bp_type_id` in legacy JSON

---

## Expected Outcomes

### After Stage 6 Debugging
- ✅ Pair (1102, 1127) in 6CAQ: `bp_type_id = 2`, quality_score matches legacy
- ✅ 6CAQ: 9 missing, 3 extra pairs (down from 10 missing, 3 extra)
- ✅ Overall match rate improvement

### After Missing Pairs Investigation
- ✅ Understanding of why 9 pairs aren't validated
- ✅ Fix or document as acceptable differences
- ✅ Further improvement in 6CAQ match rate

### After Broader Testing
- ✅ Confidence that Stage 6 works across structures
- ✅ Verified quality score adjustments
- ✅ Overall match rate improvement on 100-set test

---

## Notes

- **Frame Accessibility**: The key issue is ensuring frames are accessible when needed. Validation uses frames, so they must exist - the question is whether they're stored on residue objects or accessed differently.

- **Step Parameters**: Once frames are accessible, step parameter calculation should work correctly as it matches the legacy algorithm exactly.

- **bp_type_id Classification**: The logic is correct and matches legacy - it just needs frames and step parameters to work.

- **Quality Score Impact**: `bp_type_id = 2` causes a -2.0 adjustment, which is significant for pair selection. This is why fixing it is important.


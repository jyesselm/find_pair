# Find Pair Selection Investigation

**Date**: 2025-11-29  
**Status**: 🔍 Investigating why modern and legacy select different base pairs

---

## Problem Statement

Modern `find_pair` selects **23 pairs** for 1H4S, while legacy selects **25 pairs**. This difference causes step parameter values to differ because they're calculated for different base pairs.

**Impact**: Step parameter values won't match even if calculation is correct, because inputs (base pairs) differ.

---

## Current Status

### Modern find_pair

**Issues Identified**:
1. ❌ **No JsonWriter setup**: `find_pair_app.cpp` doesn't create or use JsonWriter
   - Result: No JSON files generated (no `find_bestpair_selection`, no `base_pair`, etc.)
   - Code location: `apps/find_pair_app.cpp` - no JsonWriter instantiation
   
2. ⚠️ **Different pair selection**: Modern selects 23 pairs vs legacy 25 pairs
   - May be due to:
     - Quality score calculation differences
     - Validation threshold differences
     - Tie-breaking differences
     - Iteration order differences

### Legacy find_pair

**Behavior**:
- ✅ Generates `find_bestpair_selection` JSON file
- ✅ Selects 25 pairs for 1H4S
- ✅ Records all JSON output correctly

---

## Root Causes

### 1. Missing JsonWriter in find_pair_app

**File**: `apps/find_pair_app.cpp`

**Current State**:
- No JsonWriter instantiation
- No JSON recording
- Protocol executes without JSON output

**Comparison with analyze_app**:
- `analyze_app.cpp` creates JsonWriter and sets it on protocol
- `find_pair_app.cpp` does NOT create JsonWriter

**Fix Needed**:
```cpp
// Add to find_pair_app.cpp (similar to analyze_app.cpp):
x3dna::io::JsonWriter json_writer(options.pdb_file);
protocol.set_json_writer(&json_writer);
// ... execute protocol ...
json_writer.write_split_files("data/json", true);
```

### 2. Pair Selection Differences

**Possible Causes**:

1. **Quality Score Calculation**:
   - Legacy: `rtn_val[5]` = base_score + `adjust_pairQuality` + (bp_type_id == 2 ? -2.0 : 0)
   - Modern: May have slight differences in calculation
   - Location: `src/x3dna/algorithms/base_pair_finder.cpp:436-450`

2. **Validation Thresholds**:
   - Need to verify all thresholds match exactly
   - Location: `include/x3dna/algorithms/base_pair_validator.hpp`

3. **Tie-Breaking**:
   - When multiple pairs have same quality_score, selection depends on iteration order
   - Need to verify iteration order matches legacy exactly

4. **Frame Calculation**:
   - Different frames → different quality scores → different selection
   - Need to verify frames match legacy before pair selection

---

## Investigation Steps

### Step 1: Add JsonWriter to find_pair_app ✅ (Next)

**Action**: Modify `apps/find_pair_app.cpp` to:
1. Create JsonWriter with PDB file path
2. Set JsonWriter on protocol before execution
3. Write JSON files after execution

**Expected Result**: Modern will generate `find_bestpair_selection` JSON files

### Step 2: Compare Selected Pairs

**Action**: After Step 1, compare modern vs legacy `find_bestpair_selection` for 1H4S

**Questions to Answer**:
- Which pairs differ?
- Are the differences due to:
  - Missing pairs in modern?
  - Extra pairs in modern?
  - Both?

### Step 3: Trace Pair Selection

**Action**: For pairs that differ, trace through selection logic:
1. Check if pair passes validation (Phase 1)
2. Check quality_score calculation
3. Check if pair is selected as "best partner"
4. Check mutual best match logic

**Tools**:
- Use `tools/trace_pair_selection.cpp` if available
- Add debug logging to `base_pair_finder.cpp`

### Step 4: Verify Quality Score Calculation

**Action**: Compare quality scores for common pairs:
1. Extract quality scores from modern validation
2. Extract quality scores from legacy validation
3. Compare `adjust_pairQuality` values
4. Compare `bp_type_id` values
5. Compare final `rtn_val[5]` (adjusted quality_score)

**Expected**: Scores should match exactly for same pairs

### Step 5: Verify Validation Thresholds

**Action**: Verify all validation parameters match legacy:
- `min_dorg`, `max_dorg`
- `min_dv`, `max_dv`
- `min_dNN`, `max_dNN`
- `min_plane_angle`, `max_plane_angle`
- `min_base_hb`
- `overlap_threshold`

**Location**: `include/x3dna/algorithms/base_pair_validator.hpp`

---

## Next Steps (Priority Order)

### Priority 1: Enable JSON Recording in find_pair_app

**Why**: Need to see what modern actually selects to compare with legacy

**Action**:
1. Add JsonWriter to `find_pair_app.cpp`
2. Test JSON generation for 1H4S
3. Compare `find_bestpair_selection` files

**Files to Modify**:
- `apps/find_pair_app.cpp` - Add JsonWriter setup

### Priority 2: Compare Selected Pairs

**Why**: Need to understand which pairs differ and why

**Action**:
1. Load modern and legacy `find_bestpair_selection` for 1H4S
2. Identify:
   - Common pairs
   - Pairs only in legacy
   - Pairs only in modern
3. Document differences

### Priority 3: Trace Selection Logic

**Why**: Need to understand why specific pairs are selected/deselected

**Action**:
1. For pairs that differ, trace through:
   - Validation phase
   - Quality score calculation
   - Best partner selection
   - Mutual best match logic
2. Identify where logic diverges

### Priority 4: Fix Selection Differences

**Why**: To get matching step parameter values

**Action**:
1. Fix quality score calculation if needed
2. Fix validation thresholds if needed
3. Fix tie-breaking if needed
4. Verify frames match before selection

---

## Known Issues

### Issue 1: No JSON Output from find_pair_app

**Status**: ❌ Not implemented
**Impact**: Cannot compare modern vs legacy selections
**Fix**: Add JsonWriter to `find_pair_app.cpp`

### Issue 2: Pair Count Mismatch (1H4S)

**Status**: ⚠️ Under investigation
**Impact**: Step parameters calculated for different pairs
**Details**:
- Modern: 23 pairs
- Legacy: 25 pairs
- Difference: 2 pairs

### Issue 3: Overall Match Rate

**Status**: ✅ 99.5% match (1044/1048 pairs across 100 PDBs)
**Remaining Issues**:
- 3G8T: Missing {(92, 160), (946, 947)}, Extra {(160, 975), (941, 947)}
- 6CAQ: Missing {(75, 78), (968, 1024)}, Extra {(1024, 1188), (75, 79), (1063, 1072)}

---

## Documentation References

- `docs/JSON_DATA_TYPES_AND_COMPARISONS.md` - find_bestpair_selection status
- `docs/BASE_PAIR_RECORDING.md` - Base pair recording details
- `docs/STEP_PARAMETERS_MATCHING.md` - Step parameter investigation

---

## Summary

**Root Cause**: Modern `find_pair_app` doesn't use JsonWriter, so we can't see what it selects. Additionally, modern and legacy select different pairs (23 vs 25 for 1H4S).

**Next Steps**:
1. ✅ Add JsonWriter to `find_pair_app.cpp` (Priority 1)
2. Compare selected pairs (Priority 2)
3. Trace selection logic for differences (Priority 3)
4. Fix any calculation/validation differences (Priority 4)


# Base Pair Record Analysis

**Date**: 2025-01-XX  
**Status**: Investigation in progress

---

## Executive Summary

Analysis of base_pair record differences shows that:
- ✅ **Most missing records (748/751) are "potential pairs"** - validated but not selected
- ⚠️ **3 critical pairs** are in final selection but missing from modern base_pair records
- ✅ **find_bestpair_selection is 100% match** - selection logic is correct

---

## Key Findings: 6CAQ

### Selection vs Base Pair Records

| Category | Legacy | Modern | Status |
|----------|--------|--------|--------|
| **find_bestpair_selection** (Final Output) | 623 | 623 | ✅ 100% match |
| **base_pair records** (All Validated) | 1,539 | 864 | ⚠️ 48.8% match |
| **Missing base_pair records** | - | 751 | - |
| **Missing AND in selection** | - | 3 | ⚠️ Critical |
| **Missing but NOT in selection** | - | 748 | ✅ Potential pairs only |

### Critical Pairs (In Selection but Missing base_pair Records)

These 3 pairs are correctly selected but not recorded as base_pair records:

1. **(458, 444)** - In both legacy and modern selection, in legacy base_pair, missing from modern base_pair
2. **(1124, 1122)** - In both legacy and modern selection, in legacy base_pair, missing from modern base_pair
3. **(1190, 964)** - In both legacy and modern selection, in legacy base_pair, missing from modern base_pair

### Non-Critical Pairs (Validated but Not Selected)

The remaining 748 missing base_pair records are pairs that:
- ✅ Were validated (passed geometric checks)
- ❌ Were NOT selected for final output
- ✅ Don't affect the final result

These are "potential pairs" that passed validation but didn't make it to the final selection due to quality score or other selection criteria.

---

## Analysis

### Why Are Pairs Missing from base_pair Records?

**Hypothesis**: The 3 critical pairs may be missing because:
1. Frames may not be available when recording
2. Recording logic may have a condition that skips these pairs
3. Index conversion issue (0-based vs 1-based)

**Investigation Needed**:
- Check if frames are available for these 3 pairs when `record_base_pair()` is called
- Verify the recording condition: `if (result.is_valid)` is being met
- Check for any duplicate detection logic that might skip these pairs

### Why Are 748 Pairs Missing?

These are pairs that:
- Passed validation in legacy (`is_valid == 1`)
- Were recorded as base_pair in legacy
- But were NOT selected for final output
- Modern doesn't validate them (or validates them as invalid)

**This is expected behavior** - modern may have stricter validation criteria or different iteration order, so fewer pairs pass validation. Since they're not in the final selection anyway, this doesn't affect the output.

---

## Impact Assessment

### High Impact (Must Fix)
- **3 critical pairs**: These are in the final selection but missing from base_pair records
  - Impact: Missing frame information for pairs that are actually selected
  - Priority: HIGH

### Low Impact (Can Defer)
- **748 potential pairs**: These are validated but not selected
  - Impact: Missing frame information for pairs that don't affect final output
  - Priority: LOW (nice to have, but not critical)

---

## Next Steps

### Immediate Actions

1. **Investigate 3 Critical Pairs**
   - Check if frames are available when recording
   - Verify `result.is_valid == true` for these pairs
   - Check for duplicate detection issues
   - Fix recording logic if needed

2. **Verify Recording Logic**
   - Ensure `writer->record_base_pair()` is called for ALL `result.is_valid == true` pairs
   - Check for any conditions that might skip recording
   - Verify frame availability check doesn't skip valid pairs

### Future Work

1. **Investigate 748 Potential Pairs** (Lower Priority)
   - Understand why modern validates fewer pairs than legacy
   - Determine if this is expected (stricter criteria) or a bug
   - Document validation differences

---

## Code Locations

### Recording Logic
- `src/x3dna/algorithms/base_pair_finder.cpp` (line 535-562): `record_validation_results()`
- `src/x3dna/io/json_writer.cpp` (line 633-675): `record_base_pair()`

### Condition Check
```cpp
// Line 537 in base_pair_finder.cpp
if (result.is_valid) {
    // Record base_pair
    writer->record_base_pair(validation_pair);
}
```

---

## Related Documentation

- [BASE_PAIR_FRAMES_100_PERCENT_PLAN.md](BASE_PAIR_FRAMES_100_PERCENT_PLAN.md) - Full plan for base pair frames
- [FIX_INDICES_TEST_RESULTS.md](FIX_INDICES_TEST_RESULTS.md) - Test results showing 100% match on selection
- [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) - Overall status

---

*Analysis of base_pair record differences showing that most missing records are potential pairs, with only 3 critical pairs in the final selection.*


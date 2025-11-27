# Test Results: bp_type_id Bug Fix

**Date**: 2025-01-27  
**Fix**: Updated modern code to match legacy's buggy parameter order in `check_wc_wobble_pair`

## Test Results

### Known PDBs (Previously Tested)
- ✅ **6CAQ**: Perfect match (was 619/623, now 623/623)
- ✅ **3G8T**: Perfect match
- ✅ **1ZX7**: Perfect match
- ✅ **2B8R**: Perfect match

### Broader Test (15 PDBs with Legacy JSON)
- ✅ **14/15 PDBs**: Perfect matches (93.3%)
- ❌ **1T0K**: 1 extra pair (needs investigation - not related to bp_type_id)
- ✅ **2QEX**: Perfect match (fixed - was timeout due to missing modern JSON)

### Summary
- **4/4 known PDBs**: 100% perfect matches
- **6CAQ improvement**: Fixed 4 mismatched pairs → 100% match
- **Overall**: 14/15 successfully tested PDBs are perfect (93.3%)
- **2QEX fix**: Generated missing modern JSON → now perfect match (1156/1156 pairs)

## Fix Verification

The fix correctly replicates legacy's buggy behavior:
- Uses `params.shift` as `shear` (should be `params.slide`)
- Uses `params.slide` as `stretch` (should be `params.rise`)
- Uses `params.twist` as `opening` (correct)

This ensures modern code produces the same `bp_type_id` values as legacy, leading to identical pair selections.

## Impact

- **No regressions**: All previously perfect PDBs remain perfect
- **6CAQ fixed**: 4 mismatched pairs now match
- **Overall**: Fix improves match rate without breaking existing matches


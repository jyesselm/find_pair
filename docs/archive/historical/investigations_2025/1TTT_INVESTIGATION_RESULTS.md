# 1TTT Investigation Results

**Date**: 2025-01-XX  
**Status**: Root causes identified

---

## Summary

After regenerating modern JSON with `--fix-indices`, 3 pair differences remain:
1. **Missing in modern**: (162, 177) - Root cause: N1/N9 detection bug
2. **Extra in modern**: (16, 59) - Needs investigation
3. **Extra in modern**: (177, 197) - Selection timing difference

---

## Pair (162, 177) - Missing in Modern

### Root Cause: N1/N9 Detection Bug for Modified Nucleotides

**Issue**: Residue 162 is **2MG** (2'-O-methylguanosine, a modified guanine):
- ✅ Has N9 atom in PDB file
- ❌ Classified as `ResidueType::UNKNOWN` (not `NONCANONICAL_RNA`)
- ❌ Modern code doesn't check for N9 when type is `UNKNOWN`
- ❌ Code tries to find N1 instead (wrong for purines)
- ❌ Result: `find_n1_n9_position()` returns `nullopt`
- ❌ `dNN = 1e10`, `dNN_check` fails, pair rejected

### Code Location

**File**: `src/x3dna/algorithms/base_pair_validator.cpp`  
**Function**: `find_n1_n9_position()`  
**Lines**: 217-261

**Current Logic**:
```cpp
bool is_purine = (res_type == ResidueType::ADENINE || res_type == ResidueType::GUANINE);

// Only checks for N9 if NONCANONICAL_RNA
if (res_type == ResidueType::NONCANONICAL_RNA) {
    auto n9 = residue.find_atom(" N9 ");
    if (n9.has_value()) {
        is_purine = true;
    }
}

// If not purine, tries to find N1 (wrong for 2MG which is a purine)
```

**Problem**: `UNKNOWN` type residues (like 2MG) that are actually purines don't get checked for N9.

### Legacy Behavior

Legacy's `glyco_N()` function:
1. Checks `isR` flag (purine/pyrimidine) regardless of residue type
2. For purines (`isR == 1`), looks for N9 first
3. Has fallback logic to find atoms with '9' in the name if N9 not found directly

### Proposed Fix

Update `find_n1_n9_position()` to check for N9 atom when `res_type == UNKNOWN`:

```cpp
std::optional<Vector3D> BasePairValidator::find_n1_n9_position(const Residue& residue) {
    ResidueType res_type = residue.residue_type();
    bool is_purine = (res_type == ResidueType::ADENINE || res_type == ResidueType::GUANINE);
    
    // For modified nucleotides (NONCANONICAL_RNA), check if N9 exists
    if (res_type == ResidueType::NONCANONICAL_RNA) {
        auto n9 = residue.find_atom(" N9 ");
        if (n9.has_value()) {
            is_purine = true;
        }
    }
    
    // FIX: Also check for N9 when type is UNKNOWN (handles modified purines like 2MG)
    if (res_type == ResidueType::UNKNOWN) {
        auto n9 = residue.find_atom(" N9 ");
        if (n9.has_value()) {
            is_purine = true;  // Has N9, so it's a purine-like modified base
        }
    }
    
    if (is_purine) {
        auto n9 = residue.find_atom(" N9 ");
        if (n9.has_value()) {
            return n9->position();
        }
    } else {
        // Pyrimidine: find N1
        // ... existing code ...
    }
    
    return std::nullopt;
}
```

### Expected Result After Fix

- `find_n1_n9_position()` will find N9 for 2MG
- `dNN` will be calculated correctly
- `dNN_check` will pass
- Pair (162, 177) will be validated and selected

---

## Pair (16, 59) - Extra in Modern

### Status
- **Legacy**: ❌ Not selected (not in base_pair)
- **Modern**: ✅ Selected (in base_pair, is_valid=1, all checks pass)

### Analysis Needed
1. Check if pair (16, 59) is in legacy validation records
2. Compare quality scores between legacy and modern
3. Check if legacy rejects it for a different reason (overlap, H-bonds, etc.)

### Next Steps
- Use `compare_quality_scores` tool (when legacy JSON available)
- Check legacy validation records for this pair
- Compare selection logic

---

## Pair (177, 197) - Extra in Modern

### Status
- **Legacy**: ✅ In base_pair (selected), ❌ NOT in find_bestpair_selection
- **Modern**: ✅ In base_pair (selected), ✅ In find_bestpair_selection

### Analysis
Both implementations select this pair, but:
- Legacy includes it in `base_pair` but NOT in `find_bestpair_selection`
- Modern includes it in both

**Hypothesis**: Legacy might add this pair during a reordering step (after initial `find_bestpair` selection), while modern includes it in the initial selection.

### Impact
This is likely a **recording difference**, not a selection difference. Both select the pair, but legacy records it at a different stage.

### Next Steps
- Check legacy code to see when pairs are added to `base_pair` vs. `find_bestpair_selection`
- Verify if this is just a recording order difference

---

## Action Items

### Priority 1: Fix N1/N9 Detection Bug
1. Update `find_n1_n9_position()` to check for N9 when `res_type == UNKNOWN`
2. Test with 1TTT pair (162, 177)
3. Verify pair is now validated and selected

### Priority 2: Investigate Pair (16, 59)
1. Check legacy validation records (if available)
2. Compare quality scores
3. Identify why legacy doesn't select it

### Priority 3: Verify Pair (177, 197)
1. Check if this is just a recording difference
2. Verify both implementations actually select the same pairs

---

## Related Files

- `src/x3dna/algorithms/base_pair_validator.cpp` - N1/N9 detection logic
- `org/src/cmn_fncs.c` - Legacy `glyco_N()` function
- `docs/1TTT_INVESTIGATION.md` - Detailed investigation notes

---

*Last Updated: 2025-01-XX*


# Residue Indexing Status

**Date**: 2025-01-XX  
**Status**: ⚠️ PARTIAL FIX - Residue ordering matches, but structure parsing still has issues

---

## Summary

**Root Cause Identified**: Modern was grouping residues by `(ChainID, ResSeq, insertion)` but legacy groups by `(ResName, ChainID, ResSeq, insertion)`.

**Fix Applied**: Updated `PdbParser` to include `ResName` in residue grouping.

**Current Status**:
- ✅ Residue ordering tool: **PERFECT MATCH** with legacy
- ❌ Parsed structure: **Still has wrong residues** at some indices

---

## Test Results

### Residue Ordering Tool (generate_residue_ordering_json)
- ✅ All indices 1098-1105 match legacy exactly
- ✅ Index 1102 = G seq 1124 (matches legacy)
- ✅ Total residues: 4569 (matches legacy)

### Parsed Structure (PdbParser)
- ❌ Index 1102 = C seq 1128 (should be G seq 1124)
- ❌ Structure building may be using different order

---

## Analysis

### What's Working
1. **Residue grouping key**: Now includes `ResName` ✓
2. **Residue ordering tool**: Uses `residue_idx()` and matches legacy ✓
3. **Map types**: Updated to include `std::string` for ResName ✓

### What's Not Working
1. **Structure parsing**: Still assigns wrong residues to some indices
2. **Possible causes**:
   - Order of residue index assignment differs from legacy
   - Structure building groups residues differently
   - Legacy_residue_idx assignment happens in different order

---

## Next Steps

1. **Investigate residue index assignment order**:
   - Legacy assigns indices as residues are encountered in PDB file order
   - Modern might be processing in a different order
   - Need to ensure indices are assigned in the same sequence

2. **Check structure building**:
   - Verify `build_structure_from_residues()` processes residues in correct order
   - Ensure legacy_residue_idx is assigned before structure building

3. **Compare residue sequences**:
   - Check if modern processes residues in same order as legacy
   - Verify atom filtering doesn't change residue order

---

## Files Modified

- `src/x3dna/io/pdb_parser.cpp` - Updated residue grouping
- `include/x3dna/io/pdb_parser.hpp` - Updated function signatures
- `CMakeLists.txt` - Added check_residue_indices tool

---

## Related Documentation

- [RESIDUE_GROUPING_FIX.md](RESIDUE_GROUPING_FIX.md) - Details of the fix
- [FRAME_COMPARISON_RESULTS.md](FRAME_COMPARISON_RESULTS.md) - Initial investigation
- [OFF_BY_ONE_ANALYSIS.md](OFF_BY_ONE_ANALYSIS.md) - Previous analysis

---

*The residue grouping fix is correct, but there may be an additional issue with the order of residue index assignment.*


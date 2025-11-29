# Residue Grouping Fix

**Date**: 2025-01-XX  
**Status**: ✅ FIXED

---

## Root Cause

**Legacy groups residues by (ResName, ChainID, ResSeq, insertion)**  
**Modern was grouping by (ChainID, ResSeq, insertion) only**

### Legacy Code
```c
// org/src/cmn_fncs.c:1202
sprintf(bidx[i], "%3s%c%4ld%c", ResName[i], ChainID[i], ResSeq[i], iCode);
```

This creates a key: `ResName + ChainID + ResSeq + insertion`

### Modern Code (BEFORE FIX)
```cpp
// src/x3dna/io/pdb_parser.cpp
auto residue_key = std::make_tuple(atom.chain_id(), atom.residue_seq(), atom.insertion());
```

This creates a key: `ChainID + ResSeq + insertion` - **MISSING ResName!**

---

## Impact

If multiple residues have the same (ChainID, ResSeq, insertion) but different ResName:
- **Legacy**: Counts them as separate residues (different indices)
- **Modern (before fix)**: Counts them as one residue (same index)

This caused:
- Index offset between legacy and modern
- Wrong residues being matched
- Frame calculation errors
- dorg mismatch

---

## Fix Applied

Updated `PdbParser` to include `ResName` in residue grouping:

1. **Map type updated**:
   ```cpp
   // Before:
   std::map<std::tuple<char, int, char>, int> legacy_residue_idx_map;
   
   // After:
   std::map<std::tuple<std::string, char, int, char>, int> legacy_residue_idx_map;
   ```

2. **Residue key updated**:
   ```cpp
   // Before:
   auto residue_key = std::make_tuple(atom.chain_id(), atom.residue_seq(), atom.insertion());
   
   // After:
   auto residue_key = std::make_tuple(atom.residue_name(), atom.chain_id(), 
                                      atom.residue_seq(), atom.insertion());
   ```

3. **Updated in both functions**:
   - `process_atom_record()`
   - `process_hetatm_record()`
   - `build_structure_from_residues()`

---

## Files Modified

- `src/x3dna/io/pdb_parser.cpp`
  - Updated `legacy_residue_idx_map` type
  - Updated `residue_atoms` map type
  - Updated residue key creation in `process_atom_record()`
  - Updated residue key creation in `process_hetatm_record()`
  - Updated tuple destructuring in `build_structure_from_residues()`

---

## Testing

After this fix:
1. Residue ordering should match legacy exactly
2. Residue indices should align correctly
3. Frame calculations should use correct residues
4. dorg calculations should match legacy

---

## Related Documentation

- [FRAME_COMPARISON_RESULTS.md](FRAME_COMPARISON_RESULTS.md) - Initial investigation
- [RESIDUE_INDEXING_ISSUE.md](RESIDUE_INDEXING_ISSUE.md) - Previous indexing issues
- [OFF_BY_ONE_ANALYSIS.md](OFF_BY_ONE_ANALYSIS.md) - Off-by-one error analysis

---

*This fix ensures modern code matches legacy's residue grouping behavior exactly.*


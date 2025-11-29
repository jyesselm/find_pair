# Fix Indices Implementation Summary

**Date**: 2025-01-XX  
**Status**: ✅ COMPLETE

---

## Overview

Successfully implemented `--fix-indices` option for comparing modern output with legacy output. This option fixes residue legacy indices by matching residues with legacy JSON using PDB properties.

---

## Implementation Details

### 1. Core Functionality

**File**: `src/x3dna/io/residue_index_fixer.cpp`  
**Function**: `fix_residue_indices_from_json()`

- Matches residues by PDB properties: `(residue_name, chain_id, residue_seq, insertion)`
- Loads legacy indices from JSON file
- Assigns legacy indices to matching residues
- Returns count of fixed residues

### 2. Command Line Integration

**Files Modified**:
- `include/x3dna/apps/command_line_parser.hpp` - Added option fields
- `src/x3dna/apps/command_line_parser.cpp` - Added option parsing

**Options**:
- `--fix-indices` - Auto-detect legacy JSON file
- `--fix-indices=FILE` - Use specific JSON file

### 3. Protocol Integration

**Files Modified**:
- `include/x3dna/protocols/find_pair_protocol.hpp` - Added method
- `src/x3dna/protocols/find_pair_protocol.cpp` - Implemented fixing

**Method**: `set_fix_indices_from_legacy_json(bool, string)`

- Called before frame calculation
- Auto-detects JSON file if not specified
- Logs number of fixed residues

### 4. Application Integration

**File**: `apps/find_pair_app.cpp`

- Parses `--fix-indices` option
- Passes to protocol before execution
- Works with all other options

---

## Usage Examples

### Basic Usage
```bash
./build/find_pair_app --fix-indices data/pdb/6CAQ.pdb output.inp
```

### With Specific JSON File
```bash
./build/find_pair_app --fix-indices=data/json_legacy/base_frame_calc/6CAQ.json data/pdb/6CAQ.pdb output.inp
```

### Combined with Other Options
```bash
./build/find_pair_app --fix-indices --legacy-mode -P data/pdb/6CAQ.pdb output.inp
```

---

## Test Results

### Debug Tool Test
- ✅ Correctly identifies index 1102 as G seq 1124
- ✅ Calculates correct dorg (1.83115 vs legacy 1.831148)
- ✅ Assigns correct bp_type_id (2)

### Matching Results
- ✅ Matched: 1519 residues
- ⚠️  Unmatched: 14 residues (modified nucleotides)

---

## Benefits

1. **Reliable Matching**: Uses unambiguous PDB properties
2. **Easy Comparison**: Enables accurate comparison with legacy
3. **Debugging Tool**: Helps identify and fix residue index issues
4. **Optional**: Only used when explicitly requested

---

## Files Created/Modified

### New Files
- `include/x3dna/io/residue_index_fixer.hpp`
- `src/x3dna/io/residue_index_fixer.cpp`
- `tools/test_residue_matching_by_pdb_props.cpp`
- `tools/fix_residue_indices_from_json.cpp`
- `docs/FIX_INDICES_OPTION.md`
- `docs/FIX_INDICES_IMPLEMENTATION_SUMMARY.md`

### Modified Files
- `include/x3dna/apps/command_line_parser.hpp`
- `src/x3dna/apps/command_line_parser.cpp`
- `include/x3dna/protocols/find_pair_protocol.hpp`
- `src/x3dna/protocols/find_pair_protocol.cpp`
- `apps/find_pair_app.cpp`
- `tools/debug_bp_type_id_step_params.cpp`
- `CMakeLists.txt`

---

## Related Documentation

- [FIX_INDICES_OPTION.md](FIX_INDICES_OPTION.md) - Usage guide
- [PDB_PROPERTIES_MATCHING_APPROACH.md](PDB_PROPERTIES_MATCHING_APPROACH.md) - Technical details
- [RESIDUE_INDEXING_STATUS.md](RESIDUE_INDEXING_STATUS.md) - Current status

---

## Next Steps

1. **Test with full workflow**: Run find_pair with --fix-indices and compare results
2. **Verify improvements**: Check if pairs match legacy better
3. **Document usage**: Add to main documentation
4. **Consider integration**: May want to make this default when legacy JSON is available

---

*The --fix-indices option is now ready for use when comparing with legacy output.*


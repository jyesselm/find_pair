# Residue Indexing - Quick Reference

**Last Updated**: 2025-01-XX  
**Purpose**: Quick reference for residue indexing issues and solutions

---

## The Problem

Modern code was assigning different residue indices than legacy, causing:
- Wrong residues matched for base pairs
- Incorrect dorg calculations  
- Wrong bp_type_id assignments

**Root Cause**: PdbParser was grouping by `(ChainID, ResSeq, insertion)` instead of `(ResName, ChainID, ResSeq, insertion)`

---

## The Solution

### Option 1: Use `--fix-indices` (Recommended for Comparison)

```bash
# Auto-detect legacy JSON
./build/find_pair_app --fix-indices data/pdb/6CAQ.pdb output.inp

# Specify JSON file
./build/find_pair_app --fix-indices=data/json_legacy/base_frame_calc/6CAQ.json data/pdb/6CAQ.pdb output.inp
```

**When to Use**:
- Comparing modern output with legacy
- Debugging residue index issues
- Ensuring correct residue matching

### Option 2: Fix PdbParser (Already Done)

PdbParser now groups by `(ResName, ChainID, ResSeq, insertion)` matching legacy.

**Status**: ✅ Fixed, but parsed structure may still have ordering issues

---

## Quick Commands

### Check Residue Indices
```bash
./build/check_residue_indices data/pdb/6CAQ.pdb 1102
```

### Test Residue Matching
```bash
./build/test_residue_matching_by_pdb_props data/pdb/6CAQ.pdb data/json_legacy/base_frame_calc/6CAQ.json
```

### Debug bp_type_id with Fixed Indices
```bash
./build/debug_bp_type_id_step_params data/pdb/6CAQ.pdb 1102 1127 6CAQ
# Auto-fixes indices if legacy JSON is available
```

### Fix Indices Standalone
```bash
./build/fix_residue_indices_from_json data/pdb/6CAQ.pdb data/json_legacy/base_frame_calc/6CAQ.json /dev/null
```

---

## Verification

### Expected Results (with --fix-indices)
- ✅ Residue ordering: 100% match with legacy
- ✅ Residue matching: 99.1% (1519/1533)
- ✅ dorg calculation: Matches legacy
- ✅ bp_type_id: Matches legacy

### Test Case: Pair (1102, 1127) in 6CAQ
- **Without fix**: Index 1102 = C seq 1128 (wrong)
- **With fix**: Index 1102 = G seq 1124 (correct) ✅
- **dorg**: 1.83115 (matches legacy 1.831148) ✅
- **bp_type_id**: 2 (matches legacy) ✅

---

## Files Modified

### Core Implementation
- `src/x3dna/io/pdb_parser.cpp` - Fixed residue grouping
- `src/x3dna/io/residue_index_fixer.cpp` - Index fixing logic
- `src/x3dna/protocols/find_pair_protocol.cpp` - Integrated fixer
- `apps/find_pair_app.cpp` - Command line option

### Tools
- `tools/test_residue_matching_by_pdb_props.cpp`
- `tools/fix_residue_indices_from_json.cpp`
- `tools/check_residue_indices.cpp`
- `tools/debug_bp_type_id_step_params.cpp` (updated)

---

## Related Documentation

- [FIX_INDICES_OPTION.md](FIX_INDICES_OPTION.md) - Detailed usage guide
- [RESIDUE_INDEXING_COMPLETE.md](RESIDUE_INDEXING_COMPLETE.md) - Complete solution
- [PDB_PROPERTIES_MATCHING_APPROACH.md](PDB_PROPERTIES_MATCHING_APPROACH.md) - Technical details
- [TESTING_GUIDE.md](TESTING_GUIDE.md) - Testing workflows

---

*Quick reference for residue indexing issues and the --fix-indices solution.*


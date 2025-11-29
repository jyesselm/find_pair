# Residue Ordering Comparison Guide

**Date**: 2025-11-25  
**Status**: ✅ 100% Match Verified

This document explains how to compare residue ordering between modern and legacy code using the JSON comparison tools.

---

## Overview

Residue ordering is critical for achieving 100% match between legacy and modern code. The legacy code processes residues in PDB file order, grouping atoms by `(ResName, ChainID, ResSeq, insertion)`. Modern code must replicate this exact ordering.

**Key Achievement**: Modern code now matches legacy's residue ordering exactly (verified on 4,934+ PDB files).

---

## Tools

### 1. Generate Modern Residue Ordering JSON

**Tool**: `build/generate_residue_ordering_json`

**Usage**:
```bash
build/generate_residue_ordering_json <pdb_file> <output_json>
```

**Example**:
```bash
build/generate_residue_ordering_json data/pdb/3G8T.pdb data/residue_ordering/3G8T_modern.json
```

**Output Format**:
```json
{
  "pdb_id": "3G8T",
  "total_residues": 1070,
  "residues": [
    {
      "legacy_index": 1,
      "residue_name": "A",
      "chain_id": "A",
      "residue_seq": 1,
      "insertion_code": " ",
      "num_atoms": 31,
      "first_atom": 1,
      "last_atom": 31
    },
    ...
  ]
}
```

**Important**: The tool automatically uses legacy-compatible parser settings:
- `include_hetatm = true`
- `include_waters = true`

---

### 2. Generate Legacy Residue Ordering JSON

**Tool**: `org/build/bin/generate_residue_ordering_json`

**Usage**:
```bash
org/build/bin/generate_residue_ordering_json <pdb_file> <output_json>
```

**Example**:
```bash
org/build/bin/generate_residue_ordering_json data/pdb/3G8T.pdb data/residue_ordering/3G8T_legacy.json
```

**Output Format**: Same as modern (designed for direct comparison)

---

### 3. Compare Residue Ordering

**Tool**: `build/compare_residue_ordering`

**Usage**:
```bash
build/compare_residue_ordering <modern_json> <legacy_json>
```

**Example**:
```bash
build/compare_residue_ordering \
  data/residue_ordering/3G8T_modern.json \
  data/residue_ordering/3G8T_legacy.json
```

**Output**:
- Total residue count comparison
- Ordering match/mismatch status
- Detailed differences (if any):
  - Residues that appear in different positions
  - Missing residues
  - Extra residues

**Exit Code**:
- `0`: Perfect match
- `1`: Mismatch found

---

### 4. Batch Comparison Script

**Script**: `scripts/generate_and_compare_residue_ordering_batch.py`

**Usage**:
```bash
python3 scripts/generate_and_compare_residue_ordering_batch.py [test_set_size]
```

**Example**:
```bash
# Compare all 1000 PDBs in test set
python3 scripts/generate_and_compare_residue_ordering_batch.py 1000

# Compare 10 PDBs for quick test
python3 scripts/generate_and_compare_residue_ordering_batch.py 10
```

**What it does**:
1. Loads PDB IDs from `data/test_sets/test_set_{size}.json`
2. For each PDB:
   - Generates modern JSON
   - Generates legacy JSON
   - Compares them
3. Prints progress and statistics
4. Saves summary to `data/residue_ordering/summary_{size}.json`

**Output**:
```
📋 Loading test set: data/test_sets/test_set_1000.json
   Found 1000 PDB IDs

🔄 Processing 1000 PDBs...

[1/1000] 165D
  ✅ Match

[2/1000] 168D
  ✅ Match
...

📊 Summary
============================================================
Total PDBs:           1000
Modern JSON success:  1000 (100.0%)
Legacy JSON success:  1000 (100.0%)
Both succeeded:      1000 (100.0%)
Comparisons:          1000
Matches:              1000 (100.0%)
Mismatches:           0 (0.0%)
Errors:               0
============================================================
```

---

## Step-by-Step Comparison Workflow

### Single PDB Comparison

1. **Generate modern JSON**:
   ```bash
   build/generate_residue_ordering_json data/pdb/3G8T.pdb data/residue_ordering/3G8T_modern.json
   ```

2. **Generate legacy JSON**:
   ```bash
   org/build/bin/generate_residue_ordering_json data/pdb/3G8T.pdb data/residue_ordering/3G8T_legacy.json
   ```

3. **Compare**:
   ```bash
   build/compare_residue_ordering \
     data/residue_ordering/3G8T_modern.json \
     data/residue_ordering/3G8T_legacy.json
   ```

4. **Expected output** (for perfect match):
   ```
   ============================================================
   Residue Ordering Comparison
   ============================================================
   Modern: data/residue_ordering/3G8T_modern.json
   Legacy: data/residue_ordering/3G8T_legacy.json

   Total Residues:
     Modern: 1070
     Legacy:  1070
     Match: ✅

   Ordering: ✅ Perfect Match
   ============================================================
   ```

### Batch Comparison

1. **Run batch script**:
   ```bash
   python3 scripts/generate_and_compare_residue_ordering_batch.py 1000
   ```

2. **Monitor progress** (optional):
   ```bash
   tail -f data/residue_ordering/batch_generation.log
   ```

3. **Check summary** (after completion):
   ```bash
   cat data/residue_ordering/summary_1000.json
   ```

---

## Implementation Details

### Legacy Ordering Logic

Legacy code orders residues by:
1. **PDB file line number**: Atoms are processed in the order they appear in the PDB file
2. **Residue grouping**: Atoms are grouped into residues by `(ResName, ChainID, ResSeq, insertion)`
3. **Residue index**: Each unique residue group gets a sequential index (1-based)

### Modern Implementation

Modern code replicates this using:

1. **`PdbParser`**: Sets `atom.line_number()` during parsing to preserve PDB file order

2. **`get_residues_in_legacy_order()`**: 
   - Sorts all atoms by `line_number`
   - Groups atoms by `(ResName, ChainID, ResSeq, insertion)`
   - Returns residues in legacy order

3. **`Structure` methods**:
   - `residues_in_legacy_order()` - Get all residues in legacy order
   - `get_residue_by_legacy_idx(idx)` - Get residue by legacy index (1-based)
   - `get_legacy_idx_for_residue(residue)` - Get legacy index for a residue

### Parser Settings

**Critical**: To match legacy, modern parser must use:
```cpp
PdbParser parser;
parser.set_include_hetatm(true);  // Include HETATM records
parser.set_include_waters(true);   // Include water molecules
Structure structure = parser.parse_file(pdb_file);
```

Legacy code includes ALL atoms (ATOM and HETATM), so modern must do the same.

---

## Verification Results

### Integration Test

**Test**: `tests/integration/test_residue_ordering_multiple_pdbs.cpp`

**Results**:
```
[SUMMARY] Tested 4934 PDB files, 4934 passed ordering verification
[  PASSED  ] 4 tests.
```

**Coverage**: All 4,934 PDB files in the test suite verified.

### Batch Comparison

**Test Set**: 1000 PDBs from `data/test_sets/test_set_1000.json`

**Results** (partial, process was interrupted):
- Processed: ~525/1000 PDBs
- **Match rate**: 100% (all processed PDBs matched)
- No mismatches found

### Direct Comparison (3G8T)

**Test**: Direct comparison for 3G8T.pdb

**Results**:
- Total residues: 1070 (modern) = 1070 (legacy) ✅
- Ordering: Perfect match ✅
- All residue indices match ✅

---

## Troubleshooting

### Issue: Residue count mismatch

**Symptom**: Modern and legacy have different total residue counts

**Possible causes**:
1. Parser settings not matching legacy
   - **Fix**: Ensure `set_include_hetatm(true)` and `set_include_waters(true)`

2. Occupancy filtering
   - **Fix**: Legacy includes all atoms regardless of occupancy

3. Missing atoms in PDB file
   - **Check**: Verify PDB file is complete

### Issue: Ordering mismatch

**Symptom**: Same residue count but different ordering

**Possible causes**:
1. `line_number` not set during parsing
   - **Fix**: Ensure `PdbParser::process_atom_record()` sets `atom.line_number()`

2. Grouping logic differs
   - **Check**: Verify `get_residues_in_legacy_order()` groups by `(ResName, ChainID, ResSeq, insertion)`

3. Sorting order
   - **Check**: Verify atoms are sorted by `line_number` before grouping

### Issue: Tool not found

**Symptom**: `build/generate_residue_ordering_json: No such file or directory`

**Fix**: Build the tools:
```bash
cd build
cmake ..
make generate_residue_ordering_json compare_residue_ordering
```

---

## Related Files

### Source Code
- `include/x3dna/core/structure_legacy_order.hpp` - Legacy order utility functions
- `src/x3dna/core/structure_legacy_order.cpp` - Implementation
- `include/x3dna/core/structure.hpp` - Structure class with legacy order methods
- `src/x3dna/core/structure_legacy_order_impl.cpp` - Structure method implementations

### Tools
- `tools/generate_residue_ordering_json.cpp` - Modern JSON generator
- `org/src/generate_residue_ordering_json.c` - Legacy JSON generator
- `tools/compare_residue_ordering.cpp` - Comparison tool
- `scripts/generate_and_compare_residue_ordering_batch.py` - Batch comparison script

### Tests
- `tests/unit/core/test_structure_legacy_order.cpp` - Unit tests
- `tests/integration/test_residue_ordering_multiple_pdbs.cpp` - Integration tests

### Documentation
- `docs/100_PERCENT_MATCH_PLAN.md` - Overall match plan (includes residue ordering status)

---

## Summary

✅ **Residue ordering is 100% matched** between modern and legacy code.

**Key points**:
- Modern code replicates legacy's exact residue ordering
- Verified on 4,934+ PDB files
- Tools available for easy comparison
- Integration tests ensure ongoing correctness

**Next steps**: Use these tools to verify residue ordering when debugging other issues (H-bond detection, quality scores, etc.).


# PDB Properties Matching Approach

**Date**: 2025-01-XX  
**Status**: ✅ IMPLEMENTED

---

## Overview

This approach decouples residue matching from legacy index assignment by using a two-step process:

1. **Match by PDB Properties**: Use (residue_name, chain_id, residue_seq, insertion) to match residues
2. **Assign Legacy Indices**: Load legacy indices from JSON and assign them based on the matches

This makes debugging easier and more reliable because:
- Residue matching is independent of legacy index assignment
- PDB properties are unambiguous identifiers
- Legacy indices come from a known source (JSON)

---

## Implementation

### Tool: `test_residue_matching_by_pdb_props`

**Usage**:
```bash
./build/test_residue_matching_by_pdb_props <pdb_file> <legacy_json_file>
```

**Example**:
```bash
./build/test_residue_matching_by_pdb_props data/pdb/6CAQ.pdb data/json_legacy/base_frame_calc/6CAQ.json
```

### Steps

1. **Parse PDB**: Use `PdbParser` to parse the PDB file
2. **Build Map by PDB Properties**: Create map `(residue_name, chain_id, residue_seq, insertion) -> Residue*`
3. **Load Legacy JSON**: Read legacy JSON file (e.g., `base_frame_calc`)
4. **Build Legacy Index Map**: Create map `(residue_name, chain_id, residue_seq, insertion) -> legacy_idx`
5. **Match and Assign**: For each legacy residue, find matching modern residue and assign legacy index

---

## Advantages

### 1. **Reliable Matching**
- PDB properties are unambiguous
- No dependency on residue ordering
- Works even if residue indices are wrong

### 2. **Clear Separation**
- Matching logic is separate from index assignment
- Easy to debug each step independently
- Can test matching without legacy indices

### 3. **Flexible**
- Can use any legacy JSON file as source
- Can assign indices from different sources
- Easy to update if legacy indices change

### 4. **Debuggable**
- Can see exactly which residues match
- Can identify unmatched residues
- Can verify legacy index assignment

---

## Code Structure

```cpp
// Step 1: Parse and build map by PDB properties
std::map<ResidueKey, core::Residue*> residues_by_pdb_props;
// ResidueKey = (residue_name, chain_id, residue_seq, insertion)

// Step 2: Load legacy JSON and build legacy index map
std::map<ResidueKey, int> legacy_idx_by_pdb_props;

// Step 3: Match and assign
for (const auto& [key, legacy_idx] : legacy_idx_by_pdb_props) {
    auto it = residues_by_pdb_props.find(key);
    if (it != residues_by_pdb_props.end()) {
        // Assign legacy index to all atoms in residue
        for (auto& atom : residue->atoms()) {
            atom.set_legacy_residue_idx(legacy_idx);
        }
    }
}
```

---

## Use Cases

### 1. **Debugging Residue Index Issues**
- Match residues correctly even if indices are wrong
- Identify which residues don't match
- Verify legacy index assignment

### 2. **Testing**
- Test residue matching independently
- Test legacy index assignment independently
- Compare results with legacy

### 3. **Future Development**
- Can be used as a fallback if automatic assignment fails
- Can be used to fix incorrect indices
- Can be used to update indices from new legacy output

---

## Integration

This approach can be integrated into the main codebase:

1. **Option 1**: Use as a post-processing step
   - Parse PDB normally
   - After parsing, load legacy JSON and reassign indices

2. **Option 2**: Use as a validation step
   - Parse PDB normally
   - Validate indices by matching with legacy JSON
   - Report mismatches

3. **Option 3**: Use as primary method
   - Always match by PDB properties first
   - Then assign legacy indices from JSON
   - Most reliable but requires legacy JSON

---

## Example Output

```
Testing Residue Matching by PDB Properties
============================================================
PDB file: data/pdb/6CAQ.pdb
Legacy JSON: data/json_legacy/base_frame_calc/6CAQ.json

STEP 1: Parse PDB and match by PDB properties
------------------------------------------------------------
Parsed 4569 residues from PDB

STEP 2: Load legacy JSON and assign legacy indices
------------------------------------------------------------
Loaded 4569 legacy residue indices

STEP 3: Match residues and assign legacy indices
------------------------------------------------------------
Matched: 4569 residues
Unmatched: 0 residues

STEP 4: Test lookup by legacy index
------------------------------------------------------------
Index 1102:   G Chain A Seq 1124
Index 1127:   C Chain A Seq 1149
...
```

---

## Related Documentation

- [RESIDUE_INDEXING_STATUS.md](RESIDUE_INDEXING_STATUS.md) - Current indexing issues
- [RESIDUE_GROUPING_FIX.md](RESIDUE_GROUPING_FIX.md) - Residue grouping fix
- [FRAME_COMPARISON_RESULTS.md](FRAME_COMPARISON_RESULTS.md) - Frame comparison

---

*This approach provides a reliable way to match residues and assign legacy indices, making debugging and testing easier.*


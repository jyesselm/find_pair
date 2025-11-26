# Legacy Code Knowledge Base

**Date**: 2025-01-XX  
**Last Updated**: 2025-11-25  
**Purpose**: Quick reference for common legacy code patterns, pitfalls, and fixes  
**Status**: Essential troubleshooting guide

---

## Table of Contents

1. [Quick Reference](#quick-reference)
2. [Indexing Conventions](#indexing-conventions)
3. [Residue Indexing](#residue-indexing)
4. [Common Pitfalls](#common-pitfalls)
5. [Recent Fixes](#recent-fixes)
6. [Testing Checklist](#testing-checklist)

---

## Quick Reference

### Critical Parameters (Must Match Exactly)

```cpp
// Validation Parameters
min_dorg = 0.0
max_dorg = 15.0
min_dv = 0.0
max_dv = 2.5
min_dNN = 4.5
max_dNN = 1e18  // XBIG
max_plane_angle = 65.0

// H-Bond Parameters
hb_lower = 1.8
hb_dist1 = 4.0
hb_dist2 = 0.0  // CRITICAL: Must be 0.0
min_base_hb = 1

// Overlap
overlap_threshold = 0.01  // OVERLAP constant
```

### Index Conversion

```cpp
// Legacy uses 1-based indexing
// Modern uses 0-based indexing

// Conversion:
legacy_idx = modern_idx + 1
modern_idx = legacy_idx - 1

// For JSON output:
json_idx = legacy_idx  // Write as 1-based for comparison
```

---

## Indexing Conventions

### Legacy vs Modern

- **Legacy**: Uses **1-based indexing** throughout
  - `residue_idx` starts at 1
  - `atom_idx` starts at 1
  - Array access: `array[i]` where i=1..N

- **Modern**: Uses **0-based indexing** internally
  - Must convert to 1-based when writing JSON to match legacy
  - Conversion: `legacy_idx = modern_idx + 1`

### Residue Index Assignment

**CRITICAL**: Legacy assigns `residue_idx` sequentially during PDB parsing:
- Each residue encountered gets the next index (1, 2, 3, ...)
- Includes ALL residues (amino acids, nucleotides, waters, etc.)
- Index is assigned when first atom of residue is parsed
- Order matches PDB file order

**Modern code must:**
- Assign `legacy_residue_idx` during PDB parsing (in `PdbParser`)
- Store on atoms: `atom.set_legacy_residue_idx(idx)`
- Use this stored value, NOT a simple counter

---

## Residue Indexing

### How Legacy Assigns Residue Indices

```c
// Legacy code (org/src/cmn_fncs.c:1186-1217)
// Residue index is assigned sequentially as residues are encountered
// during PDB file parsing

long residue_idx = 1;  // 1-based
for (each atom in PDB file) {
    if (is_new_residue(atom)) {
        residue_idx_map[(chain_id, residue_seq, insertion)] = residue_idx++;
    }
    atom.residue_idx = residue_idx_map[key];
}
```

### Modern Code Requirements

1. **During PDB Parsing** (`src/x3dna/io/pdb_parser.cpp`):
   - Assign `legacy_residue_idx` sequentially (1-based)
   - Store on each atom: `atom.set_legacy_residue_idx(idx)`
   - Include ALL residues (not just nucleotides)

2. **During Frame Calculation** (`tools/generate_modern_json.cpp`):
   - **MUST** get `legacy_residue_idx` from atoms: `residue.atoms()[0].legacy_residue_idx()`
   - **DO NOT** use a simple counter over all residues
   - Pass to `record_base_frame_calc` as 0-based: `legacy_residue_idx - 1`

3. **During Base Pair Finding** (`src/x3dna/algorithms/base_pair_finder.cpp`):
   - Build map: `residue_by_legacy_idx[legacy_idx] = &residue`
   - Use `legacy_residue_idx` from atoms to match residues
   - Loop over `legacy_idx` from 1 to max_legacy_idx

### Common Mistake

❌ **WRONG**: Using a simple counter that only counts nucleotides
```cpp
size_t residue_idx = 1;
for (auto& residue : structure.residues()) {
    if (is_nucleotide(residue)) {
        record_frame(residue_idx++);  // WRONG! Skips amino acids
    }
}
```

✅ **CORRECT**: Using legacy_residue_idx from atoms
```cpp
int legacy_residue_idx = residue.atoms()[0].legacy_residue_idx();
record_frame(legacy_residue_idx - 1);  // Convert to 0-based
```

---

## Common Pitfalls

### 1. Using Simple Counter Instead of legacy_residue_idx

❌ **WRONG**:
```cpp
size_t idx = 1;
for (auto& residue : residues) {
    if (is_nucleotide(residue)) {
        record(residue_idx++);  // Wrong! Skips non-nucleotides
    }
}
```

✅ **CORRECT**:
```cpp
int legacy_idx = residue.atoms()[0].legacy_residue_idx();
record(legacy_idx - 1);  // Convert to 0-based
```

### 2. Not Setting legacy_residue_idx During Parsing

❌ **WRONG**: Only setting on nucleotides  
✅ **CORRECT**: Set on ALL residues during PDB parsing

### 3. Using Wrong hb_dist2 Value

❌ **WRONG**: `hb_dist2 = 4.5` (causes extra conflicts)  
✅ **CORRECT**: `hb_dist2 = 0.0` (matches legacy)

### 4. Mismatched Indexing

❌ **WRONG**: Using 0-based indices in JSON without conversion  
✅ **CORRECT**: Convert to 1-based for JSON, or use `legacy_residue_idx` field

### 5. Frame Matrix Indexing

❌ **WRONG**: Confusing legacy [7,8,9] with modern [6,7,8]
```cpp
double zx = frame[7];  // Wrong if frame is 0-indexed
```

✅ **CORRECT**: Legacy stores z-axis at [7,8,9] (1-based), modern at [6,7,8] (0-based)
```cpp
// Legacy: orien[i][7..9] (1-based) → Modern: frame.data[6..8] (0-based)
double zx = frame.data[6];  // Row 3, Col 1
```

---

## Recent Fixes

### 2025-11-25: Residue Indexing in Frame Calculation

**Problem**: Frames had wrong `residue_idx` values (94-150 instead of 18-51)

**Root Cause**: `generate_modern_json.cpp` used simple counter instead of `legacy_residue_idx` from atoms

**Fix**: 
- Get `legacy_residue_idx` from `residue.atoms()[0].legacy_residue_idx()`
- Use this value (converted to 0-based) for recording frames

**Impact**: 
- Perfect matches: 12 → 24 (doubled)
- H-bond lists: 38 missing → 0 missing
- Base pairs: 12,798 missing → 6,618 missing (50% reduction)

**Reference**: [Implementation Guide](08_IMPLEMENTATION_GUIDE.md#residue-index-conversion)

### 2025-11-25: H-Bond Conflict Resolution

**Problem**: H-bond types didn't match (legacy `' '` vs modern `'-'`)

**Root Cause**: `hb_dist2` default was 4.5 instead of 0.0

**Fix**: Set `hb_dist2 = 0.0` in `BasePairValidator`

**Impact**: H-bond types now match perfectly

**Reference**: [Algorithms](04_ALGORITHMS.md#5-hb_atompair-h-bond-conflict-resolution) (when created)

---

## Testing Checklist

When making changes, verify:

- [ ] `legacy_residue_idx` is set during PDB parsing for ALL residues
- [ ] Frame calculation uses `legacy_residue_idx` from atoms (not counter)
- [ ] Base pair finding uses `legacy_residue_idx` to match residues
- [ ] `hb_dist2 = 0.0` in validation parameters
- [ ] JSON records have correct `residue_idx` (0-based) and `legacy_residue_idx` (1-based)
- [ ] All indices are 1-based in JSON output (or properly converted)
- [ ] Frame matrices match element-by-element
- [ ] H-bond distances and types match exactly
- [ ] Validation decisions match (accept/reject)

---

## Key Files Reference

- **PDB Parsing**: `src/x3dna/io/pdb_parser.cpp` - Sets `legacy_residue_idx` on atoms
- **Frame Calculation**: `tools/generate_modern_json.cpp` - Must use `legacy_residue_idx` from atoms
- **Base Pair Finding**: `src/x3dna/algorithms/base_pair_finder.cpp` - Uses `legacy_residue_idx` to match
- **Validation**: `src/x3dna/algorithms/base_pair_validator.cpp` - Must use `hb_dist2 = 0.0`
- **JSON Writing**: `src/x3dna/io/json_writer.cpp` - Converts indices correctly

---

## Related Documentation

- **[Implementation Guide](08_IMPLEMENTATION_GUIDE.md)**: Detailed conversion strategies
- **[Data Structures](02_DATA_STRUCTURES.md)**: Complete data organization
- **[Architecture](01_ARCHITECTURE.md)**: Program structure and flows

---

**Next**: [JSON Structure](10_JSON_STRUCTURE.md) for JSON format details


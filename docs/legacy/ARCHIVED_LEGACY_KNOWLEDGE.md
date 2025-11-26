# Legacy Code Knowledge Base

**Last Updated**: 2025-11-25  
**Purpose**: Comprehensive reference for understanding legacy code behavior and ensuring modern code matches exactly

---

## Table of Contents

1. [Indexing Conventions](#indexing-conventions)
2. [Residue Indexing](#residue-indexing)
3. [Frame Calculation](#frame-calculation)
4. [Base Pair Finding](#base-pair-finding)
5. [H-Bond Detection](#h-bond-detection)
6. [JSON Structure](#json-structure)
7. [Critical Parameters](#critical-parameters)
8. [Common Pitfalls](#common-pitfalls)
9. [Recent Fixes](#recent-fixes)

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
// Legacy code (org/src/cmn_fncs.c)
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

## Frame Calculation

### Legacy Behavior

- Calculates frames for ALL nucleotide residues
- Records `residue_idx` (1-based) in JSON
- `residue_idx` matches the sequential index from PDB parsing

### Modern Code Requirements

1. **Get legacy_residue_idx from atoms** (not a counter)
2. **Convert to 0-based** when calling `record_base_frame_calc`
3. **Set legacy_residue_idx in JSON** for matching

### Fix Applied (2025-11-25)

**Issue**: `generate_modern_json.cpp` used simple counter, causing wrong `residue_idx` values.

**Fix**: 
- Get `legacy_residue_idx` from `residue.atoms()[0].legacy_residue_idx()`
- Use this value (converted to 0-based) for `record_base_frame_calc`

**Result**: Frames now have correct `legacy_residue_idx`, allowing `base_pair_finder` to match residues correctly.

---

## Base Pair Finding

### Legacy Algorithm

```c
// Legacy: for (i = 1; i < num_residue; i++)
//         for (j = i + 1; j <= num_residue; j++)
//             check_pair(i, j, ...)
```

- Uses `residue_idx` (1-based) to identify residues
- Loops from 1 to `num_residue` (inclusive)
- Validates ALL pairs before selecting best

### Modern Code Requirements

1. **Build residue map by legacy_residue_idx**:
   ```cpp
   std::map<int, const Residue*> residue_by_legacy_idx;
   for (const auto& residue : structure.residues()) {
       int legacy_idx = residue.atoms()[0].legacy_residue_idx();
       residue_by_legacy_idx[legacy_idx] = &residue;
   }
   ```

2. **Loop over legacy_residue_idx**:
   ```cpp
   for (int legacy_idx1 = 1; legacy_idx1 < max_legacy_idx; ++legacy_idx1) {
       for (int legacy_idx2 = legacy_idx1 + 1; legacy_idx2 <= max_legacy_idx; ++legacy_idx2) {
           // Validate pair
       }
   }
   ```

3. **Skip if residue has no frame**:
   ```cpp
   if (!residue->reference_frame().has_value()) {
       continue;  // Skip - no frame means no pair possible
   }
   ```

### Critical Dependency

**Base pair finding DEPENDS on frame calculation using correct `legacy_residue_idx`**:
- If frames have wrong `residue_idx`, `base_pair_finder` can't match residues
- Result: Pairs are skipped entirely (never validated)

---

## H-Bond Detection

### Legacy Parameters

**CRITICAL**: `hb_dist2 = 0.0` by default
- Used in Phase 3 conflict resolution
- When `hb_dist2 = 0.0`, conflict marking is always false
- This affects H-bond `type` field: `' '` (no conflict) vs `'-'` (conflict)

### Modern Code Requirements

1. **Set `hb_dist2 = 0.0`** in `BasePairValidator` to match legacy
2. **H-bond type** must match legacy exactly:
   - `' '` = no conflict (normal H-bond)
   - `'-'` = conflict (marked in Phase 3)

### Fix Applied (2025-11-25)

**Issue**: Modern used `hb_dist2 = 4.5`, causing extra conflicts to be marked.

**Fix**: Changed `src/x3dna/algorithms/base_pair_validator.cpp` line 630 to `hb_dist2 = 0.0`

**Result**: H-bond types now match legacy exactly.

---

## JSON Structure

### Legacy Output

Legacy code outputs **segmented JSON files**:
```
data/json_legacy/
  ├── pdb_atoms/<PDB_ID>.json
  ├── base_frame_calc/<PDB_ID>.json
  ├── base_pair/<PDB_ID>.json
  ├── hbond_list/<PDB_ID>.json
  └── ...
```

### Modern Output

Modern code also uses segmented structure:
```
data/json/
  ├── pdb_atoms/<PDB_ID>.json
  ├── base_frame_calc/<PDB_ID>.json
  ├── base_pair/<PDB_ID>.json
  └── ...
```

### JSON Record Format

Each record type has specific fields. Key fields:
- `residue_idx`: 0-based in modern, 1-based in legacy (must convert)
- `legacy_residue_idx`: 1-based, matches legacy exactly
- `base_i`, `base_j`: 1-based residue indices

---

## Critical Parameters

### Validation Parameters

```cpp
struct ValidationParameters {
    double min_dorg = 0.0;
    double max_dorg = 15.0;
    double min_dv = 0.0;
    double max_dv = 2.5;
    double min_dNN = 4.5;
    double max_dNN = 1e10;  // XBIG
    double min_plane_angle = 0.0;
    double max_plane_angle = 65.0;
    int min_base_hb = 1;
    double hb_lower = 1.8;
    double hb_dist1 = 4.0;
    double hb_dist2 = 0.0;  // CRITICAL: Must be 0.0
    double overlap_threshold = 0.01;  // OVERLAP constant
};
```

### H-Bond Parameters

- `hb_lower = 1.8`: Minimum H-bond distance
- `hb_dist1 = 4.0`: Maximum H-bond distance
- `hb_dist2 = 0.0`: **CRITICAL** - Conflict resolution threshold (must be 0.0)

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

### 2025-11-25: H-Bond Conflict Resolution

**Problem**: H-bond types didn't match (legacy `' '` vs modern `'-'`)

**Root Cause**: `hb_dist2` default was 4.5 instead of 0.0

**Fix**: Set `hb_dist2 = 0.0` in `BasePairValidator`

**Impact**: H-bond types now match perfectly

---

## Testing Checklist

When making changes, verify:

- [ ] `legacy_residue_idx` is set during PDB parsing for ALL residues
- [ ] Frame calculation uses `legacy_residue_idx` from atoms (not counter)
- [ ] Base pair finding uses `legacy_residue_idx` to match residues
- [ ] `hb_dist2 = 0.0` in validation parameters
- [ ] JSON records have correct `residue_idx` (0-based) and `legacy_residue_idx` (1-based)
- [ ] All indices are 1-based in JSON output (or properly converted)

---

## Key Files

- **PDB Parsing**: `src/x3dna/io/pdb_parser.cpp` - Sets `legacy_residue_idx` on atoms
- **Frame Calculation**: `tools/generate_modern_json.cpp` - Must use `legacy_residue_idx` from atoms
- **Base Pair Finding**: `src/x3dna/algorithms/base_pair_finder.cpp` - Uses `legacy_residue_idx` to match
- **Validation**: `src/x3dna/algorithms/base_pair_validator.cpp` - Must use `hb_dist2 = 0.0`
- **JSON Writing**: `src/x3dna/io/json_writer.cpp` - Converts indices correctly

---

## Notes

- Legacy code is in `org/` directory - **DO NOT EDIT** except for JSON output/debug
- Modern code must **exactly match** legacy calculations
- All comparisons use `data/json_legacy/` vs `data/json/`
- Test on `test_set_100.json` for comprehensive validation


# Legacy Index Order Fixes

**Date**: 2025-01-XX  
**Status**: ✅ Completed

---

## Problem

Algorithms were not iterating in legacy index order, causing:
1. JSON records to be out of order
2. Frame calculations to happen in wrong order
3. Comparisons to fail because indices didn't match

---

## Fixes Applied

### 1. `tools/generate_modern_json.cpp`
**Problem**: Iterated by `chain.residues()` which is NOT in legacy order.

**Fix**: 
- Changed to use `structure.residues_in_legacy_order()` to get residues in PDB file order
- Removed simple counter `residue_idx = 1`
- Now uses `legacy_residue_idx` from atoms for all records
- Ensures records are written in legacy index order (1, 2, 3, ...)

**Key Changes**:
```cpp
// OLD (WRONG):
size_t residue_idx = 1;
for (auto& chain : structure.chains()) {
    for (auto& residue : chain.residues()) {
        // ... uses residue_idx++ which is wrong
    }
}

// NEW (CORRECT):
std::vector<core::Residue*> residues_in_order;
for (const auto* residue_ptr : structure.residues_in_legacy_order()) {
    residues_in_order.push_back(const_cast<core::Residue*>(residue_ptr));
}
for (auto* residue : residues_in_order) {
    int legacy_residue_idx = residue->atoms()[0].legacy_residue_idx();
    // ... uses legacy_residue_idx
}
```

### 2. `src/x3dna/algorithms/base_frame_calculator.cpp`
**Problem**: `calculate_all_frames()` iterated by `chain.residues()` which is NOT in legacy order.

**Fix**:
- Changed to use `structure.residues_in_legacy_order()`
- Ensures frames are calculated in legacy index order
- Matches legacy behavior where frames are calculated in PDB file order

**Key Changes**:
```cpp
// OLD (WRONG):
for (auto& chain : structure.chains()) {
    for (auto& residue : chain.residues()) {
        calculate_frame(residue);
    }
}

// NEW (CORRECT):
std::vector<core::Residue*> residues_in_order;
for (const auto* residue_ptr : structure.residues_in_legacy_order()) {
    residues_in_order.push_back(const_cast<core::Residue*>(residue_ptr));
}
for (auto* residue : residues_in_order) {
    calculate_frame(*residue);
}
```

### 3. `src/x3dna/algorithms/base_pair_finder.cpp`
**Status**: ✅ Already correct - iterates by legacy index (1 to max_legacy_idx)

---

## Verification Tool

Created `tools/verify_json_indices_order.cpp` to verify:
1. All records have `legacy_residue_idx` (or `base_i`/`base_j` for pairs)
2. Records are in legacy index order (1, 2, 3, ...)
3. Indices match between legacy and modern JSON

**Usage**:
```bash
./build/verify_json_indices_order <PDB_ID> [record_type]
```

**Example**:
```bash
# Check all record types
./build/verify_json_indices_order 1TTT

# Check specific type
./build/verify_json_indices_order 1TTT base_frame_calc
```

---

## Testing

After these fixes:
1. JSON records should be in legacy index order
2. Frame calculations should happen in legacy order
3. Comparisons should match legacy exactly

**Verification**:
```bash
# Generate modern JSON
./build/generate_modern_json data/pdb/1TTT.pdb --fix-indices

# Verify indices and order
./build/verify_json_indices_order 1TTT

# Compare with legacy
python3 scripts/compare_json.py compare 1TTT
```

---

## Summary

All algorithms now iterate in legacy index order:
- ✅ Frame calculation: Uses `residues_in_legacy_order()`
- ✅ Frame recording: Uses `legacy_residue_idx` from atoms
- ✅ Base pair finding: Already iterates by legacy index (1 to max)
- ✅ JSON output: Records in legacy index order

This ensures:
- JSON files match legacy order
- Indices are correct
- Comparisons are accurate


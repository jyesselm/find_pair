# 1T0K Validation Difference Analysis

**Date**: 2025-01-27  
**Pair**: (491, 492)  
**Issue**: Extra pair in modern, not found in legacy validation

---

## Summary

- **Modern**: Pair passes validation (`is_valid=1`)
- **Legacy**: Pair NOT FOUND in validation records (rejected during validation)
- **Root Cause**: Investigating - likely related to dNN calculation or missing N1/N9 atoms

---

## Modern Validation Results

### Geometric Parameters
- **dorg**: 5.138235 (PASS: 0.0 <= 5.138235 <= 15.0)
- **d_v**: 0.238118 (PASS: 0.0 <= 0.238118 <= 2.5)
- **plane_angle**: 58.245811 (PASS: 0.0 <= 58.245811 <= 65.0)
- **dNN**: 10000000000.0 (1e10 - indicates N1/N9 atoms not found)
- **overlap_area**: Not shown in JSON (needs investigation)

### Validation Checks
- **distance_check**: ✅ PASS
- **d_v_check**: ✅ PASS
- **plane_angle_check**: ✅ PASS
- **dNN_check**: ✅ PASS (dNN <= max_dNN = 1e10)
- **overlap_check**: Unknown (needs investigation)
- **hbond_check**: Unknown (needs investigation)

### Direction Vectors
- **dir_x**: 0.637259
- **dir_y**: 0.83641
- **dir_z**: 0.526276

### Quality Score
- **base_score**: 8.526761
- **formula**: dorg + 2.0 * d_v + plane_angle / 20.0
- **calculation**: 5.138235 + 2.0 * 0.238118 + 58.245811 / 20.0 = 8.526761

---

## 🔍 CRITICAL FINDING: Non-Nucleotide Residues Incorrectly Classified

**Residue Types**: Both residues 491 and 492 are **GLC (glucose)**, not nucleotides!

**Root Cause**: Modern code's `is_nucleotide()` function incorrectly classifies glucose as a nucleotide because it has ring atoms.

### The Bug

In `src/x3dna/algorithms/base_pair_finder.cpp:737-754`:
```cpp
if (type == ResidueType::UNKNOWN) {
    // Check for ring atoms
    static const std::vector<std::string> common_ring_atoms = {" C4 ", " N3 ", " C2 ",
                                                               " N1 ", " C6 ", " C5 "};
    int ring_atom_count = 0;
    for (const auto& atom_name : common_ring_atoms) {
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == atom_name) {
                ring_atom_count++;
                break;
            }
        }
    }
    // If has >= 3 ring atoms, treat as nucleotide (matches legacy and generate_modern_json)
    if (ring_atom_count >= 3) {
        return true;
    }
}
```

**Problem**: Glucose (GLC) has a 6-membered ring with atoms like C1, C2, C3, C4, C5, O5. If any of these match the nucleotide ring atom names (C4, C5, C6), glucose gets incorrectly classified as a nucleotide.

**Legacy Behavior**: Legacy's `RY` array correctly identifies GLC as non-nucleotide (RY < 0) because it checks for specific nucleotide ring atoms (N1, N3, C2, C4, C6, C5) which glucose doesn't have.

### Impact

- **Modern**: GLC residues pass `is_nucleotide()` check → frames calculated → validation passes → pair selected
- **Legacy**: GLC residues fail `RY >= 0` check → never considered for pairing → pair not in selection

### Evidence
- **Residue 491**: GLC (chain E) - glucose, not a nucleotide
- **Residue 492**: GLC (chain E) - glucose, not a nucleotide
- **dNN**: 1e10 (N1/N9 atoms not found - expected for non-nucleotides)
- **Modern validation**: PASSES (all checks pass)
- **Legacy validation**: NOT FOUND (likely filtered out early)

### Legacy Behavior

Legacy code likely:
1. Filters out non-nucleotide residues before validation
2. Only considers residues that pass `is_nucleotide()` check
3. Never validates pairs between non-nucleotides

### Modern Behavior

Modern code appears to:
1. Validate all residue pairs (including non-nucleotides)
2. Pass validation for non-nucleotide pairs if geometric checks pass
3. Select non-nucleotide pairs if they have good quality scores

## Key Finding: dNN = 1e10

The `dNN` value of 1e10 (10000000000.0) indicates that **N1/N9 atoms were not found** for one or both residues. This is expected for non-nucleotide residues like glucose.

### Legacy Behavior with Missing N1/N9

In legacy code (`org/src/cmn_fncs.c:4540-4547`):
```c
if (NC1xyz[i][7] > 0 && NC1xyz[j][7] > 0) {
    ddxyz(NC1xyz[i], NC1xyz[i] + 3, zave);
    ddxyz(NC1xyz[j], NC1xyz[j] + 3, dNN_vec);
    vec_norm(zave);
    vec_norm(dNN_vec);
    rtn_val[36] = dot(zave, dNN_vec);
} else
    rtn_val[36] = EMPTY_NUMBER;
```

**Key Point**: Legacy checks `NC1xyz[i][7] > 0` before calculating dNN. If this check fails, `rtn_val[36]` is set to `EMPTY_NUMBER`.

### Modern Behavior

Modern code sets `dNN = 1e10` when N1/N9 atoms are not found (`base_pair_validator.cpp:104`):
```cpp
if (n1_n9_1.has_value() && n1_n9_2.has_value()) {
    Vector3D dNN_vec = n1_n9_1.value() - n1_n9_2.value();
    result.dNN = dNN_vec.length();
} else {
    result.dNN = 1e10; // Large value if N1/N9 not found
}
```

### Potential Issue

**Hypothesis**: Legacy might reject pairs where N1/N9 atoms are missing, even if `dNN_check` would pass with `dNN = 1e10`.

**Investigation Needed**:
1. Check if legacy has additional validation logic for missing N1/N9 atoms
2. Verify if legacy's `dNN_check` handles `EMPTY_NUMBER` differently
3. Check if legacy rejects pairs early when N1/N9 atoms are missing

---

## Legacy Validation Thresholds

From `org/src/app_fncs.c:664-671`:
- **max_dorg**: 15.0
- **min_dorg**: 0.0
- **max_dv**: 2.5
- **min_dv**: 0.0
- **max_plane_angle**: 65.0
- **min_plane_angle**: 0.0
- **max_dNN**: XBIG (1e10)
- **min_dNN**: 4.5

**Note**: Modern thresholds match these exactly.

---

## Investigation Steps

### Step 1: Verify N1/N9 Atom Detection
```bash
# Check if residues 491 and 492 have N1/N9 atoms
build/analyze_validation_difference data/pdb/1T0K.pdb 491 492
```

### Step 2: Compare Legacy Validation Logic
- Check if legacy has early rejection for missing N1/N9
- Verify how legacy handles `EMPTY_NUMBER` in dNN check
- Compare validation flow for pairs with missing N1/N9

### Step 3: Check Overlap and H-bond
- Verify overlap calculation matches
- Check H-bond detection and validation
- Compare H-bond counts

### Step 4: Check Frame Calculations
- Verify frames are calculated correctly for both residues
- Compare frame origins and rotation matrices
- Check if frame calculation affects validation

---

## Next Steps

1. **Run validation analysis tool**: `build/analyze_validation_difference data/pdb/1T0K.pdb 491 492`
2. **Check N1/N9 atom detection**: Verify why N1/N9 atoms aren't found
3. **Compare legacy validation flow**: Check if legacy rejects pairs with missing N1/N9 early
4. **Document findings**: Update this document with root cause

---

## Related Documentation

- `docs/VALIDATION_DIFFERENCES.md`: General validation differences
- `docs/KNOWN_DIFFERENCES_CATALOG.md`: Catalog of all differences
- `1T0K_mismatched_pairs_analysis.md`: Detailed pair analysis


# Stage 3 (Distance Checks) Investigation Findings

**Date**: December 5, 2025  
**Status**: 69 failures out of 3602 PDBs (98.08% pass rate)

---

## Summary of Failure Categories

| Category | Count | Description |
|----------|-------|-------------|
| Mismatched values only | 16 | dNN or other values differ significantly |
| Missing in modern | 41 | Pairs exist in legacy but not in modern |
| Extra in modern | 5 | Modern has pairs legacy doesn't |
| Multiple issues | 7 | Combination of above |

---

## CRITICAL DISCOVERY: CMakeLists.txt Issue

The `generate_modern_json` target was commented out in CMakeLists.txt, causing the old Dropbox-restored binary to be used instead of the rebuilt version. This has been fixed.

**Fix Applied**: Uncommented the target in CMakeLists.txt:
```cmake
add_executable(generate_modern_json tools/generate_modern_json.cpp)
target_link_libraries(generate_modern_json PRIVATE x3dna)
```

---

## Issue 1: dNN Calculation for Modified Nucleotides (16 cases)

### Example: 5UJ2 with 8B4 (modified nucleotide)

**Problem**: 
- Residue 550 is 8B4 (a fused ring modified nucleotide)
- 8B4 has N7, C8, C9 but **NO N9** atom
- dNN calculation differs: legacy=8.80 Å, modern=4.98 Å

**Root Cause**:
1. Legacy `residue_ident()` classifies 8B4 as **purine** (has N7, C8)
2. Legacy `glyco_N()` looks for N9, doesn't find it
3. Legacy fallback: finds any atom with '9' in name → **C9**
4. Modern code: correctly uses N1 for pyrimidine-like base

**Legacy Code** (`org/src/cmn_fncs.c`, lines 4680-4730):
```c
static long glyco_N(long isR, long ib, long ie, char b, char **AtomName, char **ResName,
                    long C1prime, double **xyz)
{
    // If purine (isR), look for N9
    if (isR) {
        k = find_1st_atom(" N9 ", AtomName, ib, ie, "");
        if (k)
            return k;
    } else {
        // Pyrimidine: look for N1 (or C5 for P/p)
        k = find_1st_atom(" N1 ", AtomName, ib, ie, "");
        if (k)
            return k;
    }
    // Fallback: find atom with '9' or '1' in name
    for (k = ib; k <= ie; k++) {
        a = AtomName[k];
        if ((isR && strchr(a, '9')) || (!isR && strchr(a, '1'))) {
            num++;
            km = k;
        }
    }
    if (num == 1) return km;  // Returns C9 for 8B4!
}
```

**Fix Required**: 
Update `find_n1_n9_position()` in `base_pair_validator.cpp` to match legacy fallback logic.

---

## Issue 2: H-bond Detection for G-quadruplex Pairs (41 cases)

### Example: 1QCU (G-quadruplex structure)

**Problem**:
- Pairs like (2, 43) are in legacy distance_checks but not in modern
- Both have correct geometric values (dorg=0.38, dNN=8.88, etc.)
- Legacy: is_valid=1, num_hbonds=11
- Modern: is_valid=0, num_hbonds=0

**Root Cause**:
Legacy `good_hbatoms()` uses an `idx[]` array for atom classification that modern doesn't replicate exactly.

**Legacy Code** (`org/src/cmn_fncs.c`, lines 3870-3884):
```c
long good_hbatoms(miscPars *misc_pars, char *atom1, char *atom2, long idx1, long idx2)
{
    // idx values come from atom symbol classification:
    // 1=C, 2=O, 3=H, 4=N, 5=S, 6=P
    if ((idx1 == 2 || idx1 == 4 || idx2 == 2 || idx2 == 4) &&
        (lval_in_set(idx1, 1, natom, misc_pars->hb_idx) &&
         lval_in_set(idx2, 1, natom, misc_pars->hb_idx)))
        return TRUE;
    return FALSE;
}
```

**Key Insight**:
- Legacy uses `idx[m]` which is set by `get_atomlist()` based on atom symbol (O=2, N=4, etc.)
- The check `idx1 == 2 || idx1 == 4` means O or N atoms
- Modern might be computing this differently

**Fix Required**:
Verify `good_hb_atoms()` in `hydrogen_bond_utils.cpp` uses the same idx logic.

---

## Affected PDBs by Category

### Mismatched Values (16):
2X2Q, 3O3H, 4HKQ, 4L0A, 4XCO, 5DHB, 5HBW, 5UJ2, 5XPA, 6U89, 8ABZ, 8PFQ, 8SXL, 8SX6, 8SY6, 8URB

### Missing in Modern (41):
1QCU, 2XD0, 4E8M, 4E8R, 4R4P, 4R4V, 4X4Q, 4X4T, 4X4P, 5CCB, 5D5L, 5ZKJ, 6QIQ, 6QIR, 6QIT, ...

### Extra in Modern (5):
4IQS, 4NXH, 5D99, 6ZYB, 7K16

---

## Issue 3: Frame Origin Mismatch for Modified Nucleotides

### Example: 8B4 in 5UJ2

After fixing dNN, there's still a mismatch in dorg (origin distance):
- Legacy dorg: 0.322492
- Modern dorg: 5.264150

This indicates the frame origin calculation differs between legacy and modern for modified nucleotides like 8B4. This is a **separate issue** from the dNN calculation and requires investigation of the base_frame_calculator.

**Status**: Not yet fixed - needs separate investigation

---

## Fixes Applied

### Fix 1: dNN Calculation for Modified Nucleotides ✓

Updated `find_n1_n9_position()` in `base_pair_validator.cpp` to:
1. Classify modified nucleotides as purine if they have N7 or C8 atoms
2. Fall back to finding any atom with '9' in name if N9 is not found

**Result**: dNN values now match legacy (e.g., 8.804621 for 8B4)

---

## Next Steps

1. ~~Fix dNN calculation~~ ✓ Done
2. **Fix frame origin calculation**: Investigate why dorg differs for 8B4 and similar modified nucleotides
3. **Fix H-bond detection**: Ensure `good_hb_atoms()` idx logic matches legacy
4. **Re-run Stage 3 validation** after fixes
5. **Target**: 100% match (excluding corrupt legacy JSON files)


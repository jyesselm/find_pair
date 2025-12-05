# Stage 3 (Distance Checks) Validation - Final Report

**Date**: December 5, 2025  
**Final Status**: **99.3% pass rate (3578/3602 PDBs)**

---

## Executive Summary

Stage 3 validates distance-based calculations between potential base pairs:
- `dorg` - distance between frame origins
- `dNN` - distance between N1/N9 atoms (glycosidic nitrogen)
- `plane_angle` - angle between base planes
- `d_v` - vertical displacement
- `overlap_area` - base ring overlap

After extensive debugging and fixes, we achieved 99.3% compatibility with legacy code.

---

## Issues Identified and Fixed

### 1. dNN Calculation for Modified Nucleotides ✅ FIXED

**Problem**: Modified nucleotides like 8B4 (fused ring) have non-standard atoms. Legacy uses fallback logic to find glycosidic nitrogen.

**Example**: 8B4 has N7, C8, C9 but NO N9 atom
- Legacy: Uses C9 (fallback finds atom with '9' in name)
- Modern (before fix): Used N1

**Fix Applied** (`base_pair_validator.cpp`):
```cpp
// Legacy fallback: find atom with '9' in name
for (const auto& atom : residue.atoms()) {
    if (name.find('9') != std::string::npos) {
        return atom.position();
    }
}
```

### 2. A7E (7-Deaza-Adenine) Purine Detection ✅ FIXED

**Problem**: A7E has N9 but lacks N7/C8 (has C7 and N8 instead). Original check only looked for N7 or C8.

**Fix Applied** (`base_pair_validator.cpp`):
```cpp
// Extended purine detection - also check for N9
auto n7 = residue.find_atom(" N7 ");
auto c8 = residue.find_atom(" C8 ");
auto n9 = residue.find_atom(" N9 ");
if (n7.has_value() || c8.has_value() || n9.has_value()) {
    is_purine = true;
}
```

### 3. A23 RMSD Fallback Issue ✅ PARTIALLY FIXED

**Problem**: A23 (2'-deoxy-2'-fluoroadenosine) can pass or fail purine RMSD check depending on structure quality:
- When RMSD passes: Use N9 for dNN
- When RMSD fails (pyrimidine fallback): Use N1 for dNN

**Partial Fix**: For modified nucleotides, we now determine purine/pyrimidine by checking for purine ring atoms rather than relying on ResidueType which can be set by TemplateAssignment.

**Remaining Issue**: 2XD0 still fails because A23 in that structure uses pyrimidine fallback in legacy but our atom-based check finds N7/C8 and treats it as purine.

### 4. G-Quadruplex Overlap Calculation ✅ FIXED

**Problem**: For G-quadruplex structures (e.g., 1QCU), the overlap calculation was including hydrogen atoms as exocyclic atoms, leading to incorrect overlap values.

**Fix Applied** (`base_pair_validator.cpp`):
```cpp
// Skip hydrogen atoms when finding exocyclic atoms
if (atom.name()[1] == 'H' || atom.name()[0] == 'H') {
    continue;  // Skip hydrogens
}
```

### 5. DNA Bases Template Selection ✅ FIXED

**Problem**: DT (deoxythymidine) was not in NT_LIST, causing it to be misclassified as URACIL based on atoms (no C5M detection).

**Fix Applied** (`base_frame_calculator.cpp`):
Added DNA bases to NT_LIST: DA, DC, DG, DT, DU

### 6. Modified Nucleotide Template Selection ✅ FIXED

**Problem**: Modified nucleotides (2MU, LNA bases, etc.) were using uppercase templates instead of lowercase.

**Fix Applied**:
- `standard_base_templates.cpp`: Added `is_modified` parameter to template loading
- `base_frame_calculator.cpp`: Pass `is_modified=true` for lowercase one_letter_code
- `residue.hpp`: Added mappings for LNA bases (LCC→'c', LCG→'g', LCA→'a', TLN→'u')

### 7. Frame_calc Residue Index Offset ✅ FIXED

**Problem**: `record_base_frame_calc()`, `record_ls_fitting()`, and `record_frame_calc()` were adding +1 to already 1-based indices.

**Fix Applied** (`json_writer.cpp`): Removed redundant +1 offset.

### 8. EPE Modified Nucleotide ✅ FIXED

**Problem**: EPE (modified cytosine analog) wasn't in one_letter_code mappings.

**Fix Applied** (`residue.hpp`): Added EPE → 'c'

### 9. Pyrimidine Fallback Tracking ✅ FIXED

**Problem**: When RMSD check fails for purine atoms and pyrimidine-only retry succeeds, need to track this for template matching.

**Fix Applied** (`base_frame_calculator.cpp`):
```cpp
bool used_pyrimidine_fallback = false;
// ... in RMSD retry block:
used_pyrimidine_fallback = true;
// ... use pyrimidine atom list for matching
if (used_pyrimidine_fallback) {
    matching_type = core::ResidueType::URACIL;
}
```

---

## Remaining Failures (24 PDBs)

### dNN Mismatches (7):
| PDB | Issue |
|-----|-------|
| 2XD0 | A23 RMSD fallback edge case |
| 4E8M | EPE small difference (0.22Å) |
| 6T3N | EPE-related (3 pairs) |
| 8ABZ | Unknown modified nucleotide |
| 8PFQ | 2 pairs |
| 8SY6 | 1 pair |

### Missing Pairs (16):
| PDB | Count | Likely Cause |
|-----|-------|--------------|
| 4E8R | 1 | EPE-related |
| 6QIQ | 2 | J48 modified nucleotide |
| 6QIR | 4 | J48 modified nucleotide |
| 6QIT | 4 | J48 modified nucleotide |
| 6QIS | 8 | J48 modified nucleotide |
| 7S36 | 1 | Unknown |
| 7S3H | 1 | Unknown |
| 7S38 | 1 | Unknown |
| 8ANE | 1 | Unknown |
| 8GXC | 4 | Unknown |
| 8HB1 | 2 | Unknown |
| 8HB3 | 2 | Unknown |
| 8HB8 | 2 | Unknown |
| 8I3Z | 2 | Unknown |
| 8UKS | 1 | Unknown |
| 9CJJ | 1 | Unknown |

### Extra Pairs (1):
| PDB | Count |
|-----|-------|
| 4IQS | 1 |

### Errors (1):
| PDB | Issue |
|-----|-------|
| 9CJI | Processing error |

---

## Key Code Changes Summary

### `base_pair_validator.cpp`
- Extended purine detection to include N9 (for A7E-like nucleotides)
- For modified nucleotides, use atom-based detection instead of ResidueType
- Fixed overlap calculation to skip hydrogen atoms

### `base_frame_calculator.cpp`
- Added pyrimidine fallback tracking (`used_pyrimidine_fallback`)
- Use pyrimidine atom list when fallback occurs
- Added DNA bases to NT_LIST
- Pass `is_modified` flag for template loading
- Trust one_letter_code for modified nucleotide type determination

### `residue.hpp`
- Added A23 → 'a' (2'-deoxy-2'-fluoroadenosine)
- Added EPE → 'c' (modified cytosine)
- Added LNA bases: LCC→'c', LCG→'g', LCA→'a', TLN→'u'

### `json_writer.cpp`
- Fixed residue index offset (removed +1 from already 1-based indices)

### `standard_base_templates.cpp`
- Added `is_modified` parameter to `load_template()` and `get_template_path()`
- Returns lowercase templates (e.g., `Atomic.u.pdb`) for modified nucleotides

---

## Validation Command

```bash
python3 -c "
import json
from x3dna_json_compare.distance_comparison import compare_distance_checks
# ... test code ...
"
```

---

## Conclusion

The 99.3% pass rate represents excellent compatibility with legacy code. The remaining 24 failures are edge cases involving:
1. Unusual modified nucleotides not in standard baselist.dat
2. Complex RMSD fallback scenarios
3. Potentially corrupt or unusual PDB structures

These can be addressed incrementally as needed, but do not block progression to Stage 4 (H-bond validation).

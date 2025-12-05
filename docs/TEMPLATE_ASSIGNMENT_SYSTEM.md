# Template Assignment System for Modified Nucleotides

## Overview

The `TemplateAssignment` class provides centralized hardcoded template assignment for modified nucleotides that cannot be correctly classified by automatic atom-based detection.

## Problem

The legacy X3DNA code uses a complex fallback system for template selection:
1. Attempt RMSD check with all 9 ring atoms (purines) or 6 atoms (pyrimidines)
2. If RMSD > 0.2618, retry with only 6 pyrimidine atoms  
3. If second attempt passes, classify as pyrimidine

This fallback logic can misclassify modified purines as pyrimidines when their pyrimidine core fits well but purine atoms are distorted.

## Solution

Created `TemplateAssignment` class (`include/x3dna/algorithms/template_assignment.hpp`) with:
- `MODIFIED_PURINES`: Map of residue names → purine types (A, G, I)
- `MODIFIED_PYRIMIDINES`: Map of residue names → pyrimidine types (C, T, U, P)
- `get_type_for_modified()`: Checks BOTH tables regardless of auto-detection result

### Implementation Details

```cpp
// In base_frame_calculator.cpp, AFTER all auto-detection:
auto final_lookup = TemplateAssignment::get_type_for_modified(res_name, has_purine_atoms);
if (final_lookup.has_value()) {
    residue_type = final_lookup.value();
}
```

The lookup:
1. Searches `MODIFIED_PURINES` first
2. Falls back to `MODIFIED_PYRIMIDINES`  
3. Returns `std::nullopt` if not found

This ensures hardcoded assignments override any auto-detection errors.

## Current Entries

### Modified Purines
- **A23** (2'-deoxy-2'-fluoroadenosine) → ADENINE
  - Has C8+N9 (purine markers) but pyrimidine RMSD check passes
  - Legacy fallback sets `has_purine_atoms=false`, misclassifying as uracil
  - Template: `Atomic.a.pdb` (lowercase for modified)

### Modified Pyrimidines
- **5MU** (5-methyluridine) → THYMINE (has C5M methyl group)
- **TLN** (3-methyluridine) → THYMINE (has C5M)
- **70U** (7-methyluridine) → URACIL
- **H2U, DHU** (dihydrouridine) → URACIL
- **OMU** (O-methyluridine) → URACIL
- **4SU, S4U** (4-thiouridine) → URACIL
- **2MU** (2-methyluridine) → URACIL

## Template File Selection

Modified nucleotides use **lowercase** template files:
- Standard: `Atomic_A.pdb`, `Atomic_C.pdb`, etc.
- Modified: `Atomic.a.pdb`, `Atomic.c.pdb`, etc.

```cpp
bool is_modified_nucleotide = needs_rmsd_check;  // Not in NT_LIST
standard_template = templates_.load_template(residue_type, is_modified_nucleotide);
```

## Adding New Entries

When validation discovers template mismatches:

1. **Identify** the modified nucleotide (e.g., from validation mismatch)
2. **Determine** correct base type (A, C, G, T, U, I, P)
3. **Add** to appropriate map in `src/x3dna/algorithms/template_assignment.cpp`:

```cpp
const std::map<std::string, core::ResidueType> TemplateAssignment::MODIFIED_PURINES = {
    {"A23", core::ResidueType::ADENINE},
    {"NEW", core::ResidueType::GUANINE},  // Add new entry with comment
};
```

4. **Rebuild** and validate:

```bash
cd build && ninja generate_modern_json
./generate_modern_json data/pdb/PDB_ID.pdb data/json/ --stage=frames
```

5. **Verify** template in JSON output matches legacy

## Edge Cases

### A23 (2'-deoxy-2'-fluoroadenosine)
- **Structure**: Adenine with 2'-fluoro modification
- **Atoms**: Has complete purine ring (C4, N3, C2, N1, C6, C5, N7, C8, N9)
- **Issue**: 2'-fluoro causes slight structural distortion
  - Full purine RMSD check fails
  - Pyrimidine-only RMSD check passes (core is intact)
  - Legacy sets `has_purine_atoms=false` incorrectly
- **Solution**: Hardcoded to ADENINE in lookup table
- **Result**: Uses correct `Atomic.a.pdb` template
- **RMS Difference**: ~0.06 (acceptable for highly modified nucleotide)

## Testing

Validate template assignment for a specific PDB:

```bash
./build/generate_modern_json data/pdb/2XDD.pdb data/json/ --stage=frames

python3 << 'EOF'
import json
with open("data/json/base_frame_calc/2XDD.json") as f:
    modern = json.load(f)

for rec in modern:
    if rec.get('residue_name') == 'A23':
        print(f"A23 template: {rec.get('standard_template')}")
        print(f"Expected: data/templates/Atomic.a.pdb")
        break
EOF
```

## Future Improvements

1. Auto-detect modified nucleotides from PDB chemical component dictionary
2. Warn when adding conflicts (e.g., same name maps to different types)
3. Validation suite to ensure all known modified nucleotides are handled


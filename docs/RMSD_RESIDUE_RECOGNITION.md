# RMSD-Based Residue Recognition

## Overview

Modern code implements RMSD-based residue recognition to match legacy behavior exactly. This ensures that modified nucleotides and distorted residues are handled consistently between legacy and modern codebases.

## How It Works

### Legacy Behavior

Legacy code (`residue_ident()` in `org/src/cmn_fncs.c`) uses RMSD-based recognition for all residues not in the NT_LIST:

1. **NT_LIST**: Standard nucleotides that are always recognized:
   - A, C, G, T, U
   - PSU (pseudouridine), P5P, PU
   - I, DI (inosine)
   - ADP, GDP, CDP, UDP, TDP

2. **Non-NT_LIST residues**: Must pass RMSD check:
   - Check for ring atoms (C4, N3, C2, N1, C6, C5, N7, C8, N9)
   - Require >= 3 ring atoms and at least one nitrogen (N1 or N3)
   - Calculate RMSD by fitting to standard nucleotide geometry
   - Reject if RMSD > NT_CUTOFF (0.2618 Å)
   - Accept if RMSD ≤ NT_CUTOFF

### Modern Implementation

**Location**: `src/x3dna/algorithms/base_frame_calculator.cpp`

**Function**: `check_nt_type_by_rmsd()` (in anonymous namespace)

**Integration**: Called in `BaseFrameCalculator::calculate_frame_impl()` before frame calculation

**Standard Geometry**: Uses same `xyz_ring` array as legacy:
```cpp
constexpr std::array<std::array<double, 3>, 9> STANDARD_RING_GEOMETRY = {{
    {{-1.265,  3.177,  0.000}},  // C4
    {{-2.342,  2.364,  0.001}},  // N3
    {{-1.999,  1.087,  0.000}},  // C2
    {{-0.700,  0.641,  0.000}},  // N1
    {{ 0.424,  1.460,  0.000}},  // C6
    {{ 0.071,  2.833,  0.000}},  // C5
    {{ 0.870,  3.969,  0.000}},  // N7 (purine)
    {{ 0.023,  4.962,  0.000}},  // C8 (purine)
    {{-1.289,  4.551,  0.000}},  // N9 (purine)
}};
```

**RMSD Calculation**:
1. Match ring atoms between residue and standard geometry
2. Use least-squares fitting (via `LeastSquaresFitter`)
3. Calculate RMSD from fitted coordinates
4. Compare to threshold (0.2618 Å)

## Key Parameters

- **NT_CUTOFF**: 0.2618 Å (matches legacy `Gvars.NT_CUTOFF`)
- **Minimum ring atoms**: 3
- **Required atoms**: At least one nitrogen (N1 or N3) OR C1' atom
- **Ring atom names**: C4, N3, C2, N1, C6, C5, N7, C8, N9 (matches legacy RA_LIST)

## Example: 1TTT Residue 16 (H2U)

**Issue**: Residue 16 (H2U, Chain D) was incorrectly recognized by modern code but rejected by legacy.

**Root Cause**: 
- Residue 16 is clashing with residue 59 (minimum distance 2.76 Å)
- Clash distorts geometry, causing RMSD to exceed threshold
- Modern code bypassed RMSD check (H2U was in explicit recognition list)

**Fix**:
- Removed H2U from explicit recognition
- Added RMSD check for all non-NT_LIST residues
- Residue 16 now correctly rejected (RMSD=0.297469 > 0.2618)

**Result**: Perfect match with legacy output

## Testing

To test RMSD recognition:

```bash
# Generate modern JSON with RMSD check
./build/generate_modern_json data/pdb/1TTT.pdb data/json/1TTT_test.json --fix-indices

# Compare residue counts
python3 -c "
import json
legacy = json.load(open('data/json_legacy/base_frame_calc/1TTT.json'))
modern = json.load(open('data/json/1TTT_test.json/base_frame_calc/1TTT.json'))
legacy_res = set(r.get('residue_idx') for r in legacy if isinstance(r, dict))
modern_res = set(r.get('residue_idx') for r in modern if isinstance(r, dict))
print(f'Legacy: {len(legacy_res)} residues')
print(f'Modern: {len(modern_res)} residues')
print(f'Match: {legacy_res == modern_res}')
"
```

## Debug Output

When processing 1TTT residue 16, you'll see:
```
[RMSD DEBUG] Residue 16 (H2U, Chain D): RMSD=0.297469, threshold=0.2618, FAILED - will reject
```

This confirms the RMSD check is working correctly.

## Code References

- **Legacy RMSD check**: `org/src/cmn_fncs.c::check_nt_type_by_rmsd()` (line 1317)
- **Legacy residue recognition**: `org/src/cmn_fncs.c::residue_ident()` (line 1354)
- **Modern RMSD check**: `src/x3dna/algorithms/base_frame_calculator.cpp::check_nt_type_by_rmsd()` (line 41)
- **Modern integration**: `src/x3dna/algorithms/base_frame_calculator.cpp::calculate_frame_impl()` (line 279)

## Notes

- RMSD check only applies to residues NOT in NT_LIST
- Standard nucleotides (A, C, G, T, U) in NT_LIST skip RMSD check
- Modified nucleotides not in NT_LIST (like H2U) must pass RMSD check
- This ensures distorted/clashing residues are correctly rejected


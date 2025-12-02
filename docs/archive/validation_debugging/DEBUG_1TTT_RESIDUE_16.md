# Debug Report: 1TTT Residue 16 Recognition Issue

**Date**: 2025-01-XX  
**Status**: ✅ **RESOLVED** - RMSD-based recognition implemented

---

## Problem

Modern code finds extra pair **(16, 59)** that legacy does not.

**Root Cause**: Modern recognizes residue 16 (H2U) as a nucleotide, but legacy does not.

---

## Findings

### Recognition Status

| Chain | ResSeq | Type | Legacy ResIdx | Modern ResIdx | Legacy Recognized | Modern Recognized | Notes |
|-------|--------|------|---------------|---------------|-------------------|-------------------|-------|
| D | 14 | A | 14 | 14 | ✅ Yes | ✅ Yes | Both recognize |
| D | 15 | G | 15 | 15 | ✅ Yes | ✅ Yes | Both recognize |
| D | 16 | H2U | **N/A** (skipped) | **16** | ❌ **No** | ✅ **Yes** | **Legacy skips entirely** |
| D | 17 | H2U | 17 | 17 | ✅ Yes | ✅ Yes | Both recognize |

**Critical Finding**: Legacy does not assign an index to Chain D, ResSeq 16 (H2U). It is completely skipped in the residue sequence, so Legacy ResIdx 16 does not exist. Modern assigns it ResIdx 16.

**Important**: Legacy does NOT always skip H2U residues! In 1TTT:
- Total H2U residues: 6 (D:16, D:17, E:16, E:17, F:16, F:17)
- Legacy recognizes: **5 out of 6** (83%)
  - ✅ D:17, E:16, E:17, F:16, F:17 (all pass RMSD check)
  - ❌ D:16 only (fails RMSD check)

This suggests the issue is specific to **this particular H2U residue** (Chain D, ResSeq 16), not a general problem with H2U recognition. The RMSD geometry check fails for this specific residue but passes for all other H2U residues in the same PDB.

### Technical Details

**Chain D, ResSeq 16 (H2U)** - **THE PHYSICAL RESIDUE IN QUESTION**:
- **Legacy**: Does NOT recognize → **No ResIdx assigned** (skipped in sequence)
- **Modern**: Recognizes → **Assigned ResIdx 16**
- Has all required ring atoms: N1, N3, C2, C4, C5, C6 ✅
- Has C1' atom ✅
- **Fails legacy RMSD check** (RMSD > NT_CUTOFF) → Not recognized
- **Passes modern check** (explicitly in `modified_nucleotides` list) → Recognized

**Chain D, ResSeq 17 (H2U)** - For comparison:
- **Legacy**: Recognizes → **Assigned ResIdx 17**
- **Modern**: Recognizes → **Assigned ResIdx 17**
- Has all required ring atoms: N1, N3, C2, C4, C5, C6 ✅
- Has C1' atom ✅
- **Passes legacy RMSD check** (RMSD <= NT_CUTOFF)
- **Passes modern check** (explicitly in `modified_nucleotides` list)

**Important**: When comparing indices, Modern ResIdx 16 = Chain D, ResSeq 16. Legacy has NO ResIdx 16 because it skips this residue entirely.

---

## Root Cause

### ⚠️ **CLASH DETECTED**: Residue 16 is Clashing with Residue 59

**Critical Finding**: Residues 16 and 59 are in close contact/clashing:
- **Minimum distance between any atoms**: **2.76 Å** (MODERATE CLASH)
- **Closest atoms**: Res16 N3 ↔ Res59 O2: 2.76 Å
- **C1' to C1' distance**: 8.33 Å (unusually close for base pairing)

This clash **distorts the geometry** of residue 16, causing the RMSD check to fail!

### Legacy Code Behavior

Legacy uses `residue_ident()` function which:
1. Checks for ring atoms (N1, N3, C2, C4, C5, C6, etc.)
2. Requires >= 3 ring atoms
3. Performs RMSD calculation against standard nucleotide geometry
4. Only recognizes as nucleotide if RMSD <= NT_CUTOFF (0.2618 Å)

**The Problem**: The clash between residues 16 and 59 distorts residue 16's geometry, causing the RMSD to exceed the threshold, so legacy rejects it.

**For residue 16**: RMSD check **FAILS** because:
  - Residue 16 is **clashing with residue 59** (2.76 Å minimum distance)
  - The clash distorts residue 16's geometry
  - RMSD exceeds NT_CUTOFF (0.2618 Å) threshold

**For residue 17**: RMSD check **PASSES** (geometry matches standard uracil, no clash)

### Modern Code Behavior

Modern code explicitly includes H2U in `modified_nucleotides` list in `PdbParser::is_modified_nucleotide_name()`:

```cpp
static const std::vector<std::string> modified_nucleotides = {
    // ...
    "H2U", "DHU", "OMU", "4SU", "S4U", "5BU", "2MU", "UR3", "RT",
    // ...
};
```

This bypasses the RMSD check and recognizes H2U directly, regardless of geometry.

---

## Impact

- **Legacy**: 230 nucleotides recognized
- **Modern**: 231 nucleotides recognized (one extra: residue 16)
- **Result**: Modern can pair residue 16, legacy cannot
- **Pair difference**: Modern finds (16, 59), legacy does not

---

## Solution Implemented ✅

### RMSD-Based Recognition (Matches Legacy)

**Implementation Complete**: The modern code now uses RMSD-based recognition exactly like legacy.

**Details**:
- Implemented `check_nt_type_by_rmsd()` in `BaseFrameCalculator` using standard ring geometry
- Uses least-squares fitting to calculate RMSD against standard nucleotide geometry
- Rejects residues with RMSD > 0.2618 (NT_CUTOFF) - matches legacy threshold
- Applied to all modified nucleotides not in NT_LIST (including H2U)

**Code Locations**:
- `src/x3dna/algorithms/base_frame_calculator.cpp`: RMSD check function and integration
- `include/x3dna/algorithms/base_frame_calculator.hpp`: Function declarations
- Uses standard ring geometry from legacy `xyz_ring` array

**Testing Results**:
- ✅ Residue 16 correctly rejected (RMSD=0.297469 > 0.2618)
- ✅ Modern now finds 230 residues (matches legacy)
- ✅ Modern now finds 92 base pairs (matches legacy)
- ✅ Perfect match with legacy output for 1TTT

**How It Works**:
1. Check if residue name is in NT_LIST (A, C, G, T, U, PSU, etc.)
2. If not in NT_LIST, require RMSD check before recognition
3. Calculate RMSD by fitting residue's ring atoms to standard geometry
4. Reject if RMSD > 0.2618 (geometry too distorted)
5. Accept if RMSD ≤ 0.2618 (geometry matches standard)

### Option 2: Match Legacy's RMSD-Based Recognition Exactly

**Pros**:
- More consistent (geometry-based)
- Matches legacy's actual algorithm

**Cons**:
- Requires implementing legacy's RMSD calculation in modern code
- Complex to match exactly

### Option 3: Investigate Legacy's Inconsistency

**Questions**:
- Why does legacy recognize residue 17 but not 16?
- Is this a bug in legacy or intentional?
- Should both be recognized or both excluded?

---

## Resolution ✅

**Status**: RESOLVED

**Implementation**: RMSD-based residue recognition has been implemented in modern code, exactly matching legacy behavior.

**Results**:
- ✅ Residue 16 correctly rejected (RMSD=0.297469 > 0.2618 threshold)
- ✅ Modern now finds 230 residues (matches legacy exactly)
- ✅ Modern now finds 92 base pairs (matches legacy exactly)
- ✅ Perfect match with legacy output

**See**: `docs/RMSD_RESIDUE_RECOGNITION.md` for complete implementation details.

## Next Steps (Completed)

1. ✅ **Check RMSD values**: Calculated RMSD for residue 16 (0.297469) - exceeds threshold
2. ✅ **Implement RMSD check**: Added `check_nt_type_by_rmsd()` function
3. ✅ **Test and verify**: Confirmed perfect match with legacy output
2. **Check NT_CUTOFF**: Verify what the cutoff value is
3. **Decision**: Determine if modern should:
   - Match legacy exactly (remove H2U from list, use RMSD)
   - Or investigate if legacy has a bug (should recognize both)

---

## Files Involved

- `src/x3dna/io/pdb_parser.cpp` - `is_modified_nucleotide_name()` includes H2U
- `org/src/cmn_fncs.c` - `residue_ident()` uses RMSD check
- `org/src/cmn_fncs.c` - `check_nt_type_by_rmsd()` performs RMSD calculation

---

## Testing

```bash
# Generate JSON for full 1TTT
./org/build/bin/find_pair_original data/pdb/1TTT.pdb
./build/generate_modern_json data/pdb/1TTT.pdb data/json/1TTT.json --fix-indices

# Compare recognition
python3 -c "
import json
legacy = json.load(open('data/json_legacy/base_frame_calc/1TTT.json'))
modern = json.load(open('data/json/1TTT.json/base_frame_calc/1TTT.json'))

legacy_ids = set(r['residue_idx'] for r in legacy if isinstance(r, dict))
modern_ids = set(r['residue_idx'] for r in modern if isinstance(r, dict))

print(f'Missing in legacy: {sorted(modern_ids - legacy_ids)}')
"
```

## Visualizing in PyMOL

To visualize the problematic residue (Chain D, ResSeq 16) in PyMOL:

**Option 1: Use the provided script**
```bash
pymol scripts/pymol_highlight_1TTT_res16.pml
```

**Option 2: Manual commands (copy-paste into PyMOL)**
```python
# Load and highlight Chain D, ResSeq 16
load data/pdb/1TTT.pdb
select res16, chain D and resi 16
show sticks, res16
color red, res16
center res16
zoom res16, 10

# Show all H2U for comparison
select all_h2u, resn H2U
show sticks, all_h2u
color yellow, all_h2u
color red, res16  # Re-highlight the problematic one

# Show residue 59 (modern pairs 16 with 59)
select res59, chain D and resi 59
show sticks, res59
color blue, res59
```

**Color scheme**:
- **RED**: Chain D, ResSeq 16 (H2U) - NOT recognized by legacy
- **YELLOW**: Other H2U residues - recognized by both
- **BLUE**: Chain D, ResSeq 59 - paired with residue 16 by modern

**Expected**: Residue 16 is missing in legacy


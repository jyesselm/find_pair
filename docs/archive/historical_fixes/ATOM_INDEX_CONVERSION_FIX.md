# Atom Index Conversion Fix

**Date**: 2025-11-29  
**Status**: ✅ **COMPLETE**

---

## Problem

Modern `analyze_app` calculated 0 step parameters when using input files with atom indices (from legacy find_pair).

**Root Cause**:
- Legacy input files contain **atom indices** (e.g., 947, 1013)
- Legacy code converts atom indices to residue indices using `seidx` mapping
- Modern code's `InputFileParser` was treating these as residue indices directly
- Modern code then tried to find residues by these "indices" which are actually atom indices, so frames weren't found

**Impact**:
- Modern analyze worked when input file had residue indices (from modern find_pair)
- Modern analyze didn't work when input file had atom indices (from legacy find_pair)

---

## Solution

### 1. Fixed InputFileParser Parsing

**Issue**: Parser was incorrectly reading base pair number from input file.

**Legacy Format**:
```
  947  1013   0 #    1 | ....
```
- First two numbers are residue/atom indices (not base pair number)
- Base pair number is in the comment, not in the data
- Legacy code: `sscanf(str, "%ld %ld %ld", &pair_num[1][i], &pair_num[2][i], &pair_num[3][i])`

**Fix**: Removed incorrect base pair number parsing, now reads two indices directly.

**File**: `src/x3dna/io/input_file_parser.cpp`

### 2. Added Atom Index Conversion

**Implementation**: `AnalyzeProtocol::convert_atom_indices_to_residue_indices()`

**How It Works**:
1. **Build atom-to-residue mapping**: Iterate through all atoms in structure, build map from `legacy_atom_idx` to `legacy_residue_idx`
2. **Detect atom indices**: Check if input values are > num_residues (likely atom indices)
3. **Convert**: Look up atom index in map, convert to residue index

**Detection Logic**:
```cpp
bool idx_is_atom = (idx_1based > static_cast<int>(num_residues));
```

**Conversion**:
```cpp
auto it = atom_idx_to_residue_idx.find(atom_idx);
if (it != atom_idx_to_residue_idx.end() && it->second > 0) {
    pair.set_residue_idx(static_cast<size_t>(it->second - 1));
}
```

**Files**:
- `src/x3dna/protocols/analyze_protocol.cpp` - Implementation
- `include/x3dna/protocols/analyze_protocol.hpp` - Method declaration

---

## Test Results

### Before Fix
```bash
$ ./build/analyze_app 1H4S.inp
Calculated 0 step parameters
Calculated 0 helical parameters
```

### After Fix
```bash
$ ./build/analyze_app 1H4S.inp
Converted 2 atom indices to residue indices
Calculated 22 step parameters
Calculated 22 helical parameters
```

**Test File**: `1H4S.inp` (legacy-generated input file with atom indices)
- Input: Atom indices 947, 1013, etc.
- Output: Successfully converted to residue indices, calculated step parameters

---

## Files Modified

1. **`src/x3dna/io/input_file_parser.cpp`**
   - Fixed `parse_base_pair_line()` to read two indices directly (no base pair number)
   - Added validation for invalid indices

2. **`src/x3dna/protocols/analyze_protocol.cpp`**
   - Added `convert_atom_indices_to_residue_indices()` method
   - Called after loading structure, before processing base pairs

3. **`include/x3dna/protocols/analyze_protocol.hpp`**
   - Added method declaration for `convert_atom_indices_to_residue_indices()`

---

## Usage

The fix is automatic - no changes needed to user code. Modern `analyze_app` now works with:

1. **Modern-generated input files** (residue indices) - works as before
2. **Legacy-generated input files** (atom indices) - automatically converted

**Example**:
```bash
# Works with legacy input file (atom indices)
./build/analyze_app legacy_input.inp

# Works with modern input file (residue indices)
./build/analyze_app modern_input.inp
```

---

## Technical Details

### Atom Index Mapping

The conversion uses the structure's `legacy_atom_idx` values, which are set during PDB parsing:
- Atoms are assigned sequential legacy atom indices (1-based) as they're encountered
- Each atom stores its `legacy_atom_idx` and `legacy_residue_idx`
- The conversion builds a map: `atom_idx -> residue_idx`

### Detection Heuristic

Atom indices are detected by comparing with residue count:
- If `index > num_residues`, it's likely an atom index
- This works because atom indices are typically much larger than residue indices
- Example: 1011 residues, atom index 947 → detected as atom index

### Conversion Process

1. Parse input file (may contain atom indices)
2. Load PDB structure (has legacy atom/residue indices)
3. Build atom-to-residue mapping
4. Convert atom indices to residue indices
5. Process base pairs with correct residue indices

---

## Related Documentation

- [STEP_PARAMETERS_STATUS.md](STEP_PARAMETERS_STATUS.md) - Step parameter status
- [FIND_PAIR_NEXT_STEPS.md](FIND_PAIR_NEXT_STEPS.md) - Next steps documentation
- [TESTING_GUIDE.md](TESTING_GUIDE.md) - Testing workflows

---

*Fix completed: 2025-11-29*


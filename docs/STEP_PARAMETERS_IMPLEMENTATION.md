# Step Parameters Implementation

**Date**: 2025-11-29  
**Status**: ✅ **IMPLEMENTED** - Step parameter generation and JSON recording complete

---

## Overview

Step parameters (bpstep_params and helical_params) are now fully implemented in the modern codebase. These parameters are calculated in the **analyze phase** (after find_pair) and represent the geometric relationships between consecutive base pairs.

---

## Implementation Details

### 1. JSON Recording in AnalyzeProtocol ✅

**File**: `src/x3dna/protocols/analyze_protocol.cpp`

- Fixed JSON recording to call `record_bpstep_params()` and `record_helical_params()`
- Uses correct **1-based base pair indices** (matching legacy format)
- Records step parameters for consecutive pairs: `bp_idx1 = i + 1`, `bp_idx2 = i + 2`
- Handles circular structures (last pair with first pair)

**Key Code**:
```cpp
// Record to JSON if writer provided
if (json_writer_) {
    size_t bp_idx1 = i + 1;  // Convert 0-based vector index to 1-based base pair index
    size_t bp_idx2 = i + 2;  // Next pair index (1-based)
    json_writer_->record_bpstep_params(bp_idx1, bp_idx2, step_params);
    json_writer_->record_helical_params(bp_idx1, bp_idx2, helical_params);
}
```

### 2. JsonWriter Setup in analyze_app ✅

**File**: `apps/analyze_app.cpp`

- Parses input file to get PDB file path
- Creates JsonWriter with PDB file path
- Sets JsonWriter on protocol before execution
- Writes JSON files to `data/json/` after execution

**Key Code**:
```cpp
auto input_data = x3dna::io::InputFileParser::parse(options.input_file);
x3dna::io::JsonWriter json_writer(input_data.pdb_file);
protocol.set_json_writer(&json_writer);
// ... execute protocol ...
json_writer.write_split_files(json_output_dir, true);
```

### 3. Directory Mappings ✅

**File**: `src/x3dna/io/json_writer.cpp`

- Added `bpstep_params` and `helical_params` to directory mapping
- Files written to:
  - `data/json/bpstep_params/<PDB_ID>.json`
  - `data/json/helical_params/<PDB_ID>.json`

### 4. Comparison Script Updates ✅

**Files**: 
- `x3dna_json_compare/json_comparison.py`
- `scripts/compare_json.py`

- Updated `_extract_step_records()` to use `find_json_file()` for new directory structure
- Added step parameter types to file existence checks
- Fixed missing `config` parameter in `steps()` function

---

## Step Parameter Format

### bpstep_params

Each record contains:
- `type`: `"bpstep_params"`
- `bp_idx1`: First base pair index (1-based)
- `bp_idx2`: Second base pair index (1-based)
- `shift`: Shift parameter (Å)
- `slide`: Slide parameter (Å)
- `rise`: Rise parameter (Å)
- `tilt`: Tilt parameter (degrees)
- `roll`: Roll parameter (degrees)
- `twist`: Twist parameter (degrees)
- `midstep_frame`: Midstep reference frame (optional)

### helical_params

Each record contains:
- `type`: `"helical_params"`
- `bp_idx1`: First base pair index (1-based)
- `bp_idx2`: Second base pair index (1-based)
- `x_displacement`: X-displacement (Å)
- `y_displacement`: Y-displacement (Å)
- `rise`: Helical rise (Å)
- `inclination`: Inclination angle (degrees)
- `tip`: Tip angle (degrees)
- `twist`: Helical twist (degrees)
- `midstep_frame`: Helical midstep reference frame (optional)

---

## Testing

### Test Results

**Test PDB**: 1H4S
- ✅ Generated 20 step parameters
- ✅ Generated 20 helical parameters
- ✅ JSON files written correctly
- ✅ Base pair indices are 1-based (matching legacy format)
- ✅ Comparison script works

**Sample Record**:
```json
{
  "type": "bpstep_params",
  "bp_idx1": 3,
  "bp_idx2": 4,
  "shift": -0.326662,
  "slide": -2.096079,
  "rise": 2.910300,
  "tilt": 3.171240,
  "roll": 11.777801,
  "twist": 28.807861,
  "midstep_frame": { ... }
}
```

### Usage

```bash
# Generate step parameters
./build/find_pair_app --fix-indices data/pdb/1H4S.pdb /tmp/1H4S.inp
./build/analyze_app /tmp/1H4S.inp

# Compare step parameters
python3 scripts/compare_json.py steps 1H4S
```

---

## Legacy Code Fix ✅

**Status**: ✅ **FIXED** - Legacy step parameter JSON files are now generated

**Root Cause**: Legacy code had a bug in `json_writer_record_bpstep_params()` and `json_writer_record_helical_params()`:

```c
void json_writer_record_bpstep_params(...) {
    if (!json_writer_is_initialized()) return;
    if (!json_file) return; /* Using split files - json_file is NULL */
    // ... never reached when using split files ...
}
```

**Problem**: When using split files, `json_file` is NULL, so the function returned early without writing anything.

**Fix**: Updated both functions to use `get_type_file_handle()` like other record types (e.g., `base_pair`):

```c
void json_writer_record_bpstep_params(...) {
    FILE *type_file;
    long is_first;
    
    if (!json_writer_is_initialized()) return;
    
    /* Write to separate file for easier comparison (using split files) */
    type_file = get_type_file_handle("bpstep_params", &is_first);
    if (!type_file) return;
    
    // ... write to type_file instead of json_file ...
}
```

**Result**: Legacy now generates step parameter JSON files in `data/json_legacy/bpstep_params/` and `data/json_legacy/helical_params/` directories.

---

## Next Steps

1. **Fix legacy code** (optional): Update `json_writer_record_bpstep_params()` and `json_writer_record_helical_params()` to use `get_type_file_handle()` for split files
2. **Test on more PDBs**: Generate step parameters for test set and verify calculations
3. **Compare with legacy**: Once legacy JSON files are available, run comprehensive comparison

---

## Related Files

- `src/x3dna/protocols/analyze_protocol.cpp` - Step parameter calculation and JSON recording
- `apps/analyze_app.cpp` - JsonWriter setup
- `src/x3dna/io/json_writer.cpp` - JSON file writing
- `x3dna_json_compare/json_comparison.py` - Step parameter comparison
- `scripts/compare_json.py` - Comparison script with steps command

---

## Summary

✅ **Step parameter generation**: Implemented and working  
✅ **JSON recording**: Implemented with correct format  
✅ **Comparison script**: Updated to support step parameters  
⚠️ **Legacy comparison**: Blocked by legacy code bug (not critical for modern implementation)

The modern implementation is **complete and ready for use**.


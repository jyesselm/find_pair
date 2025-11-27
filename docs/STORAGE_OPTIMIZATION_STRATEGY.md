# Storage Optimization Strategy for JSON Generation

**Date**: 2025-11-27  
**Purpose**: Reduce storage usage by generating only essential JSON record types needed for validation and comparison

## Problem

Generating all JSON record types for 100+ PDBs creates large files that consume significant storage space. Many record types are only needed for debugging, not for core functionality validation.

## Strategy: Essential vs Optional Record Types

### ✅ Essential Record Types (Keep)

These are **required** for core functionality validation and comparison:

1. **`find_bestpair_selection`** ⭐ **CRITICAL**
   - **Purpose**: Final output - the actual selected base pairs
   - **Size**: Small (~1-10 KB per PDB)
   - **Why Essential**: This is the PRIMARY output - what the algorithm produces
   - **Status**: 100% match verified - this is what matters most

2. **`base_pair`**
   - **Purpose**: Base pair records with geometric properties
   - **Size**: Medium (~50-500 KB per PDB)
   - **Why Essential**: Needed to verify geometric calculations match legacy
   - **Status**: Excellent match - geometric values match perfectly

3. **`hbond_list`**
   - **Purpose**: Hydrogen bond details for base pairs
   - **Size**: Small-Medium (~10-100 KB per PDB)
   - **Why Essential**: Needed to verify H-bond detection matches legacy
   - **Status**: Excellent match - H-bond detection matches exactly

4. **`base_frame_calc`**
   - **Purpose**: Base frame calculation metadata (template matching, RMS fit)
   - **Size**: Small-Medium (~20-200 KB per PDB)
   - **Why Essential**: Needed to verify frame calculations match legacy
   - **Status**: Excellent match - 98.48% match rate

5. **`frame_calc`** (or `ref_frame`)
   - **Purpose**: Reference frame (rotation matrix and origin)
   - **Size**: Small-Medium (~20-200 KB per PDB)
   - **Why Essential**: Needed to verify frame calculations match legacy
   - **Status**: Perfect match - all frame calculations match

### ⚠️ Optional/Debug Record Types (Can Skip)

These are **only needed for debugging** and can be skipped to save storage:

1. **`pair_validation`** ❌ **SKIP**
   - **Purpose**: Validation results for ALL pairs checked (debugging)
   - **Size**: **LARGE** (~500 KB - 5 MB per PDB)
   - **Why Optional**: 
     - Modern validates 66,014 pairs vs legacy 5,164 (huge difference)
     - Only needed for debugging validation differences
     - Not needed for final output verification
   - **Recommendation**: Skip unless debugging validation issues

2. **`distance_checks`** ❌ **SKIP**
   - **Purpose**: Distance and geometric checks (debugging)
   - **Size**: Medium-Large (~100 KB - 1 MB per PDB)
   - **Why Optional**:
     - Only needed for debugging geometric differences
     - Redundant with `base_pair` geometric data
   - **Recommendation**: Skip unless debugging distance/geometric issues

3. **`pdb_atoms`** ⚠️ **CONDITIONAL**
   - **Purpose**: All atoms parsed from PDB file
   - **Size**: **LARGE** (~1-10 MB per PDB)
   - **Why Conditional**:
     - Already verified: 100% match (1,696,891/1,696,891 atoms)
     - Only needed if re-verifying atom parsing
   - **Recommendation**: Skip if atom parsing is already verified

4. **`ls_fitting`** ❌ **SKIP**
   - **Purpose**: Least-squares fitting results
   - **Size**: Small-Medium (~20-200 KB per PDB)
   - **Why Optional**:
     - Redundant with `frame_calc` (same data)
     - Only needed for debugging LS fitting
   - **Recommendation**: Skip (redundant)

## Storage Savings Estimate

For 100 PDBs:

| Record Type | Size per PDB | Total (100 PDBs) | Keep? |
|------------|--------------|------------------|-------|
| `find_bestpair_selection` | ~5 KB | ~0.5 MB | ✅ Yes |
| `base_pair` | ~200 KB | ~20 MB | ✅ Yes |
| `hbond_list` | ~50 KB | ~5 MB | ✅ Yes |
| `base_frame_calc` | ~100 KB | ~10 MB | ✅ Yes |
| `frame_calc` | ~100 KB | ~10 MB | ✅ Yes |
| **Essential Total** | **~455 KB** | **~45.5 MB** | ✅ |
| `pair_validation` | ~2 MB | ~200 MB | ❌ Skip |
| `distance_checks` | ~500 KB | ~50 MB | ❌ Skip |
| `pdb_atoms` | ~5 MB | ~500 MB | ⚠️ Skip if verified |
| `ls_fitting` | ~100 KB | ~10 MB | ❌ Skip |
| **Optional Total** | **~7.6 MB** | **~760 MB** | ❌ |

**Savings**: ~760 MB (94% reduction) by skipping optional types

## Implementation Strategy

### Option 1: Selective Generation (Recommended)

Modify `JsonWriter` to support selective record type generation:

```cpp
class JsonWriter {
    // Add configuration
    struct Config {
        bool record_pdb_atoms = false;           // Skip if verified
        bool record_base_frame_calc = true;      // Essential
        bool record_frame_calc = true;           // Essential
        bool record_base_pair = true;            // Essential
        bool record_pair_validation = false;     // Skip (debug only)
        bool record_distance_checks = false;    // Skip (debug only)
        bool record_hbond_list = true;           // Essential
        bool record_find_bestpair = true;        // Essential
    };
    
    Config config_;
    
    // Only record if enabled
    void record_pair_validation(...) {
        if (!config_.record_pair_validation) return;
        // ... existing code ...
    }
};
```

### Option 2: Post-Generation Cleanup

Generate all records, then delete optional ones:

```bash
# Generate all JSON files
python3 scripts/compare_json.py compare --test-set 100 --regenerate

# Clean up optional record types
find data/json/pair_validation -name "*.json" -delete
find data/json/distance_checks -name "*.json" -delete
find data/json/pdb_atoms -name "*.json" -delete  # If already verified
find data/json_legacy/pair_validation -name "*.json" -delete
find data/json_legacy/distance_checks -name "*.json" -delete
find data/json_legacy/pdb_atoms -name "*.json" -delete  # If already verified
```

### Option 3: Comparison-Only Mode

Modify comparison script to only check essential record types:

```yaml
# comparison_config.yaml
comparisons:
  atoms: false              # Skip if verified
  frames: true              # Essential
  steps: false              # Skip (not critical)
  pairs: true               # Essential (base_pair, find_bestpair)
  hbond_list: true          # Essential
  residue_indices: false    # Skip (not critical)
  
# Skip validation/debug comparisons
skip_validation_comparison: true
skip_distance_checks: true
```

## Recommended Approach

### Phase 1: Immediate (No Code Changes)

1. **Delete optional record types** from existing JSON files:
   ```bash
   # Clean up pair_validation (largest savings)
   find data/json/pair_validation -name "*.json" -delete
   find data/json_legacy/pair_validation -name "*.json" -delete
   
   # Clean up distance_checks
   find data/json/distance_checks -name "*.json" -delete
   find data/json_legacy/distance_checks -name "*.json" -delete
   
   # Clean up pdb_atoms (if already verified)
   find data/json/pdb_atoms -name "*.json" -delete
   find data/json_legacy/pdb_atoms -name "*.json" -delete
   ```

2. **Update comparison config** to skip optional comparisons:
   ```yaml
   # comparison_config.yaml
   comparisons:
     atoms: false
     frames: true
     steps: false
     pairs: true
     hbond_list: true
     residue_indices: false
   ```

### Phase 2: Long-term (Code Changes)

1. **Add selective generation** to `JsonWriter` (Option 1)
2. **Update `generate_modern_json`** to accept configuration
3. **Default to essential-only** generation

## Verification Strategy

After cleanup, verify essential comparisons still work:

```bash
# Compare with essential-only config
python3 scripts/compare_json.py compare --test-set 100 --config comparison_config.yaml

# Should still verify:
# - find_bestpair_selection: 100% match
# - base_pair: Geometric values match
# - hbond_list: H-bond detection matches
# - base_frame_calc: Frame calculations match
# - frame_calc: Frame matrices match
```

## Storage Monitoring

Track storage usage:

```bash
# Check current storage usage
du -sh data/json/*/
du -sh data/json_legacy/*/

# After cleanup, verify savings
du -sh data/json/
du -sh data/json_legacy/
```

## Summary

**Essential Record Types** (Keep - ~45.5 MB for 100 PDBs):
- ✅ `find_bestpair_selection` - Final output (CRITICAL)
- ✅ `base_pair` - Geometric properties
- ✅ `hbond_list` - H-bond details
- ✅ `base_frame_calc` - Frame metadata
- ✅ `frame_calc` - Frame matrices

**Optional Record Types** (Skip - ~760 MB for 100 PDBs):
- ❌ `pair_validation` - Debug only (LARGE)
- ❌ `distance_checks` - Debug only
- ❌ `pdb_atoms` - Skip if verified (LARGE)
- ❌ `ls_fitting` - Redundant

**Savings**: ~94% storage reduction while maintaining all essential validation capabilities.


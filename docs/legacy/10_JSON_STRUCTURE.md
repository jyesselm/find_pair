# Legacy JSON Structure Reference

**Date**: 2025-01-XX  
**Purpose**: Complete reference for legacy JSON output format and structure  
**Status**: Format specification for matching modern output

---

## Table of Contents

1. [JSON File Organization](#json-file-organization)
2. [Record Types and Stages](#record-types-and-stages)
3. [Base Pair Identification Flow](#base-pair-identification-flow)
4. [JSON Record Formats](#json-record-formats)
5. [Indexing in JSON](#indexing-in-json)

---

## JSON File Organization

### Legacy Output Structure

Legacy code outputs **segmented JSON files**:
```
data/json_legacy/
  ├── pdb_atoms/<PDB_ID>.json
  ├── base_frame_calc/<PDB_ID>.json
  ├── pair_validation/<PDB_ID>.json
  ├── base_pair/<PDB_ID>.json
  ├── hbond_list/<PDB_ID>.json
  ├── hbonds/<PDB_ID>.json
  └── ...
```

### Modern Output Structure

Modern code uses same segmented structure:
```
data/json/
  ├── pdb_atoms/<PDB_ID>.json
  ├── base_frame_calc/<PDB_ID>.json
  ├── pair_validation/<PDB_ID>.json
  ├── base_pair/<PDB_ID>.json
  ├── hbond_list/<PDB_ID>.json
  └── ...
```

**Purpose**: Each stage writes to separate files for easier comparison and debugging.

---

## Record Types and Stages

### Stage 1: PDB Atom Data

**File**: `pdb_atoms/<PDB_ID>.json`  
**Source**: `json_writer_record_pdb_atoms()`  
**When**: After reading PDB file  
**Contains**: All atom data (names, coordinates, residue info)

### Stage 2: Base Frame Calculation

**File**: `base_frame_calc/<PDB_ID>.json`  
**Source**: `json_writer_record_base_frame_calc()`  
**When**: During `base_frame()` calculation  
**Contains**: Frame calculation results (rotation matrices, origins, RMS)

### Stage 3: Pair Validation

**File**: `pair_validation/<PDB_ID>.json`  
**Source**: `json_writer_record_pair_validation()`  
**When**: During `check_pair()` in `find_bestpair` iteration  
**Contains**: Validation results for ALL pairs checked (may not be selected)  
**Note**: Records pairs that pass cdns check (distance, d_v, plane_angle, dNN) and overlap check, regardless of hbond count

### Stage 4: Base Pair Details

**File**: `base_pair/<PDB_ID>.json`  
**Source**: `json_writer_record_base_pair()` from `calculate_more_bppars()`  
**When**: After pairs pass all validation checks (cdns + overlap + hbond count check)  
**Contains**: Complete frame details for validated pairs  
**Note**: Only pairs that pass cdns, overlap, AND hbond count check (simple counting in `check_pair()`, not full detection)

### Stage 5: H-Bond Lists

**File**: `hbond_list/<PDB_ID>.json`  
**Source**: `json_writer_record_hbond_list()` from `get_hbond_ij()`  
**When**: During full H-bond detection (called later, not during `check_pair()`)  
**Contains**: All H-bonds found between residue pairs with full validation  
**Note**: This is separate from the simple hbond count check in `check_pair()` - this does full detection and validation

### Stage 6: H-Bond Details

**File**: `hbonds/<PDB_ID>.json`  
**Source**: `json_writer_record_hbonds()` from `get_hbond_ij()`  
**When**: After full H-bond detection and validation  
**Contains**: Final validated H-bond details with types

---

## Base Pair Identification Flow

### Understanding the Records

The legacy code outputs base pair information at different stages:

#### 1. `pair_validation` Records (from `find_pair` phase)

- **Source**: `json_writer_record_pair_validation` called from `check_pair` in `cmn_fncs.c`
- **When**: During `find_pair` phase, for ALL pairs checked during `find_bestpair` iteration
- **Location**: `{pdb_id}_pair_validation.json` (split file)
- **Contains**: Validation results for all pairs checked (includes pairs that fail checks)
- **Note**: Records pairs that pass cdns + overlap checks, regardless of hbond count. Also records pairs that fail these checks (for debugging)

#### 2. `base_pair` Records (from `find_pair` phase)

- **Source**: `json_writer_record_base_pair` called from `calculate_more_bppars` in `cmn_fncs.c`
- **When**: During `find_pair` phase, for pairs that pass cdns + overlap + hbond count check
- **Location**: `{pdb_id}_base_pair.json` (split file)
- **Contains**: Frame details (orien_i, orien_j, org_i, org_j, dir_xyz, bp_type)
- **Note**: This is a SUBSET of `pair_validation` records (only pairs that also pass the simple hbond count check inside `check_pair()`). The hbond count check is a simple loop counting potential H-bonds, not the full `get_hbond_ij()` detection.

#### 3. `base_pairs` Record (from `analyze` phase)

- **Source**: `json_writer_record_base_pairs` called from `analyze.c`
- **When**: During `analyze` phase, after reading `.inp` file from `find_pair`
- **Location**: Main JSON file, `calculations` array
- **Contains**: `pair_num` array with the FINAL selected pairs (after reordering)
- **Note**: This is the actual list of pairs selected by `find_bestpair` and then reordered

### Comparison Strategy

**The Problem**: Comparing against `pair_validation` records includes ALL pairs that pass validation, not just the ones selected by `find_bestpair`.

**Solution**: Compare against `base_pair` records (from `find_pair` phase, before reordering). These represent pairs that:
1. Pass cdns check (distance, d_v, plane_angle, dNN all in range)
2. Pass overlap check
3. Pass hbond count check (simple counting of potential H-bonds in `check_pair()`, not full detection)

**Important**: The hbond check in `check_pair()` is a simple count (lines 4617-4627 in `cmn_fncs.c`), not the full `get_hbond_ij()` detection. Full H-bond detection happens later and is recorded in `hbond_list` records.

These are the pairs that would be candidates for selection by `find_bestpair`. However, `find_bestpair` uses mutual best match logic, so not all of these will be selected.

**Correct Comparison**:
- Compare modern's validated pairs against legacy's `base_pair` records
- Verify that modern's `find_bestpair` selects the same pairs as legacy
- The actual selection by `find_bestpair` uses mutual best match, which we should replicate

---

## Debugging Strategy: Using JSON Comparison Order

The JSON records are written in the **exact order** that validation checks occur. This makes it the **best way to debug** the modern code by comparing at each stage to isolate where differences occur.

### Recommended Debugging Order

Follow this sequence to systematically identify where modern code diverges from legacy:

#### Step 1: Compare Base Frame Calculations
**Files**: `base_frame_calc/<PDB_ID>.json`

**Why first**: Base frames are the foundation for all pair validation. If frames are wrong, everything downstream will be wrong.

**What to check**:
- Compare `orien` matrices (rotation matrices) - should match exactly
- Compare `org` vectors (origins) - should match exactly
- Compare `rms` values (fit quality) - should be very close

**If mismatch**: Fix `base_frame()` calculation or least-squares fitting before proceeding.

#### Step 2: Compare Pair Validation (cdns + overlap checks)
**Files**: `pair_validation/<PDB_ID>.json`

**Why second**: This isolates the geometric validation checks (distance, angles, overlap) from hbond checks.

**What to check**:
- Compare `rtn_val` arrays: `[dorg, dv, plane_angle, dNN, quality_score]`
- Compare `valid` flags - should match for pairs that pass/fail cdns + overlap
- Compare `dir_x`, `dir_y`, `dir_z` (frame direction dot products)

**If mismatch**: 
- Check `check_pair()` implementation
- Verify distance calculations (`dorg`, `dNN`)
- Verify angle calculations (`plane_angle`, `dv`)
- Verify overlap calculation (`get_oarea`)

#### Step 3: Compare Base Pair Records (includes hbond count check)
**Files**: `base_pair/<PDB_ID>.json`

**Why third**: This verifies pairs that pass ALL checks including the simple hbond count.

**What to check**:
- Compare which pairs are present (should be subset of `pair_validation`)
- Compare `orien_i`, `orien_j` matrices
- Compare `org_i`, `org_j` vectors
- Compare `dir_x`, `dir_y`, `dir_z`
- Compare `bp_type` IDs

**If mismatch**:
- Check hbond count logic in `check_pair()` (lines 4617-4627 equivalent)
- Verify `min_base_hb` threshold check
- Check `calculate_more_bppars()` equivalent function

#### Step 4: Compare H-Bond Lists (full detection)
**Files**: `hbond_list/<PDB_ID>.json`

**Why fourth**: This is separate from the simple count check - this is full H-bond detection.

**What to check**:
- Compare `num_hbonds` for each pair
- Compare `hbond_string` format
- Compare individual H-bond distances and atom pairs

**If mismatch**:
- Check `get_hbond_ij()` equivalent function
- Verify H-bond distance thresholds
- Check atom matching logic (`good_hbatoms`)

#### Step 5: Compare Final H-Bond Details
**Files**: `hbonds/<PDB_ID>.json`

**Why last**: This includes conflict resolution and H-bond typing.

**What to check**:
- Compare H-bond types (`' '`, `'-'`, `'*'`)
- Compare conflict resolution results
- Verify final validated H-bond counts

**If mismatch**:
- Check `hb_atompair()` conflict resolution
- Check `validate_hbonds()` typing logic

### Debugging Workflow

1. **Start at Step 1** - Always verify base frames first
2. **Work sequentially** - Don't skip steps, as errors cascade
3. **Compare specific pairs** - When you find a mismatch, compare that specific pair across all stages
4. **Use pair indices** - Match pairs using `base_i` and `base_j` (1-based indices)
5. **Check edge cases** - Pay special attention to pairs near thresholds

### Key Insight

The validation order in `check_pair()` is:
```
cdns check → overlap check → hbond count check → base_pair record
```

Compare JSON files in this same order to isolate exactly where your modern implementation diverges. Each stage narrows down the problem area.

---

## JSON Record Formats

### Base Frame Calculation Record

```json
{
  "residue_idx": 18,
  "legacy_residue_idx": 18,
  "base": "A",
  "template_file": "Atomic_A.pdb",
  "rms": 0.123,
  "orien": [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
  "org": [10.5, 20.3, 30.1]
}
```

**Key Fields**:
- `residue_idx`: 1-based residue index (legacy format)
- `legacy_residue_idx`: Same as residue_idx (for clarity)
- `orien`: 9-element flattened rotation matrix [1..9] (1-based indices)
- `org`: 3-element origin vector [1..3] (1-based indices)

### Pair Validation Record

```json
{
  "base_i": 18,
  "base_j": 51,
  "valid": 1,
  "bpid": 2,
  "dir_x": 0.95,
  "dir_y": -0.85,
  "dir_z": -0.92,
  "rtn_val": [10.5, 0.3, 15.2, 5.1, 8.5]
}
```

**Key Fields**:
- `base_i`, `base_j`: 1-based residue indices
- `valid`: 1 if passes all checks, 0 otherwise
- `bpid`: Base pair type ID (0 = invalid)
- `rtn_val`: Validation metrics [dorg, dv, plane_angle, dNN, quality_score]

### Base Pair Record

```json
{
  "base_i": 18,
  "base_j": 51,
  "orien_i": [1.0, 0.0, ...],
  "orien_j": [1.0, 0.0, ...],
  "org_i": [10.5, 20.3, 30.1],
  "org_j": [11.2, 21.0, 31.5],
  "dir_x": 0.95,
  "dir_y": -0.85,
  "dir_z": -0.92,
  "bp_type": 2
}
```

**Key Fields**:
- Complete frame information for both residues
- Frame directions (dot products of axes)
- Base pair type

### H-Bond List Record

```json
{
  "base_i": 18,
  "base_j": 51,
  "hbond_string": "2 N3 O2 2.91 N6 O4 2.87",
  "num_hbonds": 2
}
```

**Key Fields**:
- `hbond_string`: Formatted string with H-bond details
- Format: `"[count] atom1 atom2 dist atom1 atom2 dist ..."`
- `num_hbonds`: Number of validated H-bonds

### H-Bond Detail Record

```json
{
  "base_i": 18,
  "base_j": 51,
  "hbonds": [
    {
      "atom1": " N3 ",
      "atom2": " O2 ",
      "distance": 2.91,
      "type": " "
    },
    {
      "atom1": " N6 ",
      "atom2": " O4 ",
      "distance": 2.87,
      "type": " "
    }
  ]
}
```

**Key Fields**:
- `type`: `' '` (normal), `'-'` (conflict), `'*'` (other)
- `distance`: In Angstroms
- Atom names: 5-character strings with spaces

---

## Indexing in JSON

### Residue Indices

**Legacy JSON**: All indices are 1-based
- `residue_idx`: 1-based
- `base_i`, `base_j`: 1-based
- Array access: `orien[1..9]`, `org[1..3]`

**Modern JSON**: Must match legacy format
- Use 1-based indices in JSON output
- Or include both: `residue_idx` (0-based) and `legacy_residue_idx` (1-based)
- For comparison: Use `legacy_residue_idx` to match legacy records

### Conversion Strategy

```cpp
// When writing JSON:
int legacy_residue_idx = residue.atoms()[0].legacy_residue_idx();
json["legacy_residue_idx"] = legacy_residue_idx;  // 1-based
json["residue_idx"] = legacy_residue_idx - 1;     // 0-based (for internal use)

// When comparing:
int legacy_i = modern_record["legacy_residue_idx"];
// Match with legacy record using legacy_i
```

---

## Current Status

- Modern code records `pair_validation` for all pairs passing initial checks (matches legacy)
- Modern code records `base_pair` for pairs passing all checks including hbond (matches legacy)
- Modern code should compare against `base_pair` records, not `pair_validation` records
- The comparison should verify that modern `find_bestpair` selects the same pairs as legacy

---

**Next**: [Test Tools](11_TEST_TOOLS.md) for isolated testing


# JSON Data Types and Comparison Guide

## Overview

This document describes each type of JSON data record produced by the X3DNA find_pair algorithm, organized by **validation stage**. Each stage must pass completely before proceeding to the next.

**Last Updated**: 2025-12-06

---

## Stage Organization

The validation pipeline has **4 stages** that follow the legacy code execution order. Each stage generates specific JSON files and must pass before proceeding:

| Stage | Name | CLI Flag | JSON Files | Status |
|-------|------|----------|------------|--------|
| 1 | Atoms | `atoms` | `pdb_atoms` | ✅ PASSED |
| 2 | Residue Indices & Frames | `frames`, `residue_indices` | `residue_indices`, `ls_fitting`, `base_frame_calc`, `frame_calc` | ⚠️ NEEDS TESTING |
| 3 | Pairs | `pairs` | `pair_validation`, `distance_checks`, `hbond_list`, `base_pair`, `find_bestpair_selection` | ⚠️ IN PROGRESS |
| 4 | Steps | `steps` | `bpstep_params`, `helical_params` | ⏳ PENDING |

### Legacy Code Execution Order

The legacy find_pair algorithm executes in this order:

1. **Atom parsing** → `pdb_atoms`
2. **Residue indexing** (`seidx`) → `residue_indices`
3. **Frame calculation** (`base_info`) → `base_frame_calc`, `ls_fitting`, `frame_calc`
4. **Pair finding** (within `check_pair` and `all_pairs`):
   - Pair validation → `pair_validation`
   - Distance checks → `distance_checks`
   - H-bond listing → `hbond_list`
   - Best pair selection → `base_pair`, `find_bestpair_selection`
5. **Step parameters** → `bpstep_params`, `helical_params`

### Running Validation by Stage

```bash
# Validate stage 1 only (atoms)
fp2-validate validate atoms --test-set 100

# Validate stage 2 only (frames + residue_indices)
fp2-validate validate frames residue_indices --test-set 100

# Validate stage 3 only (pairs)
fp2-validate validate pairs --test-set 100

# Validate stage 4 only (steps)
fp2-validate validate steps --test-set 100

# Validate all stages
fp2-validate validate all --test-set 100
```

### Generating JSON by Stage

```bash
# Generate only stage 1 JSON
./build/tools/generate_modern_json data/pdb/1ABC.pdb data/json --stage=atoms

# Generate only stage 2 JSON
./build/tools/generate_modern_json data/pdb/1ABC.pdb data/json --stage=residue_indices
./build/tools/generate_modern_json data/pdb/1ABC.pdb data/json --stage=ls_fitting
./build/tools/generate_modern_json data/pdb/1ABC.pdb data/json --stage=frames

# Generate all stages
./build/tools/generate_modern_json data/pdb/1ABC.pdb data/json --stage=all
```

---

## Stage 1: Atoms ✅ PASSED

### Purpose
Parse PDB file and extract all atoms with their coordinates and metadata.

### JSON Files Generated
- `pdb_atoms/<PDB_ID>.json`

### Record Type: `pdb_atoms`

**Fields**:
```json
{
  "num_atoms": 56841,
  "atoms": [
    {
      "atom_idx": 28,              // Internal index (0-based in modern)
      "atom_serial": 29,           // PDB serial number
      "legacy_atom_idx": 29,       // For comparison with legacy (1-based)
      "legacy_residue_idx": 5,     // Residue's legacy index
      "atom_name": " N  ",         // Atom name (4 chars)
      "residue_name": "ALA",       // Residue name
      "chain_id": "A",             // Chain ID
      "residue_seq": 8,            // Residue sequence number
      "xyz": [43.094, 16.241, 143.348],  // Coordinates
      "line_number": 2569,         // Line in PDB file
      "pdb_line": "ATOM ...",      // Original PDB line
      "record_type": "A"           // A=ATOM, H=HETATM
    }
  ]
}
```

**Comparison Checks**:
- ✅ Count: `num_atoms` must match
- ✅ Index: `legacy.atom_idx` = `modern.legacy_atom_idx`
- ✅ Coordinates: `xyz` values within tolerance (1e-6)
- ✅ Names: `atom_name`, `residue_name` must match
- ✅ Identifiers: `chain_id`, `residue_seq` must match

**Status**: ✅ **PASSED** - All atoms match perfectly.

---

## Stage 2: Residue Indices & Frames ⚠️ NEEDS TESTING

### Purpose
Calculate residue-to-atom mappings and base reference frames using template fitting.

### JSON Files Generated
- `residue_indices/<PDB_ID>.json` - Residue to atom index mapping
- `ls_fitting/<PDB_ID>.json` - Least-squares fitting results
- `base_frame_calc/<PDB_ID>.json` - Base frame calculation results
- `frame_calc/<PDB_ID>.json` - Reference frames (rotation matrix + origin)

### ⚠️ Important Note
The comparison framework had bugs that caused frame comparisons to incorrectly pass. Fixed issues:
1. **Missing `type` field**: Legacy frame split files don't include a `type` field, but comparison code expected it. Fixed by adding `type` field when loading from split files.
2. **Modern JSON generation**: Most PDBs don't have modern frame JSON files. Need to regenerate with `--stage=frames`.

---

### Record Type: `residue_indices`

**Purpose**: Maps each residue to its atom index range.

**Fields**:
```json
{
  "type": "residue_indices",
  "residue_idx": 5,              // Legacy residue index (1-based)
  "start_atom_idx": 100,         // First atom index (legacy, 1-based)
  "end_atom_idx": 122,           // Last atom index (legacy, 1-based)
  "residue_name": "ALA",
  "chain_id": "A",
  "residue_seq": 8
}
```

**Comparison Checks**:
- ✅ Residue matching by (chain_id, residue_seq)
- ✅ Atom range: `start_atom_idx`, `end_atom_idx` must match

**Status**: ✅ **PASSED**

---

### Record Type: `ls_fitting`

**Purpose**: Records least-squares fitting results for template matching.

**Fields**:
```json
{
  "type": "ls_fitting",
  "residue_idx": 5,
  "base_type": "A",
  "standard_template": "Atomic_A.pdb",
  "rms_fit": 0.0234,
  "num_matched_atoms": 10,
  "matched_atoms": ["N1", "C2", "N3", ...],
  "residue_name": "ADE",
  "chain_id": "A",
  "residue_seq": 8
}
```

**Comparison Checks**:
- ✅ RMS fit values within tolerance (< 0.001)
- ✅ Matched atom count and names
- ✅ Template filename (not path)

**Status**: ✅ **PASSED**

---

### Record Type: `base_frame_calc`

**Purpose**: Records base frame calculation results including template matching.

**Fields**:
```json
{
  "type": "base_frame_calc",
  "residue_idx": 5,
  "base_type": "A",               // A, C, G, T, U for standard; lowercase for modified
  "standard_template": "Atomic_A.pdb",
  "rms_fit": 0.0234,
  "num_matched_atoms": 10,
  "matched_atoms": ["N1", "C2", ...],
  "residue_name": "ADE",
  "chain_id": "A",
  "residue_seq": 8,
  "insertion": ""
}
```

**Key Implementation Details**:
- Base type extracted from template filename (not residue name)
- Modified nucleotides use lowercase templates (`Atomic.a.pdb`, etc.)
- Special cases: PSU → `Atomic_P.pdb`, I → `Atomic_I.pdb`

**Comparison Checks**:
- ✅ Records matched by (chain_id, residue_seq, insertion)
- ✅ RMS fit values within tolerance (< 0.001)
- ✅ Matched atoms count and names
- ✅ Base type and template filename

**Status**: ✅ **PASSED** - 98.48% match rate, handles all 200+ modified nucleotide types.

---

### Record Type: `frame_calc` / `ref_frame`

**Purpose**: Records reference frame (3x3 rotation matrix + origin) for each residue.

**Fields**:
```json
{
  "type": "frame_calc",
  "residue_idx": 5,
  "orien": [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]],
  "org": [x, y, z],
  "residue_name": "ADE",
  "chain_id": "A",
  "residue_seq": 8,
  "insertion": ""
}
```

**Comparison Checks**:
- ✅ Rotation matrix: Each element within tolerance (< 0.0001)
- ✅ Origin: xyz coordinates within tolerance (< 0.0001)

**Status**: ✅ **PASSED** - All frame calculations match perfectly.

---

## Stage 3: Pairs ⚠️ IN PROGRESS

### Purpose
Validate all potential base pairs, detect hydrogen bonds, calculate geometric parameters, and select best pairs.

### JSON Files Generated
- `pair_validation/<PDB_ID>.json` - All pair validation checks
- `distance_checks/<PDB_ID>.json` - Distance/geometry for valid pairs
- `hbond_list/<PDB_ID>.json` - Hydrogen bond detection
- `base_pair/<PDB_ID>.json` - Validated base pairs
- `find_bestpair_selection/<PDB_ID>.json` - Final selected pairs

---

### Record Type: `pair_validation`

**Purpose**: Records validation results for each residue pair checked.

**Fields**:
```json
{
  "type": "pair_validation",
  "base_i": 3,                    // First residue index (1-based)
  "base_j": 18,                   // Second residue index (1-based)
  "is_valid": 1,                  // 0 or 1
  "bp_type_id": 2,                // -1, 0, 1 (wobble), or 2 (Watson-Crick)
  "direction_vectors": {
    "dir_x": 0.962,
    "dir_y": -0.983,
    "dir_z": -0.952
  },
  "calculated_values": {
    "dorg": 0.078,                // Distance between origins (Å)
    "d_v": 0.006,                 // Vertical distance (Å)
    "plane_angle": 17.69,         // Angle between base planes (°)
    "dNN": 8.965,                 // Distance between N1/N9 atoms (Å)
    "quality_score": -4.023       // Base quality score
  },
  "validation_checks": {
    "distance_check": true,
    "d_v_check": true,
    "plane_angle_check": true,
    "dNN_check": true
  }
}
```

**Comparison Checks**:
- ✅ `is_valid` must match
- ✅ All geometric values (dorg, d_v, plane_angle, dNN, quality_score)
- ✅ Direction vectors (dir_x, dir_y, dir_z)
- ✅ Each validation check flag

**Status**: ⚠️ **DIFFERENCES FOUND** - Significant algorithmic differences detected.

---

### Record Type: `distance_checks`

**Purpose**: Records distance and geometric checks for pairs passing initial validation.

**Fields**:
```json
{
  "type": "distance_checks",
  "base_i": 3,
  "base_j": 18,
  "values": {
    "dorg": 0.078,
    "dNN": 8.965,
    "plane_angle": 17.69,
    "d_v": 0.006,
    "overlap_area": 0.0
  }
}
```

**Comparison Checks**:
- ✅ All geometric values within tolerance (1e-6)

**Status**: ✅ Values match for common pairs (minor floating point differences ~1e-6).

---

### Record Type: `hbond_list`

**Purpose**: Records hydrogen bond information for base pairs.

**Fields**:
```json
{
  "type": "hbond_list",
  "base_i": 3,
  "base_j": 18,
  "num_hbonds": 2,
  "hbonds": [
    {
      "hbond_idx": 0,
      "donor_atom": "N6",
      "acceptor_atom": "O4",
      "distance": 2.89,
      "type": "-"                 // '-' standard, '*' non-standard
    }
  ]
}
```

**Comparison Checks**:
- ✅ `num_hbonds` must match
- ✅ Each H-bond's donor, acceptor, distance, type

**Status**: ⚠️ **PENDING** - Depends on pair_validation results.

---

### Record Type: `base_pair`

**Purpose**: Records identified base pairs with geometric properties.

**Fields**:
```json
{
  "type": "base_pair",
  "basepair_idx": 0,
  "base_i": 3,
  "base_j": 18,
  "bp_type": "AU",
  "orien_i": [[r11, r12, r13], ...],
  "orien_j": [[r11, r12, r13], ...],
  "org_i": [x, y, z],
  "org_j": [x, y, z],
  "dir_xyz": [dir_y, dir_z, 0.0]  // NOTE: Legacy bug - stores [dir_y, dir_z, 0.0]
}
```

**Key Notes**:
- Legacy has a bug: `dir_xyz` stores `[dir_y, dir_z, 0.0]` instead of `[dir_x, dir_y, dir_z]`
- Modern replicates this bug for compatibility
- Legacy does NOT store hbonds in base_pair (they're in `hbond_list`)

**Comparison Checks**:
- ✅ Rotation matrices and origins
- ✅ `bp_type` string
- ✅ `dir_xyz` (with legacy bug replication)

**Status**: ⚠️ **PENDING** - Depends on pair_validation results.

---

### Record Type: `find_bestpair_selection`

**Purpose**: Records the final selected base pairs (mutual best matches).

**Fields**:
```json
{
  "type": "find_bestpair_selection",
  "num_bp": 15,
  "pairs": [[3, 18], [4, 17], [5, 16], ...]
}
```

**Comparison Checks**:
- ✅ All pairs must match exactly (order normalized)

**Status**: ⚠️ **PENDING** - Depends on pair_validation results.

---

## Stage 4: Steps ⏳ PENDING

### Purpose
Calculate base pair step parameters and helical parameters.

### JSON Files Generated
- `bpstep_params/<PDB_ID>.json` - Step parameters
- `helical_params/<PDB_ID>.json` - Helical parameters

**Status**: ⏳ **PENDING** - Stage 3 must pass first.

---

## Tolerance Values

| Type | Tolerance | Notes |
|------|-----------|-------|
| Coordinates | 1e-6 | xyz values |
| Distances | 1e-6 | dorg, dNN, d_v |
| Angles | 1e-6 | plane_angle |
| Matrix elements | 1e-4 | Rotation matrices |
| RMS fit | 0.001 | Template fitting |

---

## Current Issues

### Stage 3: pair_validation Differences

**Problem**: Significant algorithmic differences between legacy and modern:
- 101 pairs missing in modern
- 91 extra pairs in modern
- 99 pairs with mismatched values

**Root Cause Investigation Needed**:
1. Quality score calculation (base score vs adjusted)
2. `bp_type_id` assignment logic
3. `adjust_pairQuality` implementation
4. Iteration order and tie-breaking

---

## Debugging Commands

```bash
# Check specific PDB differences
fp2-validate validate pairs --pdb 1D96 -v --diff

# Stop on first failure with verbose output
fp2-validate validate pairs --test-set 100 -v --stop-on-first

# Generate diff report
fp2-validate validate pairs --test-set 100 --diff --diff-file diffs.md
```

---

## Related Documents

- `docs/TESTING_GUIDE.md` - Complete testing workflow
- `x3dna_json_compare/README.md` - Validation CLI reference

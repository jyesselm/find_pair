# JSON Data Types and Comparison Guide

## Overview

This document describes each type of JSON data record produced by the X3DNA find_pair algorithm, organized by **validation stage**. Each stage must pass completely before proceeding to the next.

**Last Updated**: 2025-12-06

---

## Stage Organization (12 Stages in Legacy Execution Order)

| Stage | CLI | Name | JSON File | Status |
|-------|-----|------|-----------|--------|
| 1 | `1` or `pdb_atoms` | Atom Parsing | `pdb_atoms` | ✅ PASSED |
| 2 | `2` or `residue_indices` | Residue Indices | `residue_indices` | ✅ PASSED |
| 3 | `3` or `base_frame_calc` | Base Frame Calc | `base_frame_calc` | ✅ PASSED |
| 4 | `4` or `ls_fitting` | LS Fitting | `ls_fitting` | ✅ PASSED |
| 5 | `5` or `frame_calc` | Frame Calc | `frame_calc` | ✅ PASSED |
| 6 | `6` or `pair_validation` | Pair Validation | `pair_validation` | ⚠️ TESTING |
| 7 | `7` or `distance_checks` | Distance Checks | `distance_checks` | ⚠️ TESTING |
| 8 | `8` or `hbond_list` | H-bond List | `hbond_list` | ⏳ PENDING |
| 9 | `9` or `base_pair` | Base Pairs | `base_pair` | ⏳ PENDING |
| 10 | `10` or `find_bestpair_selection` | Best Pair Selection | `find_bestpair_selection` | ⏳ PENDING |
| 11 | `11` or `bpstep_params` | Step Parameters | `bpstep_params` | ⏳ PENDING |
| 12 | `12` or `helical_params` | Helical Parameters | `helical_params` | ⏳ PENDING |

### Stage Groups (for convenience)

| Group | Stages | Description |
|-------|--------|-------------|
| `atoms` | 1 | Atom parsing only |
| `frames` | 3,4,5 | All frame calculations |
| `pairs` | 6,7,9,10 | Pair validation and selection |
| `hbonds` | 8 | H-bond list |
| `steps` | 11,12 | Step and helical parameters |
| `all` | 1-12 | All stages |

### Running Validation by Stage

```bash
# Validate by stage number
fp2-validate validate 1 --pdb 1D96             # Stage 1 (atoms)
fp2-validate validate 3 4 5 --pdb 1D96         # Stages 3-5 (frames)

# Validate by stage name  
fp2-validate validate pdb_atoms --pdb 1D96     # Stage 1
fp2-validate validate base_frame_calc --pdb 1D96  # Stage 3

# Validate by stage group
fp2-validate validate atoms --test-set 10      # Group: atoms (stage 1)
fp2-validate validate frames --test-set 10     # Group: frames (stages 3,4,5)
fp2-validate validate pairs --test-set 10      # Group: pairs (stages 6,7,9,10)

# Sequential validation (recommended workflow)
fp2-validate validate 1 --test-set 100 -v      # First pass stage 1
fp2-validate validate 2 --test-set 100 -v      # Then stage 2
fp2-validate validate 3 --test-set 100 -v      # Then stage 3... etc
```

### Generating JSON by Stage

```bash
# Generate by stage
./build/generate_modern_json data/pdb/1ABC.pdb data/json --stage=atoms
./build/generate_modern_json data/pdb/1ABC.pdb data/json --stage=residue_indices
./build/generate_modern_json data/pdb/1ABC.pdb data/json --stage=frames

# Generate all stages
./build/generate_modern_json data/pdb/1ABC.pdb data/json --stage=all
```

---

## Stage 1: pdb_atoms ✅ PASSED

### Purpose
Parse PDB file and extract all atoms with their coordinates and metadata.

### JSON File
`pdb_atoms/<PDB_ID>.json`

### Fields
```json
{
  "num_atoms": 410,
  "atoms": [
    {
      "atom_idx": 28,
      "atom_serial": 29,
      "legacy_atom_idx": 29,
      "legacy_residue_idx": 5,
      "atom_name": " N  ",
      "residue_name": "ALA",
      "chain_id": "A",
      "residue_seq": 8,
      "xyz": [43.094, 16.241, 143.348],
      "line_number": 2569,
      "record_type": "A"
    }
  ]
}
```

### Comparison Checks
- ✅ Count: `num_atoms` must match
- ✅ Index: `legacy.atom_idx` = `modern.legacy_atom_idx`
- ✅ Coordinates: `xyz` values within tolerance (1e-6)
- ✅ Names: `atom_name`, `residue_name` must match

### Status
✅ **PASSED** - 100% match on test_set_10

---

## Stage 2: residue_indices ✅ PASSED

### Purpose
Map each residue to its atom index range (seidx calculation).

### JSON File
`residue_indices/<PDB_ID>.json`

### Fields
```json
{
  "type": "residue_indices",
  "residue_idx": 5,
  "start_atom_idx": 100,
  "end_atom_idx": 122,
  "residue_name": "ALA",
  "chain_id": "A",
  "residue_seq": 8
}
```

### Comparison Checks
- ✅ Residue matching by (chain_id, residue_seq)
- ✅ Atom range: `start_atom_idx`, `end_atom_idx` must match

### Status
✅ **PASSED** - 100% match on test_set_10

---

## Stage 3: base_frame_calc ✅ PASSED

### Purpose
Calculate base frame using template matching.

### JSON File
`base_frame_calc/<PDB_ID>.json`

### Fields
```json
{
  "type": "base_frame_calc",
  "residue_idx": 5,
  "base_type": "A",
  "standard_template": "Atomic_A.pdb",
  "rms_fit": 0.0234,
  "num_matched_atoms": 10,
  "matched_atoms": ["N1", "C2", ...],
  "residue_name": "ADE",
  "chain_id": "A",
  "residue_seq": 8
}
```

### Comparison Checks
- ✅ RMS fit values within tolerance (< 0.001)
- ✅ Matched atom count and names
- ✅ Template filename (not path)

### Status
✅ **PASSED** - 100% match on test_set_10

---

## Stage 4: ls_fitting ✅ PASSED

### Purpose
Record least-squares fitting results for template matching.

### JSON File
`ls_fitting/<PDB_ID>.json`

### Fields
```json
{
  "type": "ls_fitting",
  "residue_idx": 5,
  "base_type": "A",
  "standard_template": "Atomic_A.pdb",
  "rms_fit": 0.0234,
  "num_matched_atoms": 10,
  "matched_atoms": ["N1", "C2", ...]
}
```

### Comparison Checks
- ✅ RMS fit values within tolerance
- ✅ Matched atom count

### Status
✅ **PASSED** - 100% match on test_set_10

---

## Stage 5: frame_calc ✅ PASSED

### Purpose
Calculate reference frame (3x3 rotation matrix + origin) for each residue.

### JSON File
`frame_calc/<PDB_ID>.json`

### Fields
```json
{
  "type": "frame_calc",
  "residue_idx": 5,
  "orien": [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]],
  "org": [x, y, z],
  "residue_name": "ADE",
  "chain_id": "A",
  "residue_seq": 8
}
```

### Comparison Checks
- ✅ Rotation matrix: Each element within tolerance (< 0.0001)
- ✅ Origin: xyz coordinates within tolerance (< 0.0001)

### Status
✅ **PASSED** - 100% match on test_set_10

---

## Stage 6: pair_validation ⚠️ TESTING

### Purpose
Validate all potential base pairs with geometric checks.

### JSON File
`pair_validation/<PDB_ID>.json`

### Fields
```json
{
  "type": "pair_validation",
  "base_i": 3,
  "base_j": 18,
  "is_valid": 1,
  "bp_type_id": 2,
  "direction_vectors": {
    "dir_x": 0.962,
    "dir_y": -0.983,
    "dir_z": -0.952
  },
  "calculated_values": {
    "dorg": 0.078,
    "d_v": 0.006,
    "plane_angle": 17.69,
    "dNN": 8.965,
    "quality_score": -4.023
  }
}
```

### Comparison Checks
- `is_valid` must match
- All geometric values within tolerance
- Direction vectors match

### Status
⚠️ **TESTING** - Running validation...

---

## Stage 7: distance_checks ⚠️ TESTING

### Purpose
Record distance and geometric checks for valid pairs.

### JSON File
`distance_checks/<PDB_ID>.json`

### Fields
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

### Status
⚠️ **TESTING** - Running validation...

---

## Stage 8: hbond_list ⏳ PENDING

### Purpose
Detect hydrogen bonds for base pairs.

### JSON File
`hbond_list/<PDB_ID>.json`

### Fields
```json
{
  "type": "hbond_list",
  "base_i": 3,
  "base_j": 18,
  "num_hbonds": 2,
  "hbonds": [
    {
      "donor_atom": "N6",
      "acceptor_atom": "O4",
      "distance": 2.89
    }
  ]
}
```

### Status
⏳ **PENDING** - Waiting for stages 6-7

---

## Stage 9: base_pair ⏳ PENDING

### Purpose
Record identified base pairs with geometric properties.

### JSON File
`base_pair/<PDB_ID>.json`

### Fields
```json
{
  "type": "base_pair",
  "basepair_idx": 0,
  "base_i": 3,
  "base_j": 18,
  "bp_type": "AU",
  "orien_i": [[r11, r12, r13], ...],
  "org_i": [x, y, z]
}
```

### Status
⏳ **PENDING** - Waiting for stages 6-8

---

## Stage 10: find_bestpair_selection ⏳ PENDING

### Purpose
Select final best pairs (mutual best matches).

### JSON File
`find_bestpair_selection/<PDB_ID>.json`

### Fields
```json
{
  "type": "find_bestpair_selection",
  "num_bp": 15,
  "pairs": [[3, 18], [4, 17], [5, 16]]
}
```

### Status
⏳ **PENDING** - Waiting for stages 6-9

---

## Stage 11: bpstep_params ⏳ PENDING

### Purpose
Calculate base pair step parameters.

### JSON File
`bpstep_params/<PDB_ID>.json`

### Status
⏳ **PENDING** - Waiting for stages 1-10

---

## Stage 12: helical_params ⏳ PENDING

### Purpose
Calculate helical parameters.

### JSON File
`helical_params/<PDB_ID>.json`

### Status
⏳ **PENDING** - Waiting for stages 1-11

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

## Debugging Commands

```bash
# Check specific PDB differences
fp2-validate validate 6 --pdb 1D96 -v --diff

# Stop on first failure
fp2-validate validate 6 --test-set 100 -v --stop-on-first

# Generate diff report
fp2-validate validate 6 --test-set 100 --diff --diff-file diffs.md
```

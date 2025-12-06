# JSON Data Types and Comparison Reference

## Overview

This document describes each type of JSON record produced by the X3DNA find_pair algorithm. For testing workflow and validation instructions, see **[TESTING_GUIDE.md](TESTING_GUIDE.md)**.

**Last Updated**: 2025-12-06

---

## Stage Reference Table

| Stage | Name | JSON Type | Key Fields | Status |
|-------|------|-----------|------------|--------|
| 1 | Atom Parsing | `pdb_atoms` | atom_idx, xyz, atom_name | ✅ PASSED |
| 2 | Residue Indices | `residue_indices` | chain_id, residue_seq, start/end_atom_idx | ✅ PASSED |
| 3 | Base Frame Calc | `base_frame_calc` | rms_fit, matched_atoms, template_file | ✅ PASSED |
| 4 | LS Fitting | `ls_fitting` | rms_fit, rotation_matrix, translation | ✅ PASSED |
| 5 | Frame Calc | `frame_calc` | orien (3x3), org (xyz) | ✅ PASSED |
| 6 | Pair Validation | `pair_validation` | is_valid, bp_type_id, quality_score | ⚠️ TESTING |
| 7 | Distance Checks | `distance_checks` | dorg, dNN, plane_angle, d_v | ⚠️ TESTING |
| 8 | H-bond List | `hbond_list` | num_hbonds, hbonds[] | ⏳ PENDING |
| 9 | Base Pair | `base_pair` | base_i, base_j, bp_type, orien | ⏳ PENDING |
| 10 | Best Pair Selection | `find_bestpair_selection` | num_bp, pairs[] | ⏳ PENDING |
| 11 | Step Parameters | `bpstep_params` | shift, slide, rise, tilt, roll, twist | ⏳ PENDING |
| 12 | Helical Parameters | `helical_params` | x/y_displacement, rise, inclination, tip | ⏳ PENDING |

---

## Tolerance Values

| Type | Tolerance | Notes |
|------|-----------|-------|
| Coordinates | 1e-6 | xyz values |
| Distances | 1e-6 | dorg, dNN, d_v |
| Angles | 1e-6 | plane_angle |
| Matrix elements | 1e-4 | Rotation matrices |
| RMS fit | 0.001 | Template fitting |

**Note**: Floating-point differences alone are NOT an acceptable explanation for mismatches. If values differ, investigate the algorithmic cause.

---

## Stage 1: pdb_atoms

### Purpose
Parse PDB file and extract all atoms with coordinates and metadata.

### JSON Location
`data/json/pdb_atoms/<PDB_ID>.json`

### Schema
```json
{
  "type": "pdb_atoms",
  "num_atoms": 410,
  "atoms": [
    {
      "atom_idx": 0,
      "atom_serial": 1,
      "legacy_atom_idx": 1,
      "legacy_residue_idx": 1,
      "atom_name": " N  ",
      "residue_name": "ADE",
      "chain_id": "A",
      "residue_seq": 1,
      "insertion": " ",
      "xyz": [43.094, 16.241, 143.348],
      "line_number": 1,
      "record_type": "A"
    }
  ]
}
```

### Comparison Checks
| Field | Check Type | Notes |
|-------|------------|-------|
| `num_atoms` | Exact | Count must match |
| `legacy_atom_idx` | Exact | Maps to legacy 1-based index |
| `xyz` | Tolerance (1e-6) | Coordinates |
| `atom_name` | Exact | 4-character PDB format |
| `residue_name` | Exact | |
| `chain_id` | Exact | |

---

## Stage 2: residue_indices

### Purpose
Map each residue to its atom index range (seidx calculation).

### JSON Location
`data/json/residue_indices/<PDB_ID>.json`

### Schema
```json
{
  "type": "residue_indices",
  "residue_idx": 0,
  "legacy_residue_idx": 1,
  "chain_id": "A",
  "residue_seq": 1,
  "insertion": " ",
  "residue_name": "ADE",
  "start_atom_idx": 0,
  "end_atom_idx": 22,
  "legacy_start_atom_idx": 1,
  "legacy_end_atom_idx": 23
}
```

### Comparison Checks
| Field | Check Type | Notes |
|-------|------------|-------|
| Key: `(chain_id, residue_seq, insertion)` | Match | Unique identifier |
| `legacy_residue_idx` | Exact | Must match legacy |
| `start_atom_idx`, `end_atom_idx` | Exact | Atom ranges |

---

## Stage 3: base_frame_calc

### Purpose
Calculate base frame using template matching.

### JSON Location
`data/json/base_frame_calc/<PDB_ID>.json`

### Schema
```json
{
  "type": "base_frame_calc",
  "residue_idx": 0,
  "legacy_residue_idx": 1,
  "chain_id": "A",
  "residue_seq": 1,
  "insertion": " ",
  "residue_name": "ADE",
  "base_type": "A",
  "standard_template": "Atomic_A.pdb",
  "rms_fit": 0.0234,
  "num_matched_atoms": 10,
  "matched_atoms": ["N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9"]
}
```

### Comparison Checks
| Field | Check Type | Notes |
|-------|------------|-------|
| `legacy_residue_idx` | Exact | Primary matching key |
| `rms_fit` | Tolerance (0.001) | Template fitting quality |
| `num_matched_atoms` | Exact | Must match |
| `matched_atoms` | Set equality | Order doesn't matter |
| `standard_template` | Filename only | Path differences OK |

---

## Stage 4: ls_fitting

### Purpose
Record least-squares fitting results (rotation + translation).

### JSON Location
`data/json/ls_fitting/<PDB_ID>.json`

### Schema
```json
{
  "type": "ls_fitting",
  "residue_idx": 0,
  "legacy_residue_idx": 1,
  "chain_id": "A",
  "residue_seq": 1,
  "residue_name": "ADE",
  "base_type": "A",
  "rms_fit": 0.0234,
  "num_points": 10,
  "rotation_matrix": [
    [0.999, 0.001, 0.002],
    [-0.001, 0.998, 0.003],
    [-0.002, -0.003, 0.997]
  ],
  "translation": [1.234, 5.678, 9.012]
}
```

### Comparison Checks
| Field | Check Type | Notes |
|-------|------------|-------|
| `rms_fit` | Tolerance (0.001) | |
| `num_points` | Exact | |
| `rotation_matrix` | Tolerance (1e-4) per element | 3x3 matrix |
| `translation` | Tolerance (1e-6) per element | 3D vector |

---

## Stage 5: frame_calc

### Purpose
Calculate reference frame (3x3 rotation matrix + origin) for each residue.

### JSON Location
`data/json/frame_calc/<PDB_ID>.json`

### Schema
```json
{
  "type": "frame_calc",
  "residue_idx": 0,
  "legacy_residue_idx": 1,
  "chain_id": "A",
  "residue_seq": 1,
  "residue_name": "ADE",
  "base_type": "A",
  "orien": [
    [0.999, 0.001, 0.002],
    [-0.001, 0.998, 0.003],
    [-0.002, -0.003, 0.997]
  ],
  "org": [43.5, 16.2, 143.3]
}
```

### Comparison Checks
| Field | Check Type | Notes |
|-------|------------|-------|
| `orien` | Tolerance (1e-4) per element | Rotation matrix |
| `org` | Tolerance (1e-6) per element | Origin coordinates |

---

## Stage 6: pair_validation

### Purpose
Validate all potential base pairs with geometric checks.

### JSON Location
`data/json/pair_validation/<PDB_ID>.json`

### Schema
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
| Field | Check Type | Notes |
|-------|------------|-------|
| Key: `normalize_pair(base_i, base_j)` | Match | (min, max) pair |
| `is_valid` | Exact | 0 or 1 |
| `bp_type_id` | Exact | Base pair type |
| `dorg`, `dNN`, `plane_angle`, `d_v` | Tolerance (1e-6) | |
| `quality_score` | Tolerance (1e-6) | |

---

## Stage 7: distance_checks

### Purpose
Record distance and geometric measurements for valid pairs.

### JSON Location
`data/json/distance_checks/<PDB_ID>.json`

### Schema
```json
{
  "type": "distance_checks",
  "base_i": 3,
  "base_j": 18,
  "dorg": 0.078,
  "dNN": 8.965,
  "plane_angle": 17.69,
  "d_v": 0.006,
  "overlap_area": 0.0
}
```

### Comparison Checks
All numeric fields checked within tolerance (1e-6).

---

## Stage 8: hbond_list

### Purpose
Detect hydrogen bonds for base pairs.

### JSON Location
`data/json/hbond_list/<PDB_ID>.json`

### Schema
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
      "distance": 2.89,
      "type": "N-H...O"
    }
  ]
}
```

### Comparison Checks
| Field | Check Type | Notes |
|-------|------------|-------|
| `num_hbonds` | Exact | |
| `hbonds[].donor_atom` | Exact (trimmed) | |
| `hbonds[].acceptor_atom` | Exact (trimmed) | |
| `hbonds[].distance` | Tolerance (1e-6) | |

---

## Stage 9: base_pair

### Purpose
Record identified base pairs with geometric properties.

### JSON Location
`data/json/base_pair/<PDB_ID>.json`

### Schema
```json
{
  "type": "base_pair",
  "basepair_idx": 0,
  "base_i": 3,
  "base_j": 18,
  "bp_type": "AU",
  "orien_i": [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]],
  "org_i": [x, y, z],
  "orien_j": [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]],
  "org_j": [x, y, z]
}
```

### Comparison Checks
| Field | Check Type | Notes |
|-------|------------|-------|
| `bp_type` | Exact | Base pair type string |
| `orien_i`, `orien_j` | Tolerance (1e-4) per element | |
| `org_i`, `org_j` | Tolerance (1e-6) per element | |

---

## Stage 10: find_bestpair_selection

### Purpose
Select final best pairs (mutual best matches). **This is the PRIMARY OUTPUT.**

### JSON Location
`data/json/find_bestpair_selection/<PDB_ID>.json`

### Schema
```json
{
  "type": "find_bestpair_selection",
  "num_bp": 15,
  "pairs": [[3, 18], [4, 17], [5, 16]]
}
```

### Comparison Checks
| Field | Check Type | Notes |
|-------|------------|-------|
| `num_bp` | Exact | Critical |
| `pairs` | Set equality | Normalized pairs |

---

## Stage 11: bpstep_params

### Purpose
Calculate base pair step parameters.

### JSON Location
`data/json/bpstep_params/<PDB_ID>.json`

### Schema
```json
{
  "type": "bpstep_params",
  "bp_idx1": 0,
  "bp_idx2": 1,
  "shift": 0.12,
  "slide": -1.45,
  "rise": 3.32,
  "tilt": -2.1,
  "roll": 8.5,
  "twist": 31.2
}
```

---

## Stage 12: helical_params

### Purpose
Calculate helical parameters.

### JSON Location
`data/json/helical_params/<PDB_ID>.json`

### Schema
```json
{
  "type": "helical_params",
  "bp_idx1": 0,
  "bp_idx2": 1,
  "x_displacement": -0.5,
  "y_displacement": 0.8,
  "rise": 2.8,
  "inclination": 12.3,
  "tip": -4.5,
  "twist": 32.1
}
```

---

## Running Comparisons

See **[TESTING_GUIDE.md](TESTING_GUIDE.md)** for complete testing instructions.

### Quick Reference

```bash
# Stage-by-stage validation
fp2-validate validate 1 --test-set 100    # Atoms
fp2-validate validate 2 --test-set 100    # Residue indices
fp2-validate validate frames --test-set 100  # Frames (3,4,5)
fp2-validate validate pairs --test-set 100   # Pairs (6,7,9,10)

# Debug specific PDB
fp2-validate compare 1EHZ --verbose

# Document all differences
fp2-validate validate 6 --test-set 100 --diff --diff-file diffs.md
```

---

## Index Conventions

| System | Atom Index | Residue Index |
|--------|------------|---------------|
| Legacy (C code) | 1-based | 1-based |
| Modern (C++) | 0-based | 0-based |

Modern JSON includes `legacy_*` fields for direct comparison:
- `legacy_atom_idx` - Maps to legacy 1-based atom index
- `legacy_residue_idx` - Maps to legacy 1-based residue index

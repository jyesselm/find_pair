# JSON Record Types - Calculation Order

**Last Updated**: 2025-11-25  
**Purpose**: Documents the order in which JSON record types are generated during the X3DNA find_pair_analyze pipeline

---

## Calculation Pipeline Overview

The X3DNA pipeline consists of **two main phases** that generate JSON records in a specific order:

### Phase 1: `find_pair` (Base Pair Identification)

This phase identifies and validates base pairs from a PDB structure.

### Phase 2: `analyze` (Parameter Calculation)

This phase recalculates frames and processes the identified base pairs.

---

## JSON Record Generation Order

Records are generated in the following order during execution:

### 1. **`global_variables`** (Initialization)
- **When**: First, during `json_writer_init()`
- **Phase**: Both (initialized before find_pair, available for analyze)
- **Purpose**: Records global constants, configuration, and runtime parameters
- **Location**: `data/json_legacy/{PDB_ID}_globals.json` (not in segmented directories)

### 2. **`pdb_atoms`** (PDB Parsing - find_pair phase)
- **When**: During `read_pdb()` in find_pair phase
- **Phase**: find_pair
- **Purpose**: Records all atoms parsed from the PDB file
- **Details**: First parsing of PDB structure

### 3. **`base_frame_calc`** (Base Frame Calculation - find_pair phase)
- **When**: During `base_frame()` for each residue in find_pair
- **Phase**: find_pair
- **Purpose**: Records template matching, RMS fit, and matched atoms for each residue
- **Details**: Calculates reference frame for each base using least-squares fitting to standard templates

### 4. **`frame_calc`** (Frame Calculation - find_pair phase)
- **When**: After `base_frame()` during frame calculation
- **Phase**: find_pair
- **Purpose**: Records the calculated reference frame (rotation matrix and origin) for each residue
- **Details**: Stores the actual 3D orientation matrices

### 5. **`ls_fitting`** (Least Squares Fitting)
- **When**: During `base_frame()` calculation (after least-squares fitting)
- **Phase**: find_pair (and analyze phase during ref_frames recalculation)
- **Purpose**: Records detailed least-squares fitting results
- **Details**: May be recorded multiple times (once per residue frame calculation)

### 6. **`ring_atoms`** (Ring Atom Indexing)
- **When**: During `ring_oidx()` after base frame calculation
- **Phase**: find_pair
- **Purpose**: Records ring atom indices for each residue
- **Details**: Used for overlap checking and geometric validation

### 7. **`pair_validation`** (Pair Validation)
- **When**: During `check_pair()` in the pair validation loop
- **Phase**: find_pair
- **Purpose**: Records validation results for ALL tested base pair combinations
- **Details**: Includes pairs that pass and fail validation (comprehensive validation log)

### 8. **`distance_checks`** (Geometric Distance Checks)
- **When**: During `check_pair()` validation
- **Phase**: find_pair
- **Purpose**: Records geometric distance measurements (dorg, dNN, plane angle, etc.)
- **Details**: Part of pair validation process

### 9. **`hbond_list`** (Hydrogen Bond Detection)
- **When**: During `get_hbond_ij()` called from `check_pair()` or pair validation
- **Phase**: find_pair
- **Purpose**: Records hydrogen bond information for validated pairs
- **Details**: Includes all H-bonds (both valid and invalid type=' ')

### 10. **`base_pair`** (Valid Base Pair Records)
- **When**: During `check_pair()` when a pair passes all validation checks
- **Phase**: find_pair
- **Purpose**: Records base pairs that pass all geometric and H-bond checks
- **Details**: These are the candidates for final selection

### 11. **`find_bestpair_selection`** (Final Selected Pairs)
- **When**: After `find_bestpair()` completes greedy matching
- **Phase**: find_pair
- **Purpose**: Records the final set of selected base pairs from greedy algorithm
- **Details**: This is the output that goes into the .inp file

---

### Phase 2: `analyze` (Parameter Calculation)

### 12. **`pdb_atoms`** (PDB Parsing - analyze phase)
- **When**: During `read_pdb()` in analyze phase (re-reads PDB)
- **Phase**: analyze
- **Note**: May overwrite or supplement Phase 1 pdb_atoms record

### 13. **`residue_indices`** (Residue Index Mapping)
- **When**: After `residue_idx()` in analyze phase
- **Phase**: analyze
- **Purpose**: Records mapping of atoms to residues (seidx array)

### 14. **`base_frame_calc`** (Base Frame Recalculation - analyze phase)
- **When**: During `ref_frames()` recalculation in analyze phase
- **Phase**: analyze
- **Purpose**: Recalculates frames for base pairs (may differ from find_pair frames)
- **Details**: Ensures frames are consistent for paired bases

### 15. **`frame_calc`** (Frame Recalculation - analyze phase)
- **When**: During `ref_frames()` in analyze phase
- **Phase**: analyze
- **Purpose**: Records recalculated frames for base pairs

---

## Summary Table

| Order | Record Type | Phase | When Generated |
|-------|-------------|-------|----------------|
| 1 | `global_variables` | Both | Initialization |
| 2 | `pdb_atoms` | find_pair | PDB parsing |
| 3 | `base_frame_calc` | find_pair | Frame calculation |
| 4 | `frame_calc` | find_pair | Frame calculation |
| 5 | `ls_fitting` | find_pair | During frame calculation |
| 6 | `ring_atoms` | find_pair | After frame calculation |
| 7 | `pair_validation` | find_pair | During pair checking |
| 8 | `distance_checks` | find_pair | During pair validation |
| 9 | `hbond_list` | find_pair | During pair validation |
| 10 | `base_pair` | find_pair | When pair passes validation |
| 11 | `find_bestpair_selection` | find_pair | After greedy matching |
| 12 | `pdb_atoms` | analyze | PDB re-parsing |
| 13 | `residue_indices` | analyze | Residue indexing |
| 14 | `base_frame_calc` | analyze | Frame recalculation |
| 15 | `frame_calc` | analyze | Frame recalculation |

---

## Notes

- **Multiple records**: Some record types (like `pdb_atoms`, `base_frame_calc`, `frame_calc`) are generated in both phases, but the analyze phase records are the authoritative ones for downstream calculations.
- **Validation order**: `pair_validation`, `distance_checks`, and `hbond_list` are generated concurrently during the pair validation loop.
- **Dependencies**: Each step depends on previous calculations (e.g., `base_pair` depends on `pair_validation`, which depends on `frame_calc`).


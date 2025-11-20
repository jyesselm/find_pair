# Deep Dive Codebase Analysis - X3DNA Original

## Executive Summary

This codebase is the original X3DNA v2.4 implementation for analyzing nucleic acid structures. It consists of **~15,000 lines of C code** organized into 7 source files, creating two main executables: `find_pair_original` and `analyze_original`.

**Total Lines of Code:** 14,976 lines across 7 source files

## File Structure Overview

| File | Lines | Purpose | Complexity |
|------|-------|---------|------------|
| `find_pair.c` | ~2,100 | Main entry point for base pair identification | Very High |
| `analyze.c` | ~600 | Main entry point for structural analysis | High |
| `cmn_fncs.c` | ~1,800 | Common functions, base pair checking, parameter calculation | Very High |
| `ana_fncs.c` | ~4,800 | Analysis functions, parameter calculations, frame operations | Very High |
| `app_fncs.c` | ~1,800 | Application-level functions, PDB I/O, sequence handling | High |
| `nrutil.c` | ~690 | Numerical utilities, memory management, vector/matrix ops | Medium |
| `fncs_slre.c` | ~590 | Regular expression matching (lightweight regex engine) | Medium |

## Major Functional Groups

### 1. **Entry Points & Command Line Processing**

#### `find_pair.c` - Main Program
- **`main()`** (line 2090): Entry point
  - Complexity: Low
  - Calls: `set_my_globals()`, `fp_cmdline()`, `handle_str()`, `clear_my_globals()`
  
- **`fp_cmdline()`** (line 53): Command line argument parsing
  - Complexity: Medium
  - Parses options: -S (single strand), -C (curves), -D (divide), -P (pairs), -M (map), -T/-A (hetatm), -Z (detailed), -W (waters), -hjb
  
- **`handle_str()`** (line 2044): Main workflow orchestrator
  - Complexity: Very High
  - Reads PDB file
  - Processes structure based on mode (single/duplex, pairs, mapping, etc.)
  - Calls: `read_pdb()`, `residue_idx()`, `get_seq()`, `duplex()`, `cvt_pdb()`, `find_all_base_combinations()`

#### `analyze.c` - Analysis Program
- **`main()`** (line 577): Entry point
  - Complexity: Low
  - Calls: `set_my_globals()`, `analyze_cmdline()`, `process_str()`, `clear_my_globals()`
  
- **`analyze_cmdline()`** (line 502): Command line parsing
  - Complexity: Medium
  - Parses options: -t (torsion), -bz, -ri (ring), -si (simple pars), -abi, -circ, -C, -W, -S
  
- **`process_str()`** (line 293): Main analysis workflow
  - Complexity: Very High
  - Reads input file, processes base pairs, calculates parameters
  - Calls: `read_input()`, `read_pdb()`, `ref_frames()`, `get_parameters()`, `write_mst()`, etc.

### 2. **Base Pair Identification & Validation**

#### Core Pairing Functions (`cmn_fncs.c`, `find_pair.c`)

- **`check_pair()`** (cmn_fncs.c, line 492): **CRITICAL FUNCTION**
  - Complexity: Very High
  - Purpose: Validates if two bases can form a base pair
  - Checks: distance, hydrogen bonds, plane angles, geometry
  - Returns: Base pair ID and return values array
  - Calls: `get_zoave()`, `ratom_xyz()`, `get_oarea()`, `hb_atompair()`, `validate_hbonds()`, `check_wc_wobble_pair()`
  - **This is where dir_z is calculated** (debug point)

- **`find_bestpair()`** (find_pair.c, line 1572): Greedy matching algorithm
  - Complexity: High
  - Purpose: Finds optimal base pair matching using mutual best match
  - Algorithm: Iterative matching until convergence
  - Calls: `best_pair()` for each residue

- **`best_pair()`** (find_pair.c, line 744): Finds best partner for a base
  - Complexity: Medium
  - Purpose: For base i, finds base j with minimum distance that passes check_pair
  - Calls: `check_pair()`

- **`all_pairs()`** (find_pair.c, line 625): Finds ALL possible pairs
  - Complexity: Very High (O(n²))
  - Purpose: Exhaustive search of all base combinations
  - Calls: `check_pair()` for every i,j combination

### 3. **Reference Frame & Coordinate Systems**

#### Frame Calculation (`cmn_fncs.c`, `ana_fncs.c`)

- **`base_frame()`** (cmn_fncs.c, line 382): **CRITICAL FUNCTION**
  - Complexity: High
  - Purpose: Calculates reference frame for each base using ring atoms
  - Algorithm: Least-squares fitting to standard base structure
  - Output: `orien[][]` (9 values: 3x3 rotation matrix) and `org[]` (3 values: origin)
  - Calls: `set_std_base_pdb()`, `read_pdb()`, `find_1st_atom()`, `ls_fitting()`, `mst2orien()`

- **`ref_frames()`** (ana_fncs.c): Calculates reference frames for base pairs
  - Complexity: Very High
  - Purpose: Determines orientation and origin for paired bases
  - Calls: `base_frame()`, `refs_right_left()`, `refs_i_j()`

- **`refs_right_left()`** (x3dna_fncs.h, line 262): Extracts frames for left/right bases
  - Complexity: Low
  - Purpose: Gets rotation matrices and origins for base pair i

- **`refs_i_j()`** (x3dna_fncs.h, line 264): Gets frames for bases i and j
  - Complexity: Low

### 4. **Parameter Calculation**

#### Base Pair Parameters (`ana_fncs.c`, `cmn_fncs.c`)

- **`bpstep_par()`** (x3dna_fncs.h, line 37): **CRITICAL FUNCTION**
  - Complexity: High
  - Purpose: Calculates 6 base pair step parameters (Shift, Slide, Rise, Tilt, Roll, Twist)
  - Input: Two rotation matrices (r1, r2) and two origins (org1, org2)
  - Output: 6 parameters in `pars[]` array
  - Algorithm: Standard 3DNA parameter calculation
  - **This is where final parameters are computed**

- **`helical_par()`** (x3dna_fncs.h, line 39): Calculates helical parameters
  - Complexity: High
  - Purpose: Alternative parameter set for helical structures
  - Similar to bpstep_par but with different reference

- **`calculate_more_bppars()`** (cmn_fncs.c): **DEBUG POINT**
  - Complexity: High
  - Purpose: Calculates additional base pair parameters
  - **This function has debug statements for frame extraction**

- **`get_parameters()`** (x3dna_fncs.h, line 51): Orchestrates parameter calculation
  - Complexity: Very High
  - Purpose: Calculates all parameters for all base pairs
  - Calls: `bpstep_par()`, `helical_par()`, `get_mtwist()`, `bz_check()`

### 5. **PDB File I/O**

#### File Reading (`app_fncs.c`)

- **`read_pdb()`** (x3dna_fncs.h, line 303): Reads PDB file
  - Complexity: Medium
  - Purpose: Parses PDB format, extracts atom coordinates
  - Returns: Number of atoms, populates arrays
  - Calls: `number_of_atoms()`, `normalize_resName_atomName()`

- **`number_of_atoms()`** (x3dna_fncs.h, line 302): Counts atoms in PDB
  - Complexity: Low
  - Purpose: Pre-allocates memory

- **`residue_idx()`** (x3dna_fncs.h, line 322): Creates residue index
  - Complexity: Medium
  - Purpose: Maps atoms to residues
  - Returns: `seidx[][]` matrix (start/end atom indices per residue)

- **`write_pdb()`** (x3dna_fncs.h, line 315): Writes PDB file
  - Complexity: Low

### 6. **Sequence & Base Identification**

#### Sequence Processing (`app_fncs.c`)

- **`get_seq()`** (x3dna_fncs.h, line 333): Extracts sequence
  - Complexity: Medium
  - Purpose: Determines base sequence and RY (purine/pyrimidine) classification
  - Output: `bseq[]` (one-letter codes), `RY[]` (1=purine, 0=pyrimidine, -1=not base)

- **`get_bpseq()`** (x3dna_fncs.h, line 335): Gets base pair sequence
  - Complexity: Medium
  - Purpose: Creates `bp_seq[][]` matrix with base letters for each strand

- **`residue_ident()`** (x3dna_fncs.h, line 328): Identifies residue type
  - Complexity: Medium
  - Purpose: Determines if residue is nucleotide, amino acid, or other
  - Returns: res_type (1=purine, 0=pyrimidine, -1=amino acid, -2=other)

### 7. **Hydrogen Bond Analysis**

#### H-Bond Detection (`cmn_fncs.c`, `app_fncs.c`)

- **`hb_atompair()`** (x3dna_fncs.h, line 467): Finds H-bond atom pairs
  - Complexity: Medium
  - Purpose: Identifies potential H-bond donor/acceptor pairs
  - Calls: `good_hbatoms()`, `update_hb_idx()`

- **`validate_hbonds()`** (x3dna_fncs.h, line 469): Validates H-bonds
  - Complexity: Medium
  - Purpose: Checks if H-bonds meet distance/angle criteria

- **`hbond_info()`** (x3dna_fncs.h, line 448): Generates H-bond information
  - Complexity: Medium
  - Purpose: Creates detailed H-bond report

### 8. **Geometric Calculations**

#### Vector & Matrix Operations (`nrutil.c`)

- **`dot()`** (nrutil.c, line 547): Vector dot product
  - Complexity: O(1)
  
- **`cross()`** (nrutil.c, line 556): Vector cross product
  - Complexity: O(1)
  
- **`veclen()`** (nrutil.c, line 573): Vector length
  - Complexity: O(1)
  
- **`vec_norm()`** (nrutil.c, line 578): Normalize vector
  - Complexity: O(1)
  
- **`magang()`** (nrutil.c, line 600): Angle between vectors
  - Complexity: O(1)
  
- **`multi_matrix()`** (nrutil.c, line 639): Matrix multiplication
  - Complexity: O(n³)

#### Least Squares Fitting (`app_fncs.c`)

- **`ls_fitting()`** (x3dna_fncs.h, line 347): Least squares fitting
  - Complexity: High
  - Purpose: Aligns two sets of points (standard base to experimental)
  - Algorithm: SVD or similar
  - Returns: Rotation matrix R and origin
  - **Critical for base_frame() calculation**

- **`ls_plane()`** (x3dna_fncs.h, line 349): Least squares plane fitting
  - Complexity: Medium
  - Purpose: Fits plane to set of points

### 9. **Base Pair Ordering & Helix Detection**

#### Ordering Functions (`find_pair.c`)

- **`re_ordering()`** (find_pair.c, line 1426): **CRITICAL FUNCTION**
  - Complexity: Very High
  - Purpose: Reorders base pairs to follow 5'→3' direction
  - Algorithm: Complex graph traversal
  - Calls: `bp_context()`, `locate_helix()`, `five2three()`, `check_zdna()`

- **`bp_context()`** (find_pair.c, line 797): Analyzes base pair context
  - Complexity: High
  - Purpose: Determines neighbors, helix breaks, overlaps
  - Calls: `bp_coplanar()`, `is_circular()`

- **`locate_helix()`** (find_pair.c, line 919): Identifies helical regions
  - Complexity: High
  - Purpose: Groups consecutive base pairs into helices
  - Algorithm: Graph traversal with end detection

- **`five2three()`** (find_pair.c, line 1296): Ensures 5'→3' ordering
  - Complexity: Very High
  - Purpose: Swaps base pairs to correct strand direction
  - Calls: `first_step()`, `check_direction()`, `check_strand2()`, `check_rise()`

### 10. **Memory Management**

#### Utilities (`nrutil.c`)

- **`dmatrix()`** / **`dmatrix_nr()`**: Allocate double matrix
- **`lmatrix()`** / **`lmatrix_nr()`**: Allocate long matrix
- **`cmatrix()`** / **`cmatrix_nr()`**: Allocate char matrix
- **`dvector()`** / **`dvector_nr()`**: Allocate double vector
- **`lvector()`** / **`lvector_nr()`**: Allocate long vector
- **`cvector()`** / **`cvector_nr()`**: Allocate char vector
- Corresponding `free_*()` functions

### 11. **String & Pattern Matching**

#### Regex Engine (`fncs_slre.c`)

- **`slre_match()`**: Main regex matching function
- **`lux_match()`**: Case-sensitive matching
- **`lux_ncmatch()`**: Case-insensitive matching
- **`lux_bcmatch()`**: Basic case-sensitive match
- Full regex compiler and VM implementation

## Function Call Graph (Major Paths)

### find_pair_original Workflow:
```
main()
  └─> set_my_globals()
  └─> fp_cmdline()
  └─> handle_str()
      └─> number_of_atoms()
      └─> read_pdb()
      └─> residue_idx()
      └─> get_seq()
      └─> populate_nt_info()
      └─> duplex() [if duplex mode]
          └─> base_info()
              └─> base_frame() [for each residue]
                  └─> set_std_base_pdb()
                  └─> read_pdb()
                  └─> find_1st_atom()
                  └─> ls_fitting()
                      └─> multi_matrix()
                      └─> transpose_matrix()
          └─> ring_oidx()
          └─> find_bestpair() [if not -p]
              └─> best_pair()
                  └─> check_pair()
                      └─> get_zoave()
                      └─> ratom_xyz()
                      └─> get_oarea()
                      └─> hb_atompair()
          └─> all_pairs() [if -p]
              └─> check_pair() [for all i,j]
          └─> re_ordering()
              └─> bp_context()
              └─> locate_helix()
              └─> five2three()
                  └─> check_direction()
                  └─> check_strand2()
          └─> write_bestpairs()
          └─> write_helix()
          └─> x3dna_input()
  └─> clear_my_globals()
```

### analyze_original Workflow:
```
main()
  └─> set_my_globals()
  └─> analyze_cmdline()
  └─> process_str()
      └─> read_input()
      └─> read_pdb()
      └─> residue_idx()
      └─> pair_checking()
      └─> drct_checking()
      └─> get_bpseq()
      └─> get_seq()
      └─> atom_list()
      └─> ref_frames()
          └─> base_frame()
      └─> get_nt_bb_torsion()
      └─> hb_information()
      └─> base_overlap()
      └─> get_parameters()
          └─> bpstep_par()
          └─> helical_par()
          └─> get_mtwist()
      └─> write_mst()
      └─> check_simple_parameters()
      └─> global_analysis()
      └─> backbone_torsion()
      └─> helix_radius()
      └─> get_helix_axis()
  └─> clear_my_globals()
```

## Complexity Analysis

### Time Complexity

| Function | Complexity | Notes |
|----------|------------|-------|
| `all_pairs()` | O(n²) | Checks all base combinations |
| `find_bestpair()` | O(n²) worst case | Greedy matching, may iterate |
| `check_pair()` | O(k) | k = number of ring atoms (~9) |
| `base_frame()` | O(n) | n = number of residues |
| `re_ordering()` | O(n²) | Graph operations on base pairs |
| `bp_context()` | O(n²) | Distance calculations |
| `ls_fitting()` | O(m³) | m = number of points (~9) |
| `read_pdb()` | O(a) | a = number of atoms |

### Space Complexity

| Function | Complexity | Notes |
|----------|------------|-------|
| `duplex()` | O(n) | Stores frames, coordinates for n residues |
| `all_pairs()` | O(n²) | Pair info matrix |
| `re_ordering()` | O(n) | Helix indices, ordering arrays |

## Critical Functions for Debugging

### 1. **`check_pair()`** (cmn_fncs.c:492)
- **Why Critical**: Validates base pairs, calculates dir_z
- **Debug Points**: 
  - dir_z calculation (line ~495 in original)
  - z-axis dot product
  - Return value array

### 2. **`calculate_more_bppars()`** (cmn_fncs.c)
- **Why Critical**: Calculates base pair parameters
- **Debug Points**:
  - Frame extraction from orien array
  - Matrix construction (r1, r2)
  - Parameter calculation

### 3. **`base_frame()`** (cmn_fncs.c:382)
- **Why Critical**: Creates reference frames
- **Debug Points**:
  - Ring atom matching
  - Least squares fitting
  - orien array population

### 4. **`bpstep_par()`** (ana_fncs.c)
- **Why Critical**: Final parameter calculation
- **Debug Points**:
  - Input matrices (r1, r2)
  - Output parameters (6 values)

### 5. **`refs_right_left()`** / **`refs_i_j()`**
- **Why Critical**: Extracts frames from orien array
- **Debug Points**:
  - Array indexing
  - Matrix extraction

## Task Grouping

### Task Group 1: **PDB Input & Parsing**
**Files**: `app_fncs.c`, `find_pair.c`, `analyze.c`
**Functions**:
- `read_pdb()`, `number_of_atoms()`, `residue_idx()`
- `get_seq()`, `get_bpseq()`, `residue_ident()`
- `normalize_resName_atomName()`, `deduce_misc()`
**Complexity**: Medium
**Purpose**: Read and parse PDB files, extract sequence information

### Task Group 2: **Reference Frame Calculation**
**Files**: `cmn_fncs.c`, `ana_fncs.c`
**Functions**:
- `base_frame()`, `ref_frames()`, `refs_right_left()`, `refs_i_j()`
- `ls_fitting()`, `ls_plane()`
- `define_frame_by_3atoms()`
**Complexity**: High
**Purpose**: Calculate coordinate frames for bases and base pairs

### Task Group 3: **Base Pair Detection**
**Files**: `cmn_fncs.c`, `find_pair.c`
**Functions**:
- `check_pair()`, `best_pair()`, `find_bestpair()`, `all_pairs()`
- `hb_atompair()`, `validate_hbonds()`, `get_hbond_ij()`
- `ratom_xyz()`, `get_oarea()`, `get_zoave()`
**Complexity**: Very High
**Purpose**: Identify and validate base pairs

### Task Group 4: **Parameter Calculation**
**Files**: `ana_fncs.c`, `cmn_fncs.c`
**Functions**:
- `bpstep_par()`, `helical_par()`, `get_parameters()`
- `calculate_more_bppars()`, `get_mtwist()`
- `simple_bp_pars()`, `simple_step_heli_pars()`
**Complexity**: High
**Purpose**: Calculate structural parameters (Shift, Slide, Rise, Tilt, Roll, Twist)

### Task Group 5: **Base Pair Ordering & Helix Detection**
**Files**: `find_pair.c`
**Functions**:
- `re_ordering()`, `bp_context()`, `locate_helix()`
- `five2three()`, `check_direction()`, `check_strand2()`
- `check_zdna()`, `is_circular()`
**Complexity**: Very High
**Purpose**: Order base pairs correctly and identify helical regions

### Task Group 6: **Geometric Utilities**
**Files**: `nrutil.c`, `app_fncs.c`
**Functions**:
- `dot()`, `cross()`, `veclen()`, `vec_norm()`, `magang()`
- `multi_matrix()`, `transpose_matrix()`, `identity_matrix()`
- `p1p2_dist()`, `torsion()`, `torsion2()`
**Complexity**: Low-Medium
**Purpose**: Vector/matrix operations and geometric calculations

### Task Group 7: **Output Generation**
**Files**: `find_pair.c`, `analyze.c`, `app_fncs.c`
**Functions**:
- `write_bestpairs()`, `write_helix()`, `write_mst()`
- `x3dna_input()`, `curves_input()`, `curves_plus_input()`
- `print_par()`, `output_ave_std()`, `print_pairinfo()`
**Complexity**: Medium
**Purpose**: Generate output files and reports

### Task Group 8: **Memory Management**
**Files**: `nrutil.c`
**Functions**:
- `dmatrix()`, `lmatrix()`, `cmatrix()`, `dvector()`, `lvector()`, `cvector()`
- `free_dmatrix()`, `free_lmatrix()`, `free_cmatrix()`, etc.
**Complexity**: Low
**Purpose**: Dynamic memory allocation and deallocation

### Task Group 9: **String & Pattern Matching**
**Files**: `fncs_slre.c`, `app_fncs.c`
**Functions**:
- `slre_match()`, `lux_match()`, `lux_ncmatch()`
- `case_strcmp()`, `case_strstr()`, `str_pmatch()`
**Complexity**: Medium
**Purpose**: Pattern matching and string operations

### Task Group 10: **Configuration & Initialization**
**Files**: `cmn_fncs.c`, `app_fncs.c`
**Functions**:
- `set_my_globals()`, `clear_my_globals()`
- `get_3dna_pars()`, `set_default_misc_pars()`
- `check_global_options()`, `overwrite_misc_pars()`
**Complexity**: Low-Medium
**Purpose**: Initialize global variables and read configuration

## Key Data Structures

### Global Variables (`struct_Gvars`)
- `misc_pars`: Parameter thresholds (hb_dist, max_dorg, etc.)
- `X3DNA_HOMEDIR`: Installation directory
- `ATOMLIST`, `BASELIST`: Atom and base name lists
- `CHAIN_MARKERS`: Chain identification markers

### Major Arrays
- `orien[][]`: 9 values per residue (3x3 rotation matrix, flattened)
- `org[][]`: 3 values per residue (origin coordinates)
- `xyz[][]`: 3 values per atom (x, y, z coordinates)
- `seidx[][]`: 2 values per residue (start/end atom indices)
- `pair_num[][]`: Base pair indices
- `bp_seq[][]`: Base pair sequence letters

## Dependencies Between Functions

### Critical Dependency Chain:
1. **PDB Reading** → `read_pdb()` → `residue_idx()` → `get_seq()`
2. **Frame Calculation** → `base_frame()` → `ls_fitting()` → `mst2orien()`
3. **Pair Detection** → `check_pair()` → `get_zoave()` → `ratom_xyz()` → `get_oarea()`
4. **Parameter Calculation** → `refs_right_left()` → `bpstep_par()` → output
5. **Ordering** → `bp_context()` → `locate_helix()` → `five2three()`

## Notes for Modernization

1. **Array Indexing**: Uses 1-based indexing (NR_END offset)
2. **Memory Management**: Manual allocation with custom wrappers
3. **Error Handling**: Uses `fatal()` for errors (exits program)
4. **Global State**: Heavy use of global variables (`Gvars`)
5. **Function Signatures**: Many functions take 10+ parameters
6. **Code Style**: Mixed naming conventions, minimal comments

## Summary Statistics

- **Total Functions**: ~400+ functions
- **Lines of Code**: 14,976
- **Source Files**: 7
- **Header Files**: 2
- **Executables**: 2
- **Complexity**: Very High (many O(n²) operations)
- **Main Algorithms**: Graph traversal, least squares fitting, greedy matching


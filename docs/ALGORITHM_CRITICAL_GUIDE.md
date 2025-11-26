# Critical Algorithm Guide for Rewriting X3DNA

This document outlines the essential concepts, algorithms, and implementation details needed to successfully rewrite the X3DNA base pair identification and parameter calculation system.

## Table of Contents
1. [Core Algorithm Overview](#core-algorithm-overview)
2. [Critical Mathematical Concepts](#critical-mathematical-concepts)
3. [Essential Data Structures](#essential-data-structures)
4. [Algorithm Workflow](#algorithm-workflow)
5. [Critical Functions & Dependencies](#critical-functions--dependencies)
6. [Implementation Details](#implementation-details)
7. [Common Pitfalls](#common-pitfalls)

---

## Core Algorithm Overview

The X3DNA algorithm has **two main phases**:

### Phase 1: Base Pair Identification (`find_pair`)
1. **Read PDB structure** → Extract atoms, residues, coordinates
2. **Identify bases** → Classify nucleotides (A, C, G, T/U) and determine RY (purine/pyrimidine)
3. **Calculate reference frames** → For each base, determine its 3D orientation using least-squares fitting to standard base templates
4. **Find base pairs** → Using greedy matching algorithm with geometric/hydrogen bond validation
5. **Order base pairs** → Reorder to follow 5'→3' direction and identify helical regions

### Phase 2: Parameter Calculation (`analyze`)
1. **Read base pair data** → From Phase 1 output (.inp file)
2. **Recalculate frames** → Ensure frames are consistent for paired bases
3. **Calculate step parameters** → For consecutive base pairs: Shift, Slide, Rise, Tilt, Roll, Twist
4. **Calculate helical parameters** → Alternative parameter set using helical axis reference

---

## Critical Mathematical Concepts

### 1. Reference Frame Representation

**Each base has a 3×3 rotation matrix (R) and origin point (O):**
- **R (orien[])**: 9 values stored as flattened 3×3 matrix: `[r11, r12, r13, r21, r22, r23, r31, r32, r33]`
- **O (org[])**: 3 values `[ox, oy, oz]` representing the base origin

**Frame extraction from array:**
```c
// orien[i] contains 9 values for base i
// Extract 3×3 matrix:
r[1][1] = orien[i][0]  // orien[i][(1-1)*3 + 1 - 1] = orien[i][0]
r[1][2] = orien[i][1]
r[1][3] = orien[i][2]
r[2][1] = orien[i][3]  // orien[i][(2-1)*3 + 1 - 1] = orien[i][3]
// etc.
// General formula: r[k][l] = orien[i][(k-1)*3 + l - 1]
// But uses 1-based indexing: r[k][l] = orien[i][(k-1)*3 + l]
```

**Critical**: The frame's z-axis is at indices `[6, 7, 8]` (indices 7-9 in 1-based):
- `z[1] = orien[i][6]`, `z[2] = orien[i][7]`, `z[3] = orien[i][8]`

### 2. Base Pair Frame Calculation

**For a base pair (i, j):**
- Left base frame: `(R1, O1)` from base i
- Right base frame: `(R2, O2)` from base j
- **Direction check**: `dir_z = dot(z1, z2)` where `z1 = orien[i][6:8]`, `z2 = orien[j][6:8]`
  - If `dir_z < 0`, bases point toward each other (correct for base pair)
  - If `dir_z > 0`, bases point away (may need to reverse one frame)

### 3. Least-Squares Fitting (Critical for base_frame)

**Algorithm**:
1. Load standard base template (from `Atomic_A.pdb`, `Atomic_C.pdb`, etc.)
2. Find matching ring atoms between experimental and standard structures
3. Use **least-squares fitting** to find rotation matrix R and translation t that minimize:
   ```
   Σ ||R·x_standard + t - x_experimental||²
   ```
4. Result: R becomes the base's rotation matrix, t becomes the origin

**Key function**: `ls_fitting()` - implements SVD-based or quaternion-based fitting

### 4. Base Pair Step Parameters (6 parameters)

**Given two consecutive base pair frames:**
- Frame 1: `(R1, O1)` for base pair i
- Frame 2: `(R2, O2)` for base pair i+1

**Calculate:**
1. **Midstep frame** `(R_mid, O_mid)`:
   - Rotation: Average of R1 and R2
   - Origin: Average of O1 and O2

2. **Transform to midpoint coordinate system**:
   - Express vectors in midpoint frame
   - Calculate relative displacement and rotation

3. **Extract 6 parameters**:
   - **Shift** (x-displacement): Translation along x-axis of midpoint frame
   - **Slide** (y-displacement): Translation along y-axis
   - **Rise** (z-displacement): Translation along z-axis
   - **Tilt** (rotation about x): Angle in degrees
   - **Roll** (rotation about y): Angle in degrees  
   - **Twist** (rotation about z): Angle in degrees

**Critical formula**: Standard 3DNA parameter calculation (implemented in `bpstep_par()`)

### 5. Hydrogen Bond Validation

**For base pair (i, j) to be valid:**
1. **Distance check**: N-N distance between bases must be in range `[min_dNN, max_dNN]` (typically 4.5 - 20 Å)
2. **Hydrogen bond check**: Must have at least `min_base_hb` hydrogen bonds (typically 1-2)
3. **H-bond criteria**:
   - Donor-acceptor distance: `hb_dist1` (typically 4.0 Å) or `hb_dist2` (typically 0.0, unused)
   - H-donor-acceptor angle: Must be within acceptable range
   - Atoms match pattern in `hb_atoms` (typically ".O.N" for O and N atoms)

4. **Plane angle check**: Angle between base planes must be `≤ max_plane_angle` (typically 65°)
5. **Origin distance check**: Distance between base origins must be `≤ max_dorg` (typically 15.0 Å)

---

## Essential Data Structures

### 1. Base Reference Frames
```c
// Per-residue arrays
double **orien;  // [1..num_residue][1..9] - 9 values per base (3×3 matrix flattened)
double **org;    // [1..num_residue][1..3] - 3 values per base (origin)
```

### 2. Atom Coordinates
```c
double **xyz;    // [1..num_atoms][1..3] - x, y, z for each atom
char **AtomName; // [1..num_atoms][0..4] - Atom names (" C1'", " N1 ", etc.)
char **ResName;  // [1..num_atoms][0..3] - Residue names ("  A", "  C", etc.)
```

### 3. Residue Indexing
```c
long **seidx;   // [1..num_residue][1..2] - [start_atom_idx, end_atom_idx] for each residue
```

### 4. Base Pair Data
```c
long **pair_num;  // [1..ds][1..num_bp] - Residue indices for each base pair
char **bp_seq;    // [1..ds][1..num_bp] - Base letters ('A', 'C', 'G', 'T')
long *RY;         // [1..num_residue] - 1=purine, 0=pyrimidine, -1=not base
```

### 5. Parameter Output
```c
double **bp_step_par;  // [1..num_step][1..6] - Shift, Slide, Rise, Tilt, Roll, Twist
double **bp_heli_par;  // [1..num_step][1..6] - Helical parameters
```

---

## Algorithm Workflow

### Step-by-Step: Base Pair Identification

```
1. READ PDB FILE
   ├─> number_of_atoms() - Count atoms
   ├─> read_pdb() - Parse PDB format
   └─> Extract: AtomName, ResName, ChainID, ResSeq, xyz, Miscs

2. INDEX RESIDUES
   ├─> residue_idx() - Map atoms to residues
   └─> Output: seidx[][] - [start_idx, end_idx] per residue

3. IDENTIFY BASES
   ├─> get_seq() - Extract sequence
   ├─> residue_ident() - Classify as base/amino acid/other
   └─> Output: bseq[], RY[] (purine/pyrimidine classification)

4. CALCULATE BASE FRAMES (CRITICAL)
   ├─> For each residue i:
   │   ├─> base_frame() - Calculate reference frame
   │   │   ├─> Load standard base PDB (Atomic_A.pdb, etc.)
   │   │   ├─> Find matching ring atoms (C1', N1/N9, C4, etc.)
   │   │   ├─> ls_fitting() - Least-squares fit to standard
   │   │   └─> Store: orien[i][1..9], org[i][1..3]
   │   └─> ring_oidx() - Index ring atoms for later use

5. FIND BASE PAIRS
   ├─> find_bestpair() - Greedy matching algorithm
   │   ├─> Iterate until no new pairs found:
   │   │   ├─> For each unmatched base i:
   │   │   │   ├─> best_pair(i) - Find best partner j
   │   │   │   │   ├─> Try all other bases j
   │   │   │   │   ├─> check_pair(i, j) - Validate pair
   │   │   │   │   └─> Return j with minimum distance that passes
   │   │   │   └─> Check if mutual: best_pair(j) == i
   │   │   │       └─> If yes, mark both as matched
   │   └─> Output: base_pairs[][] - Matrix of pair data

6. REORDER BASE PAIRS
   ├─> re_ordering() - Complex reordering algorithm
   │   ├─> bp_context() - Analyze neighbor relationships
   │   ├─> locate_helix() - Identify helical regions
   │   ├─> five2three() - Ensure 5'→3' ordering
   │   └─> Output: bp_idx[], helix_marker[]

7. WRITE OUTPUT
   └─> x3dna_input() - Write .inp file for analyze phase
```

### Step-by-Step: Parameter Calculation

```
1. READ INPUT FILE (.inp from find_pair)
   ├─> read_input() - Parse base pair data
   └─> Extract: pdbfile, num_bp, pair_num[][], ds, ip, hetatm

2. READ PDB FILE (again)
   └─> read_pdb() - Get atom coordinates

3. RECALCULATE REFERENCE FRAMES
   ├─> ref_frames() - Calculate frames for base pairs
   │   ├─> For each base pair (i, j):
   │   │   ├─> base_frame(i) - Calculate frame for base i
   │   │   ├─> base_frame(j) - Calculate frame for base j
   │   │   ├─> Check orientation (dir_z = dot(z_i, z_j))
   │   │   └─> Adjust if needed (reverse one frame)
   │   └─> Store: bp_orien[][], bp_org[][]

4. CALCULATE STEP PARAMETERS
   ├─> get_parameters() - Calculate for all steps
   │   └─> For each consecutive pair (i, i+1):
   │       ├─> refs_i_j(i, i+1) - Extract frames
   │       ├─> bpstep_par(r1, o1, r2, o2, pars, ...) - Calculate 6 params
   │       ├─> helical_par(r1, o1, r2, o2, pars, ...) - Calculate helical
   │       └─> Store: bp_step_par[][], bp_heli_par[][]

5. OUTPUT RESULTS
   └─> write_mst() - Write parameter file
```

---

## Critical Functions & Dependencies

### Must Implement Correctly (in order of dependency):

1. **`ls_fitting()`** - Least-squares fitting
   - **Input**: Two sets of 3D points (standard base, experimental base)
   - **Output**: Rotation matrix R, translation vector t, RMS fit
   - **Algorithm**: SVD or quaternion-based Kabsch algorithm
   - **Critical**: This determines base frame accuracy

2. **`base_frame()`** - Calculate base reference frame
   - **Depends on**: `ls_fitting()`, standard base templates
   - **Output**: `orien[i][1..9]`, `org[i][1..3]`
   - **Critical**: Frame orientation determines all subsequent calculations

3. **`get_zoave()`** - Calculate average z-axis
   - **Purpose**: Get base normal vector from ring atoms
   - **Used in**: `check_pair()` for plane angle calculation

4. **`check_pair()`** - Validate base pair
   - **Depends on**: `base_frame()`, `get_zoave()`, `ratom_xyz()`, `hb_atompair()`
   - **Checks**: Distance, H-bonds, plane angle, geometry
   - **Output**: Base pair type ID and validation status
   - **Critical**: This determines which pairs are valid

5. **`refs_i_j()`** - Extract frames for two bases
   - **Purpose**: Get rotation matrices and origins for bases i and j
   - **Input**: `orien[][]`, `org[][]`, base indices
   - **Output**: Two 3×3 matrices `r1[][]`, `r2[][]` and two origins `o1[]`, `o2[]`
   - **Critical**: Must extract frames correctly from flattened array

6. **`bpstep_par()`** - Calculate 6 step parameters
   - **Depends on**: `refs_i_j()` 
   - **Input**: Two frames `(r1, o1)` and `(r2, o2)`
   - **Output**: 6 parameters in array: `[Shift, Slide, Rise, Tilt, Roll, Twist]`
   - **Algorithm**: Standard 3DNA parameter calculation
   - **Critical**: Final output - must match 3DNA reference implementation exactly

7. **`helical_par()`** - Calculate helical parameters
   - Similar to `bpstep_par()` but uses helical axis reference

---

## Implementation Details

### Critical Array Indexing

**1-based indexing with NR_END offset:**
```c
#define NR_END 1

// Allocate: dmatrix(1, n, 1, m) creates array[n+1][m+1]
// Access: arr[i][j] where i ∈ [1, n], j ∈ [1, m]
// Actual memory: arr[0..n][0..m], but indices [0][*] and [*][0] unused
```

**Frame array indexing:**
```c
// orien[i] is 9 values: [0..8] or [1..9] depending on implementation
// Standard form (0-based): orien[i][0..8]
// Extract matrix element: r[k][l] = orien[i][(k-1)*3 + l] for 1-based k,l
// But code uses: r[k][l] = orien[i][(k-1)*3 + l - 1] for 1-based everything
```

### Critical Coordinate System

**Base frame definition:**
- **x-axis**: Points from base to sugar (C1' direction)
- **y-axis**: Perpendicular to x, in base plane
- **z-axis**: Normal to base plane (points up from base)
- **Origin**: Typically center of base ring (N1/N9 for purines, N1 for pyrimidines)

**Direction convention:**
- For paired bases, z-axes should point toward each other: `dot(z1, z2) < 0`
- If `dot(z1, z2) > 0`, one frame must be reversed (flip y and z)

### Critical Parameter Thresholds

From `misc_pars` (can be configured):
- `min_dNN = 4.5 Å` - Minimum N-N distance
- `max_dNN = XBIG` (~1e18) - Maximum N-N distance (effectively infinite)
- `max_dorg = 15.0 Å` - Maximum origin distance
- `max_plane_angle = 65°` - Maximum angle between base planes
- `min_base_hb = 1` - Minimum number of hydrogen bonds
- `hb_dist1 = 4.0 Å` - Maximum H-bond distance
- `helix_break = 7.5 Å` - Distance threshold for helix breaks

### Standard Base Templates

**Required PDB files** (loaded from X3DNA_HOMEDIR):
- `Atomic_A.pdb` - Adenine standard structure
- `Atomic_C.pdb` - Cytosine standard structure
- `Atomic_G.pdb` - Guanine standard structure
- `Atomic_T.pdb` - Thymine standard structure
- `Atomic_U.pdb` - Uracil standard structure

These contain ideal base geometries with ring atoms at standard positions.

---

## Common Pitfalls

### 1. Array Indexing Errors
- **Problem**: Mixing 0-based and 1-based indexing
- **Solution**: Be consistent - original uses 1-based with NR_END offset

### 2. Frame Extraction Errors
- **Problem**: Incorrectly extracting 3×3 matrix from flattened array
- **Solution**: Carefully map `orien[i][k]` to matrix element `r[row][col]`

### 3. Direction/Handedness Errors
- **Problem**: Not checking `dir_z` or reversing frames incorrectly
- **Solution**: Always check `dot(z1, z2) < 0` for paired bases

### 4. Least-Squares Fitting Errors
- **Problem**: Incorrect rotation matrix from fitting
- **Solution**: Use robust SVD-based algorithm, ensure proper atom matching

### 5. Parameter Calculation Errors
- **Problem**: Final parameters don't match 3DNA reference
- **Solution**: Use exact same algorithm as `bpstep_par()` - don't deviate

### 6. Hydrogen Bond Detection Errors
- **Problem**: Missing valid H-bonds or accepting invalid ones
- **Solution**: Match exact distance/angle criteria, check atom types correctly

### 7. Base Pair Matching Errors
- **Problem**: Missing valid pairs or finding spurious ones
- **Solution**: Use same greedy matching algorithm, check all validation criteria

### 8. Coordinate System Confusion
- **Problem**: Wrong frame orientation or origin
- **Solution**: Match standard base template exactly, use same reference atoms

---

## Testing Strategy

### Essential Test Cases:
1. **Single base pair** - Verify frame calculation and parameters
2. **B-DNA helix** - Verify step parameters match reference values
3. **A-DNA helix** - Different geometry, test parameter ranges
4. **Z-DNA helix** - Special handling needed (check_zdna function)
5. **Mismatched pairs** - Ensure validation rejects invalid pairs
6. **Circular structures** - Test 5'→3' ordering with circular sequences

### Validation:
- Compare outputs with original 3DNA v2.4
- Use DSSR (3DNA successor) as reference
- Test on PDB structures: 1BNA, 1EHZ, 2BNA, etc.

---

## Summary

**Critical for successful rewrite:**
1. ✅ Correct least-squares fitting implementation
2. ✅ Accurate base frame calculation and storage
3. ✅ Proper frame extraction from flattened arrays
4. ✅ Exact `bpstep_par()` algorithm matching 3DNA reference
5. ✅ Correct hydrogen bond validation
6. ✅ Proper base pair matching algorithm
7. ✅ Correct handling of frame directions/orientations
8. ✅ Matching parameter thresholds and validation criteria

**Most error-prone areas:**
- Array indexing (1-based vs 0-based)
- Frame matrix extraction
- Least-squares fitting
- Parameter calculation algorithm

Focus on getting these right first, as errors cascade through the entire pipeline.


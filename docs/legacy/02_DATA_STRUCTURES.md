# Legacy Code Data Structures

**Date**: 2025-01-XX  
**Purpose**: Complete reference for all data structures, memory allocation patterns, and indexing conventions  
**Status**: Comprehensive guide for understanding legacy data organization

---

## Table of Contents

1. [Indexing Conventions](#indexing-conventions)
2. [Atom-Level Data Structures](#atom-level-data-structures)
3. [Residue-Level Data Structures](#residue-level-data-structures)
4. [Base Frame Data Structures](#base-frame-data-structures)
5. [Base Pair Data Structures](#base-pair-data-structures)
6. [Geometric Data Structures](#geometric-data-structures)
7. [Memory Allocation Pattern](#memory-allocation-pattern)
8. [Array Access Patterns](#array-access-patterns)
9. [Data Structure Examples](#data-structure-examples)

---

## Indexing Conventions

### 1-Based Indexing

**Critical**: All arrays in legacy code use **1-based indexing**, not 0-based.

**Definition**:
```c
#define NR_END 1

// Array indices start at 1, not 0
double *vec = dvector(1, n);  // Valid indices: [1..n]
vec[1]  // First element
vec[n]  // Last element
vec[0]  // Undefined/unused (memory allocated but not accessed)
```

**Impact on Implementation**:
- Modern C++ code uses 0-based indexing
- **Conversion required**: `legacy_index = modern_index + 1`
- **Range**: Legacy uses `[1..n]`, modern uses `[0..n-1]`

**Example**:
```c
// Legacy code:
for (i = 1; i <= num_residue; i++) {
    process_residue(i);  // i ∈ [1, num_residue]
}

// Modern equivalent:
for (i = 0; i < num_residue; i++) {
    process_residue(i);  // i ∈ [0, num_residue-1]
    // Legacy index would be: legacy_i = i + 1
}
```

---

## Atom-Level Data Structures

### Core Atom Arrays

All atom arrays are allocated with `num_atoms` elements, indexed from 1 to `num_atoms`.

#### Atom Names
```c
char **AtomName;  // [1..num_atoms][0..4]
// Format: 5-character string (padded with spaces)
// Examples: " C1'", " N1 ", " O2 ", " P  "
// Note: [0..4] means indices 0-4 (inclusive), total 5 chars
// Access: AtomName[i] points to char array, AtomName[i][j] for j ∈ [0,4]
```

**Allocation**:
```c
AtomName = cmatrix(1, num_atoms, 0, 4);
// Creates: num_atoms rows, each with 5 characters (indices 0-4)
```

**Example Usage**:
```c
// Check if atom is " C1'"
if (is_equal_string(AtomName[i], " C1'")) {
    // Process C1' atom
}
```

---

#### Residue Names
```c
char **ResName;  // [1..num_atoms][0..3]
// Format: 3-character string (padded with spaces)
// Examples: "  A", "  C", "  G", "  T", "  U", " GLY", " PRO"
// Note: All atoms in a residue share the same ResName value
```

**Allocation**:
```c
ResName = cmatrix(1, num_atoms, 0, 3);
```

---

#### Chain Identifiers
```c
char *ChainID;  // [1..num_atoms]
// Format: Single character
// Examples: 'A', 'B', ' ', '1'
// Note: Single character per atom (not a string array)
```

**Allocation**:
```c
ChainID = cvector(1, num_atoms);
```

---

#### Residue Sequence Numbers
```c
long *ResSeq;  // [1..num_atoms]
// Format: Integer residue number from PDB file
// Examples: 1, 2, 3, ..., 100
// Note: PDB residue numbers (may have gaps, may not start at 1)
```

**Allocation**:
```c
ResSeq = lvector(1, num_atoms);
```

---

#### Atomic Coordinates
```c
double **xyz;  // [1..num_atoms][1..3]
// Format: 3D coordinates
// xyz[i][1] = x-coordinate
// xyz[i][2] = y-coordinate
// xyz[i][3] = z-coordinate
// Units: Angstroms (Å)
```

**Allocation**:
```c
xyz = dmatrix(1, num_atoms, 1, 3);
```

**Example Usage**:
```c
// Calculate distance between atoms i and j
double dist = sqrt(
    pow(xyz[j][1] - xyz[i][1], 2) +
    pow(xyz[j][2] - xyz[i][2], 2) +
    pow(xyz[j][3] - xyz[i][3], 2)
);
```

---

#### Miscellaneous Fields
```c
char **Miscs;  // [1..num_atoms][0..NMISC]
// Format: Array of miscellaneous PDB fields
// Miscs[i][0]: Record type ('A' for ATOM, 'H' for HETATM)
// Miscs[i][1]: Alternative location indicator
// Miscs[i][2]: Insertion code
// Miscs[i][3..NMISC]: Additional fields
// NMISC: Typically 4 or more
```

**Allocation**:
```c
Miscs = cmatrix(1, num_atoms, 0, NMISC);
```

**Example Usage**:
```c
// Check insertion code
char insertion = Miscs[i][2];  // 'A', 'B', ' ', etc.

// Check if ATOM record
if (Miscs[i][0] == 'A') {
    // Standard ATOM record
}
```

---

### Atom Type Indices

```c
long *idx;  // [1..num_atoms]
// Format: Atom type classification
// idx[i] = 1: Nitrogen (N)
// idx[i] = 2: Oxygen (O)
// idx[i] = 3: Other special atoms
// idx[i] = 4: Other atoms (e.g., Carbon)
// Note: Used for H-bond validation (good_hbatoms)
```

**Allocation**:
```c
idx = lvector(1, num_atoms);
```

**Populated by**:
```c
atom_idx(num_atoms, AtomName, NULL, idx);
// Sets idx[i] based on AtomName[i]
```

---

## Residue-Level Data Structures

### Residue Index Mapping

**Purpose**: Maps sequential residue indices to atom ranges.

```c
long **seidx;  // [1..num_residue][1..2]
// seidx[i][1] = start atom index (first atom of residue i)
// seidx[i][2] = end atom index (last atom of residue i, inclusive)
// 
// Example:
//   seidx[5][1] = 42
//   seidx[5][2] = 50
//   Means residue 5 contains atoms 42, 43, 44, ..., 50
//   
//   To iterate atoms in residue i:
//   for (atom_idx = seidx[i][1]; atom_idx <= seidx[i][2]; atom_idx++)
```

**Allocation**:
```c
seidx = lmatrix(1, num_residue, 1, 2);
// Or: seidx = residue_idx(num_atoms, ResSeq, Miscs, ChainID, ResName, &num_residue);
```

**How It's Created** (`residue_idx()`):
1. For each atom, create unique residue identifier: `"{ResName}{ChainID}{ResSeq}{InsertionCode}"`
2. Find boundaries where residue identifier changes
3. Set `seidx[i][1]` = first atom of residue i
4. Set `seidx[i][2]` = last atom of residue i

**Critical Properties**:
- **Sequential numbering**: Residues numbered 1, 2, 3, ..., `num_residue`
- **All residues included**: Amino acids, nucleotides, water, ions, etc.
- **Order preserved**: Maintains order from PDB file
- **1-based indexing**: Residue indices start at 1

---

### Residue Classification

```c
long *RY;  // [1..num_residue]
// Format: Purine/Pyrimidine classification
// RY[i] = 1:  Purine (A, G)
// RY[i] = 0:  Pyrimidine (C, T, U)
// RY[i] = -1: Not a nucleotide (amino acid, water, etc.)
// RY[i] = -2: Invalid/unrecognized
```

**Populated by**:
```c
get_seq(num_residue, seidx, AtomName, ResName, ..., bseq, RY);
// Calls residue_ident() for each residue
```

**Usage**:
```c
// Skip non-nucleotides
if (RY[i] < 0) {
    continue;  // Not a nucleotide
}

// Check if purine
if (RY[i] == 1) {
    // Purine (A or G)
    RingAtom_num = 9;  // 9 ring atoms
} else {
    // Pyrimidine (C, T, or U)
    RingAtom_num = 6;  // 6 ring atoms
}
```

---

### Base Sequence

```c
char *bseq;  // [1..num_residue]
// Format: Single character per residue
// Examples: 'A', 'C', 'G', 'T', 'U', 'I', 'X'
// 'X': Unknown/invalid base
// Note: Only nucleotides have valid base letters
```

**Allocation**:
```c
bseq = cvector(1, num_residue);
```

**Populated by**:
```c
get_seq(num_residue, seidx, AtomName, ..., bseq, RY);
```

**Usage**:
```c
// Get base type for residue i
char base = bseq[i];  // 'A', 'C', 'G', 'T', etc.

// Check base type
if (bseq[i] == 'A') {
    // Adenine
}
```

---

## Base Frame Data Structures

### Rotation Matrices

**Critical**: Rotation matrices are stored as **flattened 9-element arrays**, not 3×3 matrices.

```c
double **orien;  // [1..num_residue][1..9]
// Format: Flattened 3×3 rotation matrix
//
// Mapping from 9-element array to 3×3 matrix:
//   R[1][1] = orien[i][1]  (row 1, col 1)
//   R[1][2] = orien[i][2]  (row 1, col 2)
//   R[1][3] = orien[i][3]  (row 1, col 3)
//   R[2][1] = orien[i][4]  (row 2, col 1)
//   R[2][2] = orien[i][5]  (row 2, col 2)
//   R[2][3] = orien[i][6]  (row 2, col 3)
//   R[3][1] = orien[i][7]  (row 3, col 1)
//   R[3][2] = orien[i][8]  (row 3, col 2)
//   R[3][3] = orien[i][9]  (row 3, col 3)
//
// Column vectors (axes):
//   x-axis: [orien[i][1], orien[i][4], orien[i][7]]
//   y-axis: [orien[i][2], orien[i][5], orien[i][8]]
//   z-axis: [orien[i][3], orien[i][6], orien[i][9]]
//
// Note: In 0-based indexing (modern), use indices [0..8]
//   x-axis: [orien[i][0], orien[i][3], orien[i][6]]
//   y-axis: [orien[i][1], orien[i][4], orien[i][7]]
//   z-axis: [orien[i][2], orien[i][5], orien[i][8]]
```

**Allocation**:
```c
orien = dmatrix(1, num_residue, 1, 9);
```

**Conversion Functions**:
```c
// Convert 3×3 matrix to 9-element array
mst2orien(orien[i], 0, R);  // R is 3×3 matrix

// Convert 9-element array to 3×3 matrix
orien2mst(orien[i], 0, mst);  // mst is 3×3 matrix
```

**Critical Details**:
- **Row-major storage**: Elements stored row-by-row
- **1-based indices**: Use `[1..9]` in legacy, `[0..8]` in modern
- **z-axis**: Points perpendicular to base plane (upward)
- **Right-handed**: Matrix represents right-handed coordinate system

---

### Frame Origins

```c
double **org;  // [1..num_residue][1..3]
// Format: Translation vector (origin)
// org[i][1] = x-coordinate
// org[i][2] = y-coordinate
// org[i][3] = z-coordinate
// Units: Angstroms (Å)
// Note: Origin is center of base (typically)
```

**Allocation**:
```c
org = dmatrix(1, num_residue, 1, 3);
```

**Calculation**:
- Set during least-squares fitting in `base_frame()`
- Calculated as: `org = ave_exyz - R × ave_sxyz`
- Represents translation to align template to experimental structure

---

### N and C1' Coordinates

```c
double **NC1xyz;  // [1..num_residue][1..7]
// Format: Special coordinates for N and C1' atoms
// NC1xyz[i][1] = N atom x-coordinate (N1 or N9)
// NC1xyz[i][2] = N atom y-coordinate
// NC1xyz[i][3] = N atom z-coordinate
// NC1xyz[i][4] = C1' atom x-coordinate
// NC1xyz[i][5] = C1' atom y-coordinate
// NC1xyz[i][6] = C1' atom z-coordinate
// NC1xyz[i][7] = Flag (1.0 if C1' exists, -1.0 if missing)
```

**Allocation**:
```c
NC1xyz = dmatrix(1, num_residue, 1, 7);
```

**Usage**:
- Used in `check_pair()` for distance calculations
- N-to-N distance: `rtn_val[4] = ||NC1xyz[i] - NC1xyz[j]||`
- C1' used for reference point

---

### O3' and P Coordinates

```c
double **o3_p;  // [1..num_residue][1..8]
// Format: O3' and P atom coordinates
// o3_p[i][1..3] = O3' coordinates (or -1.0 if missing)
// o3_p[i][4..6] = P coordinates (or -1.0 if missing)
// o3_p[i][7] = Flag for O3' (1.0 if exists, -1.0 if missing)
// o3_p[i][8] = Flag for P (1.0 if exists, -1.0 if missing)
```

**Allocation**:
```c
o3_p = dmatrix(1, num_residue, 1, 8);
```

**Usage**:
- Used in `re_ordering()` to determine 5'→3' direction
- Used in step parameter calculations

---

## Base Pair Data Structures

### Base Pair Residue Indices

```c
long **pair_num;  // [1..ds][1..num_bp]
// Format: Residue indices for each base pair
// pair_num[1][k] = residue index of first strand (strand 1)
// pair_num[2][k] = residue index of second strand (strand 2)
// ds = 1: Single strand
// ds = 2: Duplex (two strands)
//
// Example:
//   pair_num[1][3] = 15  // First strand, pair 3, residue 15
//   pair_num[2][3] = 8   // Second strand, pair 3, residue 8
//   Means pair 3 is between residue 15 (strand 1) and residue 8 (strand 2)
```

**Allocation**:
```c
pair_num = lmatrix(1, ds, 1, num_bp);
```

**Storage in analyze**:
```c
// In analyze, additional info stored:
// pair_num[ds+1][k]: Helix marker (if applicable)
```

---

### Base Pair Sequence

```c
char **bp_seq;  // [1..ds][1..num_bp]
// Format: Base letters for each strand
// bp_seq[1][k] = base letter for strand 1, pair k
// bp_seq[2][k] = base letter for strand 2, pair k
// Examples: 'A', 'C', 'G', 'T', 'U'
```

**Allocation**:
```c
bp_seq = cmatrix(0, ds, 1, num_bp);  // Note: starts at 0, not 1!
// Access: bp_seq[1..ds][1..num_bp]
```

**Usage**:
```c
// Get base pair sequence
char base1 = bp_seq[1][k];  // First strand
char base2 = bp_seq[2][k];  // Second strand
```

---

### Watson-Crick Information

```c
long *WC_info;  // [1..num_bp]
// Format: Watson-Crick pairing classification
// WC_info[k] = 1: Watson-Crick pair (A-T, G-C)
// WC_info[k] = 2: Wobble pair (G-T, G-U)
// WC_info[k] = 0: Other/non-standard
```

**Populated by**:
```c
check_Watson_Crick(num_bp, bp_seq, WC_info);
```

---

## Geometric Data Structures

### Ring Atom Indices

```c
long **ring_atom;  // [1..num_residue][1..19]
// Format: Indices of ring atoms for overlap calculation
//
// Structure:
//   ring_atom[i][1..9]   = Ring atom indices (or 0 if not present)
//   ring_atom[i][10]     = Number of ring atoms found
//   ring_atom[i][11..19] = Connected atom indices (for projection)
//
// Example for Adenine (9 ring atoms):
//   ring_atom[i][1] = 42   (atom index of " C4 ")
//   ring_atom[i][2] = 43   (atom index of " N3 ")
//   ring_atom[i][3] = 44   (atom index of " C2 ")
//   ...
//   ring_atom[i][9] = 50   (atom index of " N9 ")
//   ring_atom[i][10] = 9   (count)
```

**Allocation**:
```c
ring_atom = lmatrix(1, num_residue, 1, 19);
```

**Populated by**:
```c
ring_oidx(num_atoms, num_residue, RY, seidx, AtomName, xyz, idx, ring_atom);
```

**Ring Atoms List** (`RA_LIST`):
```c
static char *RingAtom[] = {
    " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
};
// First 6: Pyrimidine ring (C, T, U)
// All 9: Purine ring (A, G)
```

**Usage**:
- Used in `get_oarea()` for overlap area calculation
- Project ring atoms onto plane, calculate polygon intersection

---

## Memory Allocation Pattern

### Allocation Functions

**From Numerical Recipes (`nrutil.c`)**:

```c
// 1D arrays (vectors)
double *dvector(long nl, long nh);      // [nl..nh]
long   *lvector(long nl, long nh);      // [nl..nh]
char   *cvector(long nl, long nh);      // [nl..nh]

// 2D arrays (matrices)
double **dmatrix(long nrl, long nrh, long ncl, long nch);  // [nrl..nrh][ncl..nch]
long   **lmatrix(long nrl, long nrh, long ncl, long nch);  // [nrl..nrh][ncl..nch]
char   **cmatrix(long nrl, long nrh, long ncl, long nch);  // [nrl..nrh][ncl..nch]
```

### NR_END Offset

**Definition**:
```c
#define NR_END 1
```

**Purpose**: Allocate extra element at index 0 (unused) to enable 1-based indexing.

**Implementation**:
```c
// Allocation for vector [1..n]
double *dvector(long nl, long nh) {
    // Allocates: (nh - nl + 1 + NR_END) elements
    // Actual memory: indices [0..nh]
    // Returns pointer adjusted so vec[1] is first valid element
    // vec[0] exists but is never accessed
}
```

**Example**:
```c
double *vec = dvector(1, 10);
// Allocates 11 elements: vec[0..10]
// Valid access: vec[1] through vec[10]
// Invalid access: vec[0] (undefined, but memory exists)
```

### Deallocation Functions

```c
// Must use corresponding free function
free_dvector(double *v, long nl, long nh);
free_lvector(long *v, long nl, long nh);
free_cvector(char *v, long nl, long nh);
free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch);
free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);
```

**Critical**: Always match allocation/deallocation ranges.

---

## Array Access Patterns

### Iterating Over Residues

```c
// Standard pattern
for (i = 1; i <= num_residue; i++) {
    // Process residue i
    ib = seidx[i][1];  // First atom
    ie = seidx[i][2];  // Last atom
    
    // Iterate atoms in residue
    for (atom_idx = ib; atom_idx <= ie; atom_idx++) {
        // Process atom atom_idx
    }
}
```

### Iterating Over Base Pairs

```c
// For duplex (ds = 2)
for (k = 1; k <= num_bp; k++) {
    long res1 = pair_num[1][k];  // Residue on strand 1
    long res2 = pair_num[2][k];  // Residue on strand 2
    
    // Access frames
    double *orien1 = orien[res1];  // Frame for residue 1
    double *org1 = org[res1];      // Origin for residue 1
    double *orien2 = orien[res2];  // Frame for residue 2
    double *org2 = org[res2];      // Origin for residue 2
}
```

### Accessing Frame Elements

```c
// Get z-axis of residue i (1-based)
double zx = orien[i][7];  // x-component of z-axis
double zy = orien[i][8];  // y-component of z-axis
double zz = orien[i][9];  // z-component of z-axis

// Calculate dot product of z-axes
double dir_z = orien[i][7]*orien[j][7] + 
               orien[i][8]*orien[j][8] + 
               orien[i][9]*orien[j][9];
```

---

## Data Structure Examples

### Example 1: Accessing Atom in Residue

```c
// Find " C1'" atom in residue i
long ib = seidx[i][1];
long ie = seidx[i][2];
long c1_prime = 0;

for (atom_idx = ib; atom_idx <= ie; atom_idx++) {
    if (is_equal_string(AtomName[atom_idx], " C1'")) {
        c1_prime = atom_idx;
        break;
    }
}

if (c1_prime > 0) {
    // Use coordinates
    double x = xyz[c1_prime][1];
    double y = xyz[c1_prime][2];
    double z = xyz[c1_prime][3];
}
```

### Example 2: Extracting Frame Matrix

```c
// Convert orien[i] to 3×3 matrix (1-based)
double R[4][4];  // [1..3][1..3] used

R[1][1] = orien[i][1]; R[1][2] = orien[i][2]; R[1][3] = orien[i][3];
R[2][1] = orien[i][4]; R[2][2] = orien[i][5]; R[2][3] = orien[i][6];
R[3][1] = orien[i][7]; R[3][2] = orien[i][8]; R[3][3] = orien[i][9];
```

### Example 3: Iterating Base Pairs in Duplex

```c
// Process all base pairs
for (k = 1; k <= num_bp; k++) {
    long res_i = pair_num[1][k];
    long res_j = pair_num[2][k];
    
    char base_i = bp_seq[1][k];
    char base_j = bp_seq[2][k];
    
    // Calculate distance between origins
    double dx = org[res_j][1] - org[res_i][1];
    double dy = org[res_j][2] - org[res_i][2];
    double dz = org[res_j][3] - org[res_i][3];
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
}
```

---

**Next**: [Core Functions](03_CORE_FUNCTIONS.md) for function-level details


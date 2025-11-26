# Legacy Helper Functions Reference

**Date**: 2025-01-XX  
**Purpose**: Reference for utility and helper functions used throughout legacy code  
**Status**: Essential utility function reference

---

## Table of Contents

1. [Memory Management](#memory-management)
2. [Vector/Matrix Operations](#vectormatrix-operations)
3. [String Operations](#string-operations)
4. [Atom/Residue Operations](#atomresidue-operations)
5. [File I/O](#file-io)
6. [Mathematical Utilities](#mathematical-utilities)

---

## Memory Management

### Allocation Functions

**Location**: `org/src/nrutil.c`

#### Vector Allocation

```c
double *dvector(long nl, long nh);      // Allocate double vector [nl..nh]
long   *lvector(long nl, long nh);      // Allocate long vector [nl..nh]
char   *cvector(long nl, long nh);      // Allocate char vector [nl..nh]
```

**Usage**:
```c
double *vec = dvector(1, n);  // Allocates [0..n], uses [1..n]
// Access: vec[1] through vec[n]
```

#### Matrix Allocation

```c
double **dmatrix(long nrl, long nrh, long ncl, long nch);  // [nrl..nrh][ncl..nch]
long   **lmatrix(long nrl, long nrh, long ncl, long nch);
char   **cmatrix(long nrl, long nrh, long ncl, long nch);
```

**Usage**:
```c
double **mat = dmatrix(1, nrows, 1, ncols);  // Allocates [0..nrows][0..ncols], uses [1..nrows][1..ncols]
// Access: mat[i][j] where i ∈ [1, nrows], j ∈ [1, ncols]
```

### Deallocation Functions

```c
free_dvector(double *v, long nl, long nh);
free_lvector(long *v, long nl, long nh);
free_cvector(char *v, long nl, long nh);
free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch);
free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);
```

**Critical**: Must match allocation ranges exactly.

---

## Vector/Matrix Operations

### Vector Operations

```c
double dot(double *a, double *b);           // Dot product
void cross(double *a, double *b, double *c);  // Cross product
double veclen(double *a);                   // Vector length
void vec_norm(double *a);                   // Normalize vector
void ddxyz(double *a, double *b, double *c);  // c = b - a
void cpxyz(double *src, double *dst);       // Copy coordinates
double p1p2_dist(double *p1, double *p2);   // Distance between points
double vec_ang(double *va, double *vb, double *vref);  // Angle between vectors
```

**Usage Examples**:
```c
double dist = p1p2_dist(xyz[i], xyz[j]);  // Distance between atoms i and j
double len = veclen(dorg);                 // Length of origin distance vector
vec_norm(zave);                            // Normalize z-axis vector
```

### Matrix Operations

```c
void mst2orien(double *orien_vec, long ioffset, double **mst);  // Matrix to array
void orien2mst(double *orien_vec, long ioffset, double **mst);  // Array to matrix
void x_y_z_2_mtx(double *x, double *y, double *z, double **mtx);  // Vectors to matrix
void mtx_2_x_y_z(double **mtx, double *x, double *y, double *z);  // Matrix to vectors
```

**Frame Conversion**:
```c
// Convert 3×3 matrix to 9-element array
mst2orien(orien[i], 0, R);  // R[1..3][1..3] → orien[i][1..9]

// Convert 9-element array to 3×3 matrix
orien2mst(orien[i], 0, mst);  // orien[i][1..9] → mst[1..3][1..3]
```

---

## String Operations

### String Manipulation

```c
char *trim(char *a);        // Remove leading/trailing whitespace
char *ltrim(char *a);       // Remove leading whitespace
char *rtrim(char *a);       // Remove trailing whitespace
long is_equal_string(const char *str1, const char *str2);  // String comparison
long is_equal_case_string(const char *str1, const char *str2);  // Case-insensitive
long is_empty_string(const char *str);  // Check if empty
char *my_strdup(const char *src);  // String duplication
```

**Usage**:
```c
if (is_equal_string(AtomName[i], " C1'")) {
    // Process C1' atom
}
```

### String Matching

```c
long strmatch_idx(char *str, char **strmat, long nb, long ne);  // Find string in array
long num_strmatch(char *str, char **strmat, long nb, long ne);   // Count matches
```

---

## Atom/Residue Operations

### Atom Finding

```c
long find_1st_atom(char *str, char **strmat, long nb, long ne, char *idmsg);
```

**Purpose**: Find first atom matching name in residue range.

**Usage**:
```c
long c1_prime = find_1st_atom(" C1'", AtomName, ib, ie, idmsg);
if (c1_prime) {
    // Atom found at index c1_prime
}
```

### Atom Type Classification

```c
void atom_idx(long num, char **AtomName, char **ResName, long *idx);
```

**Purpose**: Classify atoms by type.

**Returns**: `idx[i]`:
- `1`: Nitrogen (N)
- `2`: Oxygen (O)
- `3`: Special atoms
- `4`: Other (e.g., Carbon)

### Residue Identification

```c
long residue_ident(char **AtomName, double **xyz, char **Miscs, long ib, long ie);
```

**Purpose**: Classify residue type.

**Returns**:
- `1`: Purine (A, G)
- `0`: Pyrimidine (C, T, U)
- `-1`: Amino acid
- `-2`: Unknown/invalid

### Base Identification

```c
char base_ident(long ib, long ie, char **AtomName, double **xyz, long RY, 
                char *ResName, char *idmsg, long num_sbase, char **baselist);
```

**Purpose**: Identify base type from atoms.

**Returns**: Base letter ('A', 'C', 'G', 'T', 'U', etc.)

---

## File I/O

### File Operations

```c
FILE *open_file(char *filename, char *mode);  // Open file with error checking
long close_file(FILE *fp);                     // Close file
long exist_file(char *filename);               // Check if file exists
void remove_file(char *filename);              // Delete file
void copy_file_pointer(FILE *fpi, FILE *fpo, char *msg);  // Copy file
```

### Line Reading

```c
char *my_getline(FILE *fp);  // Read line from file (allocates memory)
```

**Usage**:
```c
char *line;
while ((line = my_getline(fp)) != NULL) {
    // Process line
    free(line);
}
```

---

## Mathematical Utilities

### Matrix Operations

```c
void cov_matrix(double **a, double **b, long nr, long nc, double **cmtx);  // Covariance
void jacobi(double **a, long n, double *d, double **v);  // Eigenvalue decomposition
void dinverse(double **a, long n, double **y);  // Matrix inverse
void dludcmp(double **a, long n, long *indx, double *d);  // LU decomposition
void dlubksb(double **a, long n, long *indx, double *b);  // LU back substitution
```

### Rotation Matrices

```c
void rotx(double ang_deg, double **rotmat);  // Rotation around x-axis
void roty(double ang_deg, double **rotmat);  // Rotation around y-axis
void rotz(double ang_deg, double **rotmat);  // Rotation around z-axis
void arb_rotation(double *va, double ang_deg, double **rot_mtx);  // Arbitrary axis
```

### Sorting

```c
void dsort(long n, double *a, long *idx);  // Sort double array
void lsort(long n, long *a, long *idx);    // Sort long array
```

---

## Distance and Angle Calculations

### Distance Functions

```c
double p1p2_dist(double *p1, double *p2);  // Euclidean distance
long within_limits(double *p1, double *p2, double lower, double upper);  // Distance check
```

**Usage**:
```c
double dist = p1p2_dist(xyz[i], xyz[j]);
if (within_limits(xyz[i], xyz[j], 1.8, 4.0)) {
    // Distance in range [1.8, 4.0]
}
```

### Angle Functions

```c
double vec_ang(double *va, double *vb, double *vref);  // Angle between vectors
double torsion(double **d);  // Torsion angle from 4 points
double torsion2(double **d);  // Alternative torsion calculation
```

---

## Coordinate Transformations

### Alignment

```c
void align2zaxis(long num, double *haxis, double **rotmat, double **xyz, double **xyzH);
```

**Purpose**: Rotate coordinates to align `haxis` with z-axis.

**Usage**: Used in overlap calculation to project bases onto xy-plane.

### Frame Operations

```c
void ref_frame_i(long bnum, double *bp_orien, double *bp_org, double **r, double *o);
```

**Purpose**: Extract frame from flattened array format.

---

## Ring Atom Operations

### Ring Atom Extraction

```c
long ratom_xyz(long *ratom_list, long only_ring, double **xyz, double *oave, double **oxyz);
```

**Purpose**: Extract ring atom coordinates relative to `oave`.

**Parameters**:
- `ratom_list[]`: Ring atom indices (ring_atom[i])
- `only_ring`: 1 = only ring atoms, 0 = include connected atoms
- `oave[]`: Reference origin
- `xyz[][]`: Atom coordinates
- `oxyz[][]`: Output coordinates (relative to oave)

**Returns**: Number of atoms extracted

### Connected Atoms

```c
void get_cntatom(long *ringlist, long **connect, long *idx);
```

**Purpose**: Find atoms connected to ring atoms (for projection).

---

## Validation Helpers

### Range Checking

```c
long dval_in_range(double val, double min, double max);  // Check if value in range
long lval_in_set(long val, long nb, long ne, long *set);  // Check if value in set
```

### Atom Validation

```c
long is_baseatom(char *atomname);  // Check if base atom
long has_atom_name(long ib, long ie, char **AtomName, char *aname);  // Check if atom exists
```

---

## Common Usage Patterns

### Iterating Over Residue Atoms

```c
long ib = seidx[i][1];
long ie = seidx[i][2];
for (long m = ib; m <= ie; m++) {
    // Process atom m
    if (is_equal_string(AtomName[m], " C1'")) {
        // Found C1' atom
    }
}
```

### Calculating Distances

```c
double dist = p1p2_dist(xyz[i], xyz[j]);
if (within_limits(xyz[i], xyz[j], lower, upper)) {
    // Distance in valid range
}
```

### Frame Access

```c
// Get z-axis of residue i
double zx = orien[i][7];  // z-axis x-component (1-based index 7)
double zy = orien[i][8];  // z-axis y-component
double zz = orien[i][9];  // z-axis z-component

// Calculate dot product
double dir_z = orien[i][7]*orien[j][7] + 
               orien[i][8]*orien[j][8] + 
               orien[i][9]*orien[j][9];
```

---

**Related**: 
- [Data Structures](02_DATA_STRUCTURES.md) for array access patterns
- [Core Functions](03_CORE_FUNCTIONS.md) for function usage
- [Algorithms](04_ALGORITHMS.md) for mathematical operations


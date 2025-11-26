# Legacy Core Functions Reference

**Date**: 2025-01-XX  
**Purpose**: Complete reference for all critical legacy functions with signatures, algorithms, and implementation details  
**Status**: Comprehensive function reference

---

## Table of Contents

1. [Frame Calculation Functions](#frame-calculation-functions)
2. [Base Pair Validation Functions](#base-pair-validation-functions)
3. [H-Bond Detection Functions](#h-bond-detection-functions)
4. [Base Pair Finding Functions](#base-pair-finding-functions)
5. [Parameter Calculation Functions](#parameter-calculation-functions)
6. [Utility Functions](#utility-functions)

---

## Frame Calculation Functions

### 1. `base_frame()` - Base Reference Frame Calculation

**Location**: `org/src/app_fncs.c:383-459`

**Signature**:
```c
void base_frame(long num_residue, char *bseq, long **seidx, long *res_type,
                char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                char **Miscs, double **xyz, char *BDIR, double **orien, double **org)
```

**Purpose**: Calculate reference frame (rotation matrix + origin) for each nucleotide base by least-squares fitting to standard templates.

**Input Parameters**:
- `num_residue`: Total number of residues (1-based count)
- `bseq[]`: Base sequence array [1..num_residue] ('A', 'C', 'G', 'T', etc.)
- `seidx[][]`: Residue start/end atom indices [1..num_residue][1..2]
- `res_type[]`: Residue type classification [1..num_residue] (RY values)
- `AtomName[][]`, `ResName[][]`, `ChainID[]`, `ResSeq[]`, `Miscs[][]`: Atom data
- `xyz[][]`: Atom coordinates [1..num_atoms][1..3]
- `BDIR`: Base directory path for template files

**Output**:
- `orien[][]`: 9-element rotation matrix per residue [1..num_residue][1..9] (flattened 3×3)
- `org[][]`: 3-element origin vector per residue [1..num_residue][1..3]

**Algorithm**:

1. **For each residue i in [1, num_residue]**:
   - Skip if `res_type[i] < 0` (not a nucleotide)
   - Get atom range: `ib = seidx[i][1]`, `ie = seidx[i][2]`

2. **Load Standard Template**:
   - Construct filename: `{BDIR}/Atomic_{bseq[i]}.pdb` (e.g., `Atomic_A.pdb`)
   - Read template PDB file using `read_pdb()`
   - Store template atoms in temporary arrays

3. **Match Ring Atoms**:
   - Ring atom list: `RA_LIST = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "}`
   - For purines (res_type[i] == 1): Use all 9 atoms
   - For pyrimidines (res_type[i] == 0): Use first 6 atoms
   - For each ring atom:
     - Find in experimental: `find_1st_atom(ring_atom, AtomName, ib, ie)`
     - Find in standard: `find_1st_atom(ring_atom, sAtomName, 1, snum)`
     - If both found: Add to matched pairs
   - Require minimum 3 matched atoms

4. **Least-Squares Fitting**:
   - Call `ls_fitting(sRing_xyz, eRing_xyz, nmatch, fitted_xyz, R, orgi)`
   - Returns rotation matrix R (3×3) and translation vector orgi
   - RMS fit value indicates quality

5. **Store Frame**:
   - Convert rotation matrix to 9-element array: `mst2orien(orien[i], 0, R)`
   - Store origin: `org[i][1..3] = orgi[1..3]`

**Key Dependencies**:
- `ls_fitting()`: Least-squares fitting algorithm
- `read_pdb()`: PDB file reading
- `find_1st_atom()`: Find atom by name in residue
- `set_std_base_pdb()`: Construct template filename
- `mst2orien()`: Convert 3×3 matrix to 9-element array

**Critical Details**:
- Uses 1-based indexing: `i = 1..num_residue`
- Frame stored as flattened 3×3 matrix: `orien[i][1..9]` maps to `R[1..3][1..3]`
- Origin stored in 1-based array: `org[i][1..3]`
- Only processes nucleotides (skips amino acids, water, etc.)
- Template files must exist in `{BDIR}/Atomic_{base}.pdb`

**Modern Equivalent**: `BaseFrameCalculator::calculate_frame()`

---

### 2. `ls_fitting()` - Least-Squares Fitting Algorithm

**Location**: `org/src/cmn_fncs.c:1707-1766`

**Signature**:
```c
double ls_fitting(double **sxyz, double **exyz, long n, double **fitted_xyz,
                  double **R, double *orgi)
```

**Purpose**: Find optimal rotation and translation to align two sets of 3D points using quaternion-based method.

**Input Parameters**:
- `sxyz[][]`: Standard/reference coordinates [1..n][1..3] (n points × 3)
- `exyz[][]`: Experimental coordinates [1..n][1..3] (n points × 3)
- `n`: Number of matched points (must be ≥ 3)

**Output**:
- `fitted_xyz[][]`: Fitted coordinates [1..n][1..3]
- `R[][]`: Rotation matrix (3×3) [1..3][1..3]
- `orgi[]`: Translation vector [1..3]
- Returns: RMS fit value

**Algorithm** (Quaternion-based):

1. **Calculate Averages**:
   ```
   ave_sxyz[k] = (1/n) Σ sxyz[i][k] for k = 1..3
   ave_exyz[k] = (1/n) Σ exyz[i][k] for k = 1..3
   ```

2. **Calculate Covariance Matrix** (3×3):
   ```
   U[k][l] = (1/n) Σ (sxyz[i][k] - ave_sxyz[k]) × (exyz[i][l] - ave_exyz[l])
   ```

3. **Build 4×4 Quaternion Matrix N**:
   ```
   N[1][1] = U[1][1] + U[2][2] + U[3][3]
   N[2][2] = U[1][1] - U[2][2] - U[3][3]
   N[3][3] = -U[1][1] + U[2][2] - U[3][3]
   N[4][4] = -U[1][1] - U[2][2] + U[3][3]
   N[1][2] = N[2][1] = U[2][3] - U[3][2]
   N[1][3] = N[3][1] = U[3][1] - U[1][3]
   N[1][4] = N[4][1] = U[1][2] - U[2][1]
   N[2][3] = N[3][2] = U[1][2] + U[2][1]
   N[2][4] = N[4][2] = U[3][1] + U[1][3]
   N[3][4] = N[4][3] = U[2][3] + U[3][2]
   ```

4. **Eigenvalue Decomposition**:
   - Use `jacobi(N, 4, D, V)` to find eigenvectors
   - Largest eigenvalue (index 4) gives optimal quaternion
   - Eigenvector: `V[1..4][4]` = quaternion components

5. **Convert Quaternion to Rotation Matrix**:
   - Let `q = [q₀, q₁, q₂, q₃] = V[1..4][4]`
   - Build 3×3 rotation matrix R

6. **Calculate Translation**:
   ```
   orgi = ave_exyz - R × ave_sxyz
   ```

7. **Calculate RMS**:
   ```
   fitted_xyz[i] = R × sxyz[i] + orgi
   RMS = sqrt((1/n) Σ ||fitted_xyz[i] - exyz[i]||²)
   ```

**Key Dependencies**:
- `cov_matrix()`: Calculate covariance matrix
- `jacobi()`: Eigenvalue decomposition
- `dot()`: Vector dot product
- `ave_dmatrix()`: Calculate average coordinates

**Critical Details**:
- Requires minimum 3 points (fatal error if `n < 3`)
- Uses quaternion method for robustness
- 1-based indexing throughout
- Rotation matrix is proper (det(R) = 1)

**See Also**: [Algorithms: ls_fitting](04_ALGORITHMS.md#11-ls_fitting) for complete mathematical derivation

---

## Base Pair Validation Functions

### 3. `check_pair()` - Base Pair Validation

**Location**: `org/src/cmn_fncs.c:4577-4646`

**Signature**:
```c
void check_pair(long i, long j, char *bseq, long **seidx, double **xyz,
                double **NC1xyz, double **orien, double **org, long *idx,
                char **AtomName, long **ring_atom, miscPars *misc_pars,
                double *rtn_val, long *bpid, double *dir_x, double *dir_y,
                double *dir_z, long network)
```

**Purpose**: Validate if two bases can form a valid base pair by checking geometric and chemical constraints.

**Input Parameters**:
- `i, j`: Residue indices (1-based)
- `bseq[]`: Base sequence
- `seidx[][]`: Residue atom ranges
- `xyz[][]`: Atom coordinates
- `NC1xyz[][]`: N and C1' coordinates per residue
- `orien[][]`, `org[][]`: Base reference frames
- `idx[]`: Atom type indices
- `AtomName[][]`: Atom names
- `ring_atom[][]`: Ring atom indices
- `misc_pars`: Validation parameters
- `network`: Network mode flag (0 = normal, 1 = conflict resolution)

**Output**:
- `rtn_val[]`: Return values [1..5]
  - `rtn_val[1]`: Origin distance (dorg)
  - `rtn_val[2]`: Vertical distance component (dv)
  - `rtn_val[3]`: Plane angle (degrees)
  - `rtn_val[4]`: N-N distance (dNN)
  - `rtn_val[5]`: Quality score
- `bpid`: Base pair type ID (0 = invalid, >0 = valid pair type)
- `dir_x`, `dir_y`, `dir_z`: Frame direction dot products

**Validation Steps**:

1. **Calculate Base Frame Geometry**:
   ```c
   get_bp_zoave(i, j, orien, org, oave, zave);  // Average origin and z-axis
   dorg = org[j] - org[i];                       // Origin distance vector
   dNN_vec = NC1xyz[j] - NC1xyz[i];              // N-N distance vector
   ```

2. **Direction Check**:
   ```c
   dir_x = dot(orien[i][1..3], orien[j][1..3]);  // x-axis alignment
   dir_y = dot(orien[i][4..6], orien[j][4..6]);  // y-axis alignment
   dir_z = dot(orien[i][7..9], orien[j][7..9]);  // z-axis alignment (CRITICAL)
   ```

3. **Calculate Validation Metrics**:
   ```c
   rtn_val[1] = veclen(dorg);                    // Origin distance
   rtn_val[2] = fabs(dot(dorg, zave));           // Vertical distance component
   rtn_val[3] = vec_ang(orien[i][7..9], orien[j][7..9], NULL);  // Plane angle
   rtn_val[4] = veclen(dNN_vec);                 // N-N distance
   ```

4. **Geometric Constraints Check** (`cdns`):
   ```c
   cdns = (rtn_val[1] >= min_dorg && rtn_val[1] <= max_dorg) &&
          (rtn_val[2] >= min_dv && rtn_val[2] <= max_dv) &&
          (rtn_val[3] <= max_plane_angle) &&
          (rtn_val[4] >= min_dNN && rtn_val[4] <= max_dNN);
   ```

5. **Overlap Check**:
   ```c
   overlap = get_oarea(i, j, ring_atom, oave, zave, xyz, 0);
   if (overlap >= OVERLAP (0.01)) {
       // REJECT
   }
   ```

6. **Hydrogen Bond Check**:
   ```c
   num_base_hb = 0;
   for (m = seidx[i][1]; m <= seidx[i][2]; m++) {
       for (n = seidx[j][1]; n <= seidx[j][2]; n++) {
           if (good_hbatoms(...) && within_limits(...)) {
               num_base_hb++;
           }
       }
   }
   // Require: num_base_hb >= min_base_hb
   ```

7. **Base Pair Type Calculation**:
   - If all checks pass: `calculate_more_bppars()` determines `bpid`
   - Quality score: `rtn_val[5] = dorg + 2*dv + plane_angle/20 + adjust_pairQuality()`

**Key Dependencies**:
- `get_bp_zoave()`: Average base origins and z-axes
- `get_oarea()`: Calculate base overlap area
- `good_hbatoms()`: Check if atom pair can form H-bond
- `within_limits()`: Check distance constraints
- `calculate_more_bppars()`: Calculate detailed pair parameters
- `adjust_pairQuality()`: Adjust quality score based on H-bonds

**Critical Parameters** (from `misc_pars`):
- `min_dorg = 0.0`, `max_dorg = 15.0`: Origin distance range
- `min_dv = 0.0`, `max_dv = 2.5`: Vertical component range
- `min_plane_angle = 0.0`, `max_plane_angle = 65.0`: Plane angle max
- `min_dNN = 4.5`, `max_dNN = XBIG (~1e18)`: N-N distance range
- `min_base_hb = 1`: Minimum H-bonds required
- `OVERLAP = 0.01`: Overlap threshold

**Network Mode** (`network=1`):
- Used in `bases_elimination()` for conflict resolution
- Only checks geometric constraints (no H-bond requirement)
- Returns `bpid = -1` if passes constraints

**Modern Equivalent**: `BasePairValidator::validate_pair()`

---

## H-Bond Detection Functions

### 4. `get_hbond_ij()` - Hydrogen Bond Detection

**Location**: `org/src/cmn_fncs.c:4021-4083`

**Signature**:
```c
void get_hbond_ij(long i, long j, char basei, char basej, miscPars *misc_pars,
                  long **seidx, long *idx, char **AtomName, double **xyz,
                  char *hb_info)
```

**Purpose**: Find all hydrogen bonds between two residues.

**Input Parameters**:
- `i, j`: Residue indices (1-based)
- `basei, basej`: Base types (A, C, G, T)
- `misc_pars`: H-bond parameters
- `seidx[][]`: Residue atom ranges
- `idx[]`: Atom type indices
- `AtomName[][]`: Atom names
- `xyz[][]`: Atom coordinates

**Output**:
- `hb_info`: String containing H-bond count and details
- Also writes to JSON via `json_writer_record_hbond_list()`

**Algorithm**:

1. **Initial Detection**:
   ```c
   num_hbonds = 0;
   for (m = seidx[i][1]; m <= seidx[i][2]; m++) {
       for (n = seidx[j][1]; n <= seidx[j][2]; n++) {
           if (good_hbatoms(misc_pars, AtomName[m], AtomName[n], idx[m], idx[n]) &&
               within_limits(xyz[n], xyz[m], hb_lower, hb_dist1)) {
               // Record H-bond
               strcpy(hb_atom1[num_hbonds], AtomName[m]);
               strcpy(hb_atom2[num_hbonds], AtomName[n]);
               hb_dist[num_hbonds] = p1p2_dist(xyz[n], xyz[m]);
               num_hbonds++;
           }
       }
   }
   ```

2. **Conflict Resolution**:
   ```c
   hb_atompair(num_hbonds, hb_atom1, hb_atom2, hb_dist, lkg_type, misc_pars);
   ```
   - Resolves conflicts when same atom forms multiple H-bonds
   - Marks conflicts by negating distances

3. **Validation**:
   ```c
   validate_hbonds(num_hbonds, hb_dist, lkg_type, hb_type, basei, basej,
                   hb_atom1, hb_atom2);
   ```
   - Assigns H-bond types: `' '` (normal), `'-'` (conflict), `'*'` (other)

4. **Output**:
   - Format: `"[count] atom1typeatom2 dist atom1typeatom2 dist ..."`
   - Only includes validated H-bonds (type != ' ')

**Key Dependencies**:
- `good_hbatoms()`: Check if atom pair can form H-bond
- `within_limits()`: Distance check
- `hb_atompair()`: Conflict resolution
- `validate_hbonds()`: H-bond validation

**Critical Parameters**:
- `hb_lower = 1.8`: Minimum H-bond distance
- `hb_dist1 = 4.0`: Maximum H-bond distance
- `hb_dist2 = 0.0`: Conflict resolution threshold (CRITICAL: must be 0.0)

**Modern Equivalent**: `BasePairValidator::find_hydrogen_bonds_detailed()`

---

### 5. `hb_atompair()` - H-Bond Conflict Resolution

**Location**: `org/src/cmn_fncs.c:3923-3987`

**Signature**:
```c
void hb_atompair(long num_hbonds, char **hb_atom1, char **hb_atom2,
                 double *hb_dist, long *lkg_type, miscPars *misc_pars)
```

**Purpose**: Resolve conflicts when atoms form multiple H-bonds (iterative algorithm).

**Input Parameters**:
- `num_hbonds`: Number of initial H-bonds
- `hb_atom1[][]`, `hb_atom2[][]`: Atom names for each H-bond
- `hb_dist[]`: H-bond distances (modified in place)
- `lkg_type[]`: Output linkage types
- `misc_pars`: Parameters

**Output**:
- Modifies `hb_dist[]`: Negates distances for conflicted H-bonds
- Sets `lkg_type[]`: Linkage type for each H-bond

**Algorithm** (Two-Phase):

**Phase 1: Conflict Detection and Marking**:
- Iterative algorithm until all H-bonds processed
- For each H-bond, find best H-bond for atom1 and atom2 (by distance)
- If same H-bond is best for both atoms: mark as conflict (negate distance)
- Continue until convergence

**Phase 2: Linkage Type Calculation**:
- For conflicted H-bonds: mark atoms as conflicted
- For non-conflicted H-bonds sharing atoms: mark appropriately
- Calculate linkage type: `lkg_type[k] = idx2[k][1] + idx2[k][2]`
  - `lkg_type = 18`: No conflicts (both atoms unique)
  - `lkg_type < 18`: Conflict detected

**Critical Details**:
- Conflicts marked by negating distances (negative = conflict)
- Iterative algorithm until convergence
- Linkage types: 18 = unique, <18 = conflicted
- `hb_dist2 = 0.0` means Phase 2 conflict marking is effectively disabled

**See Also**: [Algorithms: hb_atompair](04_ALGORITHMS.md#5-hb_atompair) for detailed algorithm

---

### 6. `validate_hbonds()` - H-Bond Validation

**Location**: `org/src/cmn_fncs.c:3989-4019`

**Signature**:
```c
long validate_hbonds(long num_hbonds, double *hb_dist, long *lkg_type,
                     char *hb_type, char basei, char basej,
                     char **hb_atom1, char **hb_atom2)
```

**Purpose**: Assign types to H-bonds and filter invalid ones.

**Input Parameters**:
- `num_hbonds`: Number of H-bonds (after conflict resolution)
- `hb_dist[]`: H-bond distances (may be negative for conflicts)
- `lkg_type[]`: Linkage types from `hb_atompair()`
- `hb_type[]`: Output types (modified in place)
- `basei, basej`: Base types
- `hb_atom1[][]`, `hb_atom2[][]`: Atom names

**Output**:
- Returns count of valid H-bonds
- Sets `hb_type[]`: `' '` (normal), `'-'` (conflict), `'*'` (other)

**Algorithm**:

1. **Initialize All Types to ' '**
2. **Assign Types for Conflicted H-bonds**:
   - Convert negative distances to positive: `hb_dist[k] = fabs(hb_dist[k])`
   - Call `donor_acceptor()` to determine type
3. **Filter Invalid H-bonds**:
   - Remove if distance > 3.6
   - Remove if wrong type and distance not in [2.6, 3.2]
4. **Count Valid H-bonds**: `count(hb_type != ' ')`

**H-Bond Types**:
- `' '` (space): Normal H-bond (or invalid/filtered)
- `'-'` (dash): Conflict H-bond (wrong donor/acceptor roles)
- `'*'` (asterisk): Other type of H-bond

**Modern Equivalent**: `BasePairValidator::validate_hbonds()`

---

## Base Pair Finding Functions

### 7. `find_bestpair()` - Greedy Base Pair Matching

**Location**: `org/src/find_pair.c:1576-1612`

**Signature**:
```c
static long find_bestpair(long nout, long **base_pairs, long num_residue,
                          char *bseq, long **seidx, long *RY, char **AtomName,
                          double **xyz, long *idx, double **orien, double **org,
                          double **NC1xyz, long **ring_atom, miscPars *misc_pars)
```

**Purpose**: Find best base pairs using greedy matching algorithm (mutual best partner).

**Input Parameters**:
- `nout`: Number of output fields
- `base_pairs[][]`: Output array for pairs
- `num_residue`: Total residues
- `bseq[]`, `seidx[][]`, `xyz[][]`, etc.: Structure data
- `misc_pars`: Validation parameters

**Output**:
- Returns number of pairs found
- Fills `base_pairs[][]` with pair data

**Algorithm** (Iterative Greedy Matching):

```c
matched_idx[] = zeros(num_residue);  // Track matched bases
num_bp = 0;

while (matched_count increases) {
    old_count = matched_count;
    
    for (i = 1; i <= num_residue; i++) {
        if (RY[i] < 0 || matched_idx[i]) continue;
        
        // Find best partner for base i
        best_pair(i, ..., pair_istat);  // Returns j, bpid, quality
        
        if (pair_istat[1] != 0) {  // Found a partner j
            // Check if j's best partner is also i (mutual)
            best_pair(pair_istat[1], ..., pair_jstat);
            
            if (i == pair_jstat[1]) {  // Mutual match!
                matched_idx[i] = 1;
                matched_idx[pair_istat[1]] = 1;
                base_pairs[++num_bp] = (i, j, bpid, ...);
            }
        }
    }
    
    new_count = count(matched_idx);
    if (new_count == old_count) break;  // No more pairs found
}
```

**Key Dependencies**:
- `best_pair()`: Find best partner for a single base
- `check_pair()`: Validate potential pairs

**Critical Details**:
- Requires mutual matching: i's best is j AND j's best is i
- Iterative: continues until no new pairs found
- Greedy: always picks best available partner
- Processes all residues before starting next iteration

**Modern Equivalent**: `BasePairFinder::find_pairs()`

---

### 8. `best_pair()` - Find Best Partner for One Base

**Location**: `org/src/find_pair.c:748-769`

**Signature**:
```c
static void best_pair(long i, long num_residue, long *RY, long **seidx,
                      double **xyz, long *idx, double **NC1xyz, long *matched_idx,
                      double **orien, double **org, long **ring_atom,
                      char **AtomName, char *bseq, miscPars *misc_pars,
                      long *pair_stat)
```

**Purpose**: Find the best (lowest quality score) valid partner for base i.

**Input Parameters**:
- `i`: Base index to find partner for
- `num_residue`, `RY[]`, etc.: Structure data
- `matched_idx[]`: Bases already paired

**Output**:
- `pair_stat[]`: Best partner data
  - `pair_stat[1]`: Partner residue index j (0 if none)
  - `pair_stat[2]`: Base pair type ID (bpid)
  - `pair_stat[3..]`: Quality scores and parameters

**Algorithm**:

```c
best_score = XBIG;  // Very large number

for (j = 1; j <= num_residue; j++) {
    if (j == i) continue;
    if (RY[j] < 0) continue;           // Skip non-nucleotides
    if (matched_idx[j]) continue;       // Skip already matched
    
    check_pair(i, j, ..., rtn_val, &bpid, ..., 0);  // network=0
    
    if (bpid != 0 && rtn_val[5] < best_score) {
        best_score = rtn_val[5];
        pair_stat[1] = j;              // Best partner
        pair_stat[2] = bpid;           // Pair type
        pair_stat[3..] = rtn_val[...];  // Quality scores
    }
}
```

**Key Dependencies**:
- `check_pair()`: Validates and scores pair
- `rtn_val[5]`: Quality score (lower = better)

**Critical Details**:
- Loops through ALL other residues
- Only considers nucleotides (RY >= 0)
- Only considers unmatched residues
- Quality score: `rtn_val[5] = dorg + 2*dv + plane_angle/20 + adjust_pairQuality()`
- Lower score = better pair

**Modern Equivalent**: `BasePairFinder::find_best_partner()`

---

## Parameter Calculation Functions

### 9. `bpstep_par()` - Base Pair Step Parameters

**Location**: `org/src/ana_fncs.c:3512-3557`

**Signature**:
```c
void bpstep_par(double **rot1, double *org1, double **rot2, double *org2,
                double *pars, double *mst_orien, double *mst_org)
```

**Purpose**: Calculate 6 base pair step parameters (Shift, Slide, Rise, Tilt, Roll, Twist).

**Input Parameters**:
- `rot1[][]`, `rot2[][]`: Rotation matrices (3×3) for two consecutive base pairs
- `org1[]`, `org2[]`: Origin vectors (3×1) for two consecutive base pairs

**Output**:
- `pars[]`: Step parameters [1..6]
  - `pars[1]`: Shift
  - `pars[2]`: Slide
  - `pars[3]`: Rise
  - `pars[4]`: Tilt
  - `pars[5]`: Roll
  - `pars[6]`: Twist
- `mst_orien[]`: Midstep frame orientation (9-element)
- `mst_org[]`: Midstep frame origin (3-element)

**Algorithm**:

1. **Extract y-axes** from both frames
2. **Calculate hinge vector** (perpendicular to both y-axes)
3. **Create parallel frames** by rotating each frame by half the angle
4. **Extract midst step z-axis** from parallel frame 2
5. **Calculate twist** (angle between y-axes around z-axis)
6. **Build midst step frame**
7. **Calculate translations** (Shift, Slide, Rise) by projecting origin difference
8. **Calculate rotations** (Tilt, Roll) from hinge vector

**See Also**: [Algorithms: bpstep_par](04_ALGORITHMS.md#12-bpstep_par) for complete mathematical derivation

---

## Utility Functions

### 10. `get_bp_zoave()` - Average Z-Axis Calculation

**Location**: `org/src/cmn_fncs.c`

**Signature**:
```c
void get_bp_zoave(long ia, long ib, double **orien, double **org,
                  double *oave, double *zave)
```

**Purpose**: Calculate average origin and z-axis for base pair validation.

**Algorithm**:
```c
// Average origin
oave[1] = (org[ia][1] + org[ib][1]) / 2;
oave[2] = (org[ia][2] + org[ib][2]) / 2;
oave[3] = (org[ia][3] + org[ib][3]) / 2;

// Average z-axis (normalized)
zave[1] = (orien[ia][7] + orien[ib][7]) / 2;
zave[2] = (orien[ia][8] + orien[ib][8]) / 2;
zave[3] = (orien[ia][9] + orien[ib][9]) / 2;
vec_norm(zave);  // Normalize
```

**Usage in `check_pair()`**:
- `rtn_val[2] = fabs(dot(dorg, zave))`: Vertical component of origin distance
- Used to check `min_dv <= rtn_val[2] <= max_dv`

**Critical Details**:
- z-axis is stored at indices `[7, 8, 9]` in `orien[][]` (1-based)
- z-axis points perpendicular to base plane (upward)
- For paired bases, z-axes should point toward each other (`dir_z < 0`)

---

**Next**: [Algorithms](04_ALGORITHMS.md) for detailed mathematical derivations


# Legacy Code Profiling and Function Mapping

**Date**: 2025-01-XX  
**Purpose**: Comprehensive profiling and mapping of legacy X3DNA code to understand critical functions and their relationships  
**Status**: Complete analysis of legacy codebase structure

---

## Table of Contents

1. [Program Architecture Overview](#program-architecture-overview)
2. [Critical Function Profiles](#critical-function-profiles)
3. [Data Structures and Indexing](#data-structures-and-indexing)
4. [Algorithm Workflows](#algorithm-workflows)
5. [Function Dependency Map](#function-dependency-map)
6. [Indexing Conventions](#indexing-conventions)
7. [Key Parameters and Constants](#key-parameters-and-constants)
8. [Comparison Points with Modern Code](#comparison-points-with-modern-code)

---

## Program Architecture Overview

### Main Entry Points

The legacy code has two main programs:

#### 1. `find_pair` (`org/src/find_pair.c`)
**Purpose**: Identifies base pairs in nucleic acid structures

**Main Flow**:
```
find_pair_main() 
  → handle_str()
    → read_pdb()              # Parse PDB file
    → residue_idx()           # Map atoms to residues
    → get_seq()               # Extract nucleotide sequence
    → base_info()             # Calculate base frames
      → base_frame()          # Core frame calculation
    → find_bestpair()         # Greedy matching algorithm
      → best_pair()           # Find best partner for each base
        → check_pair()        # Validate potential pairs
          → get_hbond_ij()    # Detect H-bonds
    → re_ordering()           # Reorder pairs by 5'→3'
    → x3dna_input()           # Write output file
```

#### 2. `analyze` (`org/src/analyze.c`)
**Purpose**: Calculates base pair and step parameters

**Main Flow**:
```
analyze_main()
  → read_input()              # Read .inp file from find_pair
  → read_pdb()                # Re-read PDB file
  → ref_frames()              # Recalculate base pair frames
    → base_frame()            # Calculate frames for each base
  → get_parameters()          # Calculate step parameters
    → bpstep_par()            # Calculate 6 base pair step parameters
    → helical_par()           # Calculate helical parameters
```

### File Organization

```
org/src/
├── find_pair.c          # Main find_pair program
├── analyze.c            # Main analyze program
├── find_pair_analyze.c  # Combined program
├── cmn_fncs.c          # Common functions (largest file, ~5000 lines)
├── ana_fncs.c          # Analyze-specific functions
├── app_fncs.c          # Application functions (frame calculation, etc.)
├── nrutil.c            # Numerical recipes utilities
└── json_writer.c       # JSON output (added for debugging)

org/include/
├── x3dna.h             # Main header with definitions
├── x3dna_fncs.h        # Function declarations
└── json_writer.h       # JSON writer interface
```

---

## Critical Function Profiles

### 1. `base_frame()` - Base Reference Frame Calculation

**Location**: `org/src/app_fncs.c:383-459`

**Purpose**: Calculate reference frame (rotation matrix + origin) for each nucleotide base by least-squares fitting to standard templates.

**Input Parameters**:
- `num_residue`: Total number of residues
- `bseq[]`: Base sequence (A, C, G, T, etc.)
- `seidx[][]`: Residue start/end atom indices
- `res_type[]`: Residue type classification
- `AtomName[][]`, `ResName[][]`, etc.: Atom data
- `xyz[][]`: Atom coordinates
- `BDIR`: Base directory for template files

**Output**:
- `orien[][]`: 9-element rotation matrix per residue (flattened 3×3)
- `org[][]`: 3-element origin vector per residue

**Algorithm**:
1. **Load Standard Template**: 
   - Construct filename: `{BDIR}/Atomic_{base}.pdb` (e.g., `Atomic_A.pdb`)
   - Read template PDB file using `read_pdb()`

2. **Match Ring Atoms**:
   - Iterate through standard ring atoms: `RA_LIST = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "}`
   - For each ring atom, find matching atom in experimental structure
   - Requires minimum 3 matched atoms

3. **Least-Squares Fitting**:
   - Call `ls_fitting()` with matched atom coordinates
   - Returns rotation matrix R and translation vector t (stored in `org`)
   - RMS fit value indicates quality

4. **Store Frame**:
   - Convert rotation matrix to 9-element array: `mst2orien(orien[i], 0, R)`
   - Store origin: `org[i] = translation vector`

**Key Dependencies**:
- `ls_fitting()`: Least-squares fitting algorithm
- `read_pdb()`: PDB file reading
- `find_1st_atom()`: Find atom by name in residue
- `set_std_base_pdb()`: Construct template filename

**Critical Details**:
- Uses 1-based indexing: `i = 1..num_residue`
- Frame stored as flattened 3×3 matrix: `orien[i][0..8]` maps to `R[1..3][1..3]`
- Origin stored in 1-based array: `org[i][1..3]`
- Only processes nucleotides (skips amino acids, water, etc.)

**Modern Equivalent**: `BaseFrameCalculator::calculate_frame()` in modern code

---

### 2. `ls_fitting()` - Least-Squares Fitting Algorithm

**Location**: `org/src/cmn_fncs.c:1707-1766`

**Purpose**: Find optimal rotation and translation to align two sets of 3D points using quaternion-based method.

**Input Parameters**:
- `sxyz[][]`: Standard/reference coordinates (n points × 3)
- `exyz[][]`: Experimental coordinates (n points × 3)
- `n`: Number of matched points
- `fitted_xyz[][]`: Output fitted coordinates
- `R[][]`: Output rotation matrix (3×3)
- `orgi[]`: Output translation vector (1×3)

**Output**:
- Returns RMS fit value
- Modifies `R` and `orgi` in place

**Algorithm** (Quaternion-based):

1. **Calculate Covariance Matrix**:
   ```
   U = cov_matrix(sxyz, exyz, n, 3)
   ```

2. **Build 4×4 Quaternion Matrix N**:
   - Construct from covariance matrix elements
   - Symmetric matrix representing rotation

3. **Eigenvalue Decomposition**:
   - Use `jacobi()` to find eigenvectors of N
   - Largest eigenvalue (index 4) gives optimal quaternion

4. **Convert Quaternion to Rotation Matrix**:
   - Convert quaternion to 3×3 rotation matrix R
   - Formula: `R[1..3][1..3] = f(V[1..4][4])` where V is eigenvector matrix

5. **Calculate Translation**:
   ```
   orgi = ave_exyz - R × ave_sxyz
   ```

6. **Apply Transformation**:
   ```
   fitted_xyz[i] = R × sxyz[i] + orgi
   ```

7. **Calculate RMS**:
   ```
   RMS = sqrt(Σ ||fitted_xyz[i] - exyz[i]||² / n)
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

**Modern Equivalent**: Should match exactly with modern implementation

---

### 3. `check_pair()` - Base Pair Validation

**Location**: `org/src/cmn_fncs.c:4577-4646`

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
- `misc_pars`: Validation parameters

**Output**:
- `rtn_val[]`: Return values (distance, angles, quality score)
- `bpid`: Base pair type ID (0 = invalid, >0 = valid pair type)
- Returns via output parameters

**Validation Steps**:

1. **Calculate Base Frame Geometry**:
   ```
   get_bp_zoave(i, j, orien, org, oave, zave)  # Average origin and z-axis
   dorg = org[j] - org[i]                      # Origin distance vector
   dNN_vec = NC1xyz[j] - NC1xyz[i]             # N-N distance vector
   ```

2. **Direction Check**:
   ```
   dir_x = dot(orien[i][0..2], orien[j][0..2])  # x-axis alignment
   dir_y = dot(orien[i][3..5], orien[j][3..5])  # y-axis alignment
   dir_z = dot(orien[i][6..8], orien[j][6..8])  # z-axis alignment (CRITICAL)
   ```

3. **Calculate Validation Metrics**:
   ```
   rtn_val[1] = |dorg|                    # Origin distance
   rtn_val[2] = |dot(dorg, zave)|         # Vertical distance component
   rtn_val[3] = angle(z1, z2)             # Plane angle (0-90°)
   rtn_val[4] = |dNN_vec|                 # N-N distance
   rtn_val[5] = quality_score             # Pair quality score
   ```

4. **Geometric Constraints Check** (`cdns`):
   ```
   cdns = (rtn_val[1] in [min_dorg, max_dorg]) &&
          (rtn_val[2] in [min_dv, max_dv]) &&
          (rtn_val[3] <= max_plane_angle) &&
          (rtn_val[4] in [min_dNN, max_dNN])
   ```

5. **Overlap Check**:
   ```
   overlap = get_oarea(i, j, ring_atom, oave, zave, xyz, 0)
   if overlap >= OVERLAP (0.01): REJECT
   ```

6. **Hydrogen Bond Check**:
   - Loop through all atom pairs between residues:
     ```
     for m in [seidx[i][1], seidx[i][2]]:
         for n in [seidx[j][1], seidx[j][2]]:
             if good_hbatoms() && within_limits():
                 num_base_hb++
     ```
   - Require: `num_base_hb >= min_base_hb` (typically 1-2)

7. **Base Pair Type Calculation**:
   - If all checks pass: `calculate_more_bppars()` determines `bpid`
   - `bpid` encodes base pair type (Watson-Crick, wobble, etc.)

**Key Dependencies**:
- `get_bp_zoave()`: Average base origins and z-axes
- `get_oarea()`: Calculate base overlap area
- `good_hbatoms()`: Check if atom pair can form H-bond
- `within_limits()`: Check distance constraints
- `calculate_more_bppars()`: Calculate detailed pair parameters

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

**Modern Equivalent**: `BasePairValidator::validate_pair()` should match logic exactly

---

### 4. `get_hbond_ij()` - Hydrogen Bond Detection

**Location**: `org/src/cmn_fncs.c:4021-4083`

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

1. **Initial Detection** (lines 4038-4048):
   ```
   for m in [seidx[i][1], seidx[i][2]]:
       for n in [seidx[j][1], seidx[j][2]]:
           if good_hbatoms(misc_pars, AtomName[m], AtomName[n], idx[m], idx[n]) &&
              within_limits(xyz[n], xyz[m], hb_lower, hb_dist1):
               # Record H-bond
               hb_atom1[num_hbonds] = AtomName[m]
               hb_atom2[num_hbonds] = AtomName[n]
               hb_dist[num_hbonds] = p1p2_dist(xyz[n], xyz[m])
               num_hbonds++
   ```

2. **Conflict Resolution**:
   ```
   hb_atompair(num_hbonds, hb_atom1, hb_atom2, hb_dist, lkg_type, misc_pars)
   ```
   - Resolves conflicts when same atom forms multiple H-bonds
   - Marks conflicts by negating distances
   - Calculates linkage types

3. **Validation**:
   ```
   validate_hbonds(num_hbonds, hb_dist, lkg_type, hb_type, basei, basej, 
                   hb_atom1, hb_atom2)
   ```
   - Assigns H-bond types: `' '` (normal), `'-'` (conflict), `'*'` (other)
   - Filters invalid H-bonds

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

**Phase 1: Conflict Detection and Marking** (lines 3932-3963):

1. Initialize `matched_idx[]` to track processed H-bonds

2. Iterate until all H-bonds processed:
   ```
   while True:
       # Find best H-bond from current atom (by distance)
       update_hb_idx(1, dtmp, ddidx, hb_dist, num_iter)  # Best for atom1
       update_hb_idx(2, dtmp, ddidx, hb_dist, num_iter)  # Best for atom2
       
       # Check all other H-bonds from same atoms
       for n in all_hbonds:
           if same atom1 and closer: update best for atom1
           if same atom2 and closer: update best for atom2
       
       # Conflict detected if same H-bond is best for both atoms
       if ddidx[1] == ddidx[2]:
           hb_dist[k] = -hb_dist[k]  # Mark conflict (negate distance)
           Mark all H-bonds sharing atoms with k as "matched"
           Reset iteration
       
       num_iter++
   ```

3. Mark conflicted H-bonds by negating distances

**Phase 2: Linkage Type Calculation** (lines 3964-3983):

1. Initialize `idx2[][]` to track atom usage

2. For each conflicted H-bond (negative distance):
   ```
   idx2[conflicted][1] = 9  # Mark atom1 as conflicted
   idx2[conflicted][2] = 9  # Mark atom2 as conflicted
   
   # Mark non-conflicted H-bonds that share atoms
   for each non-conflicted H-bond:
       if shares atom1: idx2[non_conflicted][1] = 1
       if shares atom2: idx2[non_conflicted][2] = 1
   ```

3. Calculate linkage type:
   ```
   lkg_type[k] = idx2[k][1] + idx2[k][2]
   ```
   - `lkg_type = 18`: No conflicts (both atoms unique)
   - `lkg_type < 18`: Conflict detected

4. Additional conflict marking (line 3982):
   ```
   if lkg_type != 18 && distance in [hb_lower, hb_dist2]:
       hb_dist[k] = -hb_dist[k]  # Mark as conflict
   ```
   - Note: Since `hb_dist2 = 0.0`, this condition is never true in practice

**Key Dependencies**:
- `update_hb_idx()`: Update best H-bond index for an atom
- Uses array indexing: 1-based throughout

**Critical Details**:
- Conflicts marked by negating distances (negative = conflict)
- Iterative algorithm until convergence
- Linkage types: 18 = unique, <18 = conflicted
- `hb_dist2 = 0.0` means Phase 2 conflict marking is effectively disabled

**Modern Equivalent**: `BasePairValidator::resolve_conflicts()` should match exactly

---

### 6. `validate_hbonds()` - H-Bond Validation

**Location**: `org/src/cmn_fncs.c:3989-4019`

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

1. **Initialize All Types to ' '**:
   ```
   for k in [1, num_hbonds]:
       hb_type[k] = ' '
       if hb_dist[k] > 0.0: continue  # Skip conflicted (negative)
   ```

2. **Assign Types for Conflicted H-bonds**:
   ```
   for k in [1, num_hbonds]:
       if hb_dist[k] <= 0.0:  # Conflicted
           hb_dist[k] = fabs(hb_dist[k])  # Make positive
           hb_type[k] = donor_acceptor(basei, basej, hb_atom1[k], hb_atom2[k])
   ```
   - `donor_acceptor()` returns: `' '`, `'-'`, or `'*'`

3. **Filter Invalid H-bonds** (if special conflicts found):
   ```
   if count(type == '-') > 0:
       for k in [1, num_hbonds]:
           if hb_type[k] == ' ': continue
           
           # Remove if distance too large
           if hb_dist[k] > 3.6:
               hb_type[k] = ' '
           
           # Remove if wrong type and distance
           if hb_type[k] == '*' && lkg_type != 18 && 
              distance not in [2.6, 3.2]:
               hb_type[k] = ' '
   ```

4. **Count Valid H-bonds**:
   ```
   valid_count = count(hb_type != ' ')
   ```

**Key Dependencies**:
- `donor_acceptor()`: Determine H-bond type based on atom roles
- Uses `fabs()` to convert negative distances to positive

**H-Bond Types**:
- `' '` (space): Normal H-bond (or invalid/filtered)
- `'-'` (dash): Conflict H-bond (wrong donor/acceptor roles)
- `'*'` (asterisk): Other type of H-bond

**Modern Equivalent**: `BasePairValidator::validate_hbonds()` should match logic

---

### 7. `find_bestpair()` - Greedy Base Pair Matching

**Location**: `org/src/find_pair.c:1576-1612`

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

```
matched_idx[] = zeros(num_residue)  # Track matched bases
num_bp = 0

while matched_count increases:
    old_count = matched_count
    
    for i in [1, num_residue]:
        if RY[i] < 0: continue        # Skip non-nucleotides
        if matched_idx[i]: continue    # Skip already matched
        
        # Find best partner for base i
        best_pair(i, ..., pair_istat)  # Returns j, bpid, quality
        
        if pair_istat[1] != 0:  # Found a partner j
            # Check if j's best partner is also i (mutual)
            best_pair(pair_istat[1], ..., pair_jstat)
            
            if i == pair_jstat[1]:  # Mutual match!
                matched_idx[i] = 1
                matched_idx[pair_istat[1]] = 1
                base_pairs[++num_bp] = (i, j, bpid, ...)
    
    new_count = count(matched_idx)
    
    if new_count == old_count: break  # No more pairs found
```

**Key Dependencies**:
- `best_pair()`: Find best partner for a single base
- `check_pair()`: Validate potential pairs

**Critical Details**:
- Requires mutual matching: i's best is j AND j's best is i
- Iterative: continues until no new pairs found
- Greedy: always picks best available partner
- Processes all residues before starting next iteration

**Modern Equivalent**: `BasePairFinder::find_pairs()` uses same algorithm

---

### 8. `best_pair()` - Find Best Partner for One Base

**Location**: `org/src/find_pair.c:748-769`

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

```
best_score = XBIG (very large number)

for j in [1, num_residue]:
    if j == i: continue
    if RY[j] < 0: continue           # Skip non-nucleotides
    if matched_idx[j]: continue       # Skip already matched
    
    check_pair(i, j, ..., rtn_val, &bpid, ..., 0)  # network=0
    
    if bpid != 0 && rtn_val[5] < best_score:
        best_score = rtn_val[5]
        pair_stat[1] = j              # Best partner
        pair_stat[2] = bpid           # Pair type
        pair_stat[3..] = rtn_val[...] # Quality scores
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

**Modern Equivalent**: `BasePairFinder::find_best_partner()` should match

---

### 9. `get_zoave()` / `get_bp_zoave()` - Average Z-Axis Calculation

**Location**: `org/src/cmn_fncs.c` (function declarations in x3dna_fncs.h)

**Purpose**: Calculate average origin and z-axis for base pair validation.

**Functions**:
- `get_zoave()`: For base pair steps
- `get_bp_zoave()`: For single base pair

**Algorithm** (`get_bp_zoave()`):
```
# Average origin
oave = (org[i] + org[j]) / 2

# Average z-axis (normalized)
z1 = orien[i][6..8]  # z-axis of base i
z2 = orien[j][6..8]  # z-axis of base j
zave = normalize((z1 + z2) / 2)
```

**Usage in `check_pair()`**:
- `rtn_val[2] = fabs(dot(dorg, zave))`: Vertical component of origin distance
- Used to check `min_dv <= rtn_val[2] <= max_dv`

**Critical Details**:
- z-axis is stored at indices `[6, 7, 8]` in `orien[][]` (0-indexed: `[6..8]`)
- z-axis points perpendicular to base plane (upward)
- For paired bases, z-axes should point toward each other (`dir_z < 0`)

---

### 10. `ls_fitting()` - Least-Squares Fitting (Detailed)

**Location**: `org/src/cmn_fncs.c:1707-1766`

**Mathematical Details**:

**Input**:
- Two sets of n points: `sxyz[][]` (standard), `exyz[][]` (experimental)

**Goal**: Find rotation R and translation t to minimize:
```
Σ ||R × sxyz[i] + t - exyz[i]||²
```

**Method**: Quaternion-based least-squares fitting

**Steps**:

1. **Calculate Covariance Matrix** (3×3):
   ```
   U[k][l] = Σ(sxyz[i][k] - ave_sxyz[k]) × (exyz[i][l] - ave_exyz[l]) / n
   ```

2. **Build 4×4 Quaternion Matrix N**:
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

3. **Eigenvalue Decomposition**:
   ```
   jacobi(N, 4, D, V)  # V contains eigenvectors
   # Largest eigenvalue is at index 4
   # Corresponding eigenvector is V[1..4][4]
   ```

4. **Convert Quaternion to Rotation Matrix**:
   - Quaternion components: `q = V[1..4][4]`
   - Rotation matrix:
     ```
     R[1][1] = q1² + q2² - q3² - q4²
     R[1][2] = 2(q2×q3 - q1×q4)
     R[1][3] = 2(q2×q4 + q1×q3)
     R[2][1] = 2(q3×q2 + q1×q4)
     R[2][2] = q1² - q2² + q3² - q4²
     R[2][3] = 2(q3×q4 - q1×q2)
     R[3][1] = 2(q4×q2 - q1×q3)
     R[3][2] = 2(q4×q3 + q1×q2)
     R[3][3] = q1² - q2² - q3² + q4²
     ```

5. **Calculate Translation**:
   ```
   ave_sxyz = average of sxyz
   ave_exyz = average of exyz
   t = ave_exyz - R × ave_sxyz
   ```

6. **Calculate RMS**:
   ```
   fitted_xyz[i] = R × sxyz[i] + t
   RMS = sqrt(Σ ||fitted_xyz[i] - exyz[i]||² / n)
   ```

**Key Properties**:
- Optimal (minimizes RMS error)
- Handles noise robustly
- Requires minimum 3 points
- Returns proper rotation (det(R) = 1)

---

## Data Structures and Indexing

### Core Data Structures

#### 1. Atom Arrays (1-based indexing)
```c
char **AtomName;     // [1..num_atoms][0..4] - Atom names (" C1'", " N1 ", etc.)
char **ResName;      // [1..num_atoms][0..3] - Residue names ("  A", "  C", etc.)
char *ChainID;       // [1..num_atoms] - Chain identifier
long *ResSeq;        // [1..num_atoms] - Residue sequence number
double **xyz;        // [1..num_atoms][1..3] - x, y, z coordinates
char **Miscs;        // [1..num_atoms][0..NMISC] - Misc fields (insertion code, etc.)
```

#### 2. Residue Indexing
```c
long **seidx;        // [1..num_residue][1..2] - [start_atom_idx, end_atom_idx]
                     // Example: seidx[5][1] = 42, seidx[5][2] = 50
                     //          means residue 5 has atoms 42-50 (inclusive)

long *RY;            // [1..num_residue] - Purine/Pyrimidine classification
                     // RY[i] = 1: Purine (A, G)
                     // RY[i] = 0: Pyrimidine (C, T, U)
                     // RY[i] = -1: Not a nucleotide

char *bseq;          // [1..num_residue] - Base sequence ('A', 'C', 'G', 'T', etc.)
```

#### 3. Base Reference Frames
```c
double **orien;      // [1..num_residue][1..9] - Rotation matrix (flattened 3×3)
                     // orien[i][0..8] maps to R[1..3][1..3]:
                     //   R[1][1] = orien[i][0], R[1][2] = orien[i][1], R[1][3] = orien[i][2]
                     //   R[2][1] = orien[i][3], R[2][2] = orien[i][4], R[2][3] = orien[i][5]
                     //   R[3][1] = orien[i][6], R[3][2] = orien[i][7], R[3][3] = orien[i][8]
                     // z-axis: orien[i][6..8] (indices 6, 7, 8)

double **org;        // [1..num_residue][1..3] - Origin (translation vector)
                     // org[i][1] = x, org[i][2] = y, org[i][3] = z
```

#### 4. Base Pair Data
```c
long **pair_num;     // [1..ds][1..num_bp] - Residue indices for each pair
                     // pair_num[1][k] = residue i (first strand)
                     // pair_num[2][k] = residue j (second strand)
                     // ds = 1: single strand, ds = 2: duplex

char **bp_seq;       // [1..ds][1..num_bp] - Base letters for each pair
long *WC_info;       // [1..num_bp] - Watson-Crick pairing info
```

#### 5. Ring Atom Indices
```c
long **ring_atom;    // [1..num_residue][1..19] - Indices of ring atoms
                     // Used for overlap calculations
```

#### 6. NC1 Coordinates
```c
double **NC1xyz;     // [1..num_residue][1..7] - Special coordinates
                     // NC1xyz[i][0..2]: N atom (N1 or N9)
                     // NC1xyz[i][3..5]: C1' atom
                     // NC1xyz[i][6]: Flag (1.0 if C1' exists)
```

### Memory Allocation Pattern

**1-based indexing with NR_END offset**:
```c
#define NR_END 1

// Allocate array for indices [1..n]
double *vec = dvector(1, n);  // Actually allocates [0..n], but [0] unused

// Access: vec[i] where i ∈ [1, n]
// Actual memory: vec[0..n], but indices [0] unused
```

**2D Arrays**:
```c
double **mat = dmatrix(1, nrows, 1, ncols);
// Allocates: mat[0..nrows][0..ncols]
// Access: mat[i][j] where i ∈ [1, nrows], j ∈ [1, ncols]
// Unused: mat[0][*], mat[*][0]
```

---

## Algorithm Workflows

### Workflow 1: Base Frame Calculation

```
base_frame(num_residue, bseq, seidx, ..., orien, org)
├─ For each residue i in [1, num_residue]:
│  ├─ if RY[i] < 0: continue  # Skip non-nucleotides
│  ├─ ib = seidx[i][1]        # Start atom index
│  ├─ ie = seidx[i][2]        # End atom index
│  ├─ Load template: Atomic_{bseq[i]}.pdb
│  ├─ Match ring atoms:
│  │  ├─ For each ring atom in RA_LIST:
│  │  │  ├─ Find in experimental: find_1st_atom(ring_atom, AtomName, ib, ie)
│  │  │  └─ Find in standard: find_1st_atom(ring_atom, sAtomName, 1, snum)
│  │  └─ Collect matched pairs
│  ├─ if nmatch < 3: skip  # Need minimum 3 atoms
│  ├─ ls_fitting(sRing_xyz, eRing_xyz, nmatch, fitted_xyz, R, orgi)
│  │  ├─ Calculate covariance matrix
│  │  ├─ Build quaternion matrix
│  │  ├─ Eigenvalue decomposition
│  │  ├─ Convert to rotation matrix R
│  │  └─ Calculate translation t
│  ├─ mst2orien(orien[i], 0, R)  # Store rotation matrix
│  └─ org[i] = orgi               # Store origin
└─ Return
```

### Workflow 2: Base Pair Finding

```
find_bestpair(...)
├─ Initialize matched_idx[] = all zeros
├─ num_bp = 0
├─ while matched_count increases:
│  ├─ For each residue i:
│  │  ├─ if RY[i] < 0 or matched_idx[i]: continue
│  │  ├─ best_pair(i, ..., pair_istat)  # Find best partner j
│  │  │  ├─ For each j != i:
│  │  │  │  ├─ if RY[j] < 0 or matched_idx[j]: continue
│  │  │  │  ├─ check_pair(i, j, ..., rtn_val, &bpid, ..., 0)
│  │  │  │  │  ├─ Calculate geometric metrics
│  │  │  │  │  ├─ Check constraints (cdns)
│  │  │  │  │  ├─ Check overlap
│  │  │  │  │  ├─ Count H-bonds
│  │  │  │  │  └─ if valid: calculate bpid
│  │  │  │  └─ if bpid != 0 and rtn_val[5] < best: update best
│  │  │  └─ Return j with lowest rtn_val[5]
│  │  ├─ if pair_istat[1] != 0:  # Found partner j
│  │  │  ├─ best_pair(j, ..., pair_jstat)  # Check if j's best is i
│  │  │  └─ if i == pair_jstat[1]:  # Mutual match!
│  │  │     ├─ matched_idx[i] = 1
│  │  │     ├─ matched_idx[j] = 1
│  │  │     └─ base_pairs[++num_bp] = (i, j, ...)
│  └─ Check if matched_count increased
└─ Return num_bp
```

### Workflow 3: H-Bond Detection and Validation

```
get_hbond_ij(i, j, basei, basej, ..., hb_info)
├─ Phase 1: Initial Detection
│  ├─ For m in [seidx[i][1], seidx[i][2]]:
│  │  ├─ For n in [seidx[j][1], seidx[j][2]]:
│  │  │  ├─ if good_hbatoms(...) && within_limits(...):
│  │  │  │  ├─ Record: hb_atom1[k], hb_atom2[k], hb_dist[k]
│  │  │  │  └─ k++
│
├─ Phase 2: Conflict Resolution (hb_atompair)
│  ├─ while not all H-bonds processed:
│  │  ├─ For current H-bond k:
│  │  │  ├─ Find best H-bond for atom1 (closest)
│  │  │  ├─ Find best H-bond for atom2 (closest)
│  │  │  └─ if same H-bond is best for both:
│  │  │     ├─ Mark conflict: hb_dist[k] = -hb_dist[k]
│  │  │     └─ Mark related H-bonds as processed
│  ├─ Calculate linkage types (lkg_type)
│  └─ Additional conflict marking (if hb_dist2 > 0)
│
├─ Phase 3: Validation (validate_hbonds)
│  ├─ For each H-bond k:
│  │  ├─ if hb_dist[k] <= 0:  # Conflicted
│  │  │  ├─ hb_dist[k] = fabs(hb_dist[k])
│  │  │  └─ hb_type[k] = donor_acceptor(...)
│  │  └─ else: hb_type[k] = ' '
│  ├─ Filter invalid H-bonds (distance, type checks)
│  └─ Count valid: return count(hb_type != ' ')
│
└─ Format output string: "[count] atom1typeatom2 dist ..."
```

---

## Function Dependency Map

### Level 1: Core Mathematical Functions
```
ls_fitting()
├─ cov_matrix()
├─ jacobi() (eigenvalue decomposition)
├─ dot() (vector operations)
├─ ave_dmatrix()
└─ ddxyz()

get_oarea()
├─ ratom_xyz()
├─ get_zoave()
└─ Vector operations
```

### Level 2: Geometric Calculation Functions
```
base_frame()
├─ ls_fitting()
├─ read_pdb()
├─ find_1st_atom()
├─ set_std_base_pdb()
└─ mst2orien()

check_pair()
├─ get_bp_zoave()
├─ get_oarea()
├─ good_hbatoms()
├─ within_limits()
├─ calculate_more_bppars()
└─ adjust_pairQuality()
```

### Level 3: Base Pair Finding Functions
```
find_bestpair()
├─ best_pair()
│  └─ check_pair()
│     ├─ get_hbond_ij()
│     │  ├─ hb_atompair()
│     │  │  ├─ update_hb_idx()
│     │  │  └─ Array operations
│     │  └─ validate_hbonds()
│     │     └─ donor_acceptor()
│     └─ (other geometric checks)
└─ Array management
```

### Level 4: High-Level Workflow Functions
```
handle_str() (find_pair)
├─ read_pdb()
├─ residue_idx()
├─ get_seq()
├─ base_info()
│  └─ base_frame()
├─ find_bestpair()
│  └─ (dependency chain above)
├─ re_ordering()
└─ x3dna_input()

process_str() (analyze)
├─ read_input()
├─ read_pdb()
├─ ref_frames()
│  └─ base_frame()
└─ get_parameters()
   ├─ refs_i_j()
   ├─ bpstep_par()
   └─ helical_par()
```

---

## Indexing Conventions

### Critical: 1-based Indexing Throughout

**All arrays use 1-based indexing**:
- Residue indices: `i = 1..num_residue`
- Atom indices: `m = 1..num_atoms`
- Base pair indices: `k = 1..num_bp`

**Array Access Pattern**:
```c
// Allocate
double *vec = dvector(1, n);  // Allocates [0..n], uses [1..n]

// Access
for (i = 1; i <= n; i++) {
    vec[i] = value;  // NOT vec[i-1]
}
```

**2D Array Access**:
```c
// Allocate
double **mat = dmatrix(1, nrows, 1, ncols);

// Access
for (i = 1; i <= nrows; i++) {
    for (j = 1; j <= ncols; j++) {
        mat[i][j] = value;  // NOT mat[i-1][j-1]
    }
}
```

**Frame Matrix Extraction**:
```c
// orien[i] is 9 values: [1..9] in 1-based, [0..8] in 0-based view
// Extract 3×3 matrix R[1..3][1..3]:
R[1][1] = orien[i][0]  // Actually orien[i][1] in 1-based view
R[1][2] = orien[i][1]
R[1][3] = orien[i][2]
R[2][1] = orien[i][3]
// etc.

// z-axis at indices [6, 7, 8] (0-based) or [7, 8, 9] (1-based)
z[1] = orien[i][6]  // x-component of z-axis
z[2] = orien[i][7]  // y-component
z[3] = orien[i][8]  // z-component
```

**Residue Index Assignment**:
- Assigned sequentially during PDB parsing: 1, 2, 3, ...
- Includes ALL residues (not just nucleotides)
- Order matches PDB file order
- Stored in `seidx[][]` mapping

**Critical Conversion**:
- Legacy uses 1-based: `residue_idx = 1..N`
- Modern uses 0-based internally: `residue_idx = 0..N-1`
- Conversion: `legacy_idx = modern_idx + 1`

---

## Key Parameters and Constants

### Validation Parameters (`miscPars`)

```c
double min_dorg = 0.0;           // Minimum origin distance (Å)
double max_dorg = 15.0;          // Maximum origin distance (Å)
double min_dv = 0.0;             // Minimum vertical component (Å)
double max_dv = 2.5;             // Maximum vertical component (Å)
double min_dNN = 4.5;            // Minimum N-N distance (Å)
double max_dNN = XBIG (~1e18);   // Maximum N-N distance (effectively infinite)
double min_plane_angle = 0.0;    // Minimum plane angle (degrees)
double max_plane_angle = 65.0;   // Maximum plane angle (degrees)
double min_base_hb = 1;          // Minimum H-bonds required
double overlap_threshold = 0.01; // OVERLAP constant
```

### H-Bond Parameters

```c
double hb_lower = 1.8;    // Minimum H-bond distance (Å)
double hb_dist1 = 4.0;    // Maximum H-bond distance (Å)
double hb_dist2 = 0.0;    // Conflict resolution threshold (CRITICAL: must be 0.0)
```

**CRITICAL**: `hb_dist2 = 0.0` means Phase 3 conflict marking in `hb_atompair()` is effectively disabled.

### Constants

```c
#define OVERLAP 0.01              // Overlap threshold
#define XBIG 1.0e+18              // Very large number (effectively infinite)
#define NR_END 1                  // 1-based indexing offset
#define RTNNUM 37                 // Number of return values from check_pair
#define PSTNUM 29                 // Number of pair statistics
#define BUF512 512                // Buffer size
#define RA_NUM 9                  // Number of ring atoms
```

### Ring Atoms

```c
static char *RingAtom[] = {
    " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", 
    " C5 ", " N7 ", " C8 ", " N9 "
};
```

- Used for least-squares fitting
- Minimum 3 atoms required for frame calculation
- Purines (A, G): All 9 atoms
- Pyrimidines (C, T, U): First 6 atoms (no N7, C8, N9)

---

## Comparison Points with Modern Code

### Critical Matching Requirements

1. **Residue Indexing**:
   - Legacy: 1-based throughout
   - Modern: Must use `legacy_residue_idx` from atoms, convert correctly
   - **MUST MATCH**: Frame calculation uses same residue indices

2. **Base Frame Calculation**:
   - Legacy: `base_frame()` → `ls_fitting()`
   - Modern: `BaseFrameCalculator::calculate_frame()`
   - **MUST MATCH**: Rotation matrices and origins exactly

3. **H-Bond Detection**:
   - Legacy: `get_hbond_ij()` → `hb_atompair()` → `validate_hbonds()`
   - Modern: `find_hydrogen_bonds_detailed()` → `resolve_conflicts()` → `validate_hbonds()`
   - **MUST MATCH**: All H-bonds, distances, types

4. **Base Pair Validation**:
   - Legacy: `check_pair()` with exact parameter thresholds
   - Modern: `BasePairValidator::validate_pair()`
   - **MUST MATCH**: All geometric checks, quality scores

5. **Greedy Matching**:
   - Legacy: `find_bestpair()` → `best_pair()` (mutual matching)
   - Modern: `BasePairFinder::find_pairs()` (same algorithm)
   - **MUST MATCH**: Pair selection order and results

### Known Differences to Address

1. **Atom Selection**:
   - Legacy: Uses `seidx[i][1]..seidx[i][2]` (specific atom range)
   - Modern: Uses all atoms in residue
   - **ACTION**: Verify atom selection matches exactly

2. **H-Bond Count**:
   - Legacy sometimes finds fewer H-bonds
   - Modern sometimes finds more
   - **ACTION**: Compare atom-by-atom, verify `good_hbatoms()` logic

3. **Distance Calculations**:
   - Legacy: `p1p2_dist(xyz[n], xyz[m])` (order: j to i)
   - Modern: `(atom1.position() - atom2.position()).length()` (order: 1 to 2)
   - **ACTION**: Verify atom order matches

4. **Conflict Resolution**:
   - Legacy: `hb_dist2 = 0.0` disables Phase 3 conflicts
   - Modern: Should match (already set to 0.0)
   - **ACTION**: Verify Phase 3 is effectively disabled

### Verification Checklist

- [ ] Frame calculation: Rotation matrices match within tolerance
- [ ] Frame calculation: Origins match within tolerance
- [ ] H-bond detection: Same atoms found for each pair
- [ ] H-bond detection: Distances match exactly
- [ ] H-bond detection: Types match (' ', '-', '*')
- [ ] Base pair validation: Same pairs validated/rejected
- [ ] Base pair validation: Quality scores match
- [ ] Greedy matching: Same pairs selected
- [ ] Residue indexing: legacy_residue_idx matches exactly

---

## Expanded Algorithm Details

This section provides deep mathematical and algorithmic details for critical functions, including step-by-step pseudocode, numerical examples, and edge cases.

### 11. `ls_fitting()` - Complete Mathematical Derivation

**Quaternion-Based Least-Squares Fitting**

**Problem Statement**: Given two sets of n matched 3D points:
- Standard/reference points: `sxyz[i]` for i = 1..n
- Experimental points: `exyz[i]` for i = 1..n

Find rotation matrix R and translation vector t that minimize:
```
E = Σ ||R·sxyz[i] + t - exyz[i]||²
```

**Mathematical Derivation**:

1. **Centering**:
   ```
   ave_sxyz = (1/n) Σ sxyz[i]
   ave_exyz = (1/n) Σ exyz[i]
   
   sxyz_centered[i] = sxyz[i] - ave_sxyz
   exyz_centered[i] = exyz[i] - ave_exyz
   ```

2. **Covariance Matrix** (3×3):
   ```
   U[k][l] = (1/n) Σ sxyz_centered[i][k] × exyz_centered[i][l]
   
   For k,l ∈ {1,2,3} (x, y, z components)
   ```

3. **Quaternion Construction**:
   - Represent rotation as unit quaternion q = [q₀, q₁, q₂, q₃]
   - Build 4×4 symmetric matrix N from covariance U:
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
   - Use Jacobi method to find eigenvalues and eigenvectors of N
   - Largest eigenvalue corresponds to optimal quaternion
   - Eigenvector at index 4: V[1..4][4] gives quaternion components

5. **Quaternion to Rotation Matrix**:
   Let q = [q₀, q₁, q₂, q₃] = V[1..4][4]
   ```
   R[1][1] = q₀² + q₁² - q₂² - q₃²
   R[1][2] = 2(q₁q₂ - q₀q₃)
   R[1][3] = 2(q₁q₃ + q₀q₂)
   R[2][1] = 2(q₂q₁ + q₀q₃)
   R[2][2] = q₀² - q₁² + q₂² - q₃²
   R[2][3] = 2(q₂q₃ - q₀q₁)
   R[3][1] = 2(q₃q₁ - q₀q₂)
   R[3][2] = 2(q₃q₂ + q₀q₁)
   R[3][3] = q₀² - q₁² - q₂² + q₃²
   ```

6. **Translation Calculation**:
   ```
   t = ave_exyz - R × ave_sxyz
   ```

7. **RMS Calculation**:
   ```
   for i = 1 to n:
       fitted_xyz[i] = R × sxyz[i] + t
       error[i] = ||fitted_xyz[i] - exyz[i]||
   
   RMS = sqrt((1/n) Σ error[i]²)
   ```

**Complete Pseudocode**:
```python
def ls_fitting(sxyz, exyz, n):
    # Step 1: Calculate averages
    ave_sxyz = [0, 0, 0]
    ave_exyz = [0, 0, 0]
    for i in range(1, n+1):
        for k in range(1, 4):
            ave_sxyz[k] += sxyz[i][k]
            ave_exyz[k] += exyz[i][k]
    for k in range(1, 4):
        ave_sxyz[k] /= n
        ave_exyz[k] /= n
    
    # Step 2: Build covariance matrix U
    U = [[0]*4 for _ in range(4)]  # 3×3, 1-indexed
    for i in range(1, n+1):
        for k in range(1, 4):
            for l in range(1, 4):
                U[k][l] += (sxyz[i][k] - ave_sxyz[k]) * (exyz[i][l] - ave_exyz[l])
    for k in range(1, 4):
        for l in range(1, 4):
            U[k][l] /= n
    
    # Step 3: Build quaternion matrix N
    N = [[0]*5 for _ in range(5)]  # 4×4, 1-indexed
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
    
    # Step 4: Eigenvalue decomposition
    D, V = jacobi(N, 4)  # D = eigenvalues, V = eigenvectors
    
    # Step 5: Extract quaternion (largest eigenvalue at index 4)
    q0 = V[1][4]
    q1 = V[2][4]
    q2 = V[3][4]
    q3 = V[4][4]
    
    # Step 6: Build rotation matrix R
    R[1][1] = q0*q0 + q1*q1 - q2*q2 - q3*q3
    R[1][2] = 2*(q1*q2 - q0*q3)
    R[1][3] = 2*(q1*q3 + q0*q2)
    R[2][1] = 2*(q2*q1 + q0*q3)
    R[2][2] = q0*q0 - q1*q1 + q2*q2 - q3*q3
    R[2][3] = 2*(q2*q3 - q0*q1)
    R[3][1] = 2*(q3*q1 - q0*q2)
    R[3][2] = 2*(q3*q2 + q0*q1)
    R[3][3] = q0*q0 - q1*q1 - q2*q2 + q3*q3
    
    # Step 7: Calculate translation
    for i in range(1, 4):
        t[i] = ave_exyz[i] - dot(R[i], ave_sxyz)
    
    # Step 8: Calculate RMS
    sum_sq_error = 0
    for i in range(1, n+1):
        fitted[i] = R × sxyz[i] + t
        error = fitted[i] - exyz[i]
        sum_sq_error += dot(error, error)
    RMS = sqrt(sum_sq_error / n)
    
    return RMS
```

**Example with Numbers**:
```
Given 3 matched points:
Standard:  [0, 0, 0], [1, 0, 0], [0, 1, 0]
Experimental: [0.1, 0.05, 0], [1.05, 0.1, 0], [0.1, 1.05, 0]

After centering and covariance calculation, quaternion method finds:
R ≈ [[0.998, -0.05, 0], [0.05, 0.998, 0], [0, 0, 1]]
t ≈ [0.1, 0.05, 0]
RMS ≈ 0.07 Å
```

**Edge Cases**:
- **n < 3**: Fatal error - need minimum 3 points
- **Degenerate points**: If all points collinear, quaternion may be unstable
- **Perfect match**: RMS = 0 when points align exactly
- **Large rotation**: Quaternion method handles all rotations robustly

---

### 12. `bpstep_par()` - Base Pair Step Parameters Calculation

**Location**: `org/src/ana_fncs.c:1396-1441`

**Purpose**: Calculate 6 base pair step parameters: Shift, Slide, Rise, Tilt, Roll, Twist

**Input**:
- `rot1[3][3]`, `org1[3]`: Frame 1 (rotation matrix + origin)
- `rot2[3][3]`, `org2[3]`: Frame 2 (rotation matrix + origin)

**Output**:
- `pars[6]`: [Shift, Slide, Rise, Tilt, Roll, Twist]

**Mathematical Algorithm**:

1. **Extract y-axes** (axis 2 in 0-based, 3 in 1-based):
   ```
   t1 = rot1[1..3][2]  # y-axis of frame 1
   t2 = rot2[1..3][2]  # y-axis of frame 2
   ```

2. **Calculate Hinge Vector**:
   ```
   hinge = cross(t1, t2)  # Perpendicular to both y-axes
   rolltilt = angle(t1, t2)  # Angle between y-axes (0-180°)
   ```

3. **Edge Case: Degenerate Hinge**:
   ```
   if |hinge| < XEPS (≈ 1e-7) && (rolltilt ≈ 0° or ≈ 180°):
       # y-axes are parallel - use x and y axes
       hinge = rot1[1..3][1] + rot2[1..3][1] + rot1[1..3][2] + rot2[1..3][2]
   ```

4. **Create Parallel Frames**:
   ```
   # Rotate frame 2 back by half the angle
   temp = rotation_matrix(hinge, -0.5 * rolltilt)
   para_bp2 = temp × rot2
   
   # Rotate frame 1 forward by half the angle
   temp = rotation_matrix(hinge, 0.5 * rolltilt)
   para_bp1 = temp × rot1
   ```

5. **Extract Midstep Z-axis**:
   ```
   mstz = para_bp2[1..3][3]  # z-axis from parallel frame 2
   ```

6. **Calculate Twist**:
   ```
   t1 = para_bp1[1..3][2]  # y-axis of parallel frame 1
   t2 = para_bp2[1..3][2]  # y-axis of parallel frame 2
   Twist = angle(t1, t2, mstz)  # Angle around mstz axis
   ```

7. **Calculate Midstep Y-axis**:
   ```
   msty = vector_between(t1, t2, rotated by Twist/2 around mstz)
   ```

8. **Calculate Midstep X-axis**:
   ```
   mstx = cross(msty, mstz)  # Orthogonal to y and z
   ```

9. **Build Midstep Frame**:
   ```
   mst_orien = [mstx, msty, mstz]  # 3×3 matrix
   mst_org = (org1 + org2) / 2  # Average origin
   ```

10. **Calculate Translations (Shift, Slide, Rise)**:
    ```
    dorg = org2 - org1  # Vector from frame 1 to frame 2
    
    # Project into midpoint frame
    Shift = dot(dorg, mstx)    # x-component
    Slide = dot(dorg, msty)    # y-component
    Rise = dot(dorg, mstz)     # z-component
    ```

11. **Calculate Rotations (Tilt, Roll)**:
    ```
    phi = angle(hinge, msty, mstz)  # Angle of hinge in midpoint frame
    
    Roll = rolltilt × cos(phi)   # Rotation around y-axis
    Tilt = rolltilt × sin(phi)   # Rotation around x-axis
    ```

**Complete Pseudocode**:
```python
def bpstep_par(rot1, org1, rot2, org2):
    # Extract y-axes (column 2 in 1-based indexing, column 3 in 0-based)
    t1 = [rot1[i][2] for i in range(1, 4)]  # y-axis of frame 1
    t2 = [rot2[i][2] for i in range(1, 4)]  # y-axis of frame 2
    
    # Calculate hinge vector
    hinge = cross(t1, t2)
    rolltilt = angle(t1, t2)  # Angle between y-axes
    
    # Handle degenerate case
    if length(hinge) < XEPS and (abs(rolltilt - 180) < XEPS or rolltilt < XEPS):
        hinge = rot1[1..3][1] + rot2[1..3][1] + rot1[1..3][2] + rot2[1..3][2]
    
    # Create parallel frames
    temp1 = rotation_matrix(hinge, -0.5 * rolltilt)
    para_bp2 = temp1 × rot2
    
    temp2 = rotation_matrix(hinge, 0.5 * rolltilt)
    para_bp1 = temp2 × rot1
    
    # Extract z-axis
    mstz = [para_bp2[i][3] for i in range(1, 4)]
    
    # Calculate twist
    t1_para = [para_bp1[i][2] for i in range(1, 4)]
    t2_para = [para_bp2[i][2] for i in range(1, 4)]
    Twist = angle(t1_para, t2_para, mstz)
    
    # Calculate midpoint y-axis
    msty = vector_at_half_angle(t1_para, t2_para, Twist/2, mstz)
    
    # Calculate midpoint x-axis
    mstx = cross(msty, mstz)
    
    # Build midpoint frame
    mst_orien = [mstx, msty, mstz]
    mst_org = (org1 + org2) / 2
    
    # Calculate translations
    dorg = org2 - org1
    Shift = dot(dorg, mstx)
    Slide = dot(dorg, msty)
    Rise = dot(dorg, mstz)
    
    # Calculate rotations
    phi = angle(hinge, msty, mstz)  # in radians
    Roll = rolltilt * cos(phi)
    Tilt = rolltilt * sin(phi)
    
    pars = [Shift, Slide, Rise, Tilt, Roll, Twist]
    return pars, mst_orien, mst_org
```

**Typical Values**:
- **B-DNA**: Shift ≈ 0.0, Slide ≈ -0.5, Rise ≈ 3.4, Tilt ≈ 0°, Roll ≈ 0°, Twist ≈ 36°
- **A-DNA**: Rise ≈ 2.6, Twist ≈ 33°
- **Z-DNA**: Roll ≈ 7°, Twist ≈ -10° (left-handed)

---

### 13. `good_hbatoms()` - H-Bond Atom Pair Validation

**Location**: `org/src/cmn_fncs.c:3864-3877`

**Purpose**: Determine if two atoms can potentially form a hydrogen bond based on atom types.

**Input**:
- `atom1`, `atom2`: Atom names (e.g., " N3 ", " O2 ")
- `idx1`, `idx2`: Atom type indices from `atom_idx()`
- `misc_pars`: Contains `hb_idx[]` array with valid H-bond atom indices

**Algorithm**:

1. **Check for Phosphate/Oxygen Exclusion**:
   ```
   PO_list = [" O1P", " O2P", " O3'", " O4'", " O5'", " N7 "]
   
   if atom1 in PO_list AND atom2 in PO_list:
       return FALSE  # Exclude phosphate-phosphate or similar
   ```

2. **Check Atom Type Indices**:
   ```
   # idx values:
   # idx = 1: Nitrogen (N)
   # idx = 2: Oxygen (O)  
   # idx = 4: Other (e.g., Carbon)
   
   # At least one atom must be O or N (idx = 2 or 1)
   if (idx1 == 2 OR idx1 == 4 OR idx2 == 2 OR idx2 == 4):
       # At least one is O or N
       if idx1 in hb_idx AND idx2 in hb_idx:
           return TRUE
   ```

**Complete Pseudocode**:
```python
def good_hbatoms(misc_pars, atom1, atom2, idx1, idx2):
    # Exclusion list
    PO = [" O1P", " O2P", " O3'", " O4'", " O5'", " N7 "]
    
    # Exclude if both are in PO list
    if atom1 in PO and atom2 in PO:
        return FALSE
    
    # Check if atoms are valid H-bond types
    natom = misc_pars.hb_idx[0]  # Number of valid types
    
    # At least one must be O (idx=2) or N (idx=1)
    if (idx1 == 2 or idx1 == 4 or idx2 == 2 or idx2 == 4):
        # Both must be in valid H-bond atom list
        if idx1 in misc_pars.hb_idx[1..natom] and \
           idx2 in misc_pars.hb_idx[1..natom]:
            return TRUE
    
    return FALSE
```

**Valid H-Bond Atom Types** (from `hb_idx[]`):
- Typically includes: N (idx=1), O (idx=2), and possibly others
- Stored in `misc_pars->hb_idx[1..natom]`

**Examples**:
- `good_hbatoms(" N3 ", " O2 ", 1, 2)` → TRUE (N-O pair)
- `good_hbatoms(" O1P", " O2P", 2, 2)` → FALSE (both phosphate)
- `good_hbatoms(" C1'", " C2'", 4, 4)` → FALSE (both carbon)

---

### 14. `donor_acceptor()` - H-Bond Type Determination

**Location**: `org/src/cmn_fncs.c:4085-4133`

**Purpose**: Determine if an H-bond has correct donor-acceptor roles based on base and atom types.

**Input**:
- `basei`, `basej`: Base types ('A', 'C', 'G', 'T', 'U')
- `hb_atom1`, `hb_atom2`: Atom names

**Output**:
- `' '`: Invalid or normal (no special marking)
- `'-'`: Conflict H-bond (wrong donor/acceptor roles)
- `'*'`: Other/non-standard H-bond

**Algorithm**:

1. **Base Type Mapping**:
   ```
   CB_LIST = "ACGITU"  # Valid bases
   inum = position of basei in CB_LIST  # 0=A, 1=C, 2=G, etc.
   jnum = position of basej in CB_LIST
   ```

2. **Atom Role Lookup**:
   ```
   # Backbone atoms with roles
   bb_da = [
       " O1P_A",  # O1P is Acceptor
       " O2P_A",  # O2P is Acceptor
       " O5'_A",  # O5' is Acceptor
       " O4'_A",  # O4' is Acceptor
       " O3'_A",  # O3' is Acceptor
       " O2'_X"   # O2' is special
   ]
   
   # Base atoms with roles (per base type)
   base_da[inum][0..5] = [
       {" N9 _?", " N7 _A", " N6 _D", " N1 _A", " N3 _A"},  # A
       {" N1 _?", " O2 _A", " N3 _A", " N4 _D"},             # C
       {" N9 _?", " N7 _A", " O6 _A", " N1 _D", " N2 _D", " N3 _A"}, # G
       # ... etc
   ]
   
   # D = Donor, A = Acceptor, ? = Ambiguous
   ```

3. **Role Extraction**:
   ```
   ia = role of atom1  # 'D', 'A', or '\0' if not found
   ja = role of atom2  # 'D', 'A', or '\0' if not found
   
   # Check backbone atoms first
   for each bb_da entry:
       if atom1 matches: ia = last character (A/X)
       if atom2 matches: ja = last character (A/X)
   
   # Check base atoms
   if ia not found:
       for each base_da[inum] entry:
           if atom1 matches: ia = role character
   
   if ja not found:
       for each base_da[jnum] entry:
           if atom2 matches: ja = role character
   ```

4. **Type Determination**:
   ```
   da_types = ["AD", "AX", "XD", "XX", "DA", "DX", "XA"]
   
   if ia and ja found:
       da = ia + ja  # Concatenate roles
       if da in da_types:
           return '-'  # Conflict (wrong roles)
   
   return '*'  # Non-standard
   ```

**Complete Pseudocode**:
```python
def donor_acceptor(basei, basej, hb_atom1, hb_atom2):
    # Base type indices
    CB_LIST = "ACGITU"
    inum = index of basei in CB_LIST  # -1 if not found
    jnum = index of basej in CB_LIST
    
    if inum < 0 or jnum < 0:
        return '*'  # Invalid bases
    
    # Backbone atom roles
    bb_da = [
        " O1P_A", " O2P_A", " O5'_A", 
        " O4'_A", " O3'_A", " O2'_X"
    ]
    
    # Base atom roles (example for A, C, G)
    base_da = [
        # A: index 0
        [" N9 _?", " N7 _A", " N6 _D", " N1 _A", " N3 _A"],
        # C: index 1
        [" N1 _?", " O2 _A", " N3 _A", " N4 _D"],
        # G: index 2
        [" N9 _?", " N7 _A", " O6 _A", " N1 _D", " N2 _D", " N3 _A"],
        # ... T, U, I
    ]
    
    # Find roles
    ia = None
    ja = None
    
    # Check backbone
    for bb in bb_da:
        if hb_atom1 matches bb[0:4]:
            ia = bb[5]  # 'A' or 'X'
        if hb_atom2 matches bb[0:4]:
            ja = bb[5]
    
    # Check base atoms
    if not ia:
        for atom_role in base_da[inum]:
            if hb_atom1 matches atom_role[0:4]:
                ia = atom_role[5]  # 'D', 'A', or '?'
    
    if not ja:
        for atom_role in base_da[jnum]:
            if hb_atom2 matches atom_role[0:4]:
                ja = atom_role[5]
    
    # Determine type
    if ia and ja:
        da = ia + ja
        conflict_types = ["AD", "AX", "XD", "XX", "DA", "DX", "XA"]
        if da in conflict_types:
            return '-'  # Conflict H-bond
    
    return '*'  # Non-standard
```

**Examples**:
- `donor_acceptor('A', 'T', " N6 ", " O4 ")` → 'D' + 'A' → Not in conflict list → '*'
- `donor_acceptor('C', 'G', " N3 ", " O6 ")` → 'A' + 'A' → 'AA' in conflicts → '-'
- `donor_acceptor('G', 'C', " O2'", " N2 ")` → 'X' + 'D' → 'XD' in conflicts → '-'

---

### 15. `adjust_pairQuality()` - Quality Score Adjustment

**Location**: `org/src/cmn_fncs.c:4553-4572`

**Purpose**: Adjust pair quality score based on H-bond characteristics.

**Algorithm**:

1. **Get H-Bond List**:
   ```
   hb_numlist(i, j, basei, basej, ..., &num_hb, num_list)
   # Returns num_list[k][0..2]:
   #   [0] = conflict flag
   #   [3] = distance (stored as integer × MFACTOR)
   ```

2. **Count "Good" H-Bonds**:
   ```
   num_good_hb = 0
   for k = 1 to num_hb:
       if num_list[k][0] != 0:  # Conflicted
           continue
       
       dval = num_list[k][3] / MFACTOR  # Convert to real distance
       
       if 2.5 <= dval <= 3.5:  # Good H-bond distance range
           num_good_hb++
   ```

3. **Calculate Adjustment**:
   ```
   if num_good_hb >= 2:
       return -3.0  # Strong bonus
   else:
       return -num_good_hb  # One good H-bond = -1.0, none = 0.0
   ```

**Complete Pseudocode**:
```python
def adjust_pairQuality(i, j, basei, basej, ...):
    num_list = allocate_matrix(BUF512, 4)
    num_hb = 0
    
    # Get H-bond list
    hb_numlist(i, j, basei, basej, ..., &num_hb, num_list)
    
    # Count good H-bonds
    num_good_hb = 0
    for k in range(1, num_hb + 1):
        if num_list[k][0] != 0:  # Conflicted
            continue
        
        dval = num_list[k][3] / MFACTOR  # Distance in Angstroms
        
        if 2.5 <= dval <= 3.5:  # Good range
            num_good_hb += 1
    
    # Adjust quality score
    if num_good_hb >= 2:
        return -3.0  # Strong bonus for 2+ good H-bonds
    else:
        return -num_good_hb  # Smaller bonus
```

**Effect on Quality Score**:
- Original: `rtn_val[5] = dorg + 2*dv + plane_angle/20`
- After adjustment: `rtn_val[5] += adjust_pairQuality()`
- **Lower score = better pair** (negative adjustment improves quality)
- Two good H-bonds: score reduced by 3.0
- One good H-bond: score reduced by 1.0

---

### 16. `check_wc_wobble_pair()` - Watson-Crick/Wobble Classification

**Location**: `org/src/ana_fncs.c:1122-1131`

**Purpose**: Classify base pair as Watson-Crick (WC) or wobble based on step parameters.

**Input**:
- `bp`: Base pair string (e.g., "GC", "AT")
- `shear`, `stretch`: Step parameters (from `bpstep_par`)
- `opening`: Another step parameter

**Output**:
- `bpid = 1`: Non-WC/Wobble pair
- `bpid = 2`: Watson-Crick or Wobble pair
- `bpid` unchanged: Outside parameter ranges

**Algorithm**:
```python
def check_wc_wobble_pair(bpid, bp, shear, stretch, opening):
    WC_LIST = ["XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"]
    
    # Check parameter ranges
    if abs(stretch) > 2.0 or abs(opening) > 60:
        return  # Outside valid range
    
    # Check shear magnitude
    if 1.8 <= abs(shear) <= 2.8:
        *bpid = 1  # Non-WC/Wobble (shear too large)
    
    # Check for WC/Wobble
    if abs(shear) <= 1.8 and bp in WC_LIST:
        *bpid = 2  # WC or Wobble pair
```

**Parameter Ranges**:
- **Valid range**: `|stretch| ≤ 2.0` AND `|opening| ≤ 60°`
- **WC/Wobble**: `|shear| ≤ 1.8` AND base pair in canonical list
- **Non-WC**: `1.8 < |shear| ≤ 2.8`

**Canonical Pairs** (WC_LIST):
- "AT", "AU" (Watson-Crick)
- "TA", "UA" (reverse)
- "GC" (Watson-Crick)
- "CG" (reverse)
- "IC" (Inosine-Cytosine, wobble)
- "CI" (reverse)

**Effect in `calculate_more_bppars()`**:
```c
if (dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0) {
    check_wc_wobble_pair(bpid, bpi, pars[1], pars[2], pars[6]);
    if (*bpid == 2) {
        rtn_val[5] -= 2.0;  // Bonus for WC/Wobble pairs
    }
}
```

---

### 17. `get_oarea()` - Base Overlap Area Calculation

**Location**: `org/src/ana_fncs.c:3327-3358`

**Purpose**: Calculate overlap area between two base rings projected onto a plane.

**Algorithm**:

1. **Get Ring Atom Coordinates**:
   ```
   n1 = ratom_xyz(ring_atom[r1], only_ring, xyz, oave, oxyz1)
   n2 = ratom_xyz(ring_atom[r2], only_ring, xyz, oave, oxyz2)
   # Returns number of ring atoms and their coordinates
   ```

2. **Align to Z-axis**:
   ```
   # Rotate both sets to align zave with z-axis
   align2zaxis(n1, zave, rotmat, oxyz1, oxyz1Z)
   align2zaxis(n2, zave, rotmat, oxyz2, oxyz2Z)
   ```

3. **Project to XY Plane**:
   ```
   for i = 1 to n1:
       a[i-1].x = oxyz1Z[i][1]  # x-coordinate
       a[i-1].y = oxyz1Z[i][2]  # y-coordinate
       # z-coordinate ignored (projection)
   
   for i = 1 to n2:
       b[i-1].x = oxyz2Z[i][1]
       b[i-1].y = oxyz2Z[i][2]
   ```

4. **Calculate Polygon Intersection**:
   ```
   overlap_area = pia_inter(a, n1, b, n2)
   # Uses polygon intersection algorithm
   ```

**Complete Pseudocode**:
```python
def get_oarea(r1, r2, ring_atom, oave, zave, xyz, only_ring):
    # Get ring atoms
    oxyz1 = allocate_matrix(9, 3)
    oxyz2 = allocate_matrix(9, 3)
    
    n1 = ratom_xyz(ring_atom[r1], only_ring, xyz, oave, oxyz1)
    n2 = ratom_xyz(ring_atom[r2], only_ring, xyz, oave, oxyz2)
    
    # Align to z-axis
    rotmat = allocate_matrix(3, 3)
    oxyz1Z = allocate_matrix(9, 3)
    oxyz2Z = allocate_matrix(9, 3)
    
    align2zaxis(n1, zave, rotmat, oxyz1, oxyz1Z)
    align2zaxis(n2, zave, rotmat, oxyz2, oxyz2Z)
    
    # Create 2D polygons
    polygon_a = []
    for i in range(1, n1 + 1):
        polygon_a.append(Point(oxyz1Z[i][1], oxyz1Z[i][2]))
    
    polygon_b = []
    for i in range(1, n2 + 1):
        polygon_b.append(Point(oxyz2Z[i][1], oxyz2Z[i][2]))
    
    # Calculate intersection area
    overlap = polygon_intersection_area(polygon_a, polygon_b)
    
    return overlap
```

**Polygon Intersection Algorithm** (`pia_inter`):
- Uses convex polygon intersection
- Converts to 2D by projecting onto plane perpendicular to zave
- Calculates overlapping area using computational geometry

**Threshold**:
- `OVERLAP = 0.01` Å²
- If overlap ≥ 0.01: Bases are too close (reject pair)
- Used in three places in `check_pair()`:
  1. Overlap in average plane
  2. Overlap using base 1's origin and z-axis
  3. Overlap using base 2's origin and z-axis

---

### 18. `re_ordering()` - Base Pair Reordering Algorithm

**Location**: `org/src/find_pair.c:1430-1469`

**Purpose**: Reorder base pairs to follow 5'→3' direction and identify helical regions.

**Main Steps**:

1. **Save Original Order**:
   ```
   Write base pair information BEFORE reordering to file
   ```

2. **Build Pair Context**:
   ```
   bp_context(num_bp, misc_pars, bp_xyz, bp_order, end_list, &num_ends, fp)
   # Analyzes neighbor relationships
   # Identifies chain ends
   ```

3. **Locate Helices**:
   ```
   locate_helix(num_bp, helix_idx, num_ends, num_helix, end_list, bp_order, bp_idx, helix_marker)
   # Groups pairs into helical regions
   ```

4. **Five-to-Three Ordering**:
   ```
   five2three(num_bp, num_helix, helix_idx, bp_idx, bp_xyz, base_pairs, o3_p, fp)
   # Reorders pairs to follow 5'→3' direction
   ```

5. **Z-DNA Check**:
   ```
   check_zdna(num_helix, helix_idx, bp_idx, bp_xyz, base_pairs, fp)
   # Special handling for Z-DNA
   ```

**Key Functions**:

**`bp_context()`**: Analyzes pair relationships
- Determines which pairs are neighbors
- Identifies chain ends (pairs with only one neighbor)
- Builds `bp_order[][]` with neighbor information

**`locate_helix()`**: Groups pairs into helices
- Finds connected components of pairs
- Marks helix boundaries
- Stores in `helix_idx[][]` array

**`five2three()`**: Reorders within each helix
- Uses O3'-P distances to determine direction
- Ensures 5'→3' ordering
- Handles parallel/antiparallel strands

**Complete Workflow**:
```
re_ordering():
    1. Write original order to file
    
    2. Build pair context:
       - For each pair, find neighbors
       - Identify chain ends (num_ends)
       - Build bp_order[1..num_bp][1..3]
    
    3. Locate helices:
       - Find connected components
       - Group pairs into helix regions
       - Store in helix_idx[1..num_helix][1..7]
    
    4. Five-to-three ordering:
       - For each helix:
           - Determine strand direction from O3'-P links
           - Reorder pairs to follow 5'→3'
           - Update bp_idx[] array
    
    5. Z-DNA check:
       - Detect Z-DNA regions
       - Apply special handling
```

---

## Summary

### Critical Functions for 100% Match

1. **`ls_fitting()`**: Must match exactly (quaternion method)
2. **`base_frame()`**: Frame calculation pipeline
3. **`check_pair()`**: Geometric validation logic
4. **`hb_atompair()`**: Conflict resolution algorithm
5. **`validate_hbonds()`**: H-bond type assignment
6. **`get_hbond_ij()`**: Initial H-bond detection
7. **`find_bestpair()`**: Greedy matching algorithm

### Most Error-Prone Areas

1. **Indexing conversions** (1-based vs 0-based)
2. **Frame matrix extraction** (flattened 9-element array)
3. **H-bond conflict resolution** (iterative algorithm)
4. **Atom selection** (seidx ranges vs all atoms)
5. **Parameter thresholds** (must match exactly)

### Recommended Testing Strategy

1. **Isolate each function** with test executables
2. **Compare step-by-step** with modern implementation
3. **Use small PDB subsets** for faster iteration
4. **Add detailed debug output** at each stage
5. **Verify atom-by-atom** for H-bond detection

---

**Document Status**: Complete profiling of all critical legacy functions  
**Last Updated**: 2025-01-XX  
**Next Steps**: Create test executables for each critical function


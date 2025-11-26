# Legacy Algorithms - Mathematical Details

**Date**: 2025-01-XX  
**Purpose**: Complete mathematical derivations, pseudocode, and numerical examples for all critical algorithms  
**Status**: Comprehensive algorithm reference

---

## Table of Contents

1. [Least-Squares Fitting](#1-least-squares-fitting)
2. [Base Pair Step Parameters](#2-base-pair-step-parameters)
3. [H-Bond Detection Workflow](#3-h-bond-detection-workflow)
4. [H-Bond Conflict Resolution](#4-h-bond-conflict-resolution)
5. [Overlap Area Calculation](#5-overlap-area-calculation)
6. [Quality Score Calculation](#6-quality-score-calculation)
7. [Greedy Matching Algorithm](#7-greedy-matching-algorithm)

---

## 1. Least-Squares Fitting

### Problem Statement

Given two sets of n matched 3D points:
- Standard/reference points: `sxyz[i]` for i = 1..n
- Experimental points: `exyz[i]` for i = 1..n

Find rotation matrix R and translation vector t that minimize:
```
E = Σ ||R·sxyz[i] + t - exyz[i]||²
```

### Mathematical Derivation

#### Step 1: Centering

```
ave_sxyz[k] = (1/n) Σ sxyz[i][k]  for k = 1..3
ave_exyz[k] = (1/n) Σ exyz[i][k]  for k = 1..3

sxyz_centered[i][k] = sxyz[i][k] - ave_sxyz[k]
exyz_centered[i][k] = exyz[i][k] - ave_exyz[k]
```

#### Step 2: Covariance Matrix (3×3)

```
U[k][l] = (1/n) Σ sxyz_centered[i][k] × exyz_centered[i][l]

For k,l ∈ {1,2,3} (x, y, z components)
```

**Physical Meaning**: Measures correlation between standard and experimental coordinates.

#### Step 3: Quaternion Construction

Represent rotation as unit quaternion q = [q₀, q₁, q₂, q₃]

Build 4×4 symmetric matrix N from covariance U:

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

#### Step 4: Eigenvalue Decomposition

- Use Jacobi method to find eigenvalues and eigenvectors of N
- Largest eigenvalue corresponds to optimal quaternion
- Eigenvector at index 4: `V[1..4][4]` gives quaternion components

**Why largest eigenvalue?**: Maximizes the objective function for optimal rotation.

#### Step 5: Quaternion to Rotation Matrix

Let `q = [q₀, q₁, q₂, q₃] = V[1..4][4]`

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

**Properties**:
- R is orthogonal: `R^T × R = I`
- R is proper: `det(R) = 1`
- R preserves handedness

#### Step 6: Translation Calculation

```
t[k] = ave_exyz[k] - Σ R[k][l] × ave_sxyz[l]  for k = 1..3

Or in matrix form:
t = ave_exyz - R × ave_sxyz
```

**Physical Meaning**: Translation needed to align centroids after rotation.

#### Step 7: RMS Calculation

```
for i = 1 to n:
    fitted_xyz[i] = R × sxyz[i] + t
    error[i] = fitted_xyz[i] - exyz[i]

RMS = sqrt((1/n) Σ ||error[i]||²)
```

**RMS Interpretation**:
- Lower RMS = better fit
- Typical values: 0.1-0.3 Å for good fits
- RMS > 1.0 Å indicates poor alignment

### Complete Pseudocode

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
    for k in range(1, 4):
        t[k] = ave_exyz[k]
        for l in range(1, 4):
            t[k] -= R[k][l] * ave_sxyz[l]
    
    # Step 8: Calculate RMS
    sum_sq_error = 0
    for i in range(1, n+1):
        fitted = [0, 0, 0]
        for k in range(1, 4):
            for l in range(1, 4):
                fitted[k] += R[k][l] * sxyz[i][l]
            fitted[k] += t[k]
        
        error_sq = 0
        for k in range(1, 4):
            diff = fitted[k] - exyz[i][k]
            error_sq += diff * diff
        sum_sq_error += error_sq
    
    RMS = sqrt(sum_sq_error / n)
    
    return RMS, R, t
```

### Numerical Example

```
Given 3 matched points:
Standard:    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]
Experimental: [0.1, 0.05, 0.0], [1.05, 0.1, 0.0], [0.1, 1.05, 0.0]

After centering and covariance calculation:
U ≈ [[1.0, 0.05, 0], [0.05, 1.0, 0], [0, 0, 0]]

Quaternion method finds:
R ≈ [[0.998, -0.05, 0], [0.05, 0.998, 0], [0, 0, 1]]
t ≈ [0.1, 0.05, 0]
RMS ≈ 0.07 Å
```

### Edge Cases

- **n < 3**: Fatal error - need minimum 3 points for unique solution
- **Degenerate points**: If all points collinear, quaternion may be unstable
- **Perfect match**: RMS = 0 when points align exactly
- **Large rotation**: Quaternion method handles all rotations robustly (0-360°)
- **Noise**: Method is robust to noise in experimental coordinates

---

## 2. Base Pair Step Parameters

### Problem Statement

Calculate 6 base pair step parameters (Shift, Slide, Rise, Tilt, Roll, Twist) from two consecutive base pair frames.

**Input**:
- Frame 1: `rot1[3][3]`, `org1[3]`
- Frame 2: `rot2[3][3]`, `org2[3]`

**Output**:
- `pars[6]`: [Shift, Slide, Rise, Tilt, Roll, Twist]

### Mathematical Algorithm

#### Step 1: Extract y-Axes

```
t1 = rot1[1..3][2]  # y-axis of frame 1 (column 2)
t2 = rot2[1..3][2]  # y-axis of frame 2 (column 2)
```

**Note**: In 1-based indexing, column 2 is y-axis. In 0-based, it's column 1.

#### Step 2: Calculate Hinge Vector

```
hinge = cross(t1, t2)  # Perpendicular to both y-axes
rolltilt = angle(t1, t2)  # Angle between y-axes (0-180°)
```

**Physical Meaning**: Hinge is the axis of rotation between the two frames.

#### Step 3: Edge Case - Degenerate Hinge

```
if |hinge| < XEPS (≈ 1e-7) && (rolltilt ≈ 0° or ≈ 180°):
    # y-axes are parallel - use x and y axes
    hinge = rot1[1..3][1] + rot2[1..3][1] + rot1[1..3][2] + rot2[1..3][2]
    normalize(hinge)
```

**Why needed**: When y-axes are parallel, cross product is zero. Need alternative method.

#### Step 4: Create Parallel Frames

```
# Rotate frame 2 back by half the angle
temp = rotation_matrix(hinge, -0.5 * rolltilt)
para_bp2 = temp × rot2

# Rotate frame 1 forward by half the angle
temp = rotation_matrix(hinge, 0.5 * rolltilt)
para_bp1 = temp × rot1
```

**Purpose**: Create two frames with parallel y-axes, making twist calculation easier.

#### Step 5: Extract Midstep Z-axis

```
mstz = para_bp2[1..3][3]  # z-axis from parallel frame 2
```

**Physical Meaning**: Average z-axis direction for the step.

#### Step 6: Calculate Twist

```
t1_para = para_bp1[1..3][2]  # y-axis of parallel frame 1
t2_para = para_bp2[1..3][2]  # y-axis of parallel frame 2
Twist = angle(t1_para, t2_para, mstz)  # Angle around mstz axis
```

**Physical Meaning**: Rotation around z-axis (helical twist).

#### Step 7: Calculate Midstep Y-axis

```
# Average of y-axes, rotated by Twist/2 around mstz
msty = rotate((t1_para + t2_para) / 2, mstz, Twist/2)
normalize(msty)
```

#### Step 8: Calculate Midstep X-axis

```
mstx = cross(msty, mstz)  # Orthogonal to y and z
normalize(mstx)
```

#### Step 9: Build Midstep Frame

```
mst_orien = [mstx, msty, mstz]  # 3×3 matrix (columns)
mst_org = (org1 + org2) / 2  # Average origin
```

#### Step 10: Calculate Translations

```
dorg = org2 - org1  # Vector from frame 1 to frame 2

# Project into midpoint frame
Shift = dot(dorg, mstx)    # x-component (lateral)
Slide = dot(dorg, msty)    # y-component (lateral)
Rise = dot(dorg, mstz)     # z-component (vertical)
```

**Physical Meanings**:
- **Shift**: Lateral displacement along x-axis
- **Slide**: Lateral displacement along y-axis
- **Rise**: Vertical separation along z-axis

#### Step 11: Calculate Rotations

```
phi = angle(hinge, msty, mstz)  # Angle of hinge in midpoint frame

Roll = rolltilt × cos(phi)   # Rotation around y-axis
Tilt = rolltilt × sin(phi)   # Rotation around x-axis
```

**Physical Meanings**:
- **Roll**: Rotation around y-axis (bending in y-z plane)
- **Tilt**: Rotation around x-axis (bending in x-z plane)

### Complete Pseudocode

```python
def bpstep_par(rot1, org1, rot2, org2):
    # Extract y-axes (column 2 in 1-based indexing)
    t1 = [rot1[i][2] for i in range(1, 4)]  # y-axis of frame 1
    t2 = [rot2[i][2] for i in range(1, 4)]  # y-axis of frame 2
    
    # Calculate hinge vector
    hinge = cross(t1, t2)
    rolltilt = angle(t1, t2)  # Angle between y-axes
    
    # Handle degenerate case
    if length(hinge) < XEPS and (abs(rolltilt - 180) < XEPS or rolltilt < XEPS):
        hinge = [rot1[i][1] + rot2[i][1] + rot1[i][2] + rot2[i][2] 
                 for i in range(1, 4)]
        normalize(hinge)
    
    # Create parallel frames
    temp1 = rotation_matrix(hinge, -0.5 * rolltilt)
    para_bp2 = matrix_multiply(temp1, rot2)
    
    temp2 = rotation_matrix(hinge, 0.5 * rolltilt)
    para_bp1 = matrix_multiply(temp2, rot1)
    
    # Extract z-axis
    mstz = [para_bp2[i][3] for i in range(1, 4)]
    
    # Calculate twist
    t1_para = [para_bp1[i][2] for i in range(1, 4)]
    t2_para = [para_bp2[i][2] for i in range(1, 4)]
    Twist = angle_around_axis(t1_para, t2_para, mstz)
    
    # Calculate midpoint y-axis
    msty = rotate_vector((t1_para + t2_para) / 2, mstz, Twist/2)
    normalize(msty)
    
    # Calculate midpoint x-axis
    mstx = cross(msty, mstz)
    normalize(mstx)
    
    # Build midpoint frame
    mst_orien = [mstx, msty, mstz]  # Columns
    mst_org = [(org1[i] + org2[i]) / 2 for i in range(1, 4)]
    
    # Calculate translations
    dorg = [org2[i] - org1[i] for i in range(1, 4)]
    Shift = dot(dorg, mstx)
    Slide = dot(dorg, msty)
    Rise = dot(dorg, mstz)
    
    # Calculate rotations
    phi = angle_in_plane(hinge, msty, mstz)  # in radians
    Roll = rolltilt * cos(phi)
    Tilt = rolltilt * sin(phi)
    
    pars = [Shift, Slide, Rise, Tilt, Roll, Twist]
    return pars, mst_orien, mst_org
```

### Typical Values

- **B-DNA**: 
  - Shift ≈ 0.0 Å
  - Slide ≈ -0.5 Å
  - Rise ≈ 3.4 Å
  - Tilt ≈ 0°
  - Roll ≈ 0°
  - Twist ≈ 36°

- **A-DNA**: 
  - Rise ≈ 2.6 Å
  - Twist ≈ 33°

- **Z-DNA**: 
  - Roll ≈ 7°
  - Twist ≈ -10° (left-handed, negative twist)

---

## 3. H-Bond Detection Workflow

### Three-Phase Process

#### Phase 1: Initial Detection

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

**Criteria**:
- `good_hbatoms()`: Atom types must be valid (N, O, etc.)
- `within_limits()`: Distance must be in [hb_lower, hb_dist1] = [1.8, 4.0] Å

#### Phase 2: Conflict Resolution

See [H-Bond Conflict Resolution](#4-h-bond-conflict-resolution) below.

#### Phase 3: Validation

```c
validate_hbonds(num_hbonds, hb_dist, lkg_type, hb_type, basei, basej,
                hb_atom1, hb_atom2);
```

**Process**:
1. Initialize all types to `' '`
2. For conflicted H-bonds (negative distance):
   - Convert to positive: `hb_dist[k] = fabs(hb_dist[k])`
   - Determine type: `hb_type[k] = donor_acceptor(...)`
3. Filter invalid H-bonds:
   - Remove if distance > 3.6 Å
   - Remove if wrong type and distance not in [2.6, 3.2] Å
4. Count valid: `count(hb_type != ' ')`

---

## 4. H-Bond Conflict Resolution

### Problem

When an atom forms multiple H-bonds, we need to determine which is the "best" one. The algorithm uses an iterative approach to mark conflicts.

### Algorithm

#### Phase 1: Conflict Detection

```c
matched_idx[] = zeros(num_hbonds);  // Track processed H-bonds
num_iter = 0;

while (not all H-bonds processed) {
    // Find best H-bond for atom1 (by distance)
    best_idx1 = find_best_hbond_for_atom(atom1, hb_atom1, hb_dist, matched_idx);
    
    // Find best H-bond for atom2 (by distance)
    best_idx2 = find_best_hbond_for_atom(atom2, hb_atom2, hb_dist, matched_idx);
    
    // If same H-bond is best for both atoms: conflict!
    if (best_idx1 == best_idx2 && best_idx1 > 0) {
        hb_dist[best_idx1] = -hb_dist[best_idx1];  // Mark conflict (negate)
        // Mark all H-bonds sharing atoms as "matched"
        mark_related_hbonds(best_idx1, matched_idx);
        num_iter = 0;  // Reset iteration
    }
    
    num_iter++;
    if (num_iter > MAX_ITER) break;  // Safety check
}
```

**Key Points**:
- Conflicts marked by **negating distances** (negative = conflict)
- Iterative until convergence
- Best H-bond = closest distance

#### Phase 2: Linkage Type Calculation

```c
idx2[][] = zeros(num_hbonds, 2);

// For conflicted H-bonds
for (k = 1; k <= num_hbonds; k++) {
    if (hb_dist[k] <= 0.0) {  // Conflicted
        idx2[k][1] = 9;  // Mark atom1 as conflicted
        idx2[k][2] = 9;  // Mark atom2 as conflicted
        
        // Mark non-conflicted H-bonds that share atoms
        for (n = 1; n <= num_hbonds; n++) {
            if (hb_dist[n] > 0.0) {  // Not conflicted
                if (same_atom(hb_atom1[k], hb_atom1[n])) {
                    idx2[n][1] = 1;  // Shares atom1
                }
                if (same_atom(hb_atom2[k], hb_atom2[n])) {
                    idx2[n][2] = 1;  // Shares atom2
                }
            }
        }
    }
}

// Calculate linkage type
for (k = 1; k <= num_hbonds; k++) {
    lkg_type[k] = idx2[k][1] + idx2[k][2];
}
```

**Linkage Types**:
- `lkg_type = 18`: No conflicts (both atoms unique: 9 + 9)
- `lkg_type < 18`: Conflict detected (at least one atom shared)

#### Phase 3: Additional Conflict Marking

```c
// Only if hb_dist2 > 0 (but it's always 0.0 in practice)
if (hb_dist2 > 0.0) {
    for (k = 1; k <= num_hbonds; k++) {
        if (lkg_type[k] != 18 && 
            hb_dist[k] >= hb_lower && hb_dist[k] <= hb_dist2) {
            hb_dist[k] = -hb_dist[k];  // Mark as conflict
        }
    }
}
```

**Note**: Since `hb_dist2 = 0.0`, this phase is effectively disabled.

### Complete Pseudocode

```python
def hb_atompair(num_hbonds, hb_atom1, hb_atom2, hb_dist, lkg_type, misc_pars):
    matched_idx = [0] * (num_hbonds + 1)  # 1-indexed
    num_iter = 0
    MAX_ITER = 1000
    
    # Phase 1: Conflict detection
    while True:
        all_matched = True
        for k in range(1, num_hbonds + 1):
            if matched_idx[k] == 0 and hb_dist[k] > 0:
                all_matched = False
                break
        
        if all_matched:
            break
        
        # Find best H-bond for each atom
        best_for_atom1 = {}
        best_for_atom2 = {}
        
        for k in range(1, num_hbonds + 1):
            if matched_idx[k] or hb_dist[k] <= 0:
                continue
            
            atom1 = hb_atom1[k]
            atom2 = hb_atom2[k]
            dist = hb_dist[k]
            
            # Track best for atom1
            if atom1 not in best_for_atom1 or dist < best_for_atom1[atom1][1]:
                best_for_atom1[atom1] = (k, dist)
            
            # Track best for atom2
            if atom2 not in best_for_atom2 or dist < best_for_atom2[atom2][1]:
                best_for_atom2[atom2] = (k, dist)
        
        # Check for conflicts
        conflict_found = False
        for atom, (k1, dist1) in best_for_atom1.items():
            if atom in best_for_atom2:
                k2, dist2 = best_for_atom2[atom]
                if k1 == k2:  # Same H-bond is best for both atoms
                    # Conflict!
                    hb_dist[k1] = -hb_dist[k1]  # Negate distance
                    # Mark related H-bonds
                    mark_related(k1, matched_idx, hb_atom1, hb_atom2)
                    conflict_found = True
                    num_iter = 0  # Reset
                    break
        
        if not conflict_found:
            # Mark all best H-bonds as processed
            for atom, (k, dist) in best_for_atom1.items():
                matched_idx[k] = 1
            for atom, (k, dist) in best_for_atom2.items():
                matched_idx[k] = 1
        
        num_iter += 1
        if num_iter > MAX_ITER:
            break
    
    # Phase 2: Linkage type calculation
    idx2 = [[0, 0] for _ in range(num_hbonds + 1)]  # 1-indexed
    
    for k in range(1, num_hbonds + 1):
        if hb_dist[k] <= 0.0:  # Conflicted
            idx2[k][0] = 9  # Atom1 conflicted
            idx2[k][1] = 9  # Atom2 conflicted
            
            # Mark non-conflicted H-bonds sharing atoms
            for n in range(1, num_hbonds + 1):
                if hb_dist[n] > 0.0:  # Not conflicted
                    if hb_atom1[k] == hb_atom1[n]:
                        idx2[n][0] = 1
                    if hb_atom2[k] == hb_atom2[n]:
                        idx2[n][1] = 1
    
    # Calculate linkage types
    for k in range(1, num_hbonds + 1):
        lkg_type[k] = idx2[k][0] + idx2[k][1]
    
    # Phase 3: Additional conflict marking (disabled if hb_dist2 = 0.0)
    if misc_pars.hb_dist2 > 0.0:
        for k in range(1, num_hbonds + 1):
            if (lkg_type[k] != 18 and 
                hb_dist[k] >= misc_pars.hb_lower and 
                hb_dist[k] <= misc_pars.hb_dist2):
                hb_dist[k] = -hb_dist[k]
```

### Example

```
Initial H-bonds found:
1. N3 (res1) - O2 (res2), distance = 2.91 Å
2. N6 (res1) - O4 (res2), distance = 2.87 Å
3. N3 (res1) - O4 (res2), distance = 3.2 Å  # Conflict with #1

After conflict resolution:
- H-bond #1: distance = -2.91 (conflicted, N3 has two H-bonds)
- H-bond #2: distance = 2.87 (kept, best for N6 and O4)
- H-bond #3: distance = -3.2 (conflicted)

Linkage types:
- H-bond #1: lkg_type = 9 + 9 = 18 (but conflicted, so distance negative)
- H-bond #2: lkg_type = 9 + 9 = 18 (no conflicts)
- H-bond #3: lkg_type = 9 + 9 = 18 (but conflicted)
```

---

## 5. Overlap Area Calculation

### Problem

Calculate the overlap area between two base rings when projected onto a plane perpendicular to the average z-axis.

### Algorithm

#### Step 1: Get Ring Atom Coordinates

```c
n1 = ratom_xyz(ring_atom[r1], only_ring, xyz, oave, oxyz1);
n2 = ratom_xyz(ring_atom[r2], only_ring, xyz, oave, oxyz2);
```

**Returns**: Number of ring atoms and their coordinates relative to `oave`.

#### Step 2: Align to Z-axis

```c
align2zaxis(n1, zave, rotmat, oxyz1, oxyz1Z);
align2zaxis(n2, zave, rotmat, oxyz2, oxyz2Z);
```

**Purpose**: Rotate coordinates so that `zave` aligns with the z-axis, making projection to xy-plane straightforward.

#### Step 3: Project to XY Plane

```c
for (i = 1; i <= n1; i++) {
    polygon_a[i-1].x = oxyz1Z[i][1];  // x-coordinate
    polygon_a[i-1].y = oxyz1Z[i][2];  // y-coordinate
    // z-coordinate ignored (projection)
}

for (i = 1; i <= n2; i++) {
    polygon_b[i-1].x = oxyz2Z[i][1];
    polygon_b[i-1].y = oxyz2Z[i][2];
}
```

#### Step 4: Calculate Polygon Intersection

```c
overlap_area = pia_inter(polygon_a, n1, polygon_b, n2);
```

**Algorithm** (`pia_inter`):
- Uses computational geometry for convex polygon intersection
- Calculates overlapping area in 2D
- Returns area in Å²

### Threshold

- `OVERLAP = 0.01` Å²
- If `overlap_area >= 0.01`: Bases are too close (reject pair)
- Used in three places in `check_pair()`:
  1. Overlap in average plane (`only_ring = 0`)
  2. Overlap using base 1's origin and z-axis
  3. Overlap using base 2's origin and z-axis

---

## 6. Quality Score Calculation

### Formula

```
quality_score = dorg + 2*dv + plane_angle/20 + adjust_pairQuality()
```

**Components**:
- `dorg`: Origin distance (Å)
- `dv`: Vertical distance component (Å)
- `plane_angle`: Angle between base planes (degrees)
- `adjust_pairQuality()`: H-bond bonus (negative value)

### H-Bond Adjustment

```c
adjust_pairQuality() returns:
  -3.0 if num_good_hb >= 2  (2+ good H-bonds)
  -num_good_hb otherwise     (1 good H-bond = -1.0, none = 0.0)
```

**Good H-bond criteria**:
- Not conflicted (`num_list[k][0] == 0`)
- Distance in [2.5, 3.5] Å

### Effect

- **Lower score = better pair** (negative adjustment improves quality)
- Two good H-bonds: score reduced by 3.0
- One good H-bond: score reduced by 1.0
- No good H-bonds: no adjustment

### Example

```
Pair with:
  dorg = 10.5 Å
  dv = 0.3 Å
  plane_angle = 15°
  2 good H-bonds

quality_score = 10.5 + 2*0.3 + 15/20 + (-3.0)
              = 10.5 + 0.6 + 0.75 - 3.0
              = 8.85
```

---

## 7. Greedy Matching Algorithm

### Problem

Find the best set of base pairs such that:
1. Each base pairs with at most one other base
2. Pairing is mutual (i's best is j AND j's best is i)
3. Quality is maximized (score minimized)

### Algorithm

```c
matched_idx[] = zeros(num_residue);
num_bp = 0;

while (matched_count increases) {
    old_count = count(matched_idx);
    
    for (i = 1; i <= num_residue; i++) {
        if (RY[i] < 0 || matched_idx[i]) continue;
        
        // Find best partner for base i
        best_pair(i, ..., pair_istat);
        // Returns: pair_istat[1] = j (partner), pair_istat[2] = bpid
        
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

### Properties

- **Greedy**: Always picks best available partner
- **Mutual**: Requires both bases to choose each other
- **Iterative**: Continues until no new pairs found
- **Order-dependent**: Processes residues in order (1..num_residue)

### Example

```
Residues: 1, 2, 3, 4, 5

Iteration 1:
  Residue 1's best: 3 (score 8.5)
  Residue 3's best: 1 (score 8.5) → Match! (1,3)
  
  Residue 2's best: 4 (score 9.2)
  Residue 4's best: 2 (score 9.2) → Match! (2,4)
  
  Residue 5: no valid partner

Result: 2 pairs found: (1,3), (2,4)
```

---

**Next**: [Helper Functions](05_HELPER_FUNCTIONS.md) for utility function reference


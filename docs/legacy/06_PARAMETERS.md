# Legacy Parameters and Constants Reference

**Date**: 2025-01-XX  
**Purpose**: Complete reference for all parameters, constants, and configuration values  
**Status**: Essential for exact matching

---

## Table of Contents

1. [Validation Parameters](#validation-parameters)
2. [H-Bond Parameters](#h-bond-parameters)
3. [Constants and Macros](#constants-and-macros)
4. [Configuration Loading](#configuration-loading)
5. [Parameter Validation Rules](#parameter-validation-rules)

---

## Validation Parameters

### Base Pair Validation Parameters

All parameters stored in `miscPars` structure (`org/include/x3dna.h:29-51`).

#### Distance Parameters

```c
double min_dorg = 0.0;      // Minimum origin distance (Å)
double max_dorg = 15.0;     // Maximum origin distance (Å)
double min_dv = 0.0;        // Minimum vertical distance component (Å)
double max_dv = 2.5;        // Maximum vertical distance component (Å)
double min_dNN = 4.5;       // Minimum N-N distance (Å)
double max_dNN = XBIG;      // Maximum N-N distance (Å) ≈ 1.0e+18
```

**Usage in `check_pair()`**:
- `rtn_val[1]` (dorg): Must be in [min_dorg, max_dorg]
- `rtn_val[2]` (dv): Must be in [min_dv, max_dv]
- `rtn_val[4]` (dNN): Must be in [min_dNN, max_dNN]

#### Angle Parameters

```c
double min_plane_angle = 0.0;   // Minimum plane angle (degrees)
double max_plane_angle = 65.0;  // Maximum plane angle (degrees)
```

**Usage in `check_pair()`**:
- `rtn_val[3]` (plane_angle): Must be ≤ max_plane_angle

**Physical Meaning**: Angle between base planes. Larger angles indicate bases are not coplanar.

#### H-Bond Requirements

```c
long min_base_hb = 1;  // Minimum number of base-base H-bonds required
```

**Usage in `check_pair()`**:
- Count H-bonds between base atoms (not backbone)
- Require: `num_base_hb >= min_base_hb`

#### Overlap Threshold

```c
#define OVERLAP 0.01  // Overlap area threshold (Å²)
```

**Usage in `check_pair()`**:
- Calculate overlap area: `overlap = get_oarea(...)`
- If `overlap >= OVERLAP`: Reject pair (bases too close)

---

## H-Bond Parameters

### Distance Limits

```c
double hb_lower = 1.8;   // Minimum H-bond distance (Å)
double hb_dist1 = 4.0;   // Maximum H-bond distance (Å)
double hb_dist2 = 0.0;   // Conflict resolution threshold (Å) - CRITICAL: Must be 0.0
```

**Usage**:
- `hb_lower`, `hb_dist1`: Used in `within_limits()` for initial H-bond detection
- `hb_dist2`: Used in `hb_atompair()` Phase 3 conflict marking (effectively disabled when 0.0)

**CRITICAL**: `hb_dist2 = 0.0` is essential for matching legacy H-bond types.

### Atom Type Specification

```c
char hb_atoms[BUF512] = ".O.N";  // Valid H-bond atom types
long hb_idx[BUF512];              // Indices for valid types (populated from hb_atoms)
```

**Format**: String with atom type symbols (e.g., ".O.N" means Oxygen and Nitrogen)

**Usage in `good_hbatoms()`**:
- Checks if atom type indices are in `hb_idx[]`
- Typically includes: N (idx=1), O (idx=2)

### H-Bond Quality Criteria

**Good H-bond distance range**: [2.5, 3.5] Å

**Used in `adjust_pairQuality()`**:
- Counts H-bonds with distance in [2.5, 3.5] Å
- Not conflicted (`num_list[k][0] == 0`)

---

## Constants and Macros

### Indexing Constants

```c
#define NR_END 1        // Offset for 1-based indexing
#define TRUE 1L         // Boolean true
#define FALSE 0L        // Boolean false
#define NO_MATCH -1L    // No match found
#define DUMMY -1L       // Dummy/invalid value
```

### Buffer Sizes

```c
#define BUF32 32        // Small buffer
#define BUF512 512      // Standard buffer
#define BUF1K 1024      // 1KB buffer
#define BUF2K 2048      // 2KB buffer
#define BUFBIG 8192     // Large buffer
```

### Numerical Constants

```c
#define XBIG 1.0e+18        // Very large number (effectively infinity)
#define XBIG_CUTOFF 1.0e+16 // Large number cutoff
#define XEPS 1.0e-7         // Small epsilon for floating-point comparisons
#define PI 3.141592653589793 // Pi constant
#define OVERLAP 0.01        // Overlap threshold (Å²)
#define MFACTOR 10000.0     // Multiplication factor for integer storage
```

**Usage**:
- `XBIG`: Used for `max_dNN` (effectively no upper limit)
- `XEPS`: Used for degenerate case detection (e.g., parallel vectors)
- `OVERLAP`: Base overlap threshold
- `MFACTOR`: Used to store distances as integers (distance × MFACTOR)

### Array Size Constants

```c
#define NUM_BASE_ATOMS BUF512           // Max atoms per base
#define NUM_RESIDUE_ATOMS BUF512        // Max atoms per residue
#define NUM_DINUCLEOTIDE_ATOMS BUFBIG   // Max atoms per dinucleotide
#define MAXBASE 30000                   // Maximum number of bases
#define MAXCH 100                       // Maximum number of chains
```

### Ring Atom List

```c
#define RA_LIST " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
```

**Usage**: Ring atoms used for frame calculation and overlap detection
- First 6: Pyrimidine ring (C, T, U)
- All 9: Purine ring (A, G)

### Base Lists

```c
#define CB_LIST "ACGITU"        // Valid base characters
#define CX_LIST "ACGITUX"       // Valid base characters (includes X for unknown)
#define WC_LIST "XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"  // Watson-Crick pairs
```

### Quality Score Constants

```c
#define RTNNUM 37       // Number of return values
#define PSTNUM 29       // Number of pair statistics
```

**Usage**:
- `rtn_val[RTNNUM]`: Array for validation return values
- `pair_stat[PSTNUM]`: Array for pair statistics

### Other Constants

```c
#define NMISC 34        // Number of miscellaneous fields
#define BOND_UPPER_LIMIT 2.5    // Upper limit for bond distance
#define O3P_UPPER 2.5           // Upper limit for O3'-P distance
#define helix_break = 7.5       // Helix break distance (Å)
#define std_curved = 0.6        // Standard curvature
#define o3p_dist = 4.5          // O3'-P distance threshold
```

---

## Configuration Loading

### Parameter File

**Location**: `{X3DNA_HOMEDIR}/misc_3dna.par`

**Format**: XML-like tags

**Example**:
```xml
<min_base_hb>1</min_base_hb>
<hb_lower>1.8</hb_lower>
<hb_dist1>4.0</hb_dist1>
<hb_dist2>0.0</hb_dist2>
<hb_atoms>.O.N</hb_atoms>
<alt_list>A1</alt_list>
<max_dorg>15.0</max_dorg>
<min_dorg>0.0</min_dorg>
<max_dv>2.5</max_dv>
<min_dv>0.0</min_dv>
<max_plane_angle>65.0</max_plane_angle>
<min_plane_angle>0.0</min_plane_angle>
<max_dNN>1.0e+18</max_dNN>
<min_dNN>4.5</min_dNN>
```

### Loading Process

1. **Set Defaults**: `set_default_misc_pars(misc_pars)`
2. **Read File**: `get_3dna_pars(misc_pars)` reads from `misc_3dna.par`
3. **Validate**: `check_misc_pars(misc_pars)` validates all parameters
4. **Mask Options**: `mask_off_options(misc_pars)` if DEBUG < DEBUG_LEVEL

### Command-Line Override

Parameters can be overridden via command-line:

```bash
find_pair -min_base_hb 2 -hb_lower 2.0 -max_dorg 12.0 ...
```

**Function**: `overwrite_misc_pars(option)` in `app_fncs.c:819-870`

---

## Parameter Validation Rules

### Validation Function

`check_misc_pars()` in `org/src/app_fncs.c:708-762` validates all parameters.

### Validation Rules

#### min_base_hb

```c
if (min_base_hb < 0 || min_base_hb > 3) {
    // Reset to default (1)
    min_base_hb = 1;
}
```

**Valid Range**: [0, 3]

#### Distance Parameters

```c
// All distance parameters must be >= 0.0
if (hb_lower < 0.0) fatal(...);
if (hb_dist1 < 0.0) fatal(...);
if (hb_dist2 < 0.0) fatal(...);
if (max_dorg < 0.0) fatal(...);
if (min_dorg < 0.0) fatal(...);
if (max_dv < 0.0) fatal(...);
if (min_dv < 0.0) fatal(...);
if (max_dNN < 0.0) fatal(...);
if (min_dNN < 0.0) fatal(...);
```

#### Range Validation

```c
// min must be <= max
if (min_dorg > max_dorg) fatal(...);
if (min_dv > max_dv) fatal(...);
if (min_plane_angle > max_plane_angle) fatal(...);
if (min_dNN > max_dNN) fatal(...);
```

#### Angle Validation

```c
// Plane angles must be in [0, 90] degrees
if (max_plane_angle < 0.0 || max_plane_angle > 90.0) fatal(...);
if (min_plane_angle < 0.0 || min_plane_angle > 90.0) fatal(...);
```

#### String Parameters

```c
// hb_atoms: Reset atom type indices
reset_alist_symbol_idx(misc_pars->hb_atoms, misc_pars->hb_idx);

// alt_list: Check for spaces
check_space_in_altlist(misc_pars->alt_list);
```

---

## Complete Parameter Structure

### miscPars Structure

```c
typedef struct {
    // H-Bond parameters
    long min_base_hb;           // Minimum H-bonds required
    double hb_lower;             // Minimum H-bond distance
    double hb_dist1;             // Maximum H-bond distance
    double hb_dist2;             // Conflict resolution threshold (CRITICAL: 0.0)
    char hb_atoms[BUF512];       // Valid H-bond atom types
    long hb_idx[BUF512];         // Atom type indices (populated)
    
    // Alternative location
    char alt_list[BUF512];       // Alternative location list
    
    // Distance parameters
    double max_dorg;             // Maximum origin distance
    double min_dorg;             // Minimum origin distance
    double max_dv;               // Maximum vertical distance
    double min_dv;               // Minimum vertical distance
    
    // Angle parameters
    double max_plane_angle;      // Maximum plane angle
    double min_plane_angle;      // Minimum plane angle
    
    // N-N distance
    double max_dNN;              // Maximum N-N distance
    double min_dNN;              // Minimum N-N distance
    
    // Helix parameters
    double helix_break;          // Helix break distance
    double std_curved;           // Standard curvature
    
    // Water parameters
    double water_dist;           // Water distance threshold
    double water_dlow;            // Water lower distance
    char water_atoms[BUF512];    // Water atom types
    long water_idx[BUF512];      // Water atom indices
    
    // O3'-P distance
    double o3p_dist;             // O3'-P distance threshold
} miscPars;
```

---

## Critical Parameters for Matching

### Must Match Exactly

These parameters are **critical** for matching legacy behavior:

1. **`hb_dist2 = 0.0`**: 
   - Affects H-bond conflict resolution
   - Must be exactly 0.0 (not 0.1, not 4.5)

2. **`OVERLAP = 0.01`**:
   - Base overlap threshold
   - Used in three overlap checks

3. **Distance ranges**:
   - `hb_lower = 1.8`, `hb_dist1 = 4.0`
   - `min_dorg = 0.0`, `max_dorg = 15.0`
   - `min_dv = 0.0`, `max_dv = 2.5`
   - `min_dNN = 4.5`, `max_dNN = XBIG`

4. **Angle limits**:
   - `max_plane_angle = 65.0`

5. **H-bond requirement**:
   - `min_base_hb = 1`

### Parameter File Location

**Default**: `{X3DNA_HOMEDIR}/misc_3dna.par`

**Modern code must**:
- Use same default values
- Load from same parameter file (if exists)
- Allow command-line override (same format)

---

## Default Values Summary

```cpp
// Validation Parameters
min_dorg = 0.0
max_dorg = 15.0
min_dv = 0.0
max_dv = 2.5
min_dNN = 4.5
max_dNN = 1.0e+18  // XBIG
min_plane_angle = 0.0
max_plane_angle = 65.0
min_base_hb = 1

// H-Bond Parameters
hb_lower = 1.8
hb_dist1 = 4.0
hb_dist2 = 0.0  // CRITICAL
hb_atoms = ".O.N"

// Other
alt_list = "A1"
helix_break = 7.5
std_curved = 0.6
water_dist = 3.2
water_dlow = 0.0
water_atoms = ".O.N"
o3p_dist = 4.5

// Constants
OVERLAP = 0.01
XBIG = 1.0e+18
XEPS = 1.0e-7
```

---

**Next**: [Workflows](07_WORKFLOWS.md) for execution flow diagrams


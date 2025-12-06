# Code Flow and Architecture

**Last Updated**: 2025-01-XX  
**Purpose**: Detailed explanation of how the modern code works and flows from input to output

---

## Overview

The modernized X3DNA code follows a two-phase architecture:
1. **find_pair**: Base pair identification and selection
2. **analyze**: Step parameter calculation

Both phases use protocol classes to orchestrate the workflow, with algorithm classes performing the actual calculations.

---

## Phase 1: find_pair - Base Pair Identification

### Entry Point

**Application**: `apps/find_pair_app.cpp`
- Parses command-line arguments
- Creates `FindPairProtocol`
- Executes protocol
- Writes output file

**Tool**: `tools/generate_modern_json.cpp`
- Similar to find_pair_app but focuses on JSON generation
- Supports `--fix-indices` option for comparison with legacy

### Protocol Flow

**Class**: `FindPairProtocol` (`src/x3dna/protocols/find_pair_protocol.cpp`)

```
execute(Structure& structure)
  ├─> Fix residue indices (if --fix-indices enabled)
  │   └─> fix_residue_indices_from_json()
  │       → Loads legacy JSON and matches residues by PDB properties
  │
  ├─> calculate_frames(structure)
  │   └─> BaseFrameCalculator::calculate_all_frames()
  │       → For each residue:
  │         ├─> Detect RNA (check for O2' atoms)
  │         ├─> Load standard template (Atomic_A.pdb, etc.)
  │         ├─> Match ring atoms
  │         ├─> Least-squares fitting
  │         └─> Store reference frame on residue
  │
  ├─> find_pairs(structure)
  │   └─> BasePairFinder::find_pairs() or find_pairs_with_recording()
  │       → Greedy matching algorithm:
  │         ├─> For each unmatched base i:
  │         │   ├─> Try all other bases j
  │         │   ├─> BasePairValidator::validate(i, j)
  │         │   │   ├─> Geometric checks
  │         │   │   ├─> Overlap calculation
  │         │   │   ├─> H-bond detection
  │         │   │   └─> Quality score calculation
  │         │   └─> Select best pair (minimum quality score)
  │         └─> Check mutual selection (best_pair(i) == j && best_pair(j) == i)
  │
  └─> JSON recording (if JsonWriter provided)
      → Records all intermediate results to JSON files
```

### Frame Calculation

**Class**: `BaseFrameCalculator` (`src/x3dna/algorithms/base_frame_calculator.cpp`)

```
calculate_frame(Residue& residue)
  ├─> Load standard template
  │   └─> PdbParser::parse(Atomic_A.pdb, etc.)
  │
  ├─> Match ring atoms
  │   └─> RingAtomMatcher::match()
  │       ├─> Get ring atom names for residue type
  │       ├─> Find matching atoms in experimental structure
  │       └─> Return matched atom pairs
  │
  ├─> Least-squares fitting
  │   └─> ls_fitting()
  │       ├─> Calculate centroids
  │       ├─> Center coordinates
  │       ├─> SVD-based rotation matrix calculation
  │       └─> Calculate translation vector
  │
  └─> Store reference frame
      └─> Residue::set_reference_frame()
          → Stores rotation matrix (3×3) and origin (3D point)
```

### Base Pair Finding

**Class**: `BasePairFinder` (`src/x3dna/algorithms/base_pair_finder.cpp`)

```
find_pairs(Structure& structure)
  ├─> Initialize matched flags for all residues
  │
  ├─> Iterate until no new pairs found:
  │   ├─> For each unmatched base i:
  │   │   ├─> best_pair(i)
  │   │   │   └─> Try all other bases j:
  │   │   │       ├─> BasePairValidator::validate(i, j)
  │   │   │       │   ├─> Calculate geometric parameters
  │   │   │       │   │   ├─> dorg (origin distance)
  │   │   │       │   │   ├─> dNN (N-N distance)
  │   │   │       │   │   ├─> plane_angle (angle between planes)
  │   │   │       │   │   └─> d_v (vertical distance)
  │   │   │       │   │
  │   │   │       │   ├─> Check geometric thresholds
  │   │   │       │   │   └─> dorg ≤ max_dorg, dNN in range, etc.
  │   │   │       │   │
  │   │   │       │   ├─> Calculate overlap area
  │   │   │       │   │   └─> calculate_overlap_area()
  │   │   │       │   │       ├─> Get ring atoms + exocyclic atoms
  │   │   │       │   │       ├─> Project to 2D plane
  │   │   │       │   │       └─> Polygon intersection
  │   │   │       │   │
  │   │   │       │   ├─> Check overlap threshold (< 0.01)
  │   │   │       │   │
  │   │   │       │   ├─> Detect hydrogen bonds
  │   │   │       │   │   └─> HydrogenBondFinder::find_hbonds()
  │   │   │       │   │       ├─> Check distance (≤ 4.0 Å)
  │   │   │       │   │       ├─> Check angle (H-donor-acceptor)
  │   │   │       │   │       └─> Check atom types (O, N)
  │   │   │       │   │
  │   │   │       │   ├─> Resolve H-bond conflicts
  │   │   │       │   │   └─> resolve_hbond_conflicts()
  │   │   │       │   │
  │   │   │       │   ├─> Calculate quality score
  │   │   │       │   │   └─> adjust_pairQuality()
  │   │   │       │   │       ├─> Base: dorg + 2.0 * d_v + plane_angle / 20.0
  │   │   │       │   │       ├─> H-bond adjustment: -3.0 (2+ good), -1.0 per good
  │   │   │       │   │       └─> Watson-Crick adjustment: -2.0 (if bp_type_id == 2)
  │   │   │       │   │
  │   │   │       │   └─> Calculate bp_type_id
  │   │   │       │       └─> calculate_bp_type_id()
  │   │   │       │           ├─> Check direction vectors
  │   │   │       │           ├─> Check Watson-Crick/Wobble (if step params available)
  │   │   │       │           └─> Return: -1 (unclassified), 0 (invalid), 1 (wobble), 2 (WC)
  │   │   │       │
  │   │   │       └─> Return validation result
  │   │   │
  │   │   └─> Select pair with minimum quality score
  │   │
  │   └─> Check mutual selection
  │       └─> If best_pair(i) == j && best_pair(j) == i:
  │           └─> Mark both as matched, add to base_pairs
  │
  └─> Return base_pairs vector
```

### JSON Recording

**Class**: `JsonWriter` (`src/x3dna/io/json_writer.cpp`)

During pair finding, the following records are written:

1. **pdb_atoms** - All atoms from PDB file
2. **base_frame_calc** - Frame calculation metadata (template matching, RMS fit)
3. **frame_calc** - Reference frames (rotation matrix + origin)
4. **pair_validation** - Validation results for all tested pairs
5. **distance_checks** - Geometric measurements (dorg, dNN, plane_angle, d_v, overlap_area)
6. **hbond_list** - Hydrogen bond lists
7. **find_bestpair_selection** - Final selected pairs (PRIMARY OUTPUT)
8. **base_pair** - Base pair records (only for selected pairs)

---

## Phase 2: analyze - Step Parameter Calculation

### Entry Point

**Application**: `apps/analyze_app.cpp`
- Parses command-line arguments
- Creates `AnalyzeProtocol`
- Sets up `JsonWriter` for step parameter recording
- Executes protocol
- Outputs step parameters

### Protocol Flow

**Class**: `AnalyzeProtocol` (`src/x3dna/protocols/analyze_protocol.cpp`)

```
execute(input_file)
  ├─> Parse .inp file
  │   └─> InputFileParser::parse()
  │       ├─> Read base pairs from .inp file
  │       ├─> Extract PDB file path
  │       └─> Return InputData structure
  │
  ├─> Load PDB structure
  │   └─> PdbParser::parse()
  │       → Loads atoms, residues, chains
  │
  ├─> Convert atom indices to residue indices (if needed)
  │   └─> convert_atom_indices_to_residue_indices()
  │       → Legacy input files may have atom indices
  │       → Converts to residue indices using structure mapping
  │
  └─> execute(structure)
      ├─> recalculate_frames(structure)
      │   └─> BaseFrameCalculator::calculate_all_frames()
      │       → Recalculates frames for all residues
      │       → Reuses frames from find_pair if available
      │
      └─> calculate_parameters(structure)
          └─> ParameterCalculator::calculate_step_parameters()
              → For each consecutive base pair:
                ├─> Extract frames (R1, O1), (R2, O2)
                ├─> Calculate midpoint frame
                ├─> Calculate 6 step parameters
                └─> Calculate helical parameters
```

### Step Parameter Calculation

**Class**: `ParameterCalculator` (`src/x3dna/algorithms/parameter_calculator.cpp`)

```
calculate_step_parameters(base_pairs, structure)
  ├─> For each consecutive pair (i, i+1):
  │   ├─> Extract frames
  │   │   └─> Get reference frames for residues in pairs
  │   │       → R1, O1 from base pair i
  │   │       → R2, O2 from base pair i+1
  │   │
  │   ├─> Calculate midpoint frame
  │   │   └─> bpstep_par_impl()
  │   │       ├─> Average rotation matrices
  │   │       ├─> Average origins
  │   │       └─> Create midpoint coordinate system
  │   │
  │   ├─> Calculate step parameters
  │   │   └─> Extract 6 parameters:
  │   │       ├─> Shift (x-displacement)
  │   │       ├─> Slide (y-displacement)
  │   │       ├─> Rise (z-displacement)
  │   │       ├─> Tilt (rotation about x)
  │   │       ├─> Roll (rotation about y)
  │   │       └─> Twist (rotation about z)
  │   │
  │   └─> Calculate helical parameters
  │       └─> Alternative parameter set using helical axis
  │
  └─> Store parameters
      └─> JsonWriter::record_bpstep_params() and record_helical_params()
```

### JSON Recording

**Class**: `JsonWriter` (`src/x3dna/io/json_writer.cpp`)

During analyze phase, the following records are written:

1. **bpstep_params** - Step parameters (Shift, Slide, Rise, Tilt, Roll, Twist)
2. **helical_params** - Helical parameters (x_displacement, y_displacement, rise, inclination, tip, twist)

---

## Data Flow

### Structure Representation

```
Structure
  ├─> Chains (vector)
  │   └─> Chain
  │       ├─> Residues (vector)
  │       │   └─> Residue
  │       │       ├─> Atoms (vector)
  │       │       │   └─> Atom
  │       │       │       ├─> Name, coordinates
  │       │       │       └─> Legacy index (for comparison)
  │       │       │
  │       │       └─> ReferenceFrame (optional)
  │       │           ├─> Rotation matrix (3×3)
  │       │           └─> Origin (3D point)
  │       │
  │       └─> Chain ID
  │
  └─> PDB ID
```

### Base Pair Representation

```
BasePair
  ├─> Residue indices (i, j)
  ├─> Quality score
  ├─> bp_type_id (-1, 0, 1, 2)
  ├─> Validation result
  │   ├─> is_valid
  │   ├─> Geometric parameters
  │   ├─> Overlap area
  │   └─> H-bonds
  └─> Reference frames (R1, O1), (R2, O2)
```

---

## Key Algorithms

### 1. Least-Squares Fitting

**Purpose**: Find rotation matrix and translation that best align standard template with experimental structure

**Algorithm**:
1. Calculate centroids of both point sets
2. Center coordinates
3. Calculate covariance matrix
4. SVD decomposition
5. Extract rotation matrix
6. Calculate translation vector

**Implementation**: `src/x3dna/algorithms/base_frame_calculator.cpp::ls_fitting()`

### 2. Overlap Calculation

**Purpose**: Calculate overlap area between two base planes

**Algorithm**:
1. Get ring atoms + one exocyclic atom per ring atom
2. Project atoms to 2D plane (base plane)
3. Create polygon from projected points
4. Calculate polygon intersection area
5. Return overlap area

**Implementation**: `src/x3dna/algorithms/base_pair_validator.cpp::calculate_overlap_area()`

### 3. Hydrogen Bond Detection

**Purpose**: Find hydrogen bonds between two bases

**Algorithm**:
1. Get potential donor/acceptor atoms (O, N)
2. Check distance (≤ 4.0 Å)
3. Check angle (H-donor-acceptor)
4. Resolve conflicts (one atom can only form one H-bond)

**Implementation**: `src/x3dna/algorithms/hydrogen_bond_finder.cpp`

### 4. Step Parameter Calculation

**Purpose**: Calculate 6 step parameters for consecutive base pairs

**Algorithm**:
1. Extract frames for both pairs
2. Calculate midpoint frame
3. Transform to midpoint coordinate system
4. Extract 6 parameters from transformation

**Implementation**: `src/x3dna/algorithms/parameter_calculator.cpp::bpstep_par_impl()`

---

## Indexing Conventions

### Modern Code (0-based)
- Array indices: `[0, n-1]`
- Residue indices: Start from 0
- Vector access: `vector[i]` where `i ∈ [0, size-1]`

### Legacy Code (1-based)
- Array indices: `[1, n]` (with NR_END offset)
- Residue indices: Start from 1
- Array access: `arr[i]` where `i ∈ [1, n]`

### Conversion
- Modern code stores legacy indices on atoms/residues for comparison
- `--fix-indices` option fixes indices from legacy JSON
- Input file parser converts atom indices to residue indices automatically

---

## Error Handling

### Residue Indexing
- If legacy JSON not found for `--fix-indices`, continues without fixing
- If JSON parse error, attempts auto-fix (missing closing brackets)
- Graceful degradation: continues processing even if index fixing fails

### Frame Calculation
- If template not found, skips residue
- If insufficient atoms matched, skips residue
- If fitting fails, skips residue

### Pair Validation
- If validation fails, pair is not selected
- All validation results recorded to JSON for debugging

---

## Performance Considerations

### Caching
- Frame calculations cached on residues
- Template loading cached

### Parallelization
- JSON comparison uses parallel processing
- Frame calculation could be parallelized (future work)

### Memory
- Structures loaded into memory
- JSON files written incrementally
- Large PDBs may require significant memory

---

## Testing Integration

### Unit Tests
- Test individual classes in isolation
- Location: `tests/unit/`

### Integration Tests
- Test component interactions
- Location: `tests/integration/`

### JSON Regression Tests
- Compare output with legacy JSON
- Tool: `scripts/compare_json.py`

---

*This document provides a detailed view of code flow and architecture. See [PROJECT_SUMMARY.md](PROJECT_SUMMARY.md) for high-level overview.*


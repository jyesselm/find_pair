# Legacy Code Architecture Overview

**Date**: 2025-01-XX  
**Purpose**: Complete architectural overview of legacy X3DNA codebase structure, organization, and execution flow  
**Status**: Comprehensive reference for understanding program structure

---

## Table of Contents

1. [Main Programs](#main-programs)
2. [File Organization](#file-organization)
3. [Module Structure](#module-structure)
4. [Execution Flows](#execution-flows)
5. [Build System](#build-system)
6. [Data Flow](#data-flow)
7. [Global State](#global-state)

---

## Main Programs

The legacy codebase consists of two main programs that work together:

### 1. `find_pair` - Base Pair Identification

**Purpose**: Identifies base pairs in nucleic acid structures

**Entry Point**: `find_pair_main()` in `org/src/find_pair.c`

**Complete Execution Flow**:

```
find_pair_main(argc, argv)
├─ Parse command line arguments (fp_cmdline)
├─ Initialize global variables (set_my_globals)
├─ handle_str(args)
│  ├─ 1. PDB File Reading
│  │   ├─ number_of_atoms(pdbfile)           # Count atoms
│  │   ├─ read_pdb(pdbfile, ...)             # Parse PDB format
│  │   │   ├─ Allocate arrays: AtomName, ResName, ChainID, ResSeq, xyz, Miscs
│  │   │   ├─ Parse ATOM/HETATM records
│  │   │   └─ Store atom coordinates and metadata
│  │   └─ Result: num atoms, atom arrays filled
│  │
│  ├─ 2. Residue Indexing
│  │   ├─ residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue)
│  │   │   ├─ Map atoms to residues
│  │   │   ├─ Create seidx[][] array: [residue][1]=start_atom, [residue][2]=end_atom
│  │   │   └─ Assign sequential residue indices (1, 2, 3, ...) to ALL residues
│  │   └─ Result: num_residue, seidx[][] mapping
│  │
│  ├─ 3. Sequence Extraction
│  │   ├─ get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY)
│  │   │   ├─ Extract base sequence: bseq[1..num_residue]
│  │   │   ├─ Classify residues: RY[] = 1 (purine), 0 (pyrimidine), -1 (not nucleotide)
│  │   │   └─ Identify base types: A, C, G, T, U, etc.
│  │   └─ Result: bseq[], RY[]
│  │
│  ├─ 4. Base Frame Calculation
│  │   ├─ base_info(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, orien, org, NC1xyz, o3_p)
│  │   │   ├─ base_frame(num_residue, bseq, seidx, RY, ..., xyz, BDIR, orien, org)
│  │   │   │   ├─ For each residue i in [1, num_residue]:
│  │   │   │   │   ├─ if RY[i] < 0: skip (not nucleotide)
│  │   │   │   │   ├─ Load template: Atomic_{bseq[i]}.pdb
│  │   │   │   │   ├─ Match ring atoms between template and experimental
│  │   │   │   │   ├─ ls_fitting(): Least-squares fitting
│  │   │   │   │   ├─ Store rotation matrix: mst2orien(orien[i], 0, R)
│  │   │   │   │   └─ Store origin: org[i] = translation
│  │   │   │   └─ Result: orien[][], org[][]
│  │   │   ├─ Extract N and C1' coordinates: NC1xyz[][]
│  │   │   └─ Extract O3' and P coordinates: o3_p[][]
│  │   └─ Result: orien[][], org[][], NC1xyz[][], o3_p[][]
│  │
│  ├─ 5. Ring Atom Indexing
│  │   ├─ ring_oidx(num, num_residue, RY, seidx, AtomName, xyz, idx, ring_atom)
│  │   │   ├─ For each residue: identify ring atoms
│  │   │   └─ Store indices in ring_atom[][]
│  │   └─ Result: ring_atom[][] for overlap calculations
│  │
│  ├─ 6. Base Pair Finding
│  │   ├─ find_bestpair(nout, base_pairs, num_residue, bseq, seidx, RY, AtomName, xyz, idx, orien, org, NC1xyz, ring_atom, misc_pars)
│  │   │   ├─ Initialize matched_idx[] = all zeros
│  │   │   ├─ while matched_count increases:
│  │   │   │   ├─ For each residue i:
│  │   │   │   │   ├─ if RY[i] < 0 or matched_idx[i]: continue
│  │   │   │   │   ├─ best_pair(i, ..., pair_istat)  # Find best partner
│  │   │   │   │   │   ├─ For each j != i:
│  │   │   │   │   │   │   ├─ check_pair(i, j, ..., rtn_val, &bpid, ..., 0)
│  │   │   │   │   │   │   └─ if valid and best score: update best
│  │   │   │   │   │   └─ Return: j, bpid, quality scores
│  │   │   │   │   ├─ if partner found:
│  │   │   │   │   │   ├─ best_pair(j, ..., pair_jstat)  # Check mutual
│  │   │   │   │   │   └─ if mutual: record pair
│  │   │   │   └─ Update matched_idx
│  │   │   └─ Result: num_bp, base_pairs[][]
│  │   └─ Result: base_pairs[][] with pair data
│  │
│  ├─ 7. Base Pair Reordering
│  │   ├─ re_ordering(num_bp, base_pairs, bp_idx, helix_marker, helix_idx, misc_pars, &num_helix, o3_p, ...)
│  │   │   ├─ bp_context(): Analyze neighbor relationships
│  │   │   ├─ locate_helix(): Group pairs into helical regions
│  │   │   ├─ five2three(): Reorder to 5'→3' direction
│  │   │   └─ check_zdna(): Special Z-DNA handling
│  │   └─ Result: bp_idx[], helix_idx[][], helix_marker[]
│  │
│  └─ 8. Output Generation
│     ├─ x3dna_input(outfile, ...)           # Write .inp file
│     └─ write_bestpairs(...)                # Write pair information
│
└─ Cleanup
   ├─ Free all allocated arrays
   ├─ clear_my_globals()
   └─ print_used_time(time0)
```

**Key Functions Called**:
- `read_pdb()`: PDB file parsing
- `residue_idx()`: Residue mapping
- `get_seq()`: Sequence extraction
- `base_frame()`: Frame calculation
- `find_bestpair()`: Pair finding
- `check_pair()`: Pair validation
- `re_ordering()`: Pair reordering

---

### 2. `analyze` - Parameter Calculation

**Purpose**: Calculates base pair and step parameters from identified pairs

**Entry Point**: `analyze_main()` in `org/src/analyze.c`

**Complete Execution Flow**:

```
analyze_main(argc, argv)
├─ Parse command line arguments (analyze_cmdline)
├─ Initialize global variables (set_my_globals)
├─ process_str(inpfile, args)
│  ├─ 1. Read Input File
│  │   ├─ read_input(inpfile, pdbfile, outfile, &ds, &num_bp, &ip, &hetatm)
│  │   │   ├─ Parse .inp file from find_pair
│  │   │   ├─ Extract: PDB file name, duplex status, pair numbers
│  │   │   └─ Load pair_num[][]: Residue indices for each base pair
│  │   └─ Result: pair_num[][], ds, num_bp
│  │
│  ├─ 2. Read PDB File (again)
│  │   ├─ number_of_atoms(pdbfile)
│  │   ├─ read_pdb(pdbfile, ...)
│  │   └─ Result: Atom arrays filled
│  │
│  ├─ 3. Residue Indexing (again)
│  │   ├─ residue_idx(...)
│  │   └─ Result: seidx[][]
│  │
│  ├─ 4. Base Pair Validation
│  │   ├─ pair_checking(ip, ds, num_residue, pdbfile, &num_bp, pair_num)
│  │   │   └─ Validate pair indices
│  │   └─ Result: Validated pair_num[][]
│  │
│  ├─ 5. Extract Base Sequences
│  │   ├─ get_bpseq(ds, num_bp, pair_num, ..., bp_seq, RY)
│  │   │   └─ Extract base letters for each pair
│  │   ├─ get_seq(...)  # Full sequence
│  │   └─ Result: bp_seq[][], bseq[], RY[]
│  │
│  ├─ 6. Recalculate Base Frames
│  │   ├─ ref_frames(ds, num_bp, pair_num, bp_seq, seidx, RY, ..., xyz, fp, orien, org, WC_info, ...)
│  │   │   ├─ For each base pair (i, j):
│  │   │   │   ├─ base_frame(): Calculate frame for base i
│  │   │   │   ├─ base_frame(): Calculate frame for base j
│  │   │   │   ├─ Check orientation: dir_z = dot(z_i, z_j)
│  │   │   │   ├─ If needed: reverse one frame
│  │   │   │   └─ Store: bp_orien[][], bp_org[][]
│  │   │   └─ check_Watson_Crick(): Classify pairs
│  │   └─ Result: orien[][], org[][], WC_info[]
│  │
│  ├─ 7. Calculate Step Parameters
│  │   ├─ get_parameters(ds, num_bp, bp_seq, orien, org, WC_info, fp, ...)
│  │   │   ├─ For each consecutive pair (i, i+1):
│  │   │   │   ├─ refs_i_j(i, i+1, orien, org, r1, o1, r2, o2)
│  │   │   │   │   └─ Extract frames from orien[][], org[][]
│  │   │   │   ├─ bpstep_par(r1, o1, r2, o2, pars, mst_orien, mst_org)
│  │   │   │   │   └─ Calculate: Shift, Slide, Rise, Tilt, Roll, Twist
│  │   │   │   ├─ helical_par(r1, o1, r2, o2, pars, ...)
│  │   │   │   │   └─ Calculate helical parameters
│  │   │   │   └─ Store: bp_step_par[][], bp_heli_par[][]
│  │   │   └─ Result: Step parameters for all pairs
│  │   └─ Result: bp_step_par[][], bp_heli_par[][]
│  │
│  └─ 8. Output Results
│     ├─ write_mst(...)                      # Write parameter file
│     └─ Various output files generated
│
└─ Cleanup
   ├─ Free all arrays
   ├─ clear_my_globals()
   └─ print_used_time(time0)
```

**Key Functions Called**:
- `read_input()`: Read .inp file
- `ref_frames()`: Recalculate frames
- `bpstep_par()`: Calculate step parameters
- `helical_par()`: Calculate helical parameters

---

### 3. `find_pair_analyze` - Combined Program

**Purpose**: Runs both find_pair and analyze in sequence

**Entry Point**: `main()` in `org/src/find_pair_analyze.c`

**Execution Flow**:
```
main(argc, argv)
├─ Parse arguments for both programs
├─ Step 1: Run find_pair
│   ├─ handle_str(&fp_args)
│   └─ Result: .inp file created
├─ Step 2: Run analyze
│   ├─ process_str(inpfile, &ana_args)
│   └─ Result: Parameter files created
└─ Cleanup
```

---

## File Organization

### Source Files (`org/src/`)

```
org/src/
├── find_pair.c              # Main find_pair program (~2100 lines)
│   ├── find_pair_main()     # Entry point
│   ├── handle_str()         # Main processing function
│   ├── duplex()             # Duplex structure processing
│   ├── find_bestpair()      # Greedy matching (static)
│   ├── best_pair()          # Find best partner (static)
│   ├── re_ordering()        # Reorder pairs (static)
│   └── Various helper functions
│
├── analyze.c                # Main analyze program (~600 lines)
│   ├── analyze_main()       # Entry point
│   ├── process_str()        # Main processing function
│   └── Helper functions
│
├── find_pair_analyze.c      # Combined program (~210 lines)
│   └── main()               # Combined entry point
│
├── cmn_fncs.c               # Common functions (~5000 lines)
│   ├── PDB I/O: read_pdb(), write_pdb(), etc.
│   ├── Residue: residue_idx(), get_seq(), residue_ident()
│   ├── Validation: check_pair(), get_hbond_ij(), validate_hbonds()
│   ├── H-bonds: hb_atompair(), good_hbatoms(), donor_acceptor()
│   ├── Geometry: get_oarea(), get_zoave(), etc.
│   ├── Math: ls_fitting(), cov_matrix(), jacobi(), etc.
│   ├── Vector ops: dot(), cross(), veclen(), etc.
│   └── Utilities: Various helper functions
│
├── ana_fncs.c               # Analyze-specific functions (~3600 lines)
│   ├── Parameter calc: bpstep_par(), helical_par()
│   ├── Frame calc: ref_frames(), check_Watson_Crick()
│   ├── Overlap: base_overlap()
│   └── Various analysis functions
│
├── app_fncs.c               # Application functions (~1800 lines)
│   ├── Frame calc: base_frame()
│   ├── Peptide: peptide_info(), peptide_frame()
│   ├── Configuration: set_default_misc_pars(), get_3dna_pars()
│   └── Various application helpers
│
├── nrutil.c                 # Numerical recipes utilities (~500 lines)
│   ├── Memory: dvector(), dmatrix(), lvector(), etc.
│   ├── Matrix ops: dinverse(), dludcmp(), dlubksb()
│   └── Utility functions
│
└── json_writer.c            # JSON output (~1700 lines)
    ├── json_writer_init()   # Initialize JSON writer
    ├── Various record_* functions
    └── JSON formatting helpers
```

### Header Files (`org/include/`)

```
org/include/
├── x3dna.h                  # Main header (~185 lines)
│   ├── Type definitions: miscPars, struct_Gvars
│   ├── Constants: BUF512, XBIG, OVERLAP, etc.
│   ├── Macros: NR_END, TRUE, FALSE, etc.
│   └── Includes x3dna_fncs.h
│
├── x3dna_fncs.h             # Function declarations (~610 lines)
│   └── All function prototypes
│
└── json_writer.h            # JSON writer interface (~100 lines)
    └── JSON writer function declarations
```

---

## Module Structure

### Module 1: PDB I/O (`cmn_fncs.c`)

**Responsibilities**:
- Read/write PDB files
- Parse ATOM/HETATM records
- Handle alternative locations
- Coordinate transformations

**Key Functions**:
- `read_pdb()`: Main PDB reader
- `number_of_atoms()`: Count atoms
- `write_pdb()`: Write PDB file
- `residue_idx()`: Map atoms to residues

**Dependencies**:
- File I/O utilities
- String manipulation
- Memory allocation

---

### Module 2: Sequence & Residue Identification (`cmn_fncs.c`)

**Responsibilities**:
- Extract nucleotide sequences
- Identify residue types
- Classify purines/pyrimidines
- Map residue indices

**Key Functions**:
- `get_seq()`: Extract sequence
- `residue_ident()`: Classify residues
- `base_ident()`: Identify base types
- `get_bpseq()`: Extract base pair sequences

**Dependencies**:
- PDB I/O
- Atom name matching

---

### Module 3: Frame Calculation (`app_fncs.c`, `ana_fncs.c`)

**Responsibilities**:
- Calculate base reference frames
- Load standard templates
- Least-squares fitting
- Frame storage/retrieval

**Key Functions**:
- `base_frame()`: Calculate frames
- `ls_fitting()`: Least-squares fitting
- `ref_frames()`: Recalculate for pairs
- `mst2orien()`: Convert matrix to array

**Dependencies**:
- Template file loading
- Mathematical operations
- Matrix operations

---

### Module 4: Base Pair Validation (`cmn_fncs.c`)

**Responsibilities**:
- Validate geometric constraints
- Check hydrogen bonds
- Calculate quality scores
- Determine pair types

**Key Functions**:
- `check_pair()`: Main validation
- `get_hbond_ij()`: H-bond detection
- `adjust_pairQuality()`: Quality adjustment
- `calculate_more_bppars()`: Pair parameters

**Dependencies**:
- Frame calculation
- H-bond detection
- Overlap calculation

---

### Module 5: H-Bond Detection (`cmn_fncs.c`)

**Responsibilities**:
- Find potential H-bonds
- Resolve conflicts
- Validate H-bond types
- Assign linkage types

**Key Functions**:
- `get_hbond_ij()`: Initial detection
- `hb_atompair()`: Conflict resolution
- `validate_hbonds()`: Type validation
- `good_hbatoms()`: Atom pair validation
- `donor_acceptor()`: Role determination

**Dependencies**:
- Atom type identification
- Distance calculations

---

### Module 6: Base Pair Finding (`find_pair.c`)

**Responsibilities**:
- Greedy matching algorithm
- Best partner selection
- Mutual matching verification
- Pair ordering

**Key Functions**:
- `find_bestpair()`: Main algorithm
- `best_pair()`: Find best partner
- `re_ordering()`: Reorder pairs

**Dependencies**:
- Validation module
- Frame calculation

---

### Module 7: Parameter Calculation (`ana_fncs.c`)

**Responsibilities**:
- Calculate step parameters
- Calculate helical parameters
- Base pair parameters
- WC/Wobble classification

**Key Functions**:
- `bpstep_par()`: Step parameters
- `helical_par()`: Helical parameters
- `check_wc_wobble_pair()`: Classification

**Dependencies**:
- Frame calculation
- Matrix operations

---

## Execution Flows

### Flow 1: Complete find_pair Execution

```
Command: find_pair 1EHZ.pdb

1. Parse arguments
   └─ Determine options: ds, pairs, detailed, etc.

2. Initialize
   ├─ set_my_globals("find_pair")
   ├─ get_3dna_pars(&misc_pars)  # Load parameters
   └─ Initialize JSON writer (if enabled)

3. Read PDB
   ├─ number_of_atoms("1EHZ.pdb") → num = 452
   ├─ Allocate: AtomName[1..452], ResName[1..452], etc.
   └─ read_pdb(...) → Fill atom arrays

4. Map Residues
   └─ residue_idx(...) → num_residue = 16, seidx[1..16][1..2]

5. Extract Sequence
   └─ get_seq(...) → bseq[1..16] = "ACGT...", RY[1..16]

6. Calculate Frames
   └─ base_frame(...) → orien[1..16][1..9], org[1..16][1..3]

7. Find Pairs
   └─ find_bestpair(...) → num_bp = 8, base_pairs[1..8][...]

8. Reorder
   └─ re_ordering(...) → bp_idx[1..8], helix_idx[1..num_helix][...]

9. Write Output
   └─ x3dna_input("1EHZ.inp", ...) → Write .inp file

10. Cleanup
    └─ Free all arrays, close files
```

---

### Flow 2: Complete analyze Execution

```
Command: analyze 1EHZ.inp

1. Parse arguments
   └─ Determine: istart, istep, bz, etc.

2. Initialize
   ├─ set_my_globals("analyze")
   └─ Initialize JSON writer

3. Read Input
   └─ read_input("1EHZ.inp", ...) → pair_num[1..2][1..8]

4. Read PDB (again)
   └─ read_pdb("1EHZ.pdb", ...) → Same as find_pair

5. Validate Pairs
   └─ pair_checking(...) → Validate pair_num[][]

6. Extract Sequences
   └─ get_bpseq(...) → bp_seq[1..2][1..8]

7. Recalculate Frames
   └─ ref_frames(...) → orien[1..2][1..8*9], org[1..2][1..8*3]

8. Calculate Parameters
   └─ get_parameters(...) → bp_step_par[1..7][1..6]

9. Write Output
   └─ write_mst(...) → Write parameter files

10. Cleanup
    └─ Free all arrays
```

---

## Build System

### Compilation Structure

The legacy code uses a simple C build system:

```
Makefile/CMakeLists.txt
├─ Compile each .c file to .o
├─ Link all .o files
└─ Create executables:
   ├─ find_pair
   ├─ analyze
   └─ find_pair_analyze
```

### Dependencies

**External Libraries**:
- Standard C library only (no external dependencies)
- Uses Numerical Recipes algorithms (included in nrutil.c)

**Internal Dependencies**:
```
find_pair.c
├─ Requires: cmn_fncs.c, app_fncs.c, nrutil.c
└─ Optional: json_writer.c

analyze.c
├─ Requires: cmn_fncs.c, ana_fncs.c, app_fncs.c, nrutil.c
└─ Optional: json_writer.c

find_pair_analyze.c
├─ Requires: find_pair.c, analyze.c
└─ All dependencies from both
```

---

## Data Flow

### Data Flow Through find_pair

```
Input: PDB file
  ↓
[PDB Parsing]
  ├─ AtomName[][], ResName[][], xyz[][], etc.
  ↓
[Residue Mapping]
  ├─ seidx[][] (atom ranges per residue)
  ↓
[Sequence Extraction]
  ├─ bseq[] (base sequence)
  ├─ RY[] (purine/pyrimidine classification)
  ↓
[Frame Calculation]
  ├─ orien[][] (rotation matrices)
  ├─ org[][] (origins)
  ├─ NC1xyz[][] (N and C1' coordinates)
  └─ o3_p[][] (O3' and P coordinates)
  ↓
[Pair Finding]
  ├─ base_pairs[][] (pair data)
  ↓
[Reordering]
  ├─ bp_idx[] (reordered indices)
  ├─ helix_idx[][] (helix regions)
  └─ helix_marker[] (helix boundaries)
  ↓
Output: .inp file
```

### Data Flow Through analyze

```
Input: .inp file + PDB file
  ↓
[Input Parsing]
  ├─ pair_num[][] (residue indices)
  ↓
[PDB Parsing] (same as find_pair)
  ├─ Atom arrays
  ↓
[Frame Recalculation]
  ├─ orien[][] (frames for pairs)
  └─ org[][] (origins for pairs)
  ↓
[Parameter Calculation]
  ├─ bp_step_par[][] (step parameters)
  └─ bp_heli_par[][] (helical parameters)
  ↓
Output: Parameter files
```

---

## Global State

### Global Variables (`struct_Gvars`)

**Location**: Defined in `x3dna.h`, instantiated in `cmn_fncs.c`

**Structure**:
```c
struct struct_Gvars {
    long DEBUG;
    long VERBOSE;
    long NUM_ELE;
    long CHAIN_CASE;
    long ALL_MODEL;
    long ATTACH_RESIDUE;
    long THREE_LETTER_NTS;
    long PDBV3;
    long ORIGINAL_COORDINATE;
    long OCCUPANCY;
    long HEADER;
    long mmcif;
    double NT_CUTOFF;
    char X3DNA_VER[BUF512];
    char X3DNA_HOMEDIR[BUF512];
    char CHAIN_MARKERS[BUF512];
    char REBUILD_CHAIN_IDS[BUF512];
    char* PROGNAME;
    char** ATOM_NAMES;
    long NUM_SATOM;
    char** ATOMLIST;
    long NUM_SBASE;
    char** BASELIST;
    char** AtomName0;
    char** ResName0;
    long Name0;
    long label_RC8_YC6;
    miscPars misc_pars;  // Validation parameters
};
```

**Key Global Variables**:
- `Gvars.misc_pars`: Validation parameters (min_base_hb, hb_lower, etc.)
- `Gvars.X3DNA_HOMEDIR`: Path to X3DNA data files
- `Gvars.PROGNAME`: Program name for error messages
- `Gvars.ATOMLIST`: Valid atom names
- `Gvars.BASELIST`: Valid base names

**Initialization**:
- `set_my_globals()`: Called at program start
- Sets program name, loads atom/base lists, initializes misc_pars
- `get_3dna_pars()`: Loads parameters from `misc_3dna.par` file

**Cleanup**:
- `clear_my_globals()`: Called at program end
- Frees allocated global arrays

---

## Integration Points

### Between find_pair and analyze

**File Interface**: `.inp` file
- Format: ASCII text
- Contains: PDB file name, duplex status, pair numbers
- Written by: `x3dna_input()` in find_pair
- Read by: `read_input()` in analyze

**Data Consistency**:
- Both programs read the same PDB file
- Both use same residue indexing (from `residue_idx()`)
- Frames may differ (recalculated in analyze)

### With External Files

**Template Files**:
- Location: `{X3DNA_HOMEDIR}/Atomic_{base}.pdb`
- Used by: `base_frame()` for least-squares fitting
- Format: Standard PDB files with ideal base geometries

**Parameter File**:
- Location: `{X3DNA_HOMEDIR}/misc_3dna.par`
- Used by: `get_3dna_pars()` to load validation parameters
- Format: XML-like format with parameter tags

**Atom/Base Lists**:
- Location: `{X3DNA_HOMEDIR}/atomlist.dat`, `baselist.dat`
- Used by: `get_atomlist()`, `get_baselist()` for valid names

---

## Summary

### Key Architectural Principles

1. **Sequential Processing**: Clear pipeline from PDB → residues → frames → pairs → parameters
2. **Global State**: Heavy use of global variables for configuration
3. **1-based Indexing**: All arrays use 1-based indexing throughout
4. **Monolithic Functions**: Large functions handle multiple responsibilities
5. **File-based Communication**: Programs communicate via `.inp` files

### Critical Design Decisions

1. **Frame Recalculation**: analyze recalculates frames (doesn't reuse from find_pair)
2. **PDB Re-reading**: analyze re-reads PDB file (doesn't pass data)
3. **Greedy Matching**: find_pair uses greedy algorithm (not exhaustive)
4. **Mutual Matching**: Requires both bases to choose each other

---

**Next**: [Data Structures](02_DATA_STRUCTURES.md) for detailed data organization


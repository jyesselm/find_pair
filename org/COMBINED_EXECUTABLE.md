# Combined find_pair + analyze Executable

## Overview

The `find_pair_analyze` executable combines both `find_pair_original` and `analyze_original` functionality into a single program that runs both operations in sequence.

## Building

The combined executable is built automatically when you compile:

```bash
cd org
mkdir -p build && cd build
cmake ..
make
```

This creates three executables in `build/bin/`:
- `find_pair_original` - Original find_pair program
- `analyze_original` - Original analyze program  
- `find_pair_analyze` - **Combined program** (new)

## Usage

```bash
./build/bin/find_pair_analyze [options] <pdb_file> [outfile]
```

### Basic Example

```bash
./build/bin/find_pair_analyze 1EHZ.pdb
```

This will:
1. Run `find_pair` on `1EHZ.pdb` → creates `1EHZ.inp`
2. Automatically run `analyze` on `1EHZ.inp` → creates analysis output

### Options

The combined executable supports options from both programs:

**find_pair options:**
- `-S`, `-1` : single strand mode
- `-C` : generate Curves input
- `-D` : divide helices
- `-P` : find all pairs
- `-M` : mapping mode
- `-T` / `-A` : include/exclude HETATM
- `-Z` : detailed output
- `-W` : include waters
- `-hjb` : find all base combinations

**analyze options:**
- `-t <file>` : calculate torsions
- `-bz` : BZ junction analysis
- `-ri` : ring information
- `-si` : simple parameters
- `-abi` : ABI format
- `-circ` : circular structure
- `-C` : continue from previous
- `-S=<start>,<step>` : step parameters

### Example with Options

```bash
# Run with detailed output and simple parameters
./build/bin/find_pair_analyze -Z -si 1EHZ.pdb

# Run with waters included
./build/bin/find_pair_analyze -W 1EHZ.pdb

# Run with custom output file
./build/bin/find_pair_analyze 1EHZ.pdb my_output.inp
```

## How It Works

1. **Step 1 - find_pair**: 
   - Reads the PDB file
   - Identifies base pairs
   - Creates an `.inp` file (input file for analyze)

2. **Step 2 - analyze**:
   - Reads the `.inp` file created by find_pair
   - Calculates structural parameters
   - Generates analysis output files

## Output Files

The combined program generates the same output files as running the two programs separately:

**From find_pair:**
- `<basename>.inp` - Input file for analyze
- `bestpairs.pdb` - Best base pairs
- `hel_regions.pdb` - Helical regions
- `ref_frames.dat` - Reference frames
- Other intermediate files

**From analyze:**
- `<basename>.out` - Main analysis output
- `bp_step.par` - Base pair step parameters
- `bp_helical.par` - Helical parameters
- Other analysis files

## Technical Details

- The combined executable includes code from both `find_pair.c` and `analyze.c`
- Functions `handle_str()` and `process_str()` were made non-static to allow calling from the combined program
- The `main()` functions in the original files were renamed to `find_pair_main()` and `analyze_main()` to avoid conflicts
- When building the combined executable, `BUILD_COMBINED=1` is defined to exclude the original main() functions

## Comparison

| Executable | Size | Functionality |
|------------|------|---------------|
| `find_pair_original` | 458 KB | Base pair identification only |
| `analyze_original` | 405 KB | Parameter analysis only |
| `find_pair_analyze` | 476 KB | **Both operations in sequence** |

## Notes

- The combined executable automatically determines the `.inp` filename based on the input PDB file
- If find_pair fails to create the `.inp` file, analyze will not run
- All options from both programs are supported
- The program shows progress messages indicating which step is running


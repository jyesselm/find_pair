# Original X3DNA Code with Debug Output

This directory contains the original x3dna-v2.4 source code with extensive debug statements to trace base pair parameter calculation.

## Purpose

To compare the exact execution flow and intermediate values between the original x3dna code and our modernized version.

## Files

- `src/` - Original x3dna source files with debug statements
- `include/` - Header files
- `bin/` - Compiled executables
- `CMakeLists.txt` - Build configuration

## Usage

```bash
# Build
cd org
mkdir build && cd build
cmake ..
make

# Run with debug output
./find_pair_original /path/to/1EHZ.pdb > 1ehz.inp
./analyze_original 1ehz.inp > debug_output.txt
```

## Key Debug Points

1. `calculate_more_bppars` - Frame extraction and parameter calculation
2. `check_pair` - When parameters are calculated
3. `bpstep_par` - Actual parameter calculation function
4. Frame values in `orien` array at calculation time


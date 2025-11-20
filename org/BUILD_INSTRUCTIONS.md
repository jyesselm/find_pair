# Building the Original X3DNA with Debug Output

## Quick Start

```bash
cd org
mkdir -p build && cd build
cmake ..
make
```

## Troubleshooting

If you get compilation errors:

1. **Missing headers**: Make sure `include/` has `x3dna.h` and `x3dna_fncs.h`
   ```bash
   ls include/*.h
   ```

2. **Missing source files**: Check that all `.c` files are in `src/`
   ```bash
   ls src/*.c
   ```

3. **Linker errors**: The original x3dna code may need additional libraries. Check the original Makefile for dependencies.

## Running with Debug

```bash
# From org/ directory
./run_debug.sh /path/to/1EHZ.pdb

# Or manually:
cd build/bin
./find_pair_original /path/to/1EHZ.pdb > 1ehz.inp 2> 1ehz_debug.txt
```

## What to Look For

The debug output will show:
- `[DEBUG calculate_more_bppars]` - Frame extraction and parameter calculation
- `[DEBUG check_pair]` - When dir_z is calculated
- Frame values in `orien` array
- Final parameters: Shift, Slide, Rise, Tilt, Roll, Twist

Compare these values with our modernized code to find discrepancies.


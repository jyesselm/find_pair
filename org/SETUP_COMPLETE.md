# Original X3DNA Debug Setup - COMPLETE ✓

## What Was Done

1. ✅ Copied original x3dna-v2.4 source files to `org/src/`
2. ✅ Added debug statements to key functions:
   - `calculate_more_bppars`: Prints frame values, extracted matrices, and calculated parameters
   - `check_pair`: Prints dir_z calculation and z-axis values
3. ✅ Created build system (CMakeLists.txt)
4. ✅ Created run script (run_debug.sh)
5. ✅ Created comparison script (compare_debug_output.py)

## Files Structure

```
org/
├── src/
│   ├── find_pair.c      (with debug)
│   ├── cmn_fncs.c       (with debug in calculate_more_bppars and check_pair)
│   ├── ana_fncs.c
│   ├── app_fncs.c
│   ├── nrutil.c
│   └── analyze.c
├── include/
│   ├── x3dna.h
│   └── x3dna_fncs.h
├── CMakeLists.txt
├── run_debug.sh
├── compare_debug_output.py
└── README.md
```

## Next Steps

1. **Build the original code:**
   ```bash
   cd org
   mkdir -p build && cd build
   cmake ..
   make
   ```

2. **Run with debug output:**
   ```bash
   cd ..
   ./run_debug.sh /path/to/1EHZ.pdb
   ```

3. **Compare with our code:**
   - Run our modernized code with similar debug output
   - Use `compare_debug_output.py` to analyze differences:
     ```bash
     python compare_debug_output.py 1EHZ_find_pair_debug.txt our_debug.txt
     ```

## Key Debug Points

The debug output will show:
- **Frame extraction**: How `orien[i]` and `orien[j]` are accessed
- **Matrix construction**: How r1 and r2 are built from orien arrays
- **dir_z calculation**: The z-axis dot product
- **Final parameters**: The 6 base pair parameters from bpstep_par

This will help identify exactly where our implementation differs from the original.


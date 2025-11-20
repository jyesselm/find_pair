# Quick Start Guide

## 1. Build the Original Code

```bash
cd org
mkdir -p build && cd build
cmake ..
make
```

## 2. Run with Debug Output

```bash
cd ..
./test_run.sh ../legacy_archive/fixtures/pdb/1EHZ.pdb
```

Or manually:
```bash
./build/bin/find_pair_original ../legacy_archive/fixtures/pdb/1EHZ.pdb > 1ehz.inp 2> 1ehz_debug.txt
```

## 3. What to Look For

The debug output (`1ehz_debug.txt`) will contain:
- `[DEBUG calculate_more_bppars]` - Shows frame extraction
- `[DEBUG] Building r1 and r2 matrices` - Shows how matrices are built
- `[DEBUG] Result:` - Shows final parameters (Shift, Slide, Rise, Tilt, Roll, Twist)

## 4. Compare with Our Code

Add similar debug statements to our code in `src/rnamake/base/util/find_pair.cpp` and compare:
- Frame values in `orien` array
- How r1 and r2 are extracted
- Final parameters

This will help identify the exact difference!

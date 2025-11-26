# Legacy Test Tools - Isolated Component Testing

**Date**: 2025-11-25  
**Purpose**: Break down legacy H-bond detection into testable components

---

## Overview

These tools isolate different steps of the legacy H-bond detection algorithm to make debugging easier. Each tool tests ONE step in isolation, making it easier to identify where differences occur.

---

## Tools

### 1. `test_hbond_initial`

**Purpose**: Test ONLY initial H-bond detection (before conflict resolution or validation)

**What it does**:
- Finds all potential H-bonds using `good_hbatoms()` and `within_limits()`
- Shows which atoms are checked (seidx range)
- Shows which H-bonds pass the initial checks
- Outputs JSON for easy parsing

**Usage**:
```bash
cd org/build
./bin/test_hbond_initial ../data/pdb/3G8T.pdb 946 947
```

**Output**:
- Lists all atoms in seidx range for each residue
- Shows each H-bond found with distance
- JSON output with all initial H-bonds

**Use case**: Compare with modern's `initial_hbonds` to see if atom selection or initial detection differs

---

### 2. `test_hbond_conflict`

**Purpose**: Test ONLY conflict resolution (hb_atompair)

**What it does**:
- Gets initial H-bonds (same as test_hbond_initial)
- Applies `hb_atompair` algorithm
- Shows which H-bonds are marked as conflicts (negative distance)
- Shows linkage types assigned

**Usage**:
```bash
cd org/build
./bin/test_hbond_conflict ../data/pdb/3G8T.pdb 946 947
```

**Output**:
- Initial H-bonds found
- After conflict resolution: which are kept, which are conflicts
- Linkage types for each H-bond
- JSON output

**Use case**: Compare with modern's `after_conflict_resolution` to see if conflict resolution differs

---

### 3. `test_hbond_validation`

**Purpose**: Test ONLY validation (validate_hbonds)

**What it does**:
- Gets initial H-bonds
- Applies conflict resolution
- Applies `validate_hbonds` to assign types (' ', '-', '*')
- Shows final H-bonds with types

**Usage**:
```bash
cd org/build
./bin/test_hbond_validation ../data/pdb/3G8T.pdb 946 947
```

**Output**:
- Initial H-bonds
- After conflict resolution
- After validation with types
- JSON output matching legacy JSON format

**Use case**: Compare with modern's `after_validation` to see if validation differs

---

## Comparison Workflow

### Step 1: Compare Initial Detection
```bash
# Legacy
cd org/build
./bin/test_hbond_initial ../data/pdb/3G8T.pdb 946 947 > legacy_initial.json

# Modern
cd ../..
./build/compare_hbond_stages data/pdb/3G8T.pdb 946 947 | grep "Stage 1" -A 50
```

**Questions**:
- Do both find the same number of initial H-bonds?
- Do both check the same atoms?
- Are distances the same?

### Step 2: Compare Conflict Resolution
```bash
# Legacy
cd org/build
./bin/test_hbond_conflict ../data/pdb/3G8T.pdb 946 947 > legacy_conflict.json

# Modern
cd ../..
./build/compare_hbond_stages data/pdb/3G8T.pdb 946 947 | grep "Stage 2" -A 50
```

**Questions**:
- Do both mark the same H-bonds as conflicts?
- Are linkage types the same?
- Are distances after conflict resolution the same?

### Step 3: Compare Validation
```bash
# Legacy
cd org/build
./bin/test_hbond_validation ../data/pdb/3G8T.pdb 946 947 > legacy_validation.json

# Modern
cd ../..
./build/compare_hbond_stages data/pdb/3G8T.pdb 946 947 | grep "Stage 3" -A 50
```

**Questions**:
- Do both assign the same types?
- Are final H-bonds the same?

---

## Simple Test Cases

### Test Case 1: Simple A-T Pair
```bash
# Find a simple A-T pair in a PDB
./bin/test_hbond_initial data/pdb/1BNA.pdb 1 2
./bin/test_hbond_conflict data/pdb/1BNA.pdb 1 2
./bin/test_hbond_validation data/pdb/1BNA.pdb 1 2
```

### Test Case 2: Problematic Pair (946, 947)
```bash
# The pair we're debugging
./bin/test_hbond_initial data/pdb/3G8T.pdb 946 947
./bin/test_hbond_conflict data/pdb/3G8T.pdb 946 947
./bin/test_hbond_validation data/pdb/3G8T.pdb 946 947
```

### Test Case 3: Pair with Modified Nucleotides
```bash
# Test with modified nucleotides
./bin/test_hbond_initial data/pdb/3KNC.pdb <res1> <res2>
```

---

## Benefits

1. **Isolation**: Each tool tests ONE step, making it easier to find where differences occur
2. **JSON Output**: Easy to parse and compare programmatically
3. **Detailed Output**: Shows intermediate steps and decisions
4. **No Algorithm Changes**: Uses exact legacy code, just isolated

---

## Next Steps

1. **Run test_hbond_initial** for pair (946, 947) and compare with modern
2. **Identify differences** at each stage
3. **Fix modern** to match legacy at each stage
4. **Verify** with simple test cases

---

## Files

- `org/src/test_hbond_initial.c` - Initial detection test
- `org/src/test_hbond_conflict.c` - Conflict resolution test
- `org/src/test_hbond_validation.c` - Validation test
- `org/CMakeLists.txt` - Build configuration


# Legacy Test Tools - Isolated Component Testing

**Date**: 2025-01-XX  
**Purpose**: Tools for testing individual legacy algorithm components in isolation  
**Status**: Reference for debugging and verification

---

## Table of Contents

1. [Overview](#overview)
2. [Available Tools](#available-tools)
3. [Comparison Workflow](#comparison-workflow)
4. [Usage Examples](#usage-examples)

---

## Overview

These tools isolate different steps of the legacy H-bond detection algorithm to make debugging easier. Each tool tests ONE step in isolation, making it easier to identify where differences occur between legacy and modern implementations.

**Benefits**:
1. **Isolation**: Each tool tests ONE step, making it easier to find where differences occur
2. **JSON Output**: Easy to parse and compare programmatically
3. **Detailed Output**: Shows intermediate steps and decisions
4. **No Algorithm Changes**: Uses exact legacy code, just isolated

---

## Available Tools

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

**Source**: `org/src/test_hbond_initial.c`

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

**Source**: `org/src/test_hbond_conflict_debug.c`

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

**Source**: `org/src/test_hbond_validation.c`

---

### 4. `test_donor_acceptor`

**Purpose**: Test donor/acceptor role determination

**What it does**:
- Tests `donor_acceptor()` function in isolation
- Shows role lookup for different base/atom combinations
- Outputs role determination results

**Usage**:
```bash
cd org/build
./bin/test_donor_acceptor <base_i> <base_j> <atom1> <atom2>
```

**Source**: `org/src/test_donor_acceptor.c`

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

## Usage Examples

### Test Case 1: Simple A-T Pair

```bash
# Find a simple A-T pair in a PDB
cd org/build
./bin/test_hbond_initial ../data/pdb/1BNA.pdb 1 2
./bin/test_hbond_conflict ../data/pdb/1BNA.pdb 1 2
./bin/test_hbond_validation ../data/pdb/1BNA.pdb 1 2
```

### Test Case 2: Problematic Pair (946, 947)

```bash
# The pair we're debugging
cd org/build
./bin/test_hbond_initial ../data/pdb/3G8T.pdb 946 947
./bin/test_hbond_conflict ../data/pdb/3G8T.pdb 946 947
./bin/test_hbond_validation ../data/pdb/3G8T.pdb 946 947
```

### Test Case 3: Pair with Modified Nucleotides

```bash
# Test with modified nucleotides
cd org/build
./bin/test_hbond_initial ../data/pdb/3KNC.pdb <res1> <res2>
```

---

## Building the Tools

All test tools are built with the legacy codebase:

```bash
cd org/build
cmake ..
make
# Tools are in org/build/bin/
```

---

## Integration with Modern Code

These tools help identify where modern implementation differs from legacy:

1. **Run legacy tool** for a specific pair → Get JSON output
2. **Run modern equivalent** for same pair → Get JSON output
3. **Compare JSON files** → Identify differences
4. **Fix modern code** → Repeat until matches

**Example**:
```bash
# Get legacy output
cd org/build
./bin/test_hbond_initial ../data/pdb/3G8T.pdb 946 947 > /tmp/legacy.json

# Get modern output (hypothetical tool)
cd ../..
./build/test_hbond_stages data/pdb/3G8T.pdb 946 947 > /tmp/modern.json

# Compare
diff /tmp/legacy.json /tmp/modern.json
```

---

## Files

- `org/src/test_hbond_initial.c` - Initial detection test
- `org/src/test_hbond_conflict_debug.c` - Conflict resolution test
- `org/src/test_hbond_validation.c` - Validation test
- `org/src/test_donor_acceptor.c` - Donor/acceptor test
- `org/CMakeLists.txt` - Build configuration

---

**Related**: [Implementation Guide](08_IMPLEMENTATION_GUIDE.md#debugging-strategies) for debugging strategies


# Simplified Development Workflow

**Goal**: Achieve 100% match to legacy code in modern C++  
**Principle**: Minimize file creation, maximize reuse, focus on comparison

---

## Core Workflow (3 Steps)

### 1. Make Changes to Modern Code
Edit existing files in `src/x3dna/` - **DO NOT create new files unless absolutely necessary**

### 2. Generate & Compare JSON
```bash
# Generate modern JSON
./build/generate_modern_json data/pdb/<PDB_ID>.pdb data/json

# Compare with legacy
python3 scripts/compare_json.py compare <PDB_ID>
```

### 3. Fix Until 100% Match
Repeat steps 1-2 until comparison shows 100% match.

---

## File Creation Rules

### ✅ DO Create Files When:
1. **New algorithm component** (e.g., new calculator class)
2. **New core data structure** (e.g., new base class)
3. **Critical bug fix requires isolation**

### ❌ DO NOT Create Files For:
1. **New test files** - Add to existing test files instead
2. **New comparison scripts** - Use `compare_json.py` with different options
3. **New documentation** - Update existing docs
4. **Helper utilities** - Add to existing utility files
5. **Debug scripts** - Add debug output to existing code

---

## Test Strategy: Regression-Only

**Single Test Goal**: Does modern output match legacy output?

### Consolidate Tests

Instead of many unit/integration tests, use **one regression test**:

```cpp
// tests/regression/test_legacy_match.cpp
// Single test that compares modern JSON to legacy JSON for all PDBs
```

**Benefits**:
- One test file to maintain
- Directly measures progress toward 100% match
- No need for complex test fixtures
- Easy to run: `make test`

### When to Add Unit Tests

**Only add unit tests when**:
- Debugging a specific calculation that's not matching
- Need to isolate a component for detailed investigation
- Component is complex enough to warrant isolated testing

**Otherwise**: Use the regression test + JSON comparison

---

## Development Workflow

### Daily Development Cycle

1. **Edit modern code** (`src/x3dna/`)
2. **Build**: `make release`
3. **Generate JSON**: `./build/generate_modern_json data/pdb/<PDB>.pdb data/json`
4. **Compare**: `python3 scripts/compare_json.py compare <PDB>`
5. **Fix differences** → Repeat

### When Differences Found

1. **Check comparison report** - Shows exactly what differs
2. **Add debug output** to modern code (don't create new debug files)
3. **Compare debug output** with legacy debug output
4. **Fix calculation** in modern code
5. **Re-run comparison** until match

### Debugging Strategy

**Add debug statements to existing code**:
```cpp
// In the calculation function
#ifdef DEBUG
    std::cerr << "[DEBUG] value: " << value << "\n";
#endif
```

**Don't create**:
- New debug scripts
- New debug test files
- Separate debug utilities

---

## Script Consolidation

### Primary Scripts (Keep These)

1. **`scripts/compare_json.py`** - Main comparison tool
   - Use for all comparisons
   - Supports all comparison types via options
   - Don't create new comparison scripts

2. **`scripts/rebuild_json.py`** - Regenerate JSON files
   - Use when you need fresh JSON
   - Don't create new regeneration scripts

### Scripts to Avoid Creating

- ❌ New comparison scripts
- ❌ New analysis scripts (use `compare_json.py --verbose`)
- ❌ New test scripts (use C++ tests)
- ❌ New helper scripts (add to existing scripts)

---

## Documentation Strategy

### Core Documents (Update These)

1. **`docs/100_PERCENT_MATCH_PLAN.md`** - Track progress toward 100% match
2. **`docs/TESTING_GUIDE.md`** - How to test and compare
3. **`docs/DATA_STRUCTURE.md`** - JSON structure reference

### Don't Create New Docs For

- New features (update existing docs)
- New workflows (update TESTING_GUIDE.md)
- New algorithms (add to 100_PERCENT_MATCH_PLAN.md)

---

## Code Organization

### Reuse Existing Files

**Before creating a new file, ask**:
1. Can this go in an existing file?
2. Can this be a function in an existing class?
3. Can this be a method in an existing utility?

**Example**: Instead of `new_utility.cpp`, add to `existing_utility.cpp`

### File Naming

**Stick to existing patterns**:
- `*_calculator.cpp` for calculation classes
- `*_finder.cpp` for search/find algorithms
- `*_validator.cpp` for validation logic
- `*_parser.cpp` for parsing code

**Don't create**:
- `*_helper.cpp` (add to existing file)
- `*_utils.cpp` (add to existing file)
- `*_debug.cpp` (add debug code to main file)

---

## Testing Philosophy

### Regression Test is Primary

**One test to rule them all**:
```cpp
TEST(RegressionTest, MatchLegacy) {
    // For each PDB:
    // 1. Generate modern JSON
    // 2. Load legacy JSON
    // 3. Compare all record types
    // 4. Assert 100% match
}
```

### Unit Tests are Secondary

**Only create unit tests when**:
- Regression test shows failure
- Need to isolate specific component
- Component is too complex for direct comparison

**Otherwise**: Trust the regression test

---

## Quick Reference

### Make a Change
1. Edit `src/x3dna/**/*.cpp` or `include/x3dna/**/*.hpp`
2. `make release`
3. `./build/generate_modern_json data/pdb/<PDB>.pdb data/json`
4. `python3 scripts/compare_json.py compare <PDB>`

### Debug a Difference
1. Add `#ifdef DEBUG` statements to modern code
2. Rebuild with `-DDEBUG`
3. Compare debug output with legacy
4. Fix calculation
5. Remove debug statements (or keep for future)

### Add New Feature
1. Add to existing file (don't create new file)
2. Update `docs/100_PERCENT_MATCH_PLAN.md`
3. Test with regression test
4. Compare JSON until 100% match

---

## Success Metrics

**Single metric**: **100% JSON match rate**

Track in: `docs/100_PERCENT_MATCH_PLAN.md`

**Don't track**:
- Code coverage (not relevant for matching goal)
- Unit test count (regression test is primary)
- File count (minimize files)

---

## Summary

**Golden Rules**:
1. ✅ Edit existing files, don't create new ones
2. ✅ Use regression test + JSON comparison
3. ✅ One comparison script (`compare_json.py`)
4. ✅ One main goal (100% match)
5. ✅ Minimal documentation (update existing)

**Anti-Patterns to Avoid**:
1. ❌ Creating new test files for each feature
2. ❌ Creating new comparison scripts
3. ❌ Creating new documentation files
4. ❌ Creating helper/utility files
5. ❌ Creating debug-only files

**Focus**: Match legacy output, not perfect code architecture.

---

## See Also

- **[100_PERCENT_MATCH_PLAN.md](100_PERCENT_MATCH_PLAN.md)** - Track progress toward 100% match


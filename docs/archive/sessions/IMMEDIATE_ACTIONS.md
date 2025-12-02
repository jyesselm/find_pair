# Immediate Actions - Start Here

**Date**: December 2, 2025  
**Priority**: CRITICAL

---

## The Problem

We've been stuck on pair mismatches because **we don't have solid index matching** between legacy and modern code. This is the foundation that everything else depends on.

Additionally, we have:
- **50+ documentation files** (mostly outdated)
- **36 C++ debugging tools** (mostly one-off)
- **60+ Python scripts** (overlapping functionality)
- **No centralized tracking** of what's working

**Bottom line**: Too much bloat, no solid foundation.

---

## The Solution

### 1. Fix Index Matching (Week 1)
Create a `ResidueTracker` that tracks every residue from the moment it's read from the PDB:
- Order read (read_index)
- Whether it was filtered (and why)
- Final modern index
- Matching legacy index from JSON

**Validation**: Must match 100% on all test PDBs before proceeding.

### 2. Unified Comparison Framework (Week 2)
One script: `scripts/unified_compare.py`
- Validates indices
- Compares JSON at each stage
- Tracks progress in CSV
- No caching by default
- Resume from failures

### 3. Massive Cleanup (Weeks 3-4)
- **Docs**: 50 files → 7 core files
- **Tools**: 36 tools → 1 tool + archive
- **Scripts**: 60 scripts → 3 scripts + archive

---

## Immediate Next Steps (TODAY)

### Step 1: Review Plans ✓
You're doing this now by reading this file.

**Review these files**:
1. `CLEANUP_AND_RESTRUCTURE_PLAN.md` - Overall strategy
2. `IMPLEMENTATION_PLAN.md` - Detailed implementation with code

**Decision needed**: Approve plan to proceed?

### Step 2: Set Up Branch
```bash
cd /Users/jyesselman2/Library/CloudStorage/Dropbox/2_code/cpp/find_pair_2
git checkout -b fix-index-matching
```

### Step 3: Implement ResidueTracker (First Focus)

**Create files**:
1. `include/x3dna/residue_tracker.hpp` (see IMPLEMENTATION_PLAN.md)
2. `src/x3dna/residue_tracker.cpp` (see IMPLEMENTATION_PLAN.md)

**Update CMakeLists.txt**:
```cmake
# Add to src/x3dna/CMakeLists.txt
set(X3DNA_SOURCES
    # ... existing files ...
    residue_tracker.cpp
)
```

**Test immediately**:
```bash
make release
# Should compile without errors
```

### Step 4: Integrate into PDB Parsing

**File to edit**: `src/x3dna/pdb_parser.cpp`

**Add**:
1. Create tracker at start of parsing
2. Add each residue as read
3. Mark filtered residues
4. Assign final indices
5. Load legacy indices
6. Validate and export mapping

### Step 5: Test on Single PDB (1H4S)

```bash
# Generate modern JSON with index validation
./build/generate_modern_json data/pdb/1H4S.pdb data/json/

# Should create:
# - data/json/*/1H4S.json (modern JSON)
# - data/index_mapping/1H4S.json (index mapping)
# - Should print validation result

# Check mapping
cat data/index_mapping/1H4S.json | python3 -m json.tool | less
```

**Expected output**:
```json
[
  {
    "read_index": 0,
    "legacy_index": 1,
    "modern_index": 0,
    "filtered": false,
    "filter_reason": "",
    "chain_id": "A",
    "residue_seq": 1,
    "insertion": "",
    "residue_name": "  G"
  },
  // ... more residues ...
]
```

### Step 6: Validate Mapping

**Check**:
- `num_modern == num_legacy` ✓
- All non-filtered have both indices ✓
- Indices are sequential ✓
- No gaps ✓

**If validation fails**:
- Document the issue
- Investigate root cause
- Fix before proceeding

### Step 7: Expand to test_set_10

```bash
# Get test set PDBs
cat resources/test_sets/test_set_10.json

# Generate mappings for all
for pdb in $(cat resources/test_sets/test_set_10.json | jq -r '.[]'); do
    echo "Processing $pdb..."
    ./build/generate_modern_json data/pdb/${pdb}.pdb data/json/
done

# Check results
ls data/index_mapping/
```

**Goal**: 10/10 PDBs validate successfully

---

## What NOT to Do

❌ **Don't** try to fix pair mismatches yet  
❌ **Don't** create new comparison scripts  
❌ **Don't** create new debugging tools  
❌ **Don't** update documentation yet  
❌ **Don't** try to compare outputs yet  

✅ **Do** focus ONLY on index matching  
✅ **Do** test thoroughly on small set first  
✅ **Do** validate at every step  
✅ **Do** document any issues found  

---

## Success Criteria for Week 1

By end of Week 1, you should have:

✅ ResidueTracker implemented and tested  
✅ Integrated into PDB parsing  
✅ Index mappings generated for test_set_10  
✅ 100% validation pass on test_set_10  
✅ Index mappings generated for test_set_100  
✅ Validation report for test_set_100  
✅ CSV file: `data/index_validation_status.csv`  

**If not 100% on test_set_10**: STOP and fix before continuing.

---

## Quick Reference

### File Structure After Week 1

```
find_pair_2/
├── include/x3dna/
│   └── residue_tracker.hpp          # NEW
├── src/x3dna/
│   ├── residue_tracker.cpp          # NEW
│   └── pdb_parser.cpp               # MODIFIED
├── data/
│   ├── index_mapping/               # NEW
│   │   ├── 1H4S.json
│   │   ├── 2BNA.json
│   │   └── ... (all test PDBs)
│   └── index_validation_status.csv  # NEW
└── CLEANUP_AND_RESTRUCTURE_PLAN.md  # THIS PLAN
```

### Commands

```bash
# Build
make release

# Test single PDB
./build/generate_modern_json data/pdb/1H4S.pdb data/json/

# Check mapping
cat data/index_mapping/1H4S.json | python3 -m json.tool | head -30

# Generate all test_set_10 mappings
for pdb in 1H4S 2BNA 3DNA 4DNA 5DNA 6DNA 7DNA 8DNA 9DNA 1AAA; do
    ./build/generate_modern_json data/pdb/${pdb}.pdb data/json/
done
```

### Validation Checks

For each PDB, verify:
1. Mapping file exists: `data/index_mapping/{PDB}.json`
2. `num_modern == num_legacy`
3. All non-filtered have both indices
4. No duplicate indices
5. Indices are sequential (0-based modern, 1-based legacy)

---

## Next Steps After Week 1

Once index validation is 100%:
1. Implement `scripts/unified_compare.py`
2. Run batch comparisons with CSV tracking
3. Archive old tools and scripts
4. Clean up documentation
5. Achieve 100% match on all stages

But **DO NOT** proceed to Week 2 until Week 1 is complete.

---

## Questions?

See:
- `CLEANUP_AND_RESTRUCTURE_PLAN.md` - Overall strategy and reasoning
- `IMPLEMENTATION_PLAN.md` - Detailed code examples and checklist

**Key principle**: Fix the foundation before building on it.



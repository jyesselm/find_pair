# Comparison Framework Quick Reference

**For**: Quick reference on current comparison capabilities and next steps  
**See Also**: `COMPARISON_IMPROVEMENT_PLAN.md` for comprehensive analysis

---

## Current Comparison Capabilities

### ✅ What We Have

1. **12-Stage Validation Pipeline**
   - Stage 0: Residue Indices
   - Stage 1: PDB Atoms
   - Stage 2: Base Frame Calc (98.48% match on 1000 PDBs ⭐)
   - Stage 3: Distance Checks
   - Stage 4: H-bond List
   - Stage 5: Pair Validation
   - Stage 6: Find Bestpair Selection (🎯 PRIMARY OUTPUT - CRITICAL)
   - Stage 7: Base Pair
   - Stage 8: Step Parameters
   - Stage 9: Helical Parameters

2. **Multiple Entry Points**
   ```bash
   # CLI script (most comprehensive)
   python3 scripts/compare_json.py compare 1EHZ
   
   # Pytest tests (systematic)
   pytest tests_python/integration/test_stage_validation.py -v
   
   # New unified CLI
   fp2-validate --pdb 1EHZ
   ```

3. **Comparison Features**
   - ✅ Parallel execution (multi-process)
   - ✅ Result caching
   - ✅ Test sets (10, 50, 100, 500, 1000 PDBs)
   - ✅ Per-stage comparison functions
   - ✅ Detailed error tracking
   - ✅ Configurable tolerance (global: 1e-6)

### ⚠️ What We're Missing

1. **Verbose Mode** - Can't see field-by-field details for a specific PDB
2. **Provenance Tracking** - Don't know WHERE values come from in code
3. **Systematic Validation** - Haven't run all stages on all PDBs
4. **Field-Specific Tolerances** - Only global tolerance
5. **Cross-Stage Validation** - Don't verify index consistency across stages

---

## How to Debug a Failing PDB (Current Method)

### Step 1: Identify the Failure
```bash
python3 scripts/compare_json.py compare 1EHZ --diff-only
```

### Step 2: Look at Specific Stage
```bash
# Focus on one stage
python3 scripts/compare_json.py frames 1EHZ

# Or
pytest tests_python/integration/test_stage_validation.py::TestStageValidation::test_stage3_distance_checks --pdb 1EHZ -v
```

### Step 3: Manual JSON Inspection
```bash
# Compare JSON files manually
cat data/json_legacy/distance_checks/1EHZ.json
cat data/json/distance_checks/1EHZ.json

# Use jq for pretty printing
jq '.' data/json_legacy/distance_checks/1EHZ.json | head -50
```

### Step 4: Read Comparison Code
```python
# Look at comparison function
# File: x3dna_json_compare/distance_comparison.py
def compare_distance_checks(legacy_records, modern_records, tolerance=1e-6):
    # ...logic here
```

### Step 5: Add Debug Prints (Manual)
```python
# Modify comparison code temporarily
print(f"Legacy dorg: {leg_rec['dorg']}")
print(f"Modern dorg: {mod_rec['dorg']}")
print(f"Diff: {abs(leg_rec['dorg'] - mod_rec['dorg'])}")
```

**Problems with Current Method:**
- ❌ Time-consuming
- ❌ Requires code modifications
- ❌ Hard to see full context
- ❌ No provenance info
- ❌ Manual JSON diffing is tedious

---

## How Debugging WILL Work (With Verbose Mode)

### One Command, All Details
```bash
python3 scripts/compare_json.py compare 1EHZ --verbose
```

**Output:**
```
================================================================================
VERBOSE COMPARISON: 1EHZ
================================================================================

STAGE 3: distance_checks
------------------------
✅ MATCH (base_i=1, base_j=2)
  dorg:         14.523 == 14.523  ✓
  dNN:          15.234 == 15.234  ✓
  plane_angle:  12.456 == 12.456  ✓
  d_v:          2.345  == 2.345   ✓
  overlap_area: 45.678 == 45.678  ✓

❌ MISMATCH (base_i=3, base_j=4)
  dorg:         14.523 == 14.524  ✗ (diff: 0.001 > tolerance 1e-6)
  
  Legacy calculation:
    File: org/src/pair_geometry.c:234
    Function: calculate_pair_dorg()
    Residues: A:15 + A:18
  
  Modern calculation:
    File: src/x3dna/algorithms/base_pair_validator.cpp:456
    Function: BasePairValidator::calculate_dorg()
    Residues: legacy_idx=3, legacy_idx=4
  
  Related Records:
    base_frame_calc: RMS 0.023 == 0.023 ✓
    hbond_list: num_hbonds 2 == 2 ✓
```

**Benefits:**
- ✅ See everything in one place
- ✅ Know WHERE values come from
- ✅ See related records
- ✅ Clear diff highlighting
- ✅ No code modification needed

---

## Quick Status Check

### Run Quick Validation
```bash
# Fast subset (typically ~100-200 PDBs)
pytest tests_python/integration/test_stage_validation.py -v

# Specific stage
pytest tests_python/integration/test_stage_validation.py::TestStageValidation::test_stage3_distance_checks -v

# With max PDBs limit
pytest tests_python/integration/test_stage_validation.py -v --max-pdbs=10
```

### Check Current Match Rates
```bash
# Run comparison on test set
python3 scripts/compare_json.py compare --test-set 100

# Stage-specific
python3 scripts/compare_json.py frames --test-set 100
python3 scripts/compare_json.py pairs --test-set 100
```

### Generate Test Data
```bash
# Regenerate legacy JSON
cd org
./build/bin/find_pair_analyze ../data/pdb/1EHZ.pdb

# Regenerate modern JSON
./build/generate_modern_json data/pdb/1EHZ.pdb data/json/
```

---

## Key Comparison Functions

### Per-Stage Comparison

| Stage | Function | File | Key Fields |
|-------|----------|------|------------|
| Atoms | `compare_atoms()` | `atom_comparison.py` | atom_idx, xyz, atom_name |
| Base Frame | `compare_frames()` | `frame_comparison.py` | rms_fit, matched_atoms, template |
| Distance | `compare_distance_checks()` | `distance_comparison.py` | dorg, dNN, plane_angle, d_v, overlap |
| H-bonds | `compare_hbond_list()` | `hbond_comparison.py` | num_hbonds, donor, acceptor |
| Validation | `compare_pair_validation()` | `pair_validation_comparison.py` | is_valid, bp_type_id, quality_score |
| Selection | `compare_find_bestpair_selection()` | `find_bestpair_comparison.py` | pairs[] (PRIMARY) |
| Base Pair | `compare_base_pairs()` | `base_pair_comparison.py` | orien_i, orien_j, org_i, org_j |
| Steps | `compare_step_params()` | `step_comparison.py` | 6 parameters |
| Helical | `compare_helical_params()` | `step_comparison.py` | 6 parameters |

### Common Pattern
```python
def compare_XXX(legacy_records, modern_records, tolerance=1e-6):
    """Compare XXX records."""
    result = XXXComparison()  # Result object
    
    # 1. Build dictionaries keyed by matching criteria
    legacy_dict = {}
    for rec in legacy_records:
        key = (rec['field1'], rec['field2'])  # Matching key
        legacy_dict[key] = rec
    
    # 2. Find common/missing/extra
    common_keys = set(legacy_dict.keys()) & set(modern_dict.keys())
    
    # 3. Compare common records
    for key in common_keys:
        leg = legacy_dict[key]
        mod = modern_dict[key]
        
        # Field-by-field comparison
        if abs(leg['value'] - mod['value']) > tolerance:
            result.mismatches.append(...)
    
    return result
```

---

## Configuration

### Current Config: `comparison_config.yaml`
```yaml
comparisons:
  atoms: false         # Disabled (limited test data)
  frames: true         # Stage 2,3,4
  steps: true          # Stage 8,9
  pairs: true          # Stage 5,6,7
  hbond_list: true     # Stage 4
  residue_indices: true  # Stage 0

tolerance: 1.0e-6      # Global tolerance

cache:
  enabled: true
  force_recompute: false
```

### Proposed Enhanced Config
```yaml
comparisons:
  atoms: true
  frames: true
  steps: true
  pairs: true
  hbond_list: true
  residue_indices: true

# Field-specific tolerances (PROPOSED)
tolerance:
  default: 1.0e-6
  by_field:
    rms_fit: 1.0e-4       # More lenient
    overlap_area: 1.0e-5
    rotation_matrix: 1.0e-4
    dorg: 1.0e-6          # Strict
    dNN: 1.0e-6

# Verbose mode settings (PROPOSED)
verbose:
  show_provenance: false   # Show calculation source
  show_related: true       # Show related records
  max_mismatches: 20       # Limit per stage
  
cache:
  enabled: true
  force_recompute: false
```

---

## Test Sets

### Available Test Sets
```
resources/test_sets/
├── test_set_10.json      # 10 PDBs (quick smoke test)
├── test_set_50.json      # 50 PDBs
├── test_set_100.json     # 100 PDBs (standard test set)
├── test_set_500.json     # 500 PDBs
└── test_set_1000.json    # 1000 PDBs (comprehensive)
```

### Usage
```bash
# Use test set
python3 scripts/compare_json.py compare --test-set 100

# Or in pytest
pytest tests_python/integration/test_stage_validation.py --test-set 100 -v
```

### Fast PDBs Subset
```
data/valid_pdbs_fast.json   # ~100-200 fast PDBs
```

---

## Priority Actions

### Immediate (This Week):
1. ✅ Create comprehensive improvement plan (DONE)
2. **Implement verbose mode** ⬅️ NEXT
   - Create `verbose_reporter.py`
   - Add `--verbose` flag to `compare_json.py`
   - Test on 10-PDB set

### Short-Term (Next 2 Weeks):
3. **Systematic validation**
   - Run all stages on 100-PDB test set
   - Document failures
   - Fix critical bugs

4. **Cross-stage validation**
   - Verify index consistency
   - Check reference validity

### Medium-Term (Next Month):
5. **Path to 100%**
   - Fix all identified bugs
   - Achieve >99% match rate
   - Complete documentation

---

## Critical Insights

### 🎯 Stage 6 is THE Primary Output
**Stage 6: find_bestpair_selection** contains the actual selected pairs - this is what the algorithm produces. ALL other stages are intermediate steps or derived data.

**Priority**: Stage 6 MUST be 100% match.

### 📊 Current Best Match Rate
**Stage 2: base_frame_calc** - 98.48% match on 1000 PDBs

This is our benchmark. We should aim for ≥98% on ALL stages.

### ⚠️ Testing Gaps
- **Atoms**: Only 1 PDB validated (7EH2, which has bugs)
- **Stages 3-5**: Good comparison code, limited large-scale validation
- **Stages 8-9**: Comparison exists, needs more testing

### 🔧 Key to 100% Match
Good testing = ability to debug = ability to fix bugs = 100% match

**Verbose mode is critical** - it's the debugging tool that makes everything else possible.

---

## Common Issues & Solutions

### Issue: "Tolerance too strict"
**Symptom**: Many false positives, fields differ by ~1e-7  
**Solution**: Use field-specific tolerances
```yaml
tolerance:
  rms_fit: 1.0e-4  # More lenient for floating point accumulation
```

### Issue: "Can't find why values differ"
**Symptom**: Values differ but don't know why  
**Solution**: Use verbose mode (COMING SOON)
```bash
python3 scripts/compare_json.py compare 1EHZ --verbose --show-provenance
```

### Issue: "Too many failures to debug"
**Symptom**: 100+ PDBs fail, overwhelming  
**Solution**: Use `--stop-on-first` flag
```bash
pytest tests_python/integration/test_stage_validation.py -x -v
```

### Issue: "Slow test execution"
**Symptom**: Tests take >30 minutes  
**Solution**: Use smaller test sets, parallel execution
```bash
# Use fast subset
pytest tests_python/integration/ --max-pdbs=10 -v

# Increase workers
python3 scripts/compare_json.py compare --test-set 100 --threads 20
```

---

## Resources

### Documentation
- **This File**: Quick reference
- `COMPARISON_IMPROVEMENT_PLAN.md`: Comprehensive analysis
- `TESTING_GUIDE.md`: Testing procedures
- `JSON_DATA_TYPES_AND_COMPARISONS.md`: Field definitions

### Code
- `x3dna_json_compare/`: Comparison module
- `scripts/compare_json.py`: Main CLI
- `tests_python/integration/`: Pytest tests

### Data
- `data/json_legacy/`: Legacy output
- `data/json/`: Modern output
- `data/valid_pdbs_fast.json`: Fast PDB subset
- `resources/test_sets/`: Predefined test sets

---

**Last Updated**: December 6, 2025  
**Status**: Planning Complete, Ready for Implementation


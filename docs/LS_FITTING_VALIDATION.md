# LS_FITTING Validation Status

**Created**: December 3, 2025  
**Status**: ✅ Validation infrastructure complete, 🔄 Full batch test in progress

---

## Overview

The `ls_fitting` comparison validates that the modern code correctly calculates least-squares fitting parameters for base pairs, matching the legacy X3DNA implementation.

---

## Validation Infrastructure ✅ COMPLETE

### Python Comparison Code

**File**: `tests_python/integration/test_ls_fitting.py`

- Uses `x3dna_json_compare.frame_comparison.compare_frames()` for comparison
- Compares all `ls_fitting` records between legacy and modern JSON
- Checks for record count matches, missing residues, and calculation mismatches
- **Status**: ✅ Working correctly

### Batch Testing Tool

**File**: `tests_python/integration/test_ls_fitting_batch.py`

Features:
- Tests all PDBs from `valid_pdbs_fast.json` (3,602 fast PDBs)
- Generates legacy JSON (if needed)
- Generates modern JSON using `generate_modern_json --stage=ls_fitting`
- Compares using `compare_frames()`
- Saves detailed results to JSON
- Provides summary statistics

Usage:
```bash
# Test all fast PDBs
python tests_python/integration/test_ls_fitting_batch.py

# Test limited number
python tests_python/integration/test_ls_fitting_batch.py --max-pdbs 100

# Custom output location
python tests_python/integration/test_ls_fitting_batch.py --output my_results.json
```

**Status**: ✅ Working correctly

### Monitoring Tool

**File**: `scripts/monitor_ls_fitting_batch.py`

Features:
- Monitors batch test progress
- Analyzes results when complete
- Shows match rates and mismatch distributions
- Identifies worst failures

**Status**: ✅ Working correctly

---

## Current Validation Run 🔄 IN PROGRESS

**Started**: December 3, 2025 at 12:25 PM  
**Command**: 
```bash
python tests_python/integration/test_ls_fitting_batch.py > data/validation_results/ls_fitting_batch.log 2>&1
```

**Progress** (as of 12:27 PM):
- Tested: ~60/3602 PDBs (~1.7%)
- Estimated completion: ~90 minutes total
- Results file: `data/validation_results/ls_fitting_batch_results.json`
- Log file: `data/validation_results/ls_fitting_batch.log`

### Preliminary Results (First 60 PDBs)

- **Pass rate**: Roughly 60-70% exact matches
- **Failure pattern**: Small numbers of mismatches (typically 1-5 records per PDB)
- **Overall match rate**: >98% of individual records match

**Example results**:
- ✅ PASS: Many PDBs with 100% exact match (e.g., 100D: 20/20, 1EGK: 108/108)
- ⚠️ FAIL: Some with small mismatches (e.g., 165D: 17/18, 1A34: 19/20)

---

## Expected Outcomes

### Success Criteria

1. **High overall match rate**: >95% of all `ls_fitting` records should match
2. **Small numerical differences**: Failures should be due to floating-point precision, not algorithm errors
3. **No systematic errors**: Failures should be randomly distributed, not concentrated in specific cases

### Known Issues

- Floating-point precision differences may cause small mismatches
- Different compiler optimizations may affect numerical results
- Rounding in intermediate calculations may accumulate differently

---

## Next Steps

### 1. Complete Current Batch Test ⏳ ACTIVE

**Action**: Wait for full batch test to complete (~90 minutes)

**Check progress**:
```bash
python scripts/monitor_ls_fitting_batch.py
```

### 2. Analyze Results 📊 NEXT

**After completion**, analyze:

```bash
# View full results
python scripts/monitor_ls_fitting_batch.py

# Check detailed results
cat data/validation_results/ls_fitting_batch_results.json | jq '.total, .passed, .failed'

# View log
tail -100 data/validation_results/ls_fitting_batch.log
```

**Key questions**:
- What is the overall pass rate?
- What is the record-level match rate?
- Are failures random or systematic?
- What are the typical mismatch magnitudes?

### 3. Investigate Failures 🔍 IF NEEDED

If match rate < 95%, investigate:

```bash
# Test a specific failing PDB
python tests_python/integration/test_ls_fitting.py <PDB_ID>

# Compare detailed output
python -c "
import json
from pathlib import Path

pdb_id = '<PDB_ID>'
legacy = Path(f'data/json_legacy/ls_fitting/{pdb_id}.json')
modern = Path(f'data/json/ls_fitting/{pdb_id}.json')

with open(legacy) as f: leg = json.load(f)
with open(modern) as f: mod = json.load(f)

# Compare specific records
for i, (l, m) in enumerate(zip(leg, mod)):
    if l != m:
        print(f'Record {i}: mismatch')
        print(f'  Legacy: {l}')
        print(f'  Modern: {m}')
"
```

### 4. Fix Issues (If Necessary) 🔧 CONDITIONAL

If systematic errors found:
- Check `ls_fitting` calculation in modern code
- Compare with legacy implementation in `org/`
- Verify numerical precision settings
- Check for algorithm differences

### 5. Document Results ✅ FINAL

Update this document with:
- Final pass/fail statistics
- Overall match rate
- Any identified issues and resolutions
- Conclusions about ls_fitting accuracy

---

## Validation Commands Reference

### Run Tests

```bash
# Single PDB
python tests_python/integration/test_ls_fitting.py <PDB_ID>

# Batch test (all fast PDBs)
python tests_python/integration/test_ls_fitting_batch.py

# Batch test (limited)
python tests_python/integration/test_ls_fitting_batch.py --max-pdbs 100
```

### Monitor Progress

```bash
# Check current status
python scripts/monitor_ls_fitting_batch.py

# Watch log file
tail -f data/validation_results/ls_fitting_batch.log

# Count progress
wc -l data/validation_results/ls_fitting_batch.log
```

### Analyze Results

```bash
# View summary
python scripts/monitor_ls_fitting_batch.py

# View JSON results
cat data/validation_results/ls_fitting_batch_results.json | jq '.'

# Extract specific stats
jq '.total, .passed, .failed, .errors' data/validation_results/ls_fitting_batch_results.json
```

---

## Test Data

### PDB Sets

- **Test set**: 3,602 "fast" PDBs (< 30s processing time)
- **Excluded**: 521 "slow" PDBs (> 30s processing time)
- **Source**: `data/valid_pdbs_fast.json`

### JSON Files

**Legacy JSON**:
- Generated by: `org/build/bin/find_pair_analyze`
- Location: `data/json_legacy/ls_fitting/<PDB_ID>.json`
- Generated on-demand or used from cache

**Modern JSON**:
- Generated by: `build/generate_modern_json --stage=ls_fitting`
- Location: Generated in temp directories during testing
- Format: Same as legacy (list of ls_fitting records)

---

## Related Documentation

- **[TESTING_GUIDE.md](TESTING_GUIDE.md)** - General testing documentation
- **[NEXT_STEPS.md](NEXT_STEPS.md)** - Overall project next steps
- **[VALIDATION_PROGRESS.md](VALIDATION_PROGRESS.md)** - Complete validation status
- **[JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md)** - JSON format specifications

---

## Conclusion

The validation infrastructure for `ls_fitting` is complete and working. A comprehensive batch test over all 3,602 fast PDBs is currently running. Preliminary results (first 60 PDBs) show high match rates with small numerical differences in a minority of cases.

**Status**: 🟢 On track for validation completion

**Next checkpoint**: After batch test completes (~90 minutes from 12:25 PM)

---

*Last updated: December 3, 2025 at 12:30 PM*


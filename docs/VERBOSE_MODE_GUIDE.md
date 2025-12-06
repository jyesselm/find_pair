# Verbose Mode User Guide

**Feature**: Detailed field-by-field comparison for debugging  
**Version**: 1.0  
**Date**: December 6, 2025

---

## Overview

Verbose mode provides detailed, human-readable comparison output showing **exactly what's being compared** for a specific PDB file. Instead of just seeing summary statistics, you can see field-by-field comparisons with clear indication of what matches and what differs.

### When to Use Verbose Mode

- ✅ **Debugging a specific PDB** - See exactly where values differ
- ✅ **Understanding a failure** - Get detailed context for mismatches
- ✅ **Validating a fix** - Confirm specific fields now match
- ✅ **Learning the system** - See what fields are compared at each stage

### When NOT to Use Verbose Mode

- ❌ **Batch validation** - Use summary mode for multiple PDBs
- ❌ **CI/CD pipelines** - Use `--quiet` for exit codes only
- ❌ **Quick status check** - Use standard comparison mode

---

## Basic Usage

### Command Syntax

```bash
python3 scripts/compare_json.py compare <PDB_ID> --verbose [OPTIONS]
```

### Simple Examples

```bash
# Verbose comparison of 1EHZ
python3 scripts/compare_json.py compare 1EHZ --verbose

# Save verbose output to file
python3 scripts/compare_json.py compare 1EHZ --verbose --output 1EHZ_comparison.txt

# Verbose comparison (inside venv)
source venv/bin/activate
python3 scripts/compare_json.py compare 1EHZ --verbose
```

### Important Notes

1. **Single PDB only** - Verbose mode only works for one PDB at a time
2. **Requires both legacy and modern JSON** - Files must exist for comparison
3. **Uses configured stages** - Respects `comparison_config.yaml` settings

---

## Output Format

### Header Section

```
================================================================================
VERBOSE COMPARISON: 1EHZ
================================================================================
Date: 2025-12-06 11:31:34
Tolerance: 1e-06
Stages: frames, steps, helical, hbond_list, pair_validation
Mode: All records
```

**Shows:**
- PDB ID being compared
- Timestamp of comparison
- Numerical tolerance used
- Which stages are being compared
- Comparison mode

### Stage Section

```
--------------------------------------------------------------------------------
STAGE 4: hbond_list
--------------------------------------------------------------------------------
Total legacy records: 35
Total modern records: 35
Common records: 34
⚠️  Missing in modern: 1
⚠️  Extra in modern: 0
❌ Mismatched records: 2
```

**Shows:**
- Stage name and number
- Record counts (legacy, modern, common)
- Missing/extra records (if any)
- Number of mismatches

### Record Comparison (Match)

```
✅ MATCH (base_i=1, base_j=2)
  Legacy source: data/json_legacy/distance_checks/1EHZ.json
  Modern source: data/json/distance_checks/1EHZ.json

  Fields:
    dorg:                14.523000       == 14.523000       ✓
    dNN:                 15.234000       == 15.234000       ✓
    plane_angle:         12.456000       == 12.456000       ✓
    d_v:                 2.345000        == 2.345000        ✓
    overlap_area:        45.678000       == 45.678000       ✓
```

**Shows:**
- ✅ Match indicator
- Record key (pair or residue identification)
- Source file paths
- Field-by-field comparison
- All values match (✓)

### Record Comparison (Mismatch)

```
❌ MISMATCH (base_i=3, base_j=4)
  Legacy source: data/json_legacy/distance_checks/1EHZ.json
  Modern source: data/json/distance_checks/1EHZ.json

  Fields:
    dorg:                14.523000       vs 14.524000       ✗ (diff: 1.000000e-03, tolerance: 1.000000e-06)
      ℹ️  Exceeds tolerance by 9.990000e-04
    dNN:                 15.234000       == 15.234000       ✓
    plane_angle:         12.456000       == 12.456000       ✓
    d_v:                 2.345000        == 2.345000        ✓
    overlap_area:        45.678000       == 45.678000       ✓
```

**Shows:**
- ❌ Mismatch indicator
- Record key
- Source file paths
- Field-by-field comparison
- **Mismatched field** with:
  - Actual values (legacy vs modern)
  - ✗ indicator
  - Difference amount
  - Tolerance threshold
  - How much it exceeds tolerance
- Matching fields (✓)

### Summary Section

```
--------------------------------------------------------------------------------
SUMMARY
--------------------------------------------------------------------------------
Stages compared: 7
Perfect matches: 6
Stages with differences: 1

Differences found:
  - hbond_list: 2 mismatches
  - distance_checks: 1 mismatch (dorg field)

Overall: ⚠️  DIFFERENCES FOUND
```

**Shows:**
- Total stages compared
- How many stages match perfectly
- How many stages have differences
- List of specific differences
- Overall status

---

## Understanding the Output

### Status Indicators

| Indicator | Meaning |
|-----------|---------|
| ✅ | Perfect match |
| ❌ | Mismatch found |
| ⚠️ | Warning (missing/extra records) |
| ✓ | Individual field matches |
| ✗ | Individual field differs |
| ℹ️ | Additional information |

### Field Formatting

**Floats:**
- Small values (<0.001): Scientific notation (e.g., `1.234567e-06`)
- Normal values: Fixed precision (e.g., `14.523000`)
- Very small values (<1e-10): Shown as `0.0`

**Arrays:**
- Short arrays (≤3 elements): `[1.23, 4.56, 7.89]`
- Long arrays: Truncated

**Strings:**
- Long strings: Truncated with `...`

### Difference Reporting

When a numerical field differs:
```
dorg: 14.523000 vs 14.524000 ✗ (diff: 1.000000e-03, tolerance: 1.000000e-06)
  ℹ️  Exceeds tolerance by 9.990000e-04
```

**Information provided:**
- Legacy value: `14.523000`
- Modern value: `14.524000`
- Absolute difference: `1.000000e-03`
- Configured tolerance: `1.000000e-06`
- How much it exceeds: `9.990000e-04`

---

## Advanced Usage

### Comparing Specific Stages

Use configuration to focus on specific stages:

```bash
# Create custom config (focus on distance checks only)
cat > custom_config.yaml << EOF
comparisons:
  atoms: false
  frames: false
  steps: false
  pairs: false
  hbond_list: false
  residue_indices: false
  distance_checks: true  # Only this stage

tolerance: 1.0e-6
EOF

# Run with custom config
python3 scripts/compare_json.py compare 1EHZ --verbose --config custom_config.yaml
```

### Adjusting Tolerance

Edit `comparison_config.yaml`:

```yaml
tolerance: 1.0e-5  # More lenient (e.g., for testing)
```

Or use a custom config file as shown above.

### Saving Multiple Reports

```bash
# Compare several PDBs and save each
for pdb in 1EHZ 1H4S 2BNA; do
    python3 scripts/compare_json.py compare $pdb --verbose \
        --output verbose_reports/${pdb}_verbose.txt
done
```

---

## Troubleshooting

### Problem: "Warning: Verbose mode only supported for single PDB"

**Cause**: You specified multiple PDBs

**Solution**: Only specify one PDB at a time
```bash
# ❌ Wrong
python3 scripts/compare_json.py compare 1EHZ 1H4S --verbose

# ✅ Correct
python3 scripts/compare_json.py compare 1EHZ --verbose
```

### Problem: "No comparison result for XYZ"

**Cause**: Missing JSON files (legacy or modern)

**Solution**: Generate missing files
```bash
# Generate legacy JSON
cd org
./build/bin/find_pair_analyze ../data/pdb/XYZ.pdb

# Generate modern JSON
cd ..
./build/generate_modern_json data/pdb/XYZ.pdb data/json/
```

### Problem: "Total legacy records: 0"

**Cause**: Stage not generated by legacy code for this PDB

**Solution**: This is expected for some stages/PDBs. The stage will show as matching if both are empty.

### Problem: Too many mismatches to read

**Cause**: Many records differ, output is overwhelming

**Solution**: 
1. Check summary first (bottom of output)
2. Focus on first few mismatches
3. Fix most common issue
4. Re-run verbose mode

---

## Comparison with Standard Mode

### Standard Mode Output

```bash
python3 scripts/compare_json.py compare 1EHZ
```

**Output:**
```
Total PDBs: 1
✅ Perfect matches: 0
⚠️  Files with differences: 1

H-BOND LIST STATISTICS
Total legacy H-bond lists: 35
Total modern H-bond lists: 35
Missing in modern: 1
Mismatched pairs: 2
```

**Use for:**
- Quick status check
- Batch validation
- Summary statistics

### Verbose Mode Output

```bash
python3 scripts/compare_json.py compare 1EHZ --verbose
```

**Output:**
```
❌ MISMATCH (base_i=18, base_j=55)
  Fields:
    num_hbonds: 4 vs 7 ✗ (diff: 3.0, tolerance: 1e-06)
      ℹ️  Exceeds tolerance by 2.999999
```

**Use for:**
- Debugging specific failures
- Understanding WHY values differ
- Validating fixes

---

## Integration with Development Workflow

### 1. Identify Failure

```bash
# Run batch validation
pytest tests_python/integration/test_stage_validation.py -v

# Find failing PDB
# Example output: FAILED test_stage4_hbond_list[1EHZ]
```

### 2. Debug with Verbose Mode

```bash
# Deep dive into specific PDB
python3 scripts/compare_json.py compare 1EHZ --verbose > debug/1EHZ_verbose.txt

# Review the output
less debug/1EHZ_verbose.txt
```

### 3. Locate Issue in Code

```bash
# Search for the mismatched field calculation
# Example: num_hbonds differs
grep -r "num_hbonds" src/
grep -r "num_hbonds" org/src/

# Compare implementations
```

### 4. Fix and Validate

```bash
# Make code changes
# ...

# Regenerate modern JSON
./build/generate_modern_json data/pdb/1EHZ.pdb data/json/

# Re-check with verbose mode
python3 scripts/compare_json.py compare 1EHZ --verbose

# Should now show: ✅ ALL STAGES MATCH PERFECTLY
```

### 5. Confirm at Scale

```bash
# Run full test suite
pytest tests_python/integration/test_stage_validation.py -v

# Or batch comparison
python3 scripts/compare_json.py compare --test-set 100
```

---

## Examples

### Example 1: Perfect Match

```bash
$ python3 scripts/compare_json.py compare 2BNA --verbose

================================================================================
VERBOSE COMPARISON: 2BNA
================================================================================
Date: 2025-12-06 11:45:23
Tolerance: 1e-06
Stages: frames, hbond_list, distance_checks
Mode: All records

--------------------------------------------------------------------------------
STAGE: distance_checks
--------------------------------------------------------------------------------
Total legacy records: 42
Total modern records: 42
Common records: 42
✅ All records match perfectly

--------------------------------------------------------------------------------
STAGE: hbond_list
--------------------------------------------------------------------------------
Total legacy records: 42
Total modern records: 42
Common records: 42
✅ All records match perfectly

--------------------------------------------------------------------------------
SUMMARY
--------------------------------------------------------------------------------
Stages compared: 2
Perfect matches: 2
Stages with differences: 0

Overall: ✅ ALL STAGES MATCH PERFECTLY
```

### Example 2: Single Field Mismatch

```bash
$ python3 scripts/compare_json.py compare 3G8T --verbose

❌ MISMATCH (base_i=5, base_j=10)
  Legacy source: data/json_legacy/distance_checks/3G8T.json
  Modern source: data/json/distance_checks/3G8T.json

  Fields:
    dorg:                12.345678       == 12.345678       ✓
    dNN:                 13.456789       vs 13.456790       ✗ (diff: 1.000000e-06, tolerance: 1.000000e-06)
      ℹ️  Exceeds tolerance by 0.000000e+00
    plane_angle:         8.901234        == 8.901234        ✓
    d_v:                 1.234567        == 1.234567        ✓
    overlap_area:        23.456789       == 23.456789       ✓
```

**Analysis**: `dNN` differs by exactly the tolerance (boundary case). May need to adjust tolerance or investigate numerical precision.

### Example 3: Multiple Stages with Differences

```bash
$ python3 scripts/compare_json.py compare 1FIR --verbose

--------------------------------------------------------------------------------
SUMMARY
--------------------------------------------------------------------------------
Stages compared: 4
Perfect matches: 2
Stages with differences: 2

Differences found:
  - hbond_list: 3 mismatches
  - distance_checks: 1 mismatch

Overall: ⚠️  DIFFERENCES FOUND
```

**Action**: Focus on `hbond_list` first (more mismatches), then `distance_checks`.

---

## Tips and Best Practices

### 1. Start with Summary

Always review the summary section first to understand the scope of differences.

### 2. Focus on First Failure

Fix the first mismatch, regenerate, and re-check. Often fixes cascade.

### 3. Save Output for Reference

```bash
# Save verbose output for later review
python3 scripts/compare_json.py compare 1EHZ --verbose \
    --output analysis/1EHZ_$(date +%Y%m%d).txt
```

### 4. Use with Git

```bash
# Before making changes
python3 scripts/compare_json.py compare 1EHZ --verbose > before.txt

# Make code changes, regenerate JSON

# After changes
python3 scripts/compare_json.py compare 1EHZ --verbose > after.txt

# Compare
diff before.txt after.txt
```

### 5. Combine with Standard Tools

```bash
# Generate verbose report
python3 scripts/compare_json.py compare 1EHZ --verbose > report.txt

# Search for specific fields
grep "dorg:" report.txt
grep "❌" report.txt  # Find all mismatches
grep "ℹ️" report.txt  # Find all info messages
```

---

## Future Enhancements

Planned improvements to verbose mode:

### Phase 2.1: Provenance Tracking (Planned)
Show WHERE in the code values are calculated:
```
  Modern calculation:
    File: src/x3dna/algorithms/base_pair_validator.cpp:456
    Function: BasePairValidator::calculate_dorg()
```

### Phase 2.2: Related Record Lookup (Planned)
Show related records for context:
```
  Related records:
    ├─ base_frame_calc: RMS=0.023 ✓
    └─ pair_validation: is_valid=1 ✓
```

### Phase 2.3: Field-Specific Tolerances (Planned)
Configure different tolerances per field:
```yaml
tolerance:
  default: 1e-6
  rms_fit: 1e-4
  overlap_area: 1e-5
```

---

## Summary

Verbose mode is your **primary debugging tool** for comparison failures:

- ✅ Shows exactly what's being compared
- ✅ Clear diff highlighting  
- ✅ Field-by-field details
- ✅ Easy to read format
- ✅ Saves to file for later analysis

**Remember**: Good debugging = Fast bug fixes = 100% match rate

---

**Last Updated**: December 6, 2025  
**Version**: 1.0  
**See Also**: 
- `COMPARISON_IMPROVEMENT_PLAN.md` - Overall improvement strategy
- `COMPARISON_QUICK_GUIDE.md` - Quick reference
- `TESTING_GUIDE.md` - Testing procedures


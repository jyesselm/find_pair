# JSON Comparison Guide

**Date**: 2025-11-25  
**Purpose**: Comprehensive guide for comparing JSON outputs between modern and legacy code

This document explains how to use all available JSON comparison tools to debug differences between modern and legacy implementations.

---

## Overview

The codebase provides multiple specialized tools for comparing different aspects of the algorithm:

1. **Full JSON Comparison** - Compare complete JSON outputs
2. **Residue Ordering** - Compare residue ordering and indexing
3. **Quality Scores** - Compare quality score calculations
4. **H-bond Detection** - Compare H-bond finding and validation
5. **Residue Identification** - Compare nucleotide recognition
6. **PDB Parsing** - Compare atom/residue parsing

---

## Table of Contents

1. [Full JSON Comparison](#1-full-json-comparison)
2. [Residue Ordering Comparison](#2-residue-ordering-comparison)
3. [Quality Score Comparison](#3-quality-score-comparison)
4. [H-bond Detection Comparison](#4-h-bond-detection-comparison)
5. [Residue Identification Comparison](#5-residue-identification-comparison)
6. [PDB Parsing Comparison](#6-pdb-parsing-comparison)
7. [Batch Comparison](#7-batch-comparison)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Full JSON Comparison

### Tool: `scripts/compare_json.py`

**Purpose**: Compare complete JSON outputs between modern and legacy implementations.

**Usage**:
```bash
# Compare all available files
python3 scripts/compare_json.py compare

# Compare specific PDB file(s)
python3 scripts/compare_json.py compare 3G8T
python3 scripts/compare_json.py compare 3G8T 6CAQ 3KNC

# Compare only atoms
python3 scripts/compare_json.py atoms 3G8T

# Compare only frames
python3 scripts/compare_json.py frames 3G8T

# Use legacy mode JSON files
python3 scripts/compare_json.py compare --legacy-mode

# Save report to file
python3 scripts/compare_json.py compare --output report.md

# Show only files with differences
python3 scripts/compare_json.py compare --diff-only

# Verbose output
python3 scripts/compare_json.py compare --verbose
```

**What it compares**:
- `pdb_atoms` - All atoms parsed from PDB
- `base_frame_calc` - Base frame calculations
- `pair_validation` - Pair validation results
- `hbond_list` - H-bond lists for pairs
- `base_pairs` - Final selected base pairs

**Output**: Detailed comparison report showing matches and differences.

**See**: `docs/JSON_DATA_TYPES_AND_COMPARISONS.md` for detailed field descriptions.

---

## 2. Residue Ordering Comparison

### Tool: `build/compare_residue_ordering`

**Purpose**: Compare residue ordering between modern and legacy (critical for 100% match).

**Prerequisites**: Generate JSON files first:
```bash
# Generate modern JSON
build/generate_residue_ordering_json data/pdb/3G8T.pdb data/residue_ordering/3G8T_modern.json

# Generate legacy JSON
org/build/bin/generate_residue_ordering_json data/pdb/3G8T.pdb data/residue_ordering/3G8T_legacy.json
```

**Usage**:
```bash
build/compare_residue_ordering \
  data/residue_ordering/3G8T_modern.json \
  data/residue_ordering/3G8T_legacy.json
```

**Output**:
- Total residue count comparison
- Ordering match/mismatch status
- Detailed differences (if any)

**Exit Code**:
- `0`: Perfect match
- `1`: Mismatch found

**Status**: ✅ **100% Match Verified** (4,934+ PDBs tested)

**See**: `docs/RESIDUE_ORDERING_COMPARISON.md` for detailed guide.

---

## 3. Quality Score Comparison

### Tool: `build/compare_quality_scores`

**Purpose**: Compare quality score calculations for specific pairs (useful for debugging missing pairs).

**Usage**:
```bash
build/compare_quality_scores \
  <pdb_file> \
  <residue1_idx> <residue2_idx> \
  [modern_json] [legacy_json]
```

**Example**:
```bash
# Compare quality scores for pair (946, 947) in 3G8T
build/compare_quality_scores \
  data/pdb/3G8T.pdb \
  946 947 \
  data/json/3G8T.json \
  data/json_legacy/3G8T.json
```

**What it compares**:
- Base quality score: `dorg + 2.0*d_v + plane_angle/20.0`
- `adjust_pairQuality`: H-bond adjustment
- `bp_type_id`: Base pair type ID
- Final adjusted quality score: `base_score + adjust_pairQuality - (bp_type_id == 2 ? 2.0 : 0.0)`
- Number of "good" H-bonds (distance in [2.5, 3.5])
- Total number of H-bonds
- Validation status
- Selection status

**Output**: Side-by-side comparison showing all quality score components.

**Use Case**: Debug why a pair is missing - check if quality scores match and if selection logic differs.

---

## 4. H-bond Detection Comparison

### Tool: `build/compare_hbond_detection`

**Purpose**: Compare H-bond detection for specific pairs (separate from base pair validation).

**Usage**:
```bash
build/compare_hbond_detection \
  <pdb_file> \
  <residue1_idx> <residue2_idx> \
  [modern_json] [legacy_json]
```

**Example**:
```bash
# Compare H-bonds for pair (946, 947) in 3G8T
build/compare_hbond_detection \
  data/pdb/3G8T.pdb \
  946 947 \
  data/json/3G8T.json \
  data/json_legacy/3G8T.json
```

**What it compares**:
- H-bonds found by legacy vs modern
- H-bond distances and types
- Number of "good" H-bonds (distance in [2.5, 3.5])
- `adjust_pairQuality` value from H-bonds

**Output**: Detailed H-bond comparison showing matches and differences.

**Use Case**: Debug H-bond detection differences that affect quality scores.

---

### Tool: `build/compare_hbond_stages`

**Purpose**: Compare H-bond detection at multiple stages (initial, after conflict resolution, after validation).

**Usage**:
```bash
build/compare_hbond_stages \
  <pdb_file> \
  <residue1_idx> <residue2_idx> \
  [legacy_json]
```

**Example**:
```bash
# Compare H-bond stages for pair (946, 947) in 3G8T
build/compare_hbond_stages \
  data/pdb/3G8T.pdb \
  946 947 \
  data/json_legacy/3G8T.json
```

**What it compares**:
- **Initial detection**: Raw H-bonds found before any processing
- **After conflict resolution**: H-bonds after `hb_atompair` step
- **After validation**: Final validated H-bonds (what gets used)

**Output**: Stage-by-stage comparison showing how H-bonds are filtered.

**Use Case**: Understand where H-bond differences originate in the pipeline.

---

### Tool: `build/compare_initial_hbonds`

**Purpose**: Compare initial raw H-bond detection (before conflict resolution and validation).

**Usage**:
```bash
build/compare_initial_hbonds \
  <pdb_file> \
  <residue1_idx> <residue2_idx> \
  [legacy_json]
```

**Example**:
```bash
# Compare initial H-bonds for pair (946, 947) in 3G8T
build/compare_initial_hbonds \
  data/pdb/3G8T.pdb \
  946 947 \
  data/json_legacy/3G8T.json
```

**What it compares**:
- Initial H-bonds found by modern (before any processing)
- Initial H-bonds found by legacy (from `test_hbond_initial` tool)

**Output**: Comparison of raw H-bond detection.

**Use Case**: Debug the very first step of H-bond detection.

---

## 5. Residue Identification Comparison

### Tool: `build/compare_residue_identification`

**Purpose**: Compare residue identification (nucleotide recognition) between modern and legacy.

**Usage**:
```bash
build/compare_residue_identification \
  <pdb_file> \
  [modern_json] [legacy_json]
```

**Example**:
```bash
# Compare residue identification for 3G8T
build/compare_residue_identification \
  data/pdb/3G8T.pdb \
  data/json/3G8T.json \
  data/json_legacy/3G8T.json
```

**What it compares**:
- Which residues are recognized as nucleotides
- Residue type classification
- Frame calculation success/failure

**Output**: Comparison of residue identification results.

**Use Case**: Debug why some residues are missing (e.g., 3KNC only recognizing 16/66 residues).

---

## 6. PDB Parsing Comparison

### Tool: `build/compare_pdb_parsing`

**Purpose**: Compare PDB parsing (atom/residue counts) between modern and legacy.

**Usage**:
```bash
build/compare_pdb_parsing \
  <pdb_file> \
  [modern_json] [legacy_json]
```

**Example**:
```bash
# Compare PDB parsing for 3G8T
build/compare_pdb_parsing \
  data/pdb/3G8T.pdb \
  data/json/3G8T.json \
  data/json_legacy/3G8T.json
```

**What it compares**:
- Total atom count
- Total residue count
- HETATM handling
- Water molecule handling

**Output**: Comparison of parsing results.

**Use Case**: Debug parsing differences that affect residue counts.

---

## 7. Batch Comparison

### Residue Ordering Batch Comparison

**Script**: `scripts/generate_and_compare_residue_ordering_batch.py`

**Usage**:
```bash
# Compare all 1000 PDBs in test set
python3 scripts/generate_and_compare_residue_ordering_batch.py 1000

# Compare 10 PDBs for quick test
python3 scripts/generate_and_compare_residue_ordering_batch.py 10
```

**What it does**:
1. Generates modern JSON for all PDBs
2. Generates legacy JSON for all PDBs
3. Compares each pair
4. Saves summary to `data/residue_ordering/summary_{size}.json`

**Output**: Progress updates and final statistics.

---

## 8. Troubleshooting

### Issue: Tool not found

**Symptom**: `build/compare_*: No such file or directory`

**Fix**: Build the tools:
```bash
cd build
cmake ..
make compare_quality_scores compare_hbond_detection compare_residue_ordering
```

### Issue: JSON file not found

**Symptom**: `Error: Cannot open file: data/json/3G8T.json`

**Fix**: Generate the JSON files first:
```bash
# Generate modern JSON
build/generate_modern_json data/pdb/3G8T.pdb data/json/3G8T.json

# Legacy JSON should already exist in data/json_legacy/
```

### Issue: Residue index mismatch

**Symptom**: Tool reports wrong residue for given index

**Fix**: Ensure you're using legacy residue indices (1-based). Use `Structure::get_residue_by_legacy_idx()` in modern code.

### Issue: Comparison shows differences but they're expected

**Symptom**: Some differences are acceptable (e.g., H-bond type classification)

**Note**: Not all differences are bugs. See `docs/100_PERCENT_MATCH_PLAN.md` for documented limitations.

---

## Quick Reference

### Comparison Tools by Use Case

| Use Case | Tool | Command |
|----------|------|---------|
| Full JSON comparison | `scripts/compare_json.py` | `python3 scripts/compare_json.py compare <pdb>` |
| Residue ordering | `build/compare_residue_ordering` | `build/compare_residue_ordering <modern> <legacy>` |
| Quality scores | `build/compare_quality_scores` | `build/compare_quality_scores <pdb> <idx1> <idx2>` |
| H-bond detection | `build/compare_hbond_detection` | `build/compare_hbond_detection <pdb> <idx1> <idx2>` |
| H-bond stages | `build/compare_hbond_stages` | `build/compare_hbond_stages <pdb> <idx1> <idx2>` |
| Initial H-bonds | `build/compare_initial_hbonds` | `build/compare_initial_hbonds <pdb> <idx1> <idx2>` |
| Residue identification | `build/compare_residue_identification` | `build/compare_residue_identification <pdb>` |
| PDB parsing | `build/compare_pdb_parsing` | `build/compare_pdb_parsing <pdb>` |

### Common Workflows

#### Debug Missing Pair

1. **Check quality scores**:
   ```bash
   build/compare_quality_scores data/pdb/3G8T.pdb 946 947
   ```

2. **Check H-bond detection**:
   ```bash
   build/compare_hbond_detection data/pdb/3G8T.pdb 946 947
   ```

3. **Check H-bond stages**:
   ```bash
   build/compare_hbond_stages data/pdb/3G8T.pdb 946 947
   ```

#### Verify Residue Ordering

1. **Generate JSON files**:
   ```bash
   build/generate_residue_ordering_json data/pdb/3G8T.pdb data/residue_ordering/3G8T_modern.json
   org/build/bin/generate_residue_ordering_json data/pdb/3G8T.pdb data/residue_ordering/3G8T_legacy.json
   ```

2. **Compare**:
   ```bash
   build/compare_residue_ordering \
     data/residue_ordering/3G8T_modern.json \
     data/residue_ordering/3G8T_legacy.json
   ```

#### Full Comparison

1. **Compare all JSON records**:
   ```bash
   python3 scripts/compare_json.py compare 3G8T
   ```

2. **Save report**:
   ```bash
   python3 scripts/compare_json.py compare 3G8T --output 3G8T_comparison.md
   ```

---

## Related Documentation

- `docs/JSON_DATA_TYPES_AND_COMPARISONS.md` - Detailed JSON structure and field descriptions
- `docs/RESIDUE_ORDERING_COMPARISON.md` - Detailed residue ordering comparison guide
- `docs/100_PERCENT_MATCH_PLAN.md` - Overall match plan and status
- `docs/HBOND_FIX_SUMMARY.md` - H-bond fix details

---

## Summary

This guide covers all available JSON comparison tools. Use the appropriate tool based on what you're debugging:

- **Full comparison**: `scripts/compare_json.py`
- **Residue ordering**: `build/compare_residue_ordering`
- **Quality scores**: `build/compare_quality_scores`
- **H-bond detection**: `build/compare_hbond_detection` or `build/compare_hbond_stages`
- **Residue identification**: `build/compare_residue_identification`
- **PDB parsing**: `build/compare_pdb_parsing`

For detailed information about JSON structure, see `docs/JSON_DATA_TYPES_AND_COMPARISONS.md`.


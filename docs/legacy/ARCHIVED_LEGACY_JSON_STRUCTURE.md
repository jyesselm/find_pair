# Legacy JSON Structure - Base Pair Identification

## Overview

The legacy code outputs base pair information at different stages:

### 1. `pair_validation` Records (from `find_pair` phase)
- **Source**: `json_writer_record_pair_validation` called from `check_pair` in `cmn_fncs.c`
- **When**: During `find_pair` phase, for ALL pairs that pass initial checks (cdns + overlap)
- **Location**: `{pdb_id}_pair_validation.json` (split file)
- **Contains**: Validation results for all pairs checked during `find_bestpair` iteration
- **Note**: This includes pairs that pass validation but may NOT be selected by `find_bestpair`

### 2. `base_pair` Records (from `find_pair` phase)
- **Source**: `json_writer_record_base_pair` called from `calculate_more_bppars` in `cmn_fncs.c`
- **When**: During `find_pair` phase, for pairs that pass cdns + overlap + hbond check
- **Location**: `{pdb_id}_base_pair.json` (split file)
- **Contains**: Frame details (orien_i, orien_j, org_i, org_j, dir_xyz, bp_type)
- **Note**: This is a SUBSET of `pair_validation` records (only pairs that also pass hbond check)

### 3. `base_pairs` Record (from `analyze` phase)
- **Source**: `json_writer_record_base_pairs` called from `analyze.c`
- **When**: During `analyze` phase, after reading `.inp` file from `find_pair`
- **Location**: Main JSON file, `calculations` array
- **Contains**: `pair_num` array with the FINAL selected pairs (after reordering)
- **Note**: This is the actual list of pairs selected by `find_bestpair` and then reordered

## The Problem

We've been comparing against `pair_validation` records, which include ALL pairs that pass validation during `find_pair`, not just the ones selected by `find_bestpair`.

**What we should compare against**: The pairs actually selected by `find_bestpair` during the `find_pair` phase, BEFORE reordering.

## Solution

The `base_pair` records are from `calculate_more_bppars`, which is only called for pairs that:
1. Pass cdns check (distance, d_v, plane_angle, dNN all in range)
2. Pass overlap check
3. Pass hbond check

These are the pairs that would be candidates for selection by `find_bestpair`. However, `find_bestpair` uses mutual best match logic, so not all of these will be selected.

**The correct comparison should be**:
- Compare against `base_pair` records (from `find_pair` phase, before reordering)
- These represent pairs that passed all validation checks during the original identification
- The actual selection by `find_bestpair` uses mutual best match, which we should replicate

## Current Status

- Modern code records `pair_validation` for all pairs passing initial checks (matches legacy)
- Modern code records `base_pair` for pairs passing all checks including hbond (matches legacy)
- Modern code should compare against `base_pair` records, not `pair_validation` records
- The comparison should verify that modern `find_bestpair` selects the same pairs as legacy


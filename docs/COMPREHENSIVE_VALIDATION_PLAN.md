# Comprehensive Validation Plan: 100% Legacy-Modern Matching

**Date**: December 4, 2025  
**Goal**: Achieve 100% match rate across ALL stages using `data/valid_pdbs_fast.json` (3602 PDBs)  
**Strategy**: Stop-on-first-failure, cleanup temp files, reuse existing infrastructure

---

## ✅ Infrastructure Restored (Dec 4-5, 2025)

The following files were restored from Dropbox:

### `build/generate_modern_json` ✅
- Executable restored and working
- Usage: `./build/generate_modern_json <pdb> <output_dir> [--stage=STAGE]`
- Stages: atoms, frames, distances, hbonds, validation, selection, steps, helical, all

### `x3dna_json_compare/` Module ✅
All source files restored:
- `__init__.py`, `atom_comparison.py`, `base_pair_comparison.py`
- `config.py`, `distance_comparison.py`, `find_bestpair_comparison.py`
- `frame_comparison.py`, `hbond_comparison.py`, `json_comparison.py`
- `json_file_finder.py`, `json_validator.py`, `models.py`
- `pair_comparison.py`, `pair_validation_comparison.py`, `parallel_executor.py`
- `pdb_utils.py`, `residue_indices_comparison.py`, `result_cache.py`
- `selection_comparison.py`, `step_comparison.py`

---

## Executive Summary

### Current Status (What's Complete)

| Stage | Description | Status | Test Results File |
|-------|-------------|--------|-------------------|
| **Stage 0** | Residue Indices | ✅ **100% COMPLETE** | `data/residue_indices_test_results.json` (3602/3602) |
| **Stage 1** | Atoms | ✅ **100% COMPLETE** | `data/atoms_test_results.json` (2320+ tested, 100% success) |
| **Stage 2** | LS Fitting | ✅ **99.92% COMPLETE** | 3599/3602 pass (3 edge cases) |
| **Stage 3** | Distance Checks | ⏳ **NOT VALIDATED** | Needs batch testing |
| **Stage 4** | H-bonds | ⏳ **NOT VALIDATED** | Needs batch testing |
| **Stage 5** | Pair Validation | ⏳ **NOT VALIDATED** | Needs batch testing |
| **Stage 6** | Find Bestpair Selection | ⏳ **NOT VALIDATED** | **PRIMARY OUTPUT** |
| **Stage 7** | Base Pair Records | ⏳ **NOT VALIDATED** | Needs batch testing |
| **Stage 8** | Step Parameters | ⏳ **NOT VALIDATED** | Needs batch testing |
| **Stage 9** | Helical Parameters | ⏳ **NOT VALIDATED** | Needs batch testing |

### Key Files Available

```
data/valid_pdbs_fast.json    # 3602 fast PDBs to test
data/json_legacy/            # Legacy JSON per record type
data/json/                   # Modern JSON per record type
tests_python/                # pytest infrastructure
scripts/                     # Utility scripts
```

### Legacy JSON Coverage

| Directory | Files | Notes |
|-----------|-------|-------|
| `pdb_atoms` | 3636 | ✅ Good coverage |
| `residue_indices` | 3916 | ✅ Good coverage |
| `ls_fitting` | 4128 | ✅ Good coverage |
| `frame_calc` | 4930 | ✅ Good coverage |
| `base_frame_calc` | 4930 | ✅ Good coverage |
| `distance_checks` | 4128 | ✅ Good coverage |
| `hbond_list` | 4129 | ✅ Good coverage |
| `pair_validation` | 4128 | ✅ Good coverage |
| `base_pair` | 4129 | ✅ Good coverage |
| `find_bestpair_selection` | 3921 | ✅ Good coverage |
| `bpstep_params` | 3540 | ✅ Good coverage |
| `helical_params` | 3540 | ✅ Good coverage |

---

## Validation Stages (In Order)

### Stage 0: Residue Indices ✅ COMPLETE

**What it tests**: PDB residue parsing, chain/seq/insertion code mapping to legacy indices  
**Legacy JSON**: `data/json_legacy/residue_indices/<PDB>.json`  
**Status**: 3602/3602 PDBs pass (100%)  
**Results**: `data/residue_indices_test_results.json`

**Already Complete** - No action needed.

---

### Stage 1: Atoms ✅ COMPLETE (Partial Run)

**What it tests**: PDB atom parsing (atom_idx, name, coords, chain, residue)  
**Legacy JSON**: `data/json_legacy/pdb_atoms/<PDB>.json`  
**Status**: 2320/2320 tested pass (100% success rate, but test was interrupted)

**Existing Test**: `tests_python/integration/test_atoms_batch.py`

**To Complete Full Run**:
```bash
pytest tests_python/integration/test_atoms_batch.py -v --max-pdbs=3602
```

---

### Stage 2: LS Fitting ✅ 99.92% COMPLETE

**What it tests**: Least-squares fitting for base frame calculation  
**Legacy JSON**: `data/json_legacy/ls_fitting/<PDB>.json`  
**Status**: 3599/3602 pass (99.92%)  
**Edge Cases**: 
- 4KI4: Legacy has 30 duplicate records
- 5EAO, 5EAQ: CVC uses non-standard atom names (N01/C01 vs N1/C2)

**Existing Test**: `tests_python/integration/test_ls_fitting_batch.py`

---

### Stage 3: Distance Checks ⏳ TODO

**What it tests**: Geometric measurements between base pairs (dorg, dNN, plane_angle, d_v, overlap_area)  
**Legacy JSON**: `data/json_legacy/distance_checks/<PDB>.json` (4128 files)  
**Dependencies**: Stage 2 (ls_fitting/frames must match)

**Key Fields to Compare**:
```json
{
  "base_i": 1,            // First residue index (1-based legacy)
  "base_j": 2,            // Second residue index (1-based legacy)
  "dorg": 10.5,           // Origin distance
  "dNN": 8.2,             // N-N distance (N9 for purines, N1 for pyrimidines)
  "plane_angle": 5.3,     // Angle between base planes (degrees)
  "d_v": 0.8,             // Vertical distance
  "overlap_area": 2.1     // Base overlap area
}
```

**Comparison Strategy**:
- Match by `(base_i, base_j)` pair - normalize to `(min, max)`
- Tolerance: 1e-6 for all numeric values
- Count differences are OK (legacy records more pairs)

---

### Stage 4: H-bond List ⏳ TODO

**What it tests**: Hydrogen bond detection for base pairs  
**Legacy JSON**: `data/json_legacy/hbond_list/<PDB>.json` (4129 files)  
**Dependencies**: Stage 3 (need pairs to test)

**Key Fields to Compare**:
```json
{
  "base_i": 1,
  "base_j": 2,
  "num_hbonds": 3,
  "hbonds": [
    {
      "donor_atom": " N6 ",
      "acceptor_atom": " O4 ",
      "distance": 2.85,
      "type": "-"           // '-' standard, '*' non-standard, ' ' invalid
    }
  ]
}
```

**Comparison Strategy**:
- Match by `(base_i, base_j)` pair
- Compare `num_hbonds` count
- Compare each H-bond's donor, acceptor, distance (tolerance: 1e-6), type

---

### Stage 5: Pair Validation ⏳ TODO

**What it tests**: Validation results for each residue pair checked  
**Legacy JSON**: `data/json_legacy/pair_validation/<PDB>.json` (4128 files)  
**Dependencies**: Stages 3, 4

**Key Fields to Compare**:
```json
{
  "base_i": 1,
  "base_j": 2,
  "is_valid": 1,          // 0 or 1
  "bp_type_id": 2,        // -1, 0, 1, or 2
  "direction_vectors": {
    "dir_x": 0.5,
    "dir_y": 0.3,
    "dir_z": 0.1
  },
  "calculated_values": {
    "dorg": 10.5,
    "d_v": 0.8,
    "plane_angle": 5.3,
    "dNN": 8.2,
    "quality_score": 12.3
  }
}
```

---

### Stage 6: Find Bestpair Selection ⭐ PRIMARY OUTPUT ⏳ TODO

**What it tests**: **THE FINAL SELECTED BASE PAIRS** - this is the primary output!  
**Legacy JSON**: `data/json_legacy/find_bestpair_selection/<PDB>.json` (3921 files)  
**Dependencies**: All previous stages

**Key Fields to Compare**:
```json
{
  "num_bp": 25,
  "pairs": [[1, 50], [2, 49], [3, 48], ...]  // [base_i, base_j] pairs
}
```

**Comparison Strategy**:
- **MUST BE 100% MATCH** for production
- Compare `num_bp` count
- Compare all pairs exactly - normalize to `(min(i,j), max(i,j))`

---

### Stage 7: Base Pair Records ⏳ TODO

**What it tests**: Detailed base pair information (orientations, origins)  
**Legacy JSON**: `data/json_legacy/base_pair/<PDB>.json` (4129 files)  
**Dependencies**: Stage 6

**Key Fields to Compare**:
```json
{
  "basepair_idx": 0,
  "base_i": 1,
  "base_j": 50,
  "bp_type": "CG",
  "orien_i": [[...], [...], [...]],   // 3x3 rotation matrix
  "orien_j": [[...], [...], [...]],
  "org_i": [x, y, z],
  "org_j": [x, y, z],
  "dir_xyz": [dir_y, dir_z, 0.0]      // Note: Legacy bug - stores [dir_y, dir_z, 0]
}
```

---

### Stage 8: Step Parameters ⏳ TODO

**What it tests**: 6 step parameters for consecutive base pairs  
**Legacy JSON**: `data/json_legacy/bpstep_params/<PDB>.json` (3540 files)  
**Dependencies**: Stage 6 (need selected pairs)

**Key Fields to Compare**:
```json
{
  "bp_idx1": 0,
  "bp_idx2": 1,
  "shift": -0.5,
  "slide": -1.2,
  "rise": 3.3,
  "tilt": 2.1,
  "roll": 5.5,
  "twist": 32.0
}
```

---

### Stage 9: Helical Parameters ⏳ TODO

**What it tests**: Helical axis parameters  
**Legacy JSON**: `data/json_legacy/helical_params/<PDB>.json` (3540 files)  
**Dependencies**: Stage 6, 8

**Key Fields to Compare**:
```json
{
  "bp_idx1": 0,
  "bp_idx2": 1,
  "x_displacement": 0.1,
  "y_displacement": -0.2,
  "rise": 3.3,
  "inclination": 5.0,
  "tip": -2.0,
  "twist": 32.0
}
```

---

## Implementation Plan

### Phase 1: Create Unified Test Framework

Create a single test script that can validate any stage with these features:
- Stop on first failure
- Generate temp files, cleanup after each PDB (unless mismatch)
- Use parallel processing (20 workers)
- Save only mismatches for debugging
- Reuse `scripts/test_utils.py` infrastructure

**File**: `tests_python/integration/test_stage_validation.py`

```python
#!/usr/bin/env python3
"""
Unified stage-by-stage validation for legacy vs modern JSON comparison.

Usage:
    # Test a specific stage
    pytest tests_python/integration/test_stage_validation.py -v -k "stage3"
    
    # Run all stages in order (stops on first failure)
    pytest tests_python/integration/test_stage_validation.py -v --stop-on-failure
    
    # Test single PDB
    pytest tests_python/integration/test_stage_validation.py -v --pdb=1H4S
"""
```

### Phase 2: Implement Stage-Specific Comparisons

For each stage, implement a comparison function in `x3dna_json_compare/`:

1. **distance_checks_comparison.py**
2. **hbond_comparison.py** (may exist)
3. **pair_validation_comparison.py** 
4. **find_bestpair_comparison.py** (may exist)
5. **base_pair_comparison.py** (may exist)
6. **step_comparison.py** (may exist)
7. **helical_comparison.py**

### Phase 3: Run Full Validation

Execute each stage sequentially:

```bash
# Stage 0 - Already complete
# Stage 1 - Already complete  
# Stage 2 - Already complete (99.92%)

# Stage 3 - Distance Checks
pytest tests_python/integration/test_stage_validation.py -v -k "stage3" --stop-on-failure

# Stage 4 - H-bonds
pytest tests_python/integration/test_stage_validation.py -v -k "stage4" --stop-on-failure

# Stage 5 - Pair Validation
pytest tests_python/integration/test_stage_validation.py -v -k "stage5" --stop-on-failure

# Stage 6 - Find Bestpair Selection (PRIMARY OUTPUT)
pytest tests_python/integration/test_stage_validation.py -v -k "stage6" --stop-on-failure

# Stage 7 - Base Pair Records
pytest tests_python/integration/test_stage_validation.py -v -k "stage7" --stop-on-failure

# Stage 8 - Step Parameters
pytest tests_python/integration/test_stage_validation.py -v -k "stage8" --stop-on-failure

# Stage 9 - Helical Parameters
pytest tests_python/integration/test_stage_validation.py -v -k "stage9" --stop-on-failure
```

---

## File Management Strategy

### Temp File Creation
- Create temp directory for each PDB: `tmp/validation/<PDB>/`
- Generate legacy JSON (copy from `data/json_legacy/`) 
- Generate modern JSON to temp directory
- Compare

### Cleanup Rules
1. **Match**: Delete all temp files immediately
2. **Mismatch**: Keep files in `tmp/validation/mismatches/<PDB>/` for debugging
3. **Error**: Keep files with error log

### Disk Space Management
- Each comparison: ~100KB-10MB per PDB (depending on size)
- With cleanup: negligible (only keep mismatches)
- Without cleanup: potentially 10-50GB total

---

## Test Results Output

### Per-Stage Results File

Save to `data/validation_results/stage<N>_<name>_results.json`:

```json
{
  "stage_id": "stage3",
  "stage_name": "Distance Checks",
  "test_date": "2025-12-04 10:00:00",
  "total_pdbs": 3602,
  "passed": 3600,
  "failed": 2,
  "elapsed_seconds": 120.5,
  "failed_pdbs": [
    {
      "pdb_id": "1ABC",
      "issue": "Count mismatch: legacy=100, modern=98",
      "details": {...}
    }
  ],
  "summary": {
    "pass_rate": 99.94
  }
}
```

---

## Quick Reference Commands

### Run Full Validation (All Stages)
```bash
cd /Users/jyesselman2/Dropbox/2_code/cpp/find_pair_2
source venv/bin/activate
pytest tests_python/integration/ -v --stop-on-failure
```

### Run Single Stage
```bash
pytest tests_python/integration/test_stage_validation.py -v -k "stage6"
```

### Run Single PDB
```bash
pytest tests_python/integration/test_stage_validation.py -v --pdb=1H4S
```

### Check Current Status
```bash
# See what stages are passing
cat data/validation_results/*.json | jq '.summary.pass_rate'
```

### Generate Modern JSON for Testing
```bash
./build/generate_modern_json data/pdb/1H4S.pdb tmp/ --stage=distance_checks
```

---

## Known Issues to Investigate

### Stage 6 (Find Bestpair Selection) Known Mismatches
From previous testing:
- **3G8T**: Missing {(92, 160), (946, 947)}, Extra {(160, 975), (941, 947)}
- **6CAQ**: Missing {(75, 78), (968, 1024)}, Extra {(1024, 1188), (75, 79), (1063, 1072)}

**Root Cause Hypothesis**:
1. Quality score calculation differences
2. Tie-breaking when multiple pairs have same quality_score
3. `adjust_pairQuality` or `bp_type_id` calculation differences

---

## Success Criteria

### For Production Release
1. ✅ Stage 0 (Residue Indices): 100% match
2. ✅ Stage 1 (Atoms): 100% match
3. ✅ Stage 2 (LS Fitting): 99.92% match (3 edge cases documented)
4. ⏳ Stage 3 (Distance Checks): Must be 100% match
5. ⏳ Stage 4 (H-bonds): Must be 100% match
6. ⏳ Stage 5 (Pair Validation): Must be 100% match
7. ⏳ **Stage 6 (Find Bestpair Selection): MUST BE 100% MATCH** ⭐
8. ⏳ Stage 7 (Base Pair Records): Must be 100% match
9. ⏳ Stage 8 (Step Parameters): Must be 100% match
10. ⏳ Stage 9 (Helical Parameters): Must be 100% match

### Tolerance Values
- **Coordinates**: ±1e-6 Å
- **Distances**: ±1e-6 Å
- **Angles**: ±1e-6°
- **Rotation matrices**: < 0.0001 per element

---

## Appendix: Existing Test Infrastructure

### Scripts Available
- `scripts/test_utils.py` - Common utilities (find_executables, generate_json, cleanup)
- `scripts/compare_json.py` - CLI comparison tool
- `tests_python/conftest.py` - pytest fixtures

### Pytest Integration Tests
- `tests_python/integration/test_atoms_batch.py`
- `tests_python/integration/test_frames_batch.py`
- `tests_python/integration/test_ls_fitting_batch.py`
- `tests_python/integration/test_ls_fitting.py`
- `tests_python/integration/test_residue_indices_batch.py`

### x3dna_json_compare Module ✅ RESTORED
All module files restored from Dropbox backup (Dec 5, 2025):
- `__init__.py` - Module exports
- `atom_comparison.py` - Atom-level JSON comparison
- `base_pair_comparison.py` - Base pair record comparison
- `config.py` - Configuration settings
- `distance_comparison.py` - Distance checks comparison
- `find_bestpair_comparison.py` - Final pair selection comparison
- `frame_comparison.py` - Frame calculation comparison
- `hbond_comparison.py` - Hydrogen bond comparison
- `json_comparison.py` - Core JSON comparison engine
- `json_file_finder.py` - File discovery utilities
- `json_validator.py` - JSON schema validation
- `models.py` - Data models
- `pair_comparison.py` - Pair-level comparison
- `pair_validation_comparison.py` - Pair validation comparison
- `parallel_executor.py` - Parallel execution utilities
- `pdb_utils.py` - PDB file utilities
- `residue_indices_comparison.py` - Residue index comparison
- `result_cache.py` - Result caching
- `selection_comparison.py` - Selection comparison
- `step_comparison.py` - Step parameter comparison



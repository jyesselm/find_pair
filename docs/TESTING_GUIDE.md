# X3DNA Testing Guide

## Overview

This guide explains how to validate that the modern C++ implementation produces identical results to the legacy C code. Testing is done by comparing JSON outputs at each stage of the algorithm.

**IMPORTANT**: Testing is done in stages. Each stage must pass completely before proceeding to the next.

---

## Quick Start

```bash
# Run all tests for all fast PDBs
fp2-validate validate all --test-set 100

# Run a specific stage
fp2-validate validate 1 --test-set 100  # Stage 1: atoms

# Test a single PDB with verbose output
fp2-validate validate --pdb 1EHZ -v

# Stop on first failure (recommended for debugging)
fp2-validate validate 3 --test-set 100 -v -s
```

---

## Test Environment Setup

### Prerequisites

1. **Build both executables:**
   ```bash
   # Build modern executable
   mkdir -p build && cd build
   cmake .. && make -j
   cd ..
   
   # Build legacy executable
   cd org
   mkdir -p build && cd build
   cmake .. && make -j
   cd ../..
   ```

2. **Verify executables exist:**
   ```bash
   fp2-validate info
   ```
   
   Expected output:
   ```
   X3DNA Validation Environment
   ========================================
   Project root: /path/to/find_pair_2
   Legacy executable: ✅ org/build/bin/find_pair_analyze
   Modern executable: ✅ build/generate_modern_json
   Fast PDBs: 3602 available
   ```

3. **PDB files must be present in `data/pdb/`**

4. **Legacy JSON files should be pre-generated in `data/json_legacy/`**

---

## Testing Stages

The algorithm has **12 stages** that must be validated in order:

| Stage | CLI ID | Name | JSON Type | Status |
|-------|--------|------|-----------|--------|
| 1 | `pdb_atoms` | Atom Parsing | `pdb_atoms` | ✅ PASSED (3602/3602) |
| 2 | `residue_indices` | Residue Index Mapping | `residue_indices` | ✅ PASSED (3602/3602) |
| 3 | `base_frame_calc` | Base Frame Calculation | `base_frame_calc` | ✅ PASSED (3602/3602) |
| 4 | `ls_fitting` | Least Squares Fitting | `ls_fitting` | ✅ PASSED (3602/3602) |
| 5 | `frame_calc` | Reference Frame Calculation | `frame_calc` | ✅ PASSED (3602/3602) |
| 6 | `pair_validation` | Pair Validation | `pair_validation` | ⚠️ TESTING |
| 7 | `distance_checks` | Distance Measurements | `distance_checks` | ⚠️ TESTING |
| 8 | `hbond_list` | Hydrogen Bond List | `hbond_list` | ⏳ PENDING |
| 9 | `base_pair` | Base Pair Records | `base_pair` | ⏳ PENDING |
| 10 | `find_bestpair_selection` | Final Pair Selection | `find_bestpair_selection` | ⏳ PENDING |
| 11 | `bpstep_params` | Step Parameters | `bpstep_params` | ⏳ PENDING |
| 12 | `helical_params` | Helical Parameters | `helical_params` | ⏳ PENDING |

### Stage Groups

For convenience, stages can be referenced by group:

| Group | Stages | Description |
|-------|--------|-------------|
| `atoms` | 1 | Atom parsing only |
| `residue` | 2 | Residue indices only |
| `frames` | 3, 4, 5 | All frame calculations |
| `pairs` | 6, 7, 9, 10 | Pair validation and selection |
| `hbonds` | 8 | H-bond detection |
| `steps` | 11, 12 | Step and helical parameters |
| `all` | 1-12 | All stages |

---

## Running Validation

### Method 1: fp2-validate CLI (Recommended)

The `fp2-validate` command is the primary testing interface.

```bash
# Validate by stage number
fp2-validate validate 1 --test-set 100      # Stage 1
fp2-validate validate 3 4 5 --test-set 100  # Stages 3-5

# Validate by stage name
fp2-validate validate pdb_atoms --test-set 100
fp2-validate validate base_frame_calc --pdb 1EHZ

# Validate by group
fp2-validate validate frames --test-set 100
fp2-validate validate pairs --test-set 100

# Options
fp2-validate validate 3 --test-set 100 -v         # Verbose output
fp2-validate validate 3 --test-set 100 -s         # Stop on first failure
fp2-validate validate 3 --test-set 100 --diff     # Document differences
```

### Method 2: pytest

For CI/CD integration or detailed test output:

```bash
# Run all integration tests
pytest tests_python/integration/ -v

# Run specific stage tests
pytest tests_python/integration/test_stage_validation.py -v -k "stage3"

# Run with pytest options
pytest tests_python/integration/ -v -x --max-pdbs=50
```

### Method 3: Direct Script Execution

```bash
# Run stage validation directly
python tests_python/integration/test_stage_validation.py stage3 --max-pdbs 100

# Run specific test scripts
python tests_python/integration/test_atoms_batch.py
python tests_python/integration/test_frames_batch.py
```

---

## Validation Workflow

### Recommended Testing Sequence

1. **Start with Stage 1 (Atoms)**
   ```bash
   fp2-validate validate 1 --test-set 10 -v
   ```
   This validates that PDB parsing produces identical atoms.

2. **Move to Stage 2 (Residue Indices)**
   ```bash
   fp2-validate validate 2 --test-set 10 -v
   ```
   Validates residue-to-atom index mapping.

3. **Validate Frames (Stages 3-5)**
   ```bash
   fp2-validate validate frames --test-set 10 -v
   ```
   This is critical - frame calculations affect all downstream stages.

4. **Continue through remaining stages in order**

### Debugging Failures

When a stage fails:

1. **Identify the failing PDB:**
   ```bash
   fp2-validate validate 3 --test-set 100 -v -s
   ```

2. **Get detailed comparison:**
   ```bash
   fp2-validate compare <PDB_ID> --verbose
   ```

3. **Inspect JSON files directly:**
   ```bash
   # Legacy JSON
   cat data/json_legacy/base_frame_calc/<PDB_ID>.json | jq .

   # Modern JSON
   cat data/json/base_frame_calc/<PDB_ID>.json | jq .
   ```

4. **Re-generate modern JSON:**
   ```bash
   ./build/generate_modern_json data/pdb/<PDB_ID>.pdb data/json --stage=frames
   ```

---

## What Gets Compared

### Tolerance Values

| Type | Tolerance | Notes |
|------|-----------|-------|
| Coordinates | 1e-6 | xyz values |
| Distances | 1e-6 | dorg, dNN, d_v |
| Angles | 1e-6 | plane_angle |
| Matrix elements | 1e-4 | Rotation matrices |
| RMS fit | 0.001 | Template fitting |

### Comparison Fields by Stage

#### Stage 1: pdb_atoms
- `num_atoms` (count must match exactly)
- `atom_idx` / `legacy_atom_idx` (must match)
- `xyz` coordinates (within tolerance)
- `atom_name`, `residue_name`, `chain_id` (exact match)

#### Stage 2: residue_indices
- Residue key `(chain_id, residue_seq, insertion)` matching
- `start_atom_idx`, `end_atom_idx` (atom ranges)
- `legacy_residue_idx` (must match)

#### Stages 3-5: Frame Calculations
- `rms_fit` (within tolerance)
- `num_matched_atoms` (exact match)
- `matched_atoms` list (set equality)
- `template_file` (filename only, not path)
- Rotation matrix elements (within tolerance)
- Origin/translation vectors (within tolerance)

#### Stages 6-7: Pair Validation
- `is_valid` flag (exact match)
- `bp_type_id` (exact match)
- Geometric values: `dorg`, `dNN`, `plane_angle`, `d_v` (within tolerance)

#### Stage 8: H-bond List
- `num_hbonds` (exact match)
- Individual H-bond: `donor_atom`, `acceptor_atom`, `distance`

#### Stages 9-10: Base Pairs
- `num_bp` (exact match)
- Pair indices (as normalized pairs)
- Rotation matrices and origins for each pair

#### Stages 11-12: Parameters
- Step parameters: shift, slide, rise, tilt, roll, twist
- Helical parameters: x/y displacement, rise, inclination, tip, twist

---

## Test Sets

Pre-defined test sets for consistent testing:

| Size | PDBs | Use Case |
|------|------|----------|
| 10 | 10 | Quick sanity check |
| 50 | 50 | Fast development iteration |
| 100 | 100 | Standard testing |
| 500 | 500 | Thorough validation |
| 1000 | 1000 | Extended validation |
| all | 3602 | Full validation suite |

```bash
# Use test sets
fp2-validate validate 3 --test-set 10
fp2-validate validate 3 --test-set 100
fp2-validate validate 3 --test-set 1000
```

---

## Test Output Files

### Results Location

Test results are saved to `data/validation_results/`:

```
data/validation_results/
├── stage1_pdb_atoms_pytest.json
├── stage2_residue_indices_pytest.json
├── stage3_distance_checks_pytest.json
├── ...
└── validation_summary.json
```

### Result File Format

```json
{
  "stage_id": "stage3",
  "stage_name": "Distance Checks",
  "test_date": "2025-12-06T10:30:00",
  "total_pdbs": 100,
  "passed": 98,
  "failed": 2,
  "skipped": 0,
  "pass_rate": 98.0,
  "elapsed_seconds": 45.2,
  "failed_pdbs": [
    {
      "pdb_id": "1ABC",
      "errors": ["dorg mismatch: 0.01 vs 0.02"]
    }
  ]
}
```

---

## Common Issues & Solutions

### Issue: "Legacy JSON not found"

**Solution**: Generate legacy JSON first:
```bash
cd org
./build/bin/find_pair_analyze ../data/pdb/<PDB_ID>.pdb
cd ..
```

### Issue: "Modern executable not found"

**Solution**: Build the project:
```bash
cd build
cmake .. && make -j
```

### Issue: RMS mismatch

**Possible causes**:
- Different atom matching order
- Modified nucleotide handling
- Template file differences

**Debug**:
```bash
fp2-validate compare <PDB_ID> --verbose
```

### Issue: residue_idx mismatch

**Note**: Legacy uses 1-based indexing, modern uses 0-based. The modern JSON should include `legacy_residue_idx` for direct comparison.

---

## CI/CD Integration

### GitHub Actions Example

```yaml
- name: Run validation tests
  run: |
    fp2-validate validate all --test-set 100 --quiet
```

### Exit Codes

- `0`: All tests passed
- `1`: One or more tests failed

---

## Expected Failures: Modified Residues

Some PDBs will have validation failures due to **modified nucleotides**. These are **not bugs** in modern code - they reflect inconsistencies in legacy code's classification system.

### What Are Modified Residues?

Modified nucleotides are non-standard bases with 3-letter codes like `9DG`, `A23`, `5MU`. They're identified by lowercase `base_type` (e.g., `a`, `g`, `c`, `t`, `u`).

### Why Do They Cause Failures?

Legacy code has inconsistent classification:
1. **Atom-based** (RY): Purine if N7/C8/N9 present
2. **Template-based** (base_type): From RMSD template matching

These can conflict, causing different N atoms to be used for dNN calculation:
- Legacy: Uses RY → might select wrong atom (e.g., C9 instead of N1)
- Modern: Uses base_type consistently → selects correct atom

### Known Problematic Residues

| Code | Name | Issue |
|------|------|-------|
| **9DG** | 9-Deazaguanine | Has N7/C8 but no N9; legacy uses C9 for dNN |
| **A23** | Adenosine cyclic phosphate | 9-atom RMSD fails; legacy classifies as pyrimidine |
| **EPE** | HEPES buffer | Not a nucleotide but gets processed as one |

### Affected PDBs

PDBs with known dNN mismatches due to modified residues:
- 1Q2R, 1Q2S, 7NQ4, 8OMR (contain 9DG)
- 2XD0, 2XDD (contain A23)
- 4E8R (contains EPE)

### Documenting Modified Residues

Generate a comprehensive report:
```bash
python3 tools/document_modified_residues.py --workers 20
# Output: data/modified_residues.json
```

📖 **See [docs/MODIFIED_RESIDUES.md](MODIFIED_RESIDUES.md) for complete documentation.**

---

## Adding New Comparison Logic

New stage comparisons should be added to:
- `tests_python/integration/test_stage_validation.py` - Main validation script
- `x3dna_json_compare/` - Comparison modules

### Template for New Stage Comparison

```python
def compare_new_stage(legacy_records: List[Dict], modern_records: List[Dict],
                      tolerance: float = 1e-6) -> Tuple[bool, List[str]]:
    """Compare new_stage records."""
    errors = []
    
    # Build lookup maps
    legacy_by_key = {make_key(r): r for r in legacy_records}
    modern_by_key = {make_key(r): r for r in modern_records}
    
    # Find common keys
    common_keys = set(legacy_by_key.keys()) & set(modern_by_key.keys())
    
    # Compare
    for key in common_keys:
        leg_rec = legacy_by_key[key]
        mod_rec = modern_by_key[key]
        
        # Compare fields...
        if leg_rec.get('field') != mod_rec.get('field'):
            errors.append(f"Key {key} field mismatch")
    
    return len(errors) == 0, errors
```

---

## Contact & Support

For issues with testing:
1. Check this guide first
2. Review `docs/JSON_DATA_TYPES_AND_COMPARISONS.md`
3. Check `docs/MODIFIED_RESIDUES.md` for modified nucleotide issues
4. Check existing comparison modules in `x3dna_json_compare/`

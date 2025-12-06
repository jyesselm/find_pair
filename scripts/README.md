# Scripts Directory

This directory contains Python scripts for JSON management and specialized analysis.

**For validation, use the unified CLI tool: `fp2-validate`**

## Unified Validation CLI

All validation tasks are now handled by a single tool:

```bash
# Install the CLI (from project root)
pip install -e .

# Basic usage
fp2-validate                      # All stages, all fast PDBs
fp2-validate frames               # Frames only
fp2-validate hbonds atoms         # Multiple stages
fp2-validate --pdb 1EHZ -v        # Single PDB, verbose
fp2-validate --stop-on-first      # Stop at first failure (debug mode)
fp2-validate --diff               # Document differences to file
fp2-validate --max 100            # First 100 PDBs only
fp2-validate --test-set 100       # Use saved test set
fp2-validate --quiet              # Exit code only (for CI)

# Stage-specific aliases
fp2-validate atoms                # Alias for: fp2-validate validate atoms
fp2-validate frames               # Alias for: fp2-validate validate frames
fp2-validate hbonds               # Alias for: fp2-validate validate hbonds
fp2-validate pairs                # Alias for: fp2-validate validate pairs
fp2-validate steps                # Alias for: fp2-validate validate steps

# Info commands
fp2-validate info                 # Show environment info
fp2-validate list-pdbs            # List available PDBs
```

## Core Scripts (Still Active)

### JSON Management
- **`rebuild_json.py`** - Generate/regenerate JSON files
  ```bash
  python3 scripts/rebuild_json.py regenerate --modern-only
  python3 scripts/rebuild_json.py regenerate --legacy-only
  ```

### Data Management
- **`download_pdbs.py`** - Download PDB files from RCSB
- **`create_fast_pdbs_json.py`** - Generate valid_pdbs_fast.json
- **`find_slow_pdbs.py`** - Identify slow-to-process PDBs

### Specialized Analysis
- **`compare_rmsd_calculation.py`** - RMSD comparison debugging
- **`manual_rmsd_calc.py`** - Manual RMSD calculations
- **`check_nucleotide_types.py`** - Nucleotide type analysis
- **`find_modified_nucleotides.py`** - Modified nucleotide detection

## Cluster Computing

See `cluster/README.md` for cluster-based validation:

```bash
python3 scripts/cluster/run_cluster_comparison.py
```

## Archived Scripts

**December 6, 2025**: Scripts consolidated into unified CLI.

Archived scripts are in `archive/`:
- `archive/validation/` - Old validate_*.py scripts → use `fp2-validate`
- `archive/testing/` - Old test_*.py scripts → use `fp2-validate`
- `archive/debug/` - One-off debugging scripts

## Migration from Old Scripts

| Old Script | New Command |
|------------|-------------|
| `validate_hbonds_batch.py` | `fp2-validate hbonds` |
| `validate_hbonds_200.py` | `fp2-validate hbonds --max 200` |
| `validate_frames_batch.py` | `fp2-validate frames` |
| `test_ls_fitting_stop_on_mismatch.py` | `fp2-validate frames --stop-on-first` |
| `test_single_pdb_ls_fitting.py 1EHZ` | `fp2-validate frames --pdb 1EHZ -v` |
| `test_all_stages_batch.py` | `fp2-validate` |
| `find_first_mismatch.py` | `fp2-validate --stop-on-first` |

## Directory Structure

```
scripts/
├── README.md                # This file
├── rebuild_json.py          # JSON generation
├── download_pdbs.py         # PDB download
├── create_fast_pdbs_json.py # Data prep
├── find_slow_pdbs.py        # Data categorization
├── compare_rmsd_calculation.py  # Specialized debugging
├── manual_rmsd_calc.py          # Specialized debugging
├── check_nucleotide_types.py    # Analysis
├── find_modified_nucleotides.py # Analysis
├── cluster/                 # Cluster computing
│   ├── README.md
│   └── *.py
└── archive/                 # Archived (consolidated into CLI)
    ├── validation/          # validate_*.py
    ├── testing/             # test_*.py
    └── debug/               # debug_*.py
```

# Scripts Directory

This directory contains Python scripts for comparing, analyzing, and managing JSON data from legacy and modern codebases.

## Core Tools (Active)

These are the main scripts that should be used for regular work:

### Comparison
- **`compare_json.py`** ⭐ - Main comparison tool for comparing legacy and modern JSON outputs
  ```bash
  python3 scripts/compare_json.py compare
  python3 scripts/compare_json.py compare 1H4S
  python3 scripts/compare_json.py compare --test-set 100
  ```

### JSON Management
- **`rebuild_json.py`** ⭐ - Unified tool for regenerating, validating, and cleaning JSON files
  ```bash
  # Regenerate all JSON
  python3 scripts/rebuild_json.py regenerate
  
  # Regenerate only modern JSON
  python3 scripts/rebuild_json.py regenerate --modern-only
  
  # Regenerate only legacy JSON
  python3 scripts/rebuild_json.py regenerate --legacy-only
  
  # Validate JSON files
  python3 scripts/rebuild_json.py validate
  
  # Clean invalid JSON files
  python3 scripts/rebuild_json.py clean --execute
  ```

### Data Management
- **`download_pdbs.py`** - Download PDB files from RCSB database
  ```bash
  python3 scripts/download_pdbs.py --pdb-ids 1H4S 2BNA
  python3 scripts/download_pdbs.py --list-file pdb_list.txt
  ```
  See `download_pdbs_README.md` for details.

## Cluster Computing

### `cluster/` - Cluster Computing Scripts

Scripts for running large-scale comparisons on computing clusters:
- `run_cluster_comparison.py` - Submit comparison jobs to cluster
- `run_single_batch.py` - Process a single batch
- `aggregate_results.py` - Collect and summarize results
- See `cluster/README.md` and `cluster/QUICKSTART.md` for details

## Archived Scripts

**December 2, 2025**: Major cleanup performed to reduce clutter.

The `archive/debugging/` directory contains ~40 one-off analysis, investigation, and debugging scripts that were used during development but are no longer actively maintained. These are kept for reference:

- **Analysis scripts** (`analyze_*.py`) - Frame/index/validation analysis
- **Investigation scripts** (`investigate_*.py`) - Pair difference investigations
- **Comparison scripts** (`compare_*.py`) - Specific comparison tools
- **Debug scripts** (`debug_*.py`) - Debugging utilities
- **Extract scripts** (`extract_*.py`) - Data extraction tools
- **Validate/Verify scripts** (`validate_*.py`, `verify_*.py`) - One-off validation

**Note**: If you need functionality from archived scripts, check the core tools first (`compare_json.py` and `rebuild_json.py`), as they have comprehensive features and are actively maintained.

## Directory Structure

```
scripts/
├── README.md                    # This file
├── compare_json.py ⭐           # Main comparison tool
├── rebuild_json.py ⭐           # JSON management tool
├── download_pdbs.py             # PDB download utility
├── download_pdbs_README.md      # Download script docs
├── cluster/                     # Cluster computing scripts
│   ├── README.md
│   ├── QUICKSTART.md
│   └── *.py
└── archive/debugging/           # Archived debugging scripts (~40 files)
    ├── analyze_*.py
    ├── compare_*.py
    ├── debug_*.py
    ├── extract_*.py
    ├── investigate_*.py
    ├── validate_*.py
    ├── verify_*.py
    └── *.sh
```

## Usage Examples

### Quick Comparison
```bash
# Compare all available PDBs
python3 scripts/compare_json.py compare

# Compare specific PDB
python3 scripts/compare_json.py compare 1H4S

# Compare with verbose output
python3 scripts/compare_json.py compare 1H4S --verbose
```

### Regenerate JSON
```bash
# Regenerate modern JSON for a specific PDB
python3 scripts/rebuild_json.py regenerate 1H4S --modern-only

# Regenerate for test set
python3 scripts/rebuild_json.py regenerate --test-set 100
```

For more detailed documentation, see [docs/TESTING_GUIDE.md](../docs/TESTING_GUIDE.md).

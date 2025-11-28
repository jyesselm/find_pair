# Scripts Directory

This directory contains Python scripts for comparing, analyzing, and managing JSON data from legacy and modern codebases.

## Core Tools (Active)

These are the main scripts that should be used for regular work:

### Comparison
- **`compare_json.py`** - Main comparison tool for comparing legacy and modern JSON outputs
  ```bash
  python3 scripts/compare_json.py compare
  python3 scripts/compare_json.py compare 1H4S
  python3 scripts/compare_json.py compare --test-set 100
  ```

### JSON Management
- **`rebuild_json.py`** - Unified tool for regenerating, validating, and cleaning JSON files
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

## Archived Scripts

The `archive/` directory contains one-off analysis, investigation, and debugging scripts that were used during development but are no longer actively maintained. These are kept for reference:

- Analysis scripts (`analyze_*.py`)
- Investigation scripts (`investigate_*.py`)
- Comparison scripts (other than `compare_json.py`)
- Debug scripts (`debug_*.py`)
- Utility scripts for specific tasks

If you need functionality from archived scripts, consider using the core tools first (`compare_json.py` and `rebuild_json.py`), as they have comprehensive features.

## Directory Structure

```
scripts/
├── README.md              # This file
├── compare_json.py        # Main comparison tool (ACTIVE)
├── rebuild_json.py        # JSON management tool (ACTIVE)
└── archive/               # Archived one-off scripts
    ├── analyze_*.py
    ├── investigate_*.py
    ├── compare_*.py
    └── ...
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

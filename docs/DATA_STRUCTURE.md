# Data Directory Structure

**Last Updated**: 2025-01-XX

---

## Overview

The data directory is organized to separate modern JSON files, legacy JSON files, and comparison data.

---

## Directory Structure

```
data/
├── config/              # Configuration files
│   └── special_residues.json
├── pdb/                 # PDB files (keep as-is)
│   └── *.pdb
├── templates/           # Template files (keep as-is)
│   └── Atomic_*.pdb
├── test_sets/          # Test set definitions
│   ├── test_set_10.json
│   ├── test_set_50.json
│   ├── test_set_100.json
│   ├── test_set_500.json
│   └── test_set_1000.json
├── json/               # Modern JSON files (all modern output)
│   ├── pdb_atoms/     # PDB atom records
│   │   └── *.json
│   ├── base_frame_calc/  # Base frame calculation records
│   │   └── *.json
│   ├── frame_calc/    # Frame calculation records
│   │   └── *.json
│   ├── base_pair/     # Base pair records
│   │   └── *.json
│   ├── pair_validation/  # Pair validation records
│   │   └── *.json
│   ├── distance_checks/  # Distance check records
│   │   └── *.json
│   ├── hbond_list/    # Hydrogen bond list records
│   │   └── *.json
│   ├── find_bestpair_selection/  # Best pair selection records
│   │   └── *.json
│   └── atoms/         # Atom comparison data (modern)
│       └── *.json
└── json_legacy/       # Legacy JSON files (reference for comparison)
    ├── pdb_atoms/     # PDB atom records
    │   └── *.json
    ├── base_frame_calc/  # Base frame calculation records
    │   └── *.json
    ├── frame_calc/    # Frame calculation records
    │   └── *.json
    ├── base_pair/     # Base pair records
    │   └── *.json
    ├── pair_validation/  # Pair validation records
    │   └── *.json
    ├── distance_checks/  # Distance check records
    │   └── *.json
    ├── hbond_list/    # Hydrogen bond list records
    │   └── *.json
    ├── find_bestpair_selection/  # Best pair selection records
    │   └── *.json
    └── atoms/         # Atom comparison data (legacy)
        └── *.json
```

---

## File Locations

### Modern JSON Files

**Segmented JSON files**: Each record type in its own directory
- `data/json/pdb_atoms/<PDB_ID>.json` - PDB atom records
- `data/json/base_frame_calc/<PDB_ID>.json` - Base frame calculation records
- `data/json/frame_calc/<PDB_ID>.json` - Frame calculation records
- `data/json/base_pair/<PDB_ID>.json` - Base pair records
- `data/json/pair_validation/<PDB_ID>.json` - Pair validation records
- `data/json/distance_checks/<PDB_ID>.json` - Distance check records
- `data/json/hbond_list/<PDB_ID>.json` - Hydrogen bond list records
- `data/json/find_bestpair_selection/<PDB_ID>.json` - Best pair selection records

**Atom comparison data**: `data/json/atoms/<PDB_ID>.json`
- Atom-specific comparison data (if needed for detailed analysis)

### Legacy JSON Files

**Segmented JSON files**: Each record type in its own directory (same structure as modern)
- `data/json_legacy/pdb_atoms/<PDB_ID>.json` - PDB atom records
- `data/json_legacy/base_frame_calc/<PDB_ID>.json` - Base frame calculation records
- `data/json_legacy/frame_calc/<PDB_ID>.json` - Frame calculation records
- `data/json_legacy/base_pair/<PDB_ID>.json` - Base pair records
- `data/json_legacy/pair_validation/<PDB_ID>.json` - Pair validation records
- `data/json_legacy/distance_checks/<PDB_ID>.json` - Distance check records
- `data/json_legacy/hbond_list/<PDB_ID>.json` - Hydrogen bond list records
- `data/json_legacy/find_bestpair_selection/<PDB_ID>.json` - Best pair selection records

**Atom comparison data**: `data/json_legacy/atoms/<PDB_ID>.json`
- Legacy atom-specific comparison data (if needed for detailed analysis)

---

## Usage

### Generate Modern JSON

```bash
# Generate modern JSON for a PDB (creates segmented files in separate directories)
./build/generate_modern_json data/pdb/1H4S.pdb data/json/
# Output: 
#   data/json/pdb_atoms/1H4S.json
#   data/json/base_pair/1H4S.json
#   data/json/pair_validation/1H4S.json
#   etc.
```

### Generate Legacy JSON

```bash
# Generate legacy JSON for a PDB (creates segmented files in separate directories)
cd org
./build/bin/find_pair_analyze ../data/pdb/1H4S.pdb
# Output: 
#   data/json_legacy/pdb_atoms/1H4S.json
#   data/json_legacy/base_pair/1H4S.json
#   data/json_legacy/pair_validation/1H4S.json
#   etc.
```

### Compare JSON Files

```bash
# Compare all PDBs
python3 scripts/compare_json.py compare

# Compare specific PDB
python3 scripts/compare_json.py compare 1H4S
```

---


**structure**:
- All modern JSON files in `data/json/` organized by record type in separate directories
- All legacy JSON files in `data/json_legacy/` organized by record type in separate directories
- No main/combined JSON files - only segmented files in type-specific directories
- Each record type has its own directory (e.g., `data/json/pdb_atoms/`, `data/json/base_pair/`)
- Atom comparison data in `atoms/` subdirectory if needed

---

*This structure keeps modern and legacy JSON files clearly separated while maintaining a clean organization.*


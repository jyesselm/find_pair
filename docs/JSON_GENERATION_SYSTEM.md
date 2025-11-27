# JSON Generation System

**Purpose**: Automatically detect and generate missing modern JSON files to prevent timeouts and missing data issues.

## Problem

When testing PDBs, scripts may timeout or fail if modern JSON files are missing. Large PDBs like 2QEX (8.3 MB, 107k lines) require modern JSON to be generated before comparison scripts can run efficiently.

## Solution

Two scripts provide comprehensive JSON file management:

### 1. `scripts/check_and_generate_json.py` (Recommended)

**Comprehensive system** for checking status and generating missing JSON files.

**Usage**:
```bash
# Check status of all JSON files
python3 scripts/check_and_generate_json.py --check

# Generate all missing modern JSON files
python3 scripts/check_and_generate_json.py --generate

# Generate for specific PDBs
python3 scripts/check_and_generate_json.py --generate 2QEX 1T0K

# Generate with limit (for testing)
python3 scripts/check_and_generate_json.py --generate --limit 10
```

**Features**:
- Checks all record types (find_bestpair_selection, base_pair, pair_validation, etc.)
- Reports missing files by type
- Handles large files with configurable timeouts
- Progress reporting and error handling

### 2. `scripts/generate_missing_modern_json.py`

**Focused script** for generating missing modern JSON files.

**Usage**:
```bash
# Dry run (list what would be generated)
python3 scripts/generate_missing_modern_json.py --dry-run

# Generate all missing
python3 scripts/generate_missing_modern_json.py

# Generate with limit
python3 scripts/generate_missing_modern_json.py --limit 5

# Generate specific PDBs
python3 scripts/generate_missing_modern_json.py 2QEX 1T0K
```

## Workflow

### Before Running Tests

1. **Check status**:
   ```bash
   python3 scripts/check_and_generate_json.py --check
   ```

2. **Generate missing files**:
   ```bash
   python3 scripts/check_and_generate_json.py --generate
   ```

3. **Run tests**:
   ```bash
   python3 scripts/analyze_mismatched_pairs.py 2QEX
   ```

### For Large PDBs

Large PDBs (like 2QEX) may take several minutes to generate. The scripts:
- Use configurable timeouts (default: 600s)
- Process sequentially to avoid memory issues
- Provide progress reporting

## Example Output

```
JSON File Status Report
======================================================================
PDBs with legacy JSON: 100
PDBs with modern JSON: 101
Available PDB files: 4934
Missing modern JSON: 0
  - With PDB file available: 0
  - Missing PDB file: 0

✅ No missing modern JSON files to generate!
```

## Integration with Test Scripts

Test scripts should check for missing JSON files before running:

```python
# In test scripts
from pathlib import Path

modern_json = Path(f'data/json/find_bestpair_selection/{pdb_id}.json')
if not modern_json.exists():
    print(f"Generating modern JSON for {pdb_id}...")
    subprocess.run(['build/generate_modern_json', 
                    f'data/pdb/{pdb_id}.pdb', 'data/json'])
```

## Benefits

1. **Prevents timeouts**: Ensures JSON files exist before running comparisons
2. **Comprehensive coverage**: Checks all record types, not just find_bestpair_selection
3. **Error handling**: Reports failures and continues processing
4. **Progress tracking**: Shows progress for batch operations
5. **Flexible**: Can generate all missing or specific PDBs

## Maintenance

Run periodically to ensure all PDBs with legacy JSON also have modern JSON:

```bash
# Weekly maintenance
python3 scripts/check_and_generate_json.py --check
python3 scripts/check_and_generate_json.py --generate
```


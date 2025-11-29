# Fix Indices from Legacy JSON Option

**Date**: 2025-01-XX  
**Status**: ✅ IMPLEMENTED

---

## Overview

The `--fix-indices` option allows you to fix residue legacy indices by matching with legacy JSON output. This is useful when comparing modern output with legacy output, as it ensures residues are matched correctly even if the initial parsing assigns different indices.

---

## Usage

### Basic Usage (Auto-detect JSON file)

```bash
./build/find_pair_app --fix-indices <pdb_file> [output_file]
```

The tool will automatically look for the legacy JSON file at:
```
data/json_legacy/base_frame_calc/<PDB_ID>.json
```

### Specify JSON File

```bash
./build/find_pair_app --fix-indices=data/json_legacy/base_frame_calc/6CAQ.json <pdb_file> [output_file]
```

### Combined with Other Options

```bash
./build/find_pair_app --fix-indices --legacy-mode 6CAQ.pdb
```

---

## How It Works

1. **Parse PDB**: Structure is parsed normally
2. **Match by PDB Properties**: Residues are matched by `(residue_name, chain_id, residue_seq, insertion)`
3. **Load Legacy JSON**: Legacy indices are loaded from the JSON file
4. **Assign Indices**: Legacy indices are assigned to matching residues
5. **Continue Processing**: Frame calculation and pair finding proceed with correct indices

---

## When to Use

### ✅ Use When:
- Comparing modern output with legacy output
- Debugging residue index mismatches
- Ensuring correct residue matching for validation
- Testing against legacy reference data

### ❌ Don't Use When:
- Running production code (should fix root cause instead)
- Legacy JSON is not available
- You want to test the actual parsing behavior

---

## Example

```bash
# Fix indices from legacy JSON and find pairs
./build/find_pair_app --fix-indices data/pdb/6CAQ.pdb output.inp

# Output:
# [INFO] Fixed 1519 residue indices from legacy JSON: data/json_legacy/base_frame_calc/6CAQ.json
# Parsing PDB file: data/pdb/6CAQ.pdb
# Finding base pairs...
# Found 1234 base pairs
```

---

## Implementation Details

### Command Line Parser
- `--fix-indices`: Enable fixing (auto-detect JSON file)
- `--fix-indices=FILE`: Enable fixing with specific JSON file

### FindPairProtocol
- `set_fix_indices_from_legacy_json(bool, string)`: Set the option
- Automatically fixes indices before frame calculation
- Uses `io::fix_residue_indices_from_json()` helper function

### Residue Index Fixer
- Matches residues by PDB properties
- Assigns legacy indices from JSON
- Returns count of fixed residues

---

## Related Documentation

- [PDB_PROPERTIES_MATCHING_APPROACH.md](PDB_PROPERTIES_MATCHING_APPROACH.md) - Detailed approach
- [RESIDUE_INDEXING_STATUS.md](RESIDUE_INDEXING_STATUS.md) - Current indexing issues
- [RESIDUE_GROUPING_FIX.md](RESIDUE_GROUPING_FIX.md) - Residue grouping fix

---

*This option provides a reliable way to match residues when comparing with legacy output.*


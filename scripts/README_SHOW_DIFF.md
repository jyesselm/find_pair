# show_json_diff.py - JSON Comparison Tool

A tool to show detailed differences between generated and legacy JSON files.

## Usage

```bash
# Basic comparison
python3 scripts/show_json_diff.py <pdb_id>

# Verbose mode (shows debug fields and field differences)
python3 scripts/show_json_diff.py <pdb_id> --verbose

# Show only missing atoms
python3 scripts/show_json_diff.py <pdb_id> --missing-only

# Show only extra atoms
python3 scripts/show_json_diff.py <pdb_id> --extra-only
```

## Examples

```bash
# Compare 1LNT
python3 scripts/show_json_diff.py 1LNT

# Verbose comparison showing debug info
python3 scripts/show_json_diff.py 1NUJ --verbose

# Show only missing atoms
python3 scripts/show_json_diff.py 1ZZN --missing-only
```

## Output

The tool shows:
- **Summary**: Total atom counts and differences
- **Missing atoms**: Atoms in legacy but not in generated JSON
- **Extra atoms**: Atoms in generated but not in legacy JSON
- **Field differences**: Differences in common atoms (with --verbose)

## Debug Information

When verbose mode is enabled, the tool shows:
- `line_number`: Line number in PDB file
- `atom_serial`: Atom serial number from PDB
- `model_number`: Model number (for multi-model PDBs)
- `occupancy`: Occupancy value
- `original_atom_name`: Original atom name before normalization
- `original_residue_name`: Original residue name before normalization

These fields are automatically included in the generated JSON for debugging purposes.


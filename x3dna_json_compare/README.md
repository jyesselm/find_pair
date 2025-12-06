# X3DNA JSON Compare Package

A Python package for comparing legacy and modern JSON outputs from X3DNA.

## CLI Tool (`fp2-validate`)

The primary way to use this package is through the `fp2-validate` CLI:

```bash
# Install (from project root)
pip install -e .

# Basic validation
fp2-validate                      # All stages, all 3602 fast PDBs
fp2-validate --pdb 1EHZ -v        # Single PDB, verbose
fp2-validate --test-set 100       # 100-PDB test set
fp2-validate frames --max 50      # Frames only, first 50 PDBs

# Debugging
fp2-validate --stop-on-first      # Stop at first failure
fp2-validate --pdb 1EHZ -v -s     # Single PDB, verbose, stop on first

# CI/Scripts
fp2-validate --quiet              # Exit code only (0=pass, 1=fail)
fp2-validate --test-set 100 -q    # Quick CI check

# Generate difference reports
fp2-validate --diff               # Save to data/validation_results/
fp2-validate --diff-file my.json  # Custom output file

# Environment info
fp2-validate info                 # Shows executables, PDB counts
fp2-validate list-pdbs            # List all fast PDBs
fp2-validate list-pdbs --test-set 50  # List PDBs in test set
```

### Command Reference

| Command | Description |
|---------|-------------|
| `fp2-validate` | Validate all stages, all fast PDBs |
| `fp2-validate validate [STAGES]` | Validate specific stages |
| `fp2-validate atoms` | Validate atom records |
| `fp2-validate frames` | Validate frame calculations |
| `fp2-validate hbonds` | Validate H-bond lists |
| `fp2-validate pairs` | Validate base pairs |
| `fp2-validate steps` | Validate step parameters |
| `fp2-validate info` | Show environment info |
| `fp2-validate list-pdbs` | List available PDBs |

### Options Reference

| Option | Description |
|--------|-------------|
| `--pdb, -p TEXT` | Specific PDB(s) to validate |
| `--max, -n INT` | Maximum number of PDBs |
| `--test-set [10\|50\|100\|500\|1000]` | Use predefined test set |
| `--workers, -w INT` | Parallel workers (default: 10) |
| `--quiet, -q` | Exit code only |
| `--verbose, -v` | Show per-PDB results |
| `--stop-on-first, -s` | Stop at first failure |
| `--diff` | Document differences |
| `--diff-file PATH` | Custom output file |

---

## Python API

For programmatic use, you can import the package directly:

### Basic Comparison

```python
from pathlib import Path
from x3dna_json_compare import JsonComparator

# Default: compare all stages
comparator = JsonComparator()

# Or selectively enable/disable comparisons
comparator = JsonComparator(
    compare_atoms=True,
    compare_frames=True,
    compare_hbond_list=True,
    compare_pairs=True,
    compare_steps=True,
    enable_cache=False  # Disable caching
)

# Compare single PDB
legacy_file = Path("data/json_legacy/1H4S.json")
modern_file = Path("data/json/1H4S.json")
pdb_file = Path("data/pdb/1H4S.pdb")

result = comparator.compare_files(legacy_file, modern_file, pdb_file, "1H4S")

if result.status == 'error':
    print(f"Errors: {result.errors}")
elif result.has_differences():
    if result.frame_comparison:
        print(f"Frame differences: {len(result.frame_comparison.mismatched_calculations)}")
    if result.atom_comparison:
        print(f"Missing atoms: {len(result.atom_comparison.missing_in_modern)}")
else:
    print("✅ No differences found")
```

### Batch Validation

```python
from x3dna_json_compare.runner import ValidationRunner

# Create runner
runner = ValidationRunner(
    stages=['frames'],           # Or ['all'] for everything
    workers=10,
    verbose=True,
    stop_on_first=False
)

# Get PDB list
pdb_ids = runner.get_pdb_list(
    test_set=100,               # Use test set 100
    # Or: specific=['1EHZ', '1BNA'],
    # Or: max_count=50
)

# Run validation
summary = runner.validate(pdb_ids)

print(f"Passed: {summary.passed}/{summary.total}")
print(f"Failed: {summary.failed}")
if summary.all_passed:
    print("✅ All validations passed!")
```

### PDB File Reading

```python
from x3dna_json_compare import PdbFileReader

reader = PdbFileReader(Path("data/pdb/1H4S.pdb"))

# Get specific atom line
line_num, line = reader.find_atom_line("A", 1, " ", " C4 ")

# Get all atoms for a residue
atoms = reader.get_residue_atoms("A", 1, " ")
for atom_name, line_num, line in atoms:
    print(f"{atom_name}: {line}")
```

---

## Package Structure

```
x3dna_json_compare/
├── __init__.py          # Package exports
├── cli.py               # CLI entry point (fp2-validate)
├── runner.py            # Validation runner
├── output.py            # Output formatting
├── pdb_list.py          # PDB list management
├── executables.py       # Executable finding
├── json_comparison.py   # Main JsonComparator class
├── models.py            # Data models (ComparisonResult, etc.)
├── atom_comparison.py   # Atom comparison
├── frame_comparison.py  # Frame comparison
├── hbond_comparison.py  # H-bond comparison
├── pdb_utils.py         # PDB file utilities
├── result_cache.py      # Comparison caching
└── README.md            # This file
```

## Main Classes

- **`JsonComparator`**: Main comparison engine
- **`ValidationRunner`**: CLI validation runner  
- **`OutputFormatter`**: Output formatting
- **`PdbFileReader`**: Efficient PDB file reader

## Models

- **`ComparisonResult`**: Complete comparison result for a PDB
- **`ValidationSummary`**: Summary of batch validation
- **`FrameComparison`**: Frame calculation comparison
- **`AtomComparison`**: Atom-level comparison
- **`HbondComparison`**: H-bond comparison

---

## Migration from Old Scripts

| Old Script | New CLI Command |
|------------|-----------------|
| `validate_hbonds_batch.py` | `fp2-validate hbonds` |
| `validate_frames_batch.py` | `fp2-validate frames` |
| `test_single_pdb_ls_fitting.py 1EHZ` | `fp2-validate frames --pdb 1EHZ -v` |
| `find_first_mismatch.py` | `fp2-validate --stop-on-first` |
| `test_all_stages_batch.py` | `fp2-validate` |

See `scripts/README.md` for complete migration guide.

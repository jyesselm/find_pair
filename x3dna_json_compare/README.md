# X3DNA Comparison Library

A Python package for comparing legacy and modern JSON outputs from X3DNA.

## Installation

### Install as package (recommended)

```bash
# From project root
pip install -e .
```

This installs the package in editable mode, so changes are immediately available.

### Use without installation

The package can also be imported directly if the project root is in your Python path:

```python
from x3dna_json_compare import JsonComparator, ComparisonResult
```

## Usage

### Basic Comparison

```python
from pathlib import Path
from x3dna_json_compare import JsonComparator

# Default: compare both atoms and frames
comparator = JsonComparator()

# Or selectively enable/disable comparisons
comparator = JsonComparator(
    compare_atoms=True,   # Compare atom records (pdb_atoms)
    compare_frames=True    # Compare frame calculations (base_frame_calc, frame_calc, ls_fitting)
)

# Compare single PDB file
legacy_file = Path("data/json_legacy/1H4S.json")
modern_file = Path("data/json/1H4S.json")
pdb_file = Path("data/pdb/1H4S.pdb")

result = comparator.compare_files(legacy_file, modern_file, pdb_file, "1H4S")

if result.status == 'error':
    print(f"Errors: {result.errors}")
elif result.has_differences():
    if result.frame_comparison:
        print(f"Found {len(result.frame_comparison.mismatched_calculations)} frame differences")
    if result.atom_comparison:
        print(f"Found {len(result.atom_comparison.missing_in_modern)} missing atoms")
else:
    print("No differences found")
```

### Batch Comparison

```python
from x3dna_json_compare import JsonComparator

comparator = JsonComparator()
pdb_ids = ["1H4S", "1A34", "1BNA"]

results = comparator.batch_compare(
    pdb_ids,
    project_root=Path("."),
    max_workers=4
)

for pdb_id, result in results.items():
    print(f"{pdb_id}: {result.status}")
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

### Direct Comparison Functions

```python
from x3dna_json_compare import compare_atoms, compare_frames, PdbFileReader

# Load JSON data
legacy_calcs = legacy_json.get('calculations', [])
modern_calcs = modern_json.get('calculations', [])

# Extract records
legacy_frames = [r for r in legacy_calcs if r.get('type') == 'frame_calc']
modern_frames = [r for r in modern_calcs if r.get('type') == 'frame_calc']

# Compare
pdb_reader = PdbFileReader(Path("data/pdb/1H4S.pdb"))
frame_comparison = compare_frames(
    legacy_frames, modern_frames, 
    Path("data/pdb/1H4S.pdb"), 
    pdb_reader
)
```

## API Reference

### Main Classes

- **`JsonComparator`**: Main comparison engine with caching support
- **`PdbFileReader`**: Efficient PDB file reader with line indexing
- **`JsonValidator`**: JSON file validation utilities
- **`ComparisonCache`**: Result caching for performance

### Models

- **`ComparisonResult`**: Complete comparison result for a PDB
- **`FrameComparison`**: Frame calculation comparison results
- **`AtomComparison`**: Atom-level comparison results
- **`FrameMismatch`**: Details of a mismatched frame calculation
- **`AtomInfo`**: Information about a single atom

### Functions

- **`compare_frames()`**: Compare frame calculations directly
- **`compare_atoms()`**: Compare atom records directly
- **`get_pdb_line()`**: Convenience function to get a PDB line

## Examples

See the scripts in `scripts/` for complete examples:
- `compare_all_json_executed.py`: Compare all JSON files
- `compare_problematic_pdbs.py`: Detailed comparison with reports
- `analyze_ring_fitting_differences.py`: Ring fitting analysis


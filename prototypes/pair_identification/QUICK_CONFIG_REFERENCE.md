# Quick Config Reference

## TL;DR

```python
from pair_identification.config import get_config

# Get config (loads defaults + environment vars automatically)
config = get_config()

# Use it
pdb_files = list(config.pdb_dir.glob("*.pdb"))
```

## Common Patterns

### Pattern 1: Use Defaults

```python
from pair_identification.config import get_config

config = get_config()  # Done!
```

### Pattern 2: Load from File

```python
from pair_identification.config import load_config
from pathlib import Path

config = load_config(Path("my_config.json"))
```

### Pattern 3: Custom Config

```python
from pair_identification.config import Config

config = Config(
    default_workers=20,
    verbose=True,
    rmsd_threshold=2.0
)
```

### Pattern 4: CLI Integration

```python
from pair_identification.config import load_config, get_config
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--config', type=Path)
args = parser.parse_args()

config = load_config(args.config) if args.config else get_config()
```

### Pattern 5: Environment Overrides

```bash
export PAIR_ID_DEFAULT_WORKERS=20
export PAIR_ID_VERBOSE=true
python my_script.py  # Automatically uses env vars
```

## CLI Commands

```bash
# Show current config
python config.py show

# Create sample config
python config.py sample config.json
python config.py sample config.yaml

# Validate config
python config.py validate config.json
```

## Key Directories

| Field | Default | Exists? |
|-------|---------|---------|
| `pdb_dir` | `data/pdb` | Yes (4940 PDBs) |
| `dssr_dir` | `data/json_dssr` | Check |
| `slot_hbonds_dir` | `data/json/slot_hbonds` | Check |
| `idealized_templates_dir` | `basepair-idealized` | Check |
| `exemplar_templates_dir` | `basepair-exemplars` | Check |
| `output_dir` | `results` | Created auto |
| `viz_output_dir` | `viz_output` | Created auto |

## Key Thresholds

| Parameter | Default | Units |
|-----------|---------|-------|
| `max_hbond_distance` | 3.5 | Angstroms |
| `min_hbond_distance` | 2.4 | Angstroms |
| `rmsd_threshold` | 1.5 | Angstroms |
| `angle_threshold` | 15.0 | Degrees |
| `n1n9_range` | (8.0, 9.5) | Angstroms |

## Environment Variables

```bash
PAIR_ID_PDB_DIR               # Override pdb_dir
PAIR_ID_OUTPUT_DIR            # Override output_dir
PAIR_ID_DEFAULT_WORKERS       # Override default_workers
PAIR_ID_VERBOSE               # Override verbose (1/true/yes)
PAIR_ID_MAX_HBOND_DISTANCE    # Override max_hbond_distance
PAIR_ID_RMSD_THRESHOLD        # Override rmsd_threshold
# ... see CONFIG_README.md for full list
```

## Example Config File

```json
{
  "pdb_dir": "data/pdb",
  "output_dir": "results",
  "default_workers": 10,
  "verbose": false,
  "max_hbond_distance": 3.5,
  "min_hbond_distance": 2.4,
  "rmsd_threshold": 1.5,
  "angle_threshold": 15.0,
  "n1n9_range": [8.0, 9.5]
}
```

## Full Documentation

See `CONFIG_README.md` for complete documentation.

## Examples

See `example_config_usage.py` for 7 working examples.

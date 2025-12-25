# Configuration System

Central configuration module for the pair_identification package.

## Features

- Dataclass-based configuration with type hints
- JSON and YAML file support
- Environment variable overrides
- Global singleton pattern
- Automatic validation and directory creation
- CLI utilities for config management

## Quick Start

### Using Default Configuration

```python
from pair_identification.config import get_config

# Get default config (loads automatically)
config = get_config()

# Access configuration
pdb_dir = config.pdb_dir
workers = config.default_workers
```

### Loading from File

```python
from pair_identification.config import load_config, set_config
from pathlib import Path

# Load from JSON or YAML
config = load_config(Path("my_config.json"))

# Set as global config
set_config(config)
```

### Creating a Config File

```bash
# Create sample JSON config
python config.py sample config.json

# Create sample YAML config
python config.py sample config.yaml

# Validate existing config
python config.py validate config.json

# Show current config
python config.py show
```

## Configuration Options

### Input Directories

| Field | Default | Description |
|-------|---------|-------------|
| `pdb_dir` | `data/pdb` | Directory containing PDB files |
| `dssr_dir` | `data/json_dssr` | Directory containing DSSR JSON output |
| `slot_hbonds_dir` | `data/json/slot_hbonds` | Directory containing slot-based H-bond JSON |
| `legacy_json_dir` | `data/json_legacy` | Directory containing legacy X3DNA JSON output |
| `modern_json_dir` | `data/json` | Directory containing modern X3DNA JSON output |

### Template Directories

| Field | Default | Description |
|-------|---------|-------------|
| `idealized_templates_dir` | `basepair-idealized` | Idealized base pair templates |
| `exemplar_templates_dir` | `basepair-exemplars` | Exemplar base pair templates |

### Output Directories

| Field | Default | Description |
|-------|---------|-------------|
| `output_dir` | `results` | Root directory for analysis output |
| `viz_output_dir` | `viz_output` | Directory for visualization output |

### Analysis Thresholds

| Field | Default | Description |
|-------|---------|-------------|
| `max_hbond_distance` | `3.5` | Maximum H-bond distance in Angstroms |
| `min_hbond_distance` | `2.4` | Minimum H-bond distance in Angstroms |
| `rmsd_threshold` | `1.5` | RMSD threshold for template matching |
| `angle_threshold` | `15.0` | Angle threshold in degrees |
| `n1n9_range` | `(8.0, 9.5)` | Valid N1-N9 distance range in Angstroms |

### Processing Options

| Field | Default | Description |
|-------|---------|-------------|
| `default_workers` | `10` | Default number of parallel workers |
| `verbose` | `false` | Enable verbose output |

## Environment Variable Overrides

All configuration options can be overridden via environment variables:

| Environment Variable | Config Field |
|---------------------|--------------|
| `PAIR_ID_PDB_DIR` | `pdb_dir` |
| `PAIR_ID_DSSR_DIR` | `dssr_dir` |
| `PAIR_ID_SLOT_HBONDS_DIR` | `slot_hbonds_dir` |
| `PAIR_ID_LEGACY_JSON_DIR` | `legacy_json_dir` |
| `PAIR_ID_MODERN_JSON_DIR` | `modern_json_dir` |
| `PAIR_ID_IDEALIZED_TEMPLATES_DIR` | `idealized_templates_dir` |
| `PAIR_ID_EXEMPLAR_TEMPLATES_DIR` | `exemplar_templates_dir` |
| `PAIR_ID_OUTPUT_DIR` | `output_dir` |
| `PAIR_ID_VIZ_OUTPUT_DIR` | `viz_output_dir` |
| `PAIR_ID_MAX_HBOND_DISTANCE` | `max_hbond_distance` |
| `PAIR_ID_MIN_HBOND_DISTANCE` | `min_hbond_distance` |
| `PAIR_ID_RMSD_THRESHOLD` | `rmsd_threshold` |
| `PAIR_ID_ANGLE_THRESHOLD` | `angle_threshold` |
| `PAIR_ID_DEFAULT_WORKERS` | `default_workers` |
| `PAIR_ID_VERBOSE` | `verbose` |

### Example

```bash
# Override output directory
export PAIR_ID_OUTPUT_DIR=/tmp/results

# Override number of workers
export PAIR_ID_DEFAULT_WORKERS=20

# Enable verbose mode
export PAIR_ID_VERBOSE=true

# Run your script
python my_analysis.py
```

## Programmatic Usage

### Creating a Custom Config

```python
from pair_identification.config import Config, set_config
from pathlib import Path

# Create custom config
config = Config(
    pdb_dir=Path("/path/to/pdbs"),
    default_workers=20,
    verbose=True,
    rmsd_threshold=2.0
)

# Validate it
config.validate()

# Set as global
set_config(config)
```

### Saving Config to File

```python
from pair_identification.config import Config
from pathlib import Path

config = Config(default_workers=20)

# Save as JSON
config.save(Path("my_config.json"))

# Save as YAML
config.save(Path("my_config.yaml"))
```

### Loading and Modifying Config

```python
from pair_identification.config import load_config
from pathlib import Path

# Load existing config
config = load_config(Path("config.json"))

# Modify
config.default_workers = 20
config.verbose = True

# Validate changes
config.validate()

# Save modified config
config.save(Path("config_modified.json"))
```

### Disabling Environment Overrides

```python
from pair_identification.config import load_config
from pathlib import Path

# Load without applying environment variables
config = load_config(
    Path("config.json"),
    apply_env_overrides=False
)
```

## Validation

The configuration system automatically validates:

1. **Input directories**: Warns if missing (doesn't error)
2. **Output directories**: Creates if missing
3. **Numeric thresholds**: Validates ranges and relationships
   - `max_hbond_distance > min_hbond_distance`
   - `rmsd_threshold > 0`
   - `angle_threshold` in `(0, 180)`
   - `n1n9_range[0] < n1n9_range[1]`
   - `default_workers >= 1`

### Example Validation

```python
from pair_identification.config import Config

# This will raise ValueError
config = Config(
    max_hbond_distance=2.0,
    min_hbond_distance=3.0  # ERROR: max < min
)
config.validate()  # Raises ValueError
```

## CLI Utilities

The config module includes a command-line interface:

```bash
# Show current configuration
python config.py show

# Create sample config file
python config.py sample my_config.json
python config.py sample my_config.yaml

# Validate configuration file
python config.py validate my_config.json
```

## Integration Example

Here's how to use the config system in your analysis scripts:

```python
#!/usr/bin/env python3
"""Example analysis script using config system."""

from pair_identification.config import get_config, load_config, print_config
from pathlib import Path
import sys

def main():
    # Option 1: Use global config (respects env vars)
    config = get_config()

    # Option 2: Load specific config file
    if len(sys.argv) > 1:
        config_path = Path(sys.argv[1])
        config = load_config(config_path)

    # Print loaded config
    print_config(config)

    # Use config in your analysis
    pdb_files = list(config.pdb_dir.glob("*.pdb"))
    print(f"\nFound {len(pdb_files)} PDB files in {config.pdb_dir}")

    # Access thresholds
    print(f"Using RMSD threshold: {config.rmsd_threshold} Å")
    print(f"Using {config.default_workers} workers")

    # ... your analysis code here ...

if __name__ == '__main__':
    main()
```

Run with:
```bash
# Use default/environment config
python my_analysis.py

# Use specific config file
python my_analysis.py my_config.json
```

## Sample Configuration Files

Two sample configuration files are provided:

1. `config.sample.json` - JSON format
2. `config.sample.yaml` - YAML format

Both contain the same default values and can be used as templates for creating your own configuration files.

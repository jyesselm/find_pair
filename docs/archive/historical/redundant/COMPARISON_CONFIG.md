# Comparison Configuration Guide

The JSON comparison tool now supports YAML-based configuration to control which comparisons are performed.

## Configuration File

The default configuration file is `comparison_config.yaml` in the project root. You can specify a different file using the `--config` option.

## Configuration Options

### Comparison Types

Enable or disable specific comparison types:

```yaml
comparisons:
  atoms: false              # Compare pdb_atoms records (atom-by-atom comparison)
  frames: true              # Compare frame calculations (base_frame_calc, frame_calc, ls_fitting)
  steps: true               # Compare step parameters (bpstep_params, helical_params)
  pairs: true               # Compare pair records (pair_validation, distance_checks, base_pair, find_bestpair_selection)
  hbond_list: true          # Compare hbond_list records (hydrogen bond details for base pairs)
  residue_indices: true     # Compare residue_indices records (residue-to-atom mappings)
```

### Tolerance

Set the floating point comparison tolerance:

```yaml
tolerance: 1.0e-6
```

### Cache Settings

Control result caching:

```yaml
cache:
  enabled: true             # Enable result caching
  force_recompute: false    # Force recomputation even if cache exists
```

## Default Configuration

The default configuration (when no file is found) has:
- `atoms: false` (disabled by default - usually not needed for regular comparisons)
- `frames: true`
- `steps: true`
- `pairs: true`
- `hbond_list: true`
- `residue_indices: true`

## Usage

### Using Default Config

```bash
# Uses comparison_config.yaml from project root
python3 scripts/compare_json.py compare --test-set 10
```

### Using Custom Config

```bash
# Specify a different config file
python3 scripts/compare_json.py compare --config my_config.yaml --test-set 10
```

### Command-Line Overrides

The `--no-cache` flag still works to disable caching regardless of config:

```bash
python3 scripts/compare_json.py compare --no-cache --test-set 10
```

## Example: Disable Atom Comparison

Since atom comparisons are usually not needed (they're typically already verified), you can disable them:

```yaml
comparisons:
  atoms: false
  frames: true
  steps: true
  pairs: true
  hbond_list: true
  residue_indices: true

tolerance: 1.0e-6

cache:
  enabled: true
  force_recompute: false
```

## What Each Comparison Checks

- **atoms**: Individual atom-by-atom comparison (coordinates, names, metadata)
- **frames**: Base frame calculations, least squares fitting, frame calculations
- **steps**: Base pair step parameters and helical parameters
- **pairs**: Pair validation results, distance checks, base pairs, and find_bestpair selections
- **hbond_list**: Detailed hydrogen bond information for each base pair
- **residue_indices**: Residue-to-atom index mappings (seidx)

## Notes

- The configuration file is automatically searched for in the project root
- If no config file is found, sensible defaults are used
- Command-line flags (like `--no-cache`) can override config settings
- Configuration is loaded per process, so each worker in parallel execution will load the same config


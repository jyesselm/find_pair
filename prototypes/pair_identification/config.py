"""Central configuration module for pair_identification package.

Provides dataclass-based configuration with support for:
- JSON/YAML file loading
- Environment variable overrides
- Global singleton pattern
- Automatic validation and directory creation
"""

import json
import os
import warnings
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Tuple, Dict, Any

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


@dataclass
class Config:
    """Configuration for pair identification analysis.

    Attributes:
        pdb_dir: Directory containing PDB files.
        dssr_dir: Directory containing DSSR JSON output.
        slot_hbonds_dir: Directory containing slot-based H-bond JSON.
        legacy_json_dir: Directory containing legacy X3DNA JSON output.
        modern_json_dir: Directory containing modern X3DNA JSON output.
        idealized_templates_dir: Directory containing idealized base pair templates.
        exemplar_templates_dir: Directory containing exemplar base pair templates.
        output_dir: Root directory for analysis output.
        viz_output_dir: Directory for visualization output (PyMOL scripts, etc).
        max_hbond_distance: Maximum H-bond distance in Angstroms.
        min_hbond_distance: Minimum H-bond distance in Angstroms.
        rmsd_threshold: RMSD threshold for template matching.
        angle_threshold: Angle threshold in degrees for geometric validation.
        n1n9_range: Valid range for N1-N9 distance in Angstroms.
        default_workers: Default number of parallel workers.
        verbose: Enable verbose output.
    """

    # Input data directories
    pdb_dir: Path = Path("data/pdb")
    dssr_dir: Path = Path("data/json_dssr")
    slot_hbonds_dir: Path = Path("data/json/slot_hbonds")
    legacy_json_dir: Path = Path("data/json_legacy")
    modern_json_dir: Path = Path("data/json")

    # Template directories
    idealized_templates_dir: Path = Path("basepair-idealized")
    exemplar_templates_dir: Path = Path("basepair-exemplars")

    # Output directories
    output_dir: Path = Path("results")
    viz_output_dir: Path = Path("viz_output")

    # Analysis thresholds
    max_hbond_distance: float = 3.5
    min_hbond_distance: float = 2.4
    rmsd_threshold: float = 1.5
    angle_threshold: float = 15.0
    n1n9_range: Tuple[float, float] = (8.0, 9.5)

    # Processing options
    default_workers: int = 10
    verbose: bool = False

    def __post_init__(self):
        """Convert string paths to Path objects."""
        # Convert all path fields to Path objects
        path_fields = [
            'pdb_dir', 'dssr_dir', 'slot_hbonds_dir',
            'legacy_json_dir', 'modern_json_dir',
            'idealized_templates_dir', 'exemplar_templates_dir',
            'output_dir', 'viz_output_dir'
        ]

        for field_name in path_fields:
            value = getattr(self, field_name)
            if isinstance(value, str):
                setattr(self, field_name, Path(value))

    def validate(self, create_output_dirs: bool = True) -> None:
        """Validate configuration.

        Args:
            create_output_dirs: If True, create output directories if missing.

        Warns:
            UserWarning: If required input directories don't exist.
        """
        # Check input directories (warn if missing)
        input_dirs = [
            ('pdb_dir', self.pdb_dir),
            ('dssr_dir', self.dssr_dir),
            ('slot_hbonds_dir', self.slot_hbonds_dir),
            ('legacy_json_dir', self.legacy_json_dir),
            ('modern_json_dir', self.modern_json_dir),
            ('idealized_templates_dir', self.idealized_templates_dir),
            ('exemplar_templates_dir', self.exemplar_templates_dir),
        ]

        for name, path in input_dirs:
            if not path.exists():
                warnings.warn(
                    f"Input directory '{name}' does not exist: {path}",
                    UserWarning
                )

        # Create output directories if requested
        if create_output_dirs:
            self.output_dir.mkdir(parents=True, exist_ok=True)
            self.viz_output_dir.mkdir(parents=True, exist_ok=True)

        # Validate numeric thresholds
        if self.max_hbond_distance <= self.min_hbond_distance:
            raise ValueError(
                f"max_hbond_distance ({self.max_hbond_distance}) must be "
                f"greater than min_hbond_distance ({self.min_hbond_distance})"
            )

        if self.rmsd_threshold <= 0:
            raise ValueError(f"rmsd_threshold must be positive: {self.rmsd_threshold}")

        if self.angle_threshold <= 0 or self.angle_threshold >= 180:
            raise ValueError(
                f"angle_threshold must be in (0, 180): {self.angle_threshold}"
            )

        n1n9_min, n1n9_max = self.n1n9_range
        if n1n9_min >= n1n9_max:
            raise ValueError(
                f"n1n9_range min ({n1n9_min}) must be less than max ({n1n9_max})"
            )

        if self.default_workers < 1:
            raise ValueError(
                f"default_workers must be at least 1: {self.default_workers}"
            )

    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary with paths as strings.

        Returns:
            Dictionary representation suitable for JSON serialization.
        """
        data = asdict(self)

        # Convert Path objects to strings
        path_fields = [
            'pdb_dir', 'dssr_dir', 'slot_hbonds_dir',
            'legacy_json_dir', 'modern_json_dir',
            'idealized_templates_dir', 'exemplar_templates_dir',
            'output_dir', 'viz_output_dir'
        ]

        for field_name in path_fields:
            data[field_name] = str(data[field_name])

        return data

    def save(self, path: Path, format: str = 'auto') -> None:
        """Save configuration to file.

        Args:
            path: Path to save configuration file.
            format: Format to use ('json', 'yaml', or 'auto' to infer from extension).

        Raises:
            ValueError: If format is yaml but PyYAML is not installed.
        """
        if format == 'auto':
            format = 'yaml' if path.suffix in ['.yml', '.yaml'] else 'json'

        data = self.to_dict()

        if format == 'json':
            with open(path, 'w') as f:
                json.dump(data, f, indent=2)
        elif format == 'yaml':
            if not HAS_YAML:
                raise ValueError("PyYAML not installed, cannot save as YAML")
            with open(path, 'w') as f:
                yaml.dump(data, f, default_flow_style=False)
        else:
            raise ValueError(f"Unknown format: {format}")


def load_config(
    config_path: Optional[Path] = None,
    apply_env_overrides: bool = True
) -> Config:
    """Load configuration from file with optional environment overrides.

    Args:
        config_path: Path to config file (JSON or YAML). If None, use defaults.
        apply_env_overrides: If True, apply environment variable overrides.

    Returns:
        Loaded and validated configuration.

    Raises:
        FileNotFoundError: If config_path specified but doesn't exist.
        ValueError: If config file format is invalid.
    """
    if config_path is None:
        config = Config()
    else:
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")

        # Determine format from extension
        if config_path.suffix in ['.yml', '.yaml']:
            if not HAS_YAML:
                raise ValueError(
                    "PyYAML not installed, cannot load YAML config. "
                    "Install with: pip install pyyaml"
                )
            with open(config_path) as f:
                data = yaml.safe_load(f)
        elif config_path.suffix == '.json':
            with open(config_path) as f:
                data = json.load(f)
        else:
            raise ValueError(
                f"Unknown config file format: {config_path.suffix}. "
                "Expected .json, .yml, or .yaml"
            )

        # Create config from loaded data
        config = Config(**data)

    # Apply environment variable overrides
    if apply_env_overrides:
        config = _apply_env_overrides(config)

    # Validate configuration
    config.validate()

    return config


def _apply_env_overrides(config: Config) -> Config:
    """Apply environment variable overrides to config.

    Environment variables:
        PAIR_ID_PDB_DIR: Override pdb_dir
        PAIR_ID_DSSR_DIR: Override dssr_dir
        PAIR_ID_SLOT_HBONDS_DIR: Override slot_hbonds_dir
        PAIR_ID_LEGACY_JSON_DIR: Override legacy_json_dir
        PAIR_ID_MODERN_JSON_DIR: Override modern_json_dir
        PAIR_ID_IDEALIZED_TEMPLATES_DIR: Override idealized_templates_dir
        PAIR_ID_EXEMPLAR_TEMPLATES_DIR: Override exemplar_templates_dir
        PAIR_ID_OUTPUT_DIR: Override output_dir
        PAIR_ID_VIZ_OUTPUT_DIR: Override viz_output_dir
        PAIR_ID_MAX_HBOND_DISTANCE: Override max_hbond_distance
        PAIR_ID_MIN_HBOND_DISTANCE: Override min_hbond_distance
        PAIR_ID_RMSD_THRESHOLD: Override rmsd_threshold
        PAIR_ID_ANGLE_THRESHOLD: Override angle_threshold
        PAIR_ID_DEFAULT_WORKERS: Override default_workers
        PAIR_ID_VERBOSE: Override verbose (1/true/yes for True)

    Args:
        config: Configuration to override.

    Returns:
        Configuration with environment overrides applied.
    """
    # Path overrides
    path_overrides = {
        'PAIR_ID_PDB_DIR': 'pdb_dir',
        'PAIR_ID_DSSR_DIR': 'dssr_dir',
        'PAIR_ID_SLOT_HBONDS_DIR': 'slot_hbonds_dir',
        'PAIR_ID_LEGACY_JSON_DIR': 'legacy_json_dir',
        'PAIR_ID_MODERN_JSON_DIR': 'modern_json_dir',
        'PAIR_ID_IDEALIZED_TEMPLATES_DIR': 'idealized_templates_dir',
        'PAIR_ID_EXEMPLAR_TEMPLATES_DIR': 'exemplar_templates_dir',
        'PAIR_ID_OUTPUT_DIR': 'output_dir',
        'PAIR_ID_VIZ_OUTPUT_DIR': 'viz_output_dir',
    }

    for env_var, field_name in path_overrides.items():
        if env_var in os.environ:
            setattr(config, field_name, Path(os.environ[env_var]))

    # Float overrides
    float_overrides = {
        'PAIR_ID_MAX_HBOND_DISTANCE': 'max_hbond_distance',
        'PAIR_ID_MIN_HBOND_DISTANCE': 'min_hbond_distance',
        'PAIR_ID_RMSD_THRESHOLD': 'rmsd_threshold',
        'PAIR_ID_ANGLE_THRESHOLD': 'angle_threshold',
    }

    for env_var, field_name in float_overrides.items():
        if env_var in os.environ:
            setattr(config, field_name, float(os.environ[env_var]))

    # Integer overrides
    if 'PAIR_ID_DEFAULT_WORKERS' in os.environ:
        config.default_workers = int(os.environ['PAIR_ID_DEFAULT_WORKERS'])

    # Boolean overrides
    if 'PAIR_ID_VERBOSE' in os.environ:
        value = os.environ['PAIR_ID_VERBOSE'].lower()
        config.verbose = value in ['1', 'true', 'yes', 'on']

    return config


# Global config singleton
_config: Optional[Config] = None


def get_config() -> Config:
    """Get global configuration instance.

    Returns:
        Global configuration, loading defaults if not yet set.
    """
    global _config
    if _config is None:
        _config = load_config()
    return _config


def set_config(config: Config) -> None:
    """Set global configuration instance.

    Args:
        config: Configuration to set as global.
    """
    global _config
    config.validate()
    _config = config


def reset_config() -> None:
    """Reset global configuration to None.

    Next call to get_config() will load defaults.
    """
    global _config
    _config = None


def print_config(config: Optional[Config] = None) -> None:
    """Print configuration to stdout in readable format.

    Args:
        config: Configuration to print. If None, use global config.
    """
    if config is None:
        config = get_config()

    print("=== Pair Identification Configuration ===\n")

    print("Input Directories:")
    print(f"  pdb_dir:                   {config.pdb_dir}")
    print(f"  dssr_dir:                  {config.dssr_dir}")
    print(f"  slot_hbonds_dir:           {config.slot_hbonds_dir}")
    print(f"  legacy_json_dir:           {config.legacy_json_dir}")
    print(f"  modern_json_dir:           {config.modern_json_dir}")

    print("\nTemplate Directories:")
    print(f"  idealized_templates_dir:   {config.idealized_templates_dir}")
    print(f"  exemplar_templates_dir:    {config.exemplar_templates_dir}")

    print("\nOutput Directories:")
    print(f"  output_dir:                {config.output_dir}")
    print(f"  viz_output_dir:            {config.viz_output_dir}")

    print("\nAnalysis Thresholds:")
    print(f"  max_hbond_distance:        {config.max_hbond_distance} Å")
    print(f"  min_hbond_distance:        {config.min_hbond_distance} Å")
    print(f"  rmsd_threshold:            {config.rmsd_threshold} Å")
    print(f"  angle_threshold:           {config.angle_threshold}°")
    print(f"  n1n9_range:                {config.n1n9_range[0]}-{config.n1n9_range[1]} Å")

    print("\nProcessing Options:")
    print(f"  default_workers:           {config.default_workers}")
    print(f"  verbose:                   {config.verbose}")

    print("\n==========================================")


def create_sample_config(path: Path, format: str = 'json') -> None:
    """Create sample configuration file.

    Args:
        path: Path to save sample config.
        format: Format to use ('json' or 'yaml').
    """
    config = Config()
    config.save(path, format=format)
    print(f"Sample configuration saved to: {path}")


if __name__ == '__main__':
    # CLI for config utilities
    import sys

    if len(sys.argv) < 2:
        print("Usage:")
        print("  python config.py show              - Show current config")
        print("  python config.py sample <path>     - Create sample config file")
        print("  python config.py validate <path>   - Validate config file")
        sys.exit(1)

    command = sys.argv[1]

    if command == 'show':
        print_config()

    elif command == 'sample':
        if len(sys.argv) < 3:
            print("Error: path required")
            print("Usage: python config.py sample <path>")
            sys.exit(1)

        path = Path(sys.argv[2])
        format_str = 'yaml' if path.suffix in ['.yml', '.yaml'] else 'json'
        create_sample_config(path, format=format_str)

    elif command == 'validate':
        if len(sys.argv) < 3:
            print("Error: path required")
            print("Usage: python config.py validate <path>")
            sys.exit(1)

        path = Path(sys.argv[2])
        try:
            config = load_config(path)
            print(f"Configuration is valid: {path}")
            print("\nLoaded configuration:")
            print_config(config)
        except Exception as e:
            print(f"Configuration is invalid: {e}")
            sys.exit(1)

    else:
        print(f"Unknown command: {command}")
        sys.exit(1)

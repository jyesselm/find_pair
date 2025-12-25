#!/usr/bin/env python3
"""Example script demonstrating config system usage in analysis workflows.

This shows best practices for:
1. Loading configuration
2. Using config in analysis
3. Command-line integration
4. Creating custom configs programmatically
"""

import sys
from pathlib import Path
from typing import Optional

from config import (
    Config,
    get_config,
    load_config,
    set_config,
    print_config
)


def example_1_simple_usage():
    """Example 1: Simplest usage - use global config."""
    print("=== Example 1: Simple Global Config ===\n")

    # Just get the config - it loads defaults automatically
    config = get_config()

    # Use it
    print(f"PDB directory: {config.pdb_dir}")
    print(f"Workers: {config.default_workers}")
    print(f"RMSD threshold: {config.rmsd_threshold} Å\n")


def example_2_load_from_file():
    """Example 2: Load config from file."""
    print("=== Example 2: Load From File ===\n")

    config_path = Path("config.sample.json")
    if not config_path.exists():
        print(f"Config file not found: {config_path}")
        return

    # Load from file
    config = load_config(config_path)

    # Set as global (optional)
    set_config(config)

    print(f"Loaded config from: {config_path}")
    print(f"Workers: {config.default_workers}\n")


def example_3_custom_config():
    """Example 3: Create custom config programmatically."""
    print("=== Example 3: Custom Config ===\n")

    # Create custom config for specific analysis
    config = Config(
        pdb_dir=Path("/custom/path/to/pdbs"),
        default_workers=20,
        verbose=True,
        rmsd_threshold=2.0,
        max_hbond_distance=4.0
    )

    # Validate before use
    config.validate(create_output_dirs=False)

    print("Created custom config:")
    print(f"  PDB dir: {config.pdb_dir}")
    print(f"  Workers: {config.default_workers}")
    print(f"  Verbose: {config.verbose}")
    print(f"  RMSD threshold: {config.rmsd_threshold}\n")


def example_4_environment_vars():
    """Example 4: Environment variable overrides."""
    print("=== Example 4: Environment Variables ===\n")

    import os

    # Set environment variables
    os.environ['PAIR_ID_DEFAULT_WORKERS'] = '15'
    os.environ['PAIR_ID_VERBOSE'] = 'true'

    # Load config - will apply env overrides
    config = load_config()

    print("Config with environment overrides:")
    print(f"  Workers: {config.default_workers} (from PAIR_ID_DEFAULT_WORKERS)")
    print(f"  Verbose: {config.verbose} (from PAIR_ID_VERBOSE)\n")

    # Clean up
    del os.environ['PAIR_ID_DEFAULT_WORKERS']
    del os.environ['PAIR_ID_VERBOSE']


def example_5_analysis_workflow(config: Optional[Config] = None):
    """Example 5: Realistic analysis workflow."""
    print("=== Example 5: Analysis Workflow ===\n")

    # Use provided config or load default
    if config is None:
        config = get_config()

    # Check input data
    pdb_files = list(config.pdb_dir.glob("*.pdb")) if config.pdb_dir.exists() else []
    print(f"Found {len(pdb_files)} PDB files in {config.pdb_dir}")

    # Use config for analysis parameters
    print(f"\nAnalysis parameters:")
    print(f"  Max H-bond distance: {config.max_hbond_distance} Å")
    print(f"  Min H-bond distance: {config.min_hbond_distance} Å")
    print(f"  RMSD threshold: {config.rmsd_threshold} Å")
    print(f"  N1-N9 range: {config.n1n9_range[0]}-{config.n1n9_range[1]} Å")

    # Use config for processing
    print(f"\nProcessing with {config.default_workers} workers")

    # Output setup
    config.output_dir.mkdir(parents=True, exist_ok=True)
    output_file = config.output_dir / "results.json"
    print(f"Output will be saved to: {output_file}\n")


def example_6_cli_integration():
    """Example 6: CLI integration with optional config file."""
    print("=== Example 6: CLI Integration ===\n")

    import argparse

    parser = argparse.ArgumentParser(
        description="Example analysis with config support"
    )
    parser.add_argument(
        '--config',
        type=Path,
        help='Path to configuration file (JSON or YAML)'
    )
    parser.add_argument(
        '--workers',
        type=int,
        help='Number of workers (overrides config)'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose output (overrides config)'
    )

    # For demo, use empty args
    args = parser.parse_args([])

    # Load config
    if args.config:
        config = load_config(args.config)
        print(f"Loaded config from: {args.config}")
    else:
        config = get_config()
        print("Using default/environment config")

    # Apply command-line overrides
    if args.workers:
        config.default_workers = args.workers
        print(f"Overriding workers: {args.workers}")

    if args.verbose:
        config.verbose = True
        print("Enabling verbose mode")

    print(f"\nFinal configuration:")
    print(f"  Workers: {config.default_workers}")
    print(f"  Verbose: {config.verbose}\n")

    return config


def example_7_save_config():
    """Example 7: Save modified config."""
    print("=== Example 7: Save Modified Config ===\n")

    # Load or create config
    config = Config(default_workers=25, verbose=True)

    # Save to file
    output_path = Path("my_custom_config.json")
    config.save(output_path)
    print(f"Saved custom config to: {output_path}")

    # Load it back
    loaded = load_config(output_path, apply_env_overrides=False)
    print(f"Loaded back - workers: {loaded.default_workers}\n")

    # Clean up
    output_path.unlink(missing_ok=True)


def main():
    """Run all examples."""
    print("=" * 60)
    print("Config System Usage Examples")
    print("=" * 60 + "\n")

    # Run examples
    example_1_simple_usage()
    example_2_load_from_file()
    example_3_custom_config()
    example_4_environment_vars()
    example_5_analysis_workflow()
    example_6_cli_integration()
    example_7_save_config()

    print("=" * 60)
    print("All examples completed!")
    print("=" * 60)


if __name__ == '__main__':
    main()

"""
Configuration loader for JSON comparison.

Supports YAML-based configuration to control which comparisons are performed.
"""

import yaml
from pathlib import Path
from typing import Dict, Any, Optional


DEFAULT_CONFIG = {
    "comparisons": {
        "atoms": False,  # pdb_atoms records
        "frames": True,  # base_frame_calc, frame_calc, ls_fitting
        "steps": True,   # bpstep_params, helical_params
        "pairs": True,   # pair_validation, distance_checks, base_pair, find_bestpair_selection
        "hbond_list": True,  # hbond_list records
        "residue_indices": True,  # residue_indices records
    },
    "tolerance": 2e-5,  # 2e-5 handles normal FP variations in trigonometric calculations
}


def load_config(config_path: Optional[Path] = None) -> Dict[str, Any]:
    """
    Load comparison configuration from YAML file.
    
    Args:
        config_path: Path to YAML config file. If None, looks for comparison_config.yaml
                    in project root, or returns default config.
    
    Returns:
        Configuration dictionary
    """
    # Try to find config file
    if config_path is None:
        # Look for comparison_config.yaml in current directory and parent directories
        current_dir = Path.cwd()
        for directory in [current_dir, current_dir.parent, current_dir.parent.parent]:
            candidate = directory / "comparison_config.yaml"
            if candidate.exists():
                config_path = candidate
                break
    
    if config_path is None or not config_path.exists():
        # Return default config
        return DEFAULT_CONFIG.copy()
    
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        # Merge with defaults to ensure all keys exist
        merged_config = DEFAULT_CONFIG.copy()
        if config:
            if "comparisons" in config:
                merged_config["comparisons"].update(config["comparisons"])
            if "tolerance" in config:
                merged_config["tolerance"] = config["tolerance"]
        
        return merged_config
    except Exception as e:
        print(f"Warning: Could not load config from {config_path}: {e}")
        print("Using default configuration.")
        return DEFAULT_CONFIG.copy()


def get_comparison_flags(config: Dict[str, Any]) -> Dict[str, bool]:
    """
    Extract comparison flags from config.
    
    Args:
        config: Configuration dictionary
    
    Returns:
        Dictionary with comparison flags (compare_atoms, compare_frames, etc.)
    """
    comparisons = config.get("comparisons", {})
    return {
        "compare_atoms": comparisons.get("atoms", False),
        "compare_frames": comparisons.get("frames", True),
        "compare_steps": comparisons.get("steps", True),
        "compare_pairs": comparisons.get("pairs", True),
        "compare_hbond_list": comparisons.get("hbond_list", True),
        "compare_residue_indices": comparisons.get("residue_indices", True),
    }


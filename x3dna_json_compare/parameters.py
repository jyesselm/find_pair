"""
X3DNA Parameters - Single source of truth for all algorithm constants.

This module loads parameters from resources/config/parameters.json and provides
easy access from Python. The same JSON file is used by C++ code.

Usage:
    from x3dna_json_compare.parameters import params

    # Access nested values
    print(params.validation.distance.max_dorg)  # 15.0
    print(params.hydrogen_bond.thresholds.good_min)  # 2.5

    # Or as dict
    print(params['validation']['distance']['max_dorg'])
"""

import json
from pathlib import Path
from typing import Any, Dict


class ParameterNamespace:
    """Allows dot-notation access to nested dict values."""

    def __init__(self, data: Dict[str, Any]):
        for key, value in data.items():
            if key.startswith('_'):  # Skip description fields
                continue
            if isinstance(value, dict):
                setattr(self, key, ParameterNamespace(value))
            else:
                setattr(self, key, value)

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def to_dict(self) -> Dict[str, Any]:
        """Convert back to dict."""
        result = {}
        for key, value in self.__dict__.items():
            if isinstance(value, ParameterNamespace):
                result[key] = value.to_dict()
            else:
                result[key] = value
        return result


class Parameters:
    """
    Singleton class for accessing X3DNA parameters.

    Parameters are loaded from resources/config/parameters.json on first access.
    """

    _instance = None
    _params = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._load()
        return cls._instance

    def _load(self) -> None:
        """Load parameters from JSON file."""
        # Find the parameters.json file
        module_dir = Path(__file__).parent.parent
        config_paths = [
            module_dir / "resources" / "config" / "parameters.json",
            module_dir.parent / "resources" / "config" / "parameters.json",
            Path(__file__).parent / "resources" / "config" / "parameters.json",
        ]

        config_path = None
        for path in config_paths:
            if path.exists():
                config_path = path
                break

        if config_path is None:
            raise FileNotFoundError(
                f"Could not find parameters.json. Searched: {config_paths}"
            )

        with open(config_path) as f:
            data = json.load(f)

        self._params = ParameterNamespace(data)
        self._raw = data

    def __getattr__(self, name: str) -> Any:
        if name.startswith('_'):
            return super().__getattribute__(name)
        return getattr(self._params, name)

    def __getitem__(self, key: str) -> Any:
        return self._raw[key]

    @property
    def raw(self) -> Dict[str, Any]:
        """Get raw dict of all parameters."""
        return self._raw

    def reload(self) -> None:
        """Reload parameters from file."""
        self._load()


# Singleton instance for easy import
params = Parameters()


# Convenience exports for common parameters
def get_validation_params() -> Dict[str, Any]:
    """Get validation parameters as dict."""
    return params.raw['validation']


def get_hbond_params() -> Dict[str, Any]:
    """Get hydrogen bond parameters as dict."""
    return params.raw['hydrogen_bond']


def get_quality_params() -> Dict[str, Any]:
    """Get quality score parameters as dict."""
    return params.raw['quality_score']


if __name__ == "__main__":
    # Quick test
    print("X3DNA Parameters")
    print("=" * 40)
    print(f"max_dorg: {params.validation.distance.max_dorg}")
    print(f"max_dv: {params.validation.distance.max_dv}")
    print(f"hb_good_min: {params.hydrogen_bond.thresholds.good_min}")
    print(f"hb_good_max: {params.hydrogen_bond.thresholds.good_max}")
    print(f"rmsd_cutoff: {params.nucleotide.rmsd_cutoff}")
    print(f"wc_bonus: {params.quality_score.wc_bonus}")
    print(f"helix_break: {params.helix.helix_break}")

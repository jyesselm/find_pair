"""
Modified Nucleotide Registry.

Maps modified residue codes (5MC, 7MG, H2U, PSU, etc.) to their parent base types
for H-bond donor/acceptor determination.
"""

import json
from pathlib import Path
from typing import Dict, Optional, Tuple

# Default path to modified nucleotides config
DEFAULT_CONFIG_PATH = Path(__file__).parent.parent.parent / "resources" / "config" / "modified_nucleotides.json"


class ModifiedResidueRegistry:
    """
    Registry for mapping modified residue codes to parent base properties.

    Uses the modified_nucleotides.json config file which contains 500+ modified
    residue mappings with their parent base types.

    Usage:
        registry = ModifiedResidueRegistry()
        parent = registry.get_parent_base("5MC")  # Returns "C"
        parent = registry.get_parent_base("PSU")  # Returns "P" (pseudouridine)
    """

    # Singleton instance
    _instance: Optional['ModifiedResidueRegistry'] = None

    def __new__(cls, config_path: Optional[Path] = None):
        """Singleton pattern - only one registry instance."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self, config_path: Optional[Path] = None):
        """
        Load registry from modified_nucleotides.json.

        Args:
            config_path: Path to config file (uses default if None)
        """
        if self._initialized:
            return

        self._parent_map: Dict[str, str] = {}
        self._type_map: Dict[str, str] = {}
        self._is_purine_map: Dict[str, bool] = {}

        path = config_path or DEFAULT_CONFIG_PATH
        self._load_config(path)
        self._initialized = True

    def _load_config(self, config_path: Path):
        """Parse modified_nucleotides.json and build lookup tables."""
        if not config_path.exists():
            print(f"Warning: Modified nucleotides config not found: {config_path}")
            return

        with open(config_path) as f:
            data = json.load(f)

        modified_data = data.get("modified_nucleotides", {})

        # Process all sections (standard_nucleotides, modified_adenines, etc.)
        for section_name, section in modified_data.items():
            if not isinstance(section, dict):
                continue

            for residue_name, props in section.items():
                if not isinstance(props, dict):
                    continue

                code = props.get("code", "")
                res_type = props.get("type", "")
                is_purine = props.get("is_purine", False)

                # Normalize residue name
                residue_name = residue_name.strip().upper()

                # Determine parent base from code
                # Lowercase = modified, uppercase = standard
                # 'c' -> 'C', 'g' -> 'G', etc.
                # 'P' -> 'P' (pseudouridine is its own type)
                # 'I' -> 'I' (inosine)
                if code:
                    parent = code.upper()
                    self._parent_map[residue_name] = parent
                    self._type_map[residue_name] = res_type
                    self._is_purine_map[residue_name] = is_purine

    def get_parent_base(self, residue_code: str) -> Optional[str]:
        """
        Get the parent base type for a residue code.

        Args:
            residue_code: 3-letter residue code (e.g., "5MC", "H2U", "PSU")

        Returns:
            Single letter parent base (A, G, C, U, T, I, P) or None if unknown
        """
        code = residue_code.strip().upper()
        return self._parent_map.get(code)

    def get_base_type(self, residue_code: str) -> Optional[str]:
        """
        Get the base type name (ADENINE, GUANINE, etc.)

        Args:
            residue_code: 3-letter residue code

        Returns:
            Type string or None
        """
        code = residue_code.strip().upper()
        return self._type_map.get(code)

    def is_purine(self, residue_code: str) -> Optional[bool]:
        """
        Check if residue is a purine (adenine/guanine-like).

        Args:
            residue_code: 3-letter residue code

        Returns:
            True for purines, False for pyrimidines, None if unknown
        """
        code = residue_code.strip().upper()
        return self._is_purine_map.get(code)

    def is_known(self, residue_code: str) -> bool:
        """Check if residue code is in the registry."""
        code = residue_code.strip().upper()
        return code in self._parent_map

    def is_modified(self, residue_code: str) -> bool:
        """
        Check if residue is a modified nucleotide (not standard A, G, C, U, T).

        Standard residues: A, C, G, U, T, DA, DC, DG, DT, DU
        """
        code = residue_code.strip().upper()
        standard = {'A', 'C', 'G', 'U', 'T', 'DA', 'DC', 'DG', 'DT', 'DU',
                    'ADE', 'CYT', 'GUA', 'URA', 'THY'}
        return code not in standard and code in self._parent_map

    @property
    def count(self) -> int:
        """Number of residues in registry."""
        return len(self._parent_map)

    def __contains__(self, residue_code: str) -> bool:
        """Allow 'in' operator: if '5MC' in registry"""
        return self.is_known(residue_code)


# Global singleton instance for convenience
_global_registry: Optional[ModifiedResidueRegistry] = None


def get_registry() -> ModifiedResidueRegistry:
    """Get the global modified residue registry instance."""
    global _global_registry
    if _global_registry is None:
        _global_registry = ModifiedResidueRegistry()
    return _global_registry


def get_parent_base(residue_code: str) -> Optional[str]:
    """Convenience function to get parent base from global registry."""
    return get_registry().get_parent_base(residue_code)

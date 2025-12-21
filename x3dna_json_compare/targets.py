"""
Comparison target management.

Defines comparison targets (legacy, baseline, etc.) that modern JSON output
can be compared against. Each target corresponds to a directory under data/.

Directory convention: data/json_<target>/
  - data/json_legacy/   -> legacy target
  - data/json_baseline/ -> baseline target
  - data/json_dssr/     -> dssr target (future)
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
import shutil


@dataclass
class ComparisonTarget:
    """Represents a comparison target (legacy, baseline, etc.)."""

    name: str
    directory: str  # e.g., "json_legacy"
    description: str

    def get_path(self, project_root: Path) -> Path:
        """Get full path to target directory."""
        return project_root / "data" / self.directory

    def get_json_file(
        self, project_root: Path, pdb_id: str, record_type: str
    ) -> Path:
        """Get path to a specific JSON file for this target."""
        return self.get_path(project_root) / record_type / f"{pdb_id}.json"

    def exists(self, project_root: Path) -> bool:
        """Check if target directory exists."""
        return self.get_path(project_root).exists()

    def has_data(self, project_root: Path) -> bool:
        """Check if target directory has any JSON files."""
        target_path = self.get_path(project_root)
        if not target_path.exists():
            return False
        # Check for any .json files in subdirectories
        return any(target_path.rglob("*.json"))


# Built-in targets
BUILTIN_TARGETS: Dict[str, ComparisonTarget] = {
    "legacy": ComparisonTarget(
        name="legacy",
        directory="json_legacy",
        description="Legacy C code output (X3DNA v2.4 reference)",
    ),
    "baseline": ComparisonTarget(
        name="baseline",
        directory="json_baseline",
        description="Baseline snapshot of modern output for regression testing",
    ),
}


def get_target(name: str) -> ComparisonTarget:
    """
    Get a comparison target by name.

    Args:
        name: Target name (e.g., "legacy", "baseline")

    Returns:
        ComparisonTarget instance

    Raises:
        ValueError: If target name is unknown
    """
    if name in BUILTIN_TARGETS:
        return BUILTIN_TARGETS[name]
    raise ValueError(
        f"Unknown target: {name}. Available: {list(BUILTIN_TARGETS.keys())}"
    )


def list_targets() -> List[ComparisonTarget]:
    """List all available comparison targets."""
    return list(BUILTIN_TARGETS.values())


def target_exists(name: str, project_root: Path) -> bool:
    """Check if target directory exists with data."""
    try:
        target = get_target(name)
        return target.has_data(project_root)
    except ValueError:
        return False


def create_baseline(
    project_root: Path,
    force: bool = False,
    record_types: Optional[List[str]] = None,
) -> bool:
    """
    Create baseline from current modern output.

    Copies data/json/* to data/json_baseline/*.

    Args:
        project_root: Project root directory
        force: If True, overwrite existing baseline
        record_types: Optional list of record types to copy (default: all)

    Returns:
        True if baseline was created successfully
    """
    modern_path = project_root / "data" / "json"
    baseline_target = get_target("baseline")
    baseline_path = baseline_target.get_path(project_root)

    if not modern_path.exists():
        raise ValueError(f"Modern JSON directory not found: {modern_path}")

    if baseline_path.exists() and not force:
        raise ValueError(
            f"Baseline already exists: {baseline_path}. Use force=True to overwrite."
        )

    # Default record types to copy
    if record_types is None:
        record_types = [
            "pdb_atoms",
            "residue_indices",
            "base_frame_calc",
            "ls_fitting",
            "frame_calc",
            "pair_validation",
            "distance_checks",
            "hbond_list",
            "base_pair",
            "find_bestpair_selection",
            "bpstep_params",
            "helical_params",
        ]

    # Create baseline directory
    baseline_path.mkdir(parents=True, exist_ok=True)

    # Copy each record type directory
    copied_count = 0
    for record_type in record_types:
        src_dir = modern_path / record_type
        dst_dir = baseline_path / record_type

        if src_dir.exists():
            if dst_dir.exists() and force:
                shutil.rmtree(dst_dir)
            if not dst_dir.exists():
                shutil.copytree(src_dir, dst_dir)
                copied_count += 1

    return copied_count > 0


def get_available_targets(project_root: Path) -> List[ComparisonTarget]:
    """
    Get list of targets that have data available.

    Args:
        project_root: Project root directory

    Returns:
        List of targets with data
    """
    return [t for t in list_targets() if t.has_data(project_root)]

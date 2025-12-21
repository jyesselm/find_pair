"""
Shared pytest fixtures and utilities for all tests.

This file is automatically loaded by pytest for all tests in tests_python/.

Fixtures provided:
- project_root_path: Path to project root
- executables: Dict with legacy and modern executable paths
- valid_pdbs_fast: List of fast-running PDB IDs from valid_pdbs_fast.json
- temp_output_dir: Temporary directory cleaned up after each test

Environment Variables:
- FORCE_MODERN_REGEN: Set to '0' to allow caching modern JSON (default: '1' = always regenerate)
  Modern JSON should ALWAYS be regenerated to test current build. Legacy JSON can be cached.

See docs/TESTING_GUIDE.md for testing workflow documentation.
"""
import pytest
import sys
import json
import shutil
import os
from pathlib import Path
from typing import Tuple, Optional, List

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Import test utilities
from scripts.test_utils import (
    find_executables,
    generate_legacy_json,
    generate_modern_json,
    load_valid_pdbs_fast,
    cleanup_output_dir
)


# ============================================================================
# Path Fixtures
# ============================================================================

@pytest.fixture(scope="session")
def project_root_path() -> Path:
    """Return the project root directory."""
    return project_root


@pytest.fixture(scope="session")
def data_path(project_root_path) -> Path:
    """Return the data directory."""
    return project_root_path / "data"


@pytest.fixture(scope="session")
def pdb_dir(data_path) -> Path:
    """Return the PDB files directory."""
    return data_path / "pdb"


@pytest.fixture(scope="session")
def legacy_json_dir(data_path) -> Path:
    """Return the legacy JSON directory."""
    return data_path / "json_legacy"


@pytest.fixture(scope="session")
def modern_json_dir(data_path) -> Path:
    """Return the modern JSON directory."""
    return data_path / "json"


# ============================================================================
# Executable Fixtures
# ============================================================================

@pytest.fixture(scope="session")
def executables(project_root_path) -> dict:
    """Find and return legacy and modern executables.
    
    Returns:
        dict with 'legacy' and 'modern' keys containing Path objects or None
    """
    legacy_exe, modern_exe = find_executables(project_root_path)
    
    if legacy_exe is None:
        pytest.skip("Legacy executable not found (org/build/bin/find_pair_analyze)")
    if modern_exe is None:
        pytest.skip("Modern executable not found (build/generate_modern_json)")
    
    return {
        "legacy": legacy_exe,
        "modern": modern_exe
    }


@pytest.fixture(scope="session")
def legacy_exe(executables) -> Path:
    """Return the legacy executable path."""
    return executables["legacy"]


@pytest.fixture(scope="session")
def modern_exe(executables) -> Path:
    """Return the modern executable path."""
    return executables["modern"]


# ============================================================================
# PDB List Fixtures
# ============================================================================

@pytest.fixture(scope="session")
def valid_pdbs_fast(project_root_path) -> List[str]:
    """Load valid_pdbs_fast.json and return list of PDB IDs.
    
    These are PDBs that have been validated to work correctly and run
    fast enough for routine testing.
    """
    try:
        return load_valid_pdbs_fast(project_root_path)
    except (FileNotFoundError, ValueError) as e:
        pytest.skip(f"Could not load valid_pdbs_fast.json: {e}")


@pytest.fixture(scope="session")
def test_set_10(valid_pdbs_fast) -> List[str]:
    """Return a test set of 10 PDBs for quick tests."""
    return valid_pdbs_fast[:10]


@pytest.fixture(scope="session")
def test_set_50(valid_pdbs_fast) -> List[str]:
    """Return a test set of 50 PDBs for moderate tests."""
    return valid_pdbs_fast[:50]


@pytest.fixture(scope="session")
def test_set_100(valid_pdbs_fast) -> List[str]:
    """Return a test set of 100 PDBs for standard tests."""
    return valid_pdbs_fast[:100]


# ============================================================================
# Temporary Directory Fixtures
# ============================================================================

@pytest.fixture(scope="function")
def temp_output_dir(tmp_path):
    """Create a temporary output directory for test files.
    
    Automatically cleaned up after each test.
    """
    output_dir = tmp_path / "test_output"
    output_dir.mkdir(parents=True, exist_ok=True)
    yield output_dir
    # Cleanup after test
    if output_dir.exists():
        shutil.rmtree(output_dir)


@pytest.fixture(scope="session")
def shared_temp_dir(tmp_path_factory):
    """Create a shared temporary directory for session-scoped tests.
    
    Useful for tests that need to share generated files.
    """
    return tmp_path_factory.mktemp("shared_output")


# ============================================================================
# PDB File Fixtures
# ============================================================================

@pytest.fixture
def pdb_file_path(project_root_path, pdb_id):
    """Get path to PDB file for a given PDB ID.
    
    Requires pdb_id to be provided by test parametrization.
    """
    pdb_file = project_root_path / "data" / "pdb" / f"{pdb_id}.pdb"
    if not pdb_file.exists():
        pytest.skip(f"PDB file not found: {pdb_file}")
    return pdb_file


def get_pdb_path(project_root: Path, pdb_id: str) -> Path:
    """Helper to get PDB file path."""
    return project_root / "data" / "pdb" / f"{pdb_id}.pdb"


# ============================================================================
# JSON Loading Fixtures
# ============================================================================

@pytest.fixture(scope="session")
def load_legacy_json(legacy_json_dir):
    """Factory fixture to load legacy JSON for a given type and PDB."""
    def _load(json_type: str, pdb_id: str) -> dict:
        json_file = legacy_json_dir / json_type / f"{pdb_id}.json"
        if not json_file.exists():
            return None
        with open(json_file) as f:
            return json.load(f)
    return _load


@pytest.fixture(scope="session")
def load_modern_json(modern_json_dir):
    """Factory fixture to load modern JSON for a given type and PDB."""
    def _load(json_type: str, pdb_id: str) -> dict:
        json_file = modern_json_dir / json_type / f"{pdb_id}.json"
        if not json_file.exists():
            return None
        with open(json_file) as f:
            return json.load(f)
    return _load


# ============================================================================
# Pytest Configuration
# ============================================================================

def pytest_addoption(parser):
    """Add custom pytest command-line options."""
    parser.addoption(
        "--max-pdbs",
        action="store",
        default=None,
        type=int,
        help="Maximum number of PDBs to test"
    )
    parser.addoption(
        "--pdb-id",
        action="store",
        default=None,
        help="Test a specific PDB ID"
    )
    parser.addoption(
        "--test-set",
        action="store",
        default=None,
        choices=["10", "50", "100", "500", "1000"],
        help="Use a saved test set of specified size"
    )


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line("markers", "integration: Integration tests requiring executables")
    config.addinivalue_line("markers", "slow: Tests that take a long time")
    config.addinivalue_line("markers", "requires_legacy: Requires legacy executable")
    config.addinivalue_line("markers", "requires_modern: Requires modern executable")

    # Ensure FORCE_MODERN_REGEN is set (default to '1' = always regenerate)
    # This ensures tests always use the current build, not cached files
    if 'FORCE_MODERN_REGEN' not in os.environ:
        os.environ['FORCE_MODERN_REGEN'] = '1'


def pytest_collection_modifyitems(config, items):
    """Modify test collection based on command-line options."""
    # If --pdb is specified, filter tests to that PDB
    pdb_id = config.getoption("--pdb")
    if pdb_id:
        # This would need to be implemented per test file
        pass
    
    # Apply markers based on test file location
    for item in items:
        if "integration" in str(item.fspath):
            item.add_marker(pytest.mark.integration)


# ============================================================================
# Comparison Utilities (for use in fixtures)
# ============================================================================

TOLERANCE_DEFAULT = 1e-6
TOLERANCE_RMS = 0.001
TOLERANCE_MATRIX = 1e-4


def values_equal(a, b, tolerance: float = TOLERANCE_DEFAULT) -> bool:
    """Compare two values with tolerance for floats."""
    if a is None and b is None:
        return True
    if a is None or b is None:
        return False
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return abs(a - b) <= tolerance
    if isinstance(a, str) and isinstance(b, str):
        return a.strip() == b.strip()
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return False
        return all(values_equal(x, y, tolerance) for x, y in zip(a, b))
    return a == b


def matrix_equal(a: list, b: list, tolerance: float = TOLERANCE_MATRIX) -> bool:
    """Compare two matrices (lists of lists) with tolerance."""
    if len(a) != len(b):
        return False
    for row_a, row_b in zip(a, b):
        if len(row_a) != len(row_b):
            return False
        for val_a, val_b in zip(row_a, row_b):
            if abs(val_a - val_b) > tolerance:
                return False
    return True

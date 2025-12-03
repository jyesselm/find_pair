"""
Shared pytest fixtures and utilities for all tests.
"""
import pytest
import sys
from pathlib import Path
import shutil
import subprocess
from typing import Tuple, Optional

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


@pytest.fixture(scope="session")
def project_root_path():
    """Return the project root directory."""
    return project_root


@pytest.fixture(scope="session")
def executables(project_root_path):
    """Find and return legacy and modern executables."""
    legacy_exe, modern_exe = find_executables(project_root_path)
    
    if legacy_exe is None:
        pytest.skip("Legacy executable not found")
    if modern_exe is None:
        pytest.skip("Modern executable not found")
    
    return {
        "legacy": legacy_exe,
        "modern": modern_exe
    }


@pytest.fixture(scope="function")
def temp_output_dir(tmp_path):
    """Create a temporary output directory for test files."""
    output_dir = tmp_path / "test_output"
    output_dir.mkdir(parents=True, exist_ok=True)
    yield output_dir
    # Cleanup after test
    if output_dir.exists():
        shutil.rmtree(output_dir)


@pytest.fixture(scope="session")
def valid_pdbs_fast(project_root_path):
    """Load valid_pdbs_fast.json and return list of PDB IDs."""
    try:
        return load_valid_pdbs_fast(project_root_path)
    except (FileNotFoundError, ValueError) as e:
        pytest.skip(f"Could not load valid_pdbs_fast.json: {e}")


@pytest.fixture
def pdb_file_path(project_root_path, pdb_id):
    """Get path to PDB file for a given PDB ID."""
    pdb_file = project_root_path / "data" / "pdb" / f"{pdb_id}.pdb"
    if not pdb_file.exists():
        pytest.skip(f"PDB file not found: {pdb_file}")
    return pdb_file


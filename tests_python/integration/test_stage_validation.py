#!/usr/bin/env python3
"""
Integration tests for stage validation.

This module provides pytest tests and a CLI interface for validating
legacy vs modern JSON outputs across all 12 stages.

Usage:
    # Via pytest
    pytest tests_python/integration/test_stage_validation.py -v
    
    # Direct CLI
    python tests_python/integration/test_stage_validation.py 3 4 5
    python tests_python/integration/test_stage_validation.py frames --max-pdbs 50
"""

import sys
from pathlib import Path
import pytest

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from validation import validate_stage
from validation.cli import main as cli_main


# Test configuration
PROJECT_ROOT = Path(__file__).parent.parent.parent
JSON_DIR = PROJECT_ROOT / "data" / "json"


def get_test_pdbs(max_pdbs: int = 50) -> list:
    """Get list of PDB IDs for testing."""
    pdb_atoms_dir = JSON_DIR / "pdb_atoms"
    if not pdb_atoms_dir.exists():
        return ["100D", "157D", "1A34", "1EFW", "1F8V"]
    
    pdbs = sorted(
        f.stem for f in pdb_atoms_dir.iterdir()
        if f.suffix == ".json" and len(f.stem) == 4
    )
    return pdbs[:max_pdbs] if pdbs else ["100D"]


# --- Pytest Tests ---

class TestStageValidation:
    """Pytest tests for each validation stage."""
    
    def test_stage_1_atoms(self):
        """Test Stage 1: Atom parsing."""
        result = validate_stage(1, get_test_pdbs(10), JSON_DIR)
        assert result.failed == 0, f"Stage 1 failed for: {result.failed_pdbs[:5]}"
    
    def test_stage_2_residues(self):
        """Test Stage 2: Residue indices."""
        result = validate_stage(2, get_test_pdbs(10), JSON_DIR)
        assert result.failed == 0, f"Stage 2 failed for: {result.failed_pdbs[:5]}"
    
    def test_stage_3_base_frame_calc(self):
        """Test Stage 3: Base frame calculation."""
        result = validate_stage(3, get_test_pdbs(10), JSON_DIR)
        assert result.failed == 0, f"Stage 3 failed for: {result.failed_pdbs[:5]}"
    
    def test_stage_4_ls_fitting(self):
        """Test Stage 4: Least squares fitting."""
        result = validate_stage(4, get_test_pdbs(10), JSON_DIR)
        assert result.failed == 0, f"Stage 4 failed for: {result.failed_pdbs[:5]}"
    
    def test_stage_5_frame_calc(self):
        """Test Stage 5: Frame calculation."""
        result = validate_stage(5, get_test_pdbs(10), JSON_DIR)
        assert result.failed == 0, f"Stage 5 failed for: {result.failed_pdbs[:5]}"
    
    def test_stage_6_pair_validation(self):
        """Test Stage 6: Pair validation."""
        result = validate_stage(6, get_test_pdbs(10), JSON_DIR)
        assert result.failed == 0, f"Stage 6 failed for: {result.failed_pdbs[:5]}"
    
    def test_stage_7_distance_checks(self):
        """Test Stage 7: Distance checks."""
        result = validate_stage(7, get_test_pdbs(10), JSON_DIR)
        assert result.failed == 0, f"Stage 7 failed for: {result.failed_pdbs[:5]}"
    
    @pytest.mark.xfail(reason="Known issue: modern uses stricter H-bond criteria than legacy")
    def test_stage_8_hbonds(self):
        """Test Stage 8: Hydrogen bonds."""
        result = validate_stage(8, get_test_pdbs(10), JSON_DIR)
        assert result.failed == 0, f"Stage 8 failed for: {result.failed_pdbs[:5]}"
    
    def test_stage_9_base_pair(self):
        """Test Stage 9: Base pairs."""
        result = validate_stage(9, get_test_pdbs(10), JSON_DIR)
        assert result.failed == 0, f"Stage 9 failed for: {result.failed_pdbs[:5]}"
    
    def test_stage_10_best_pair(self):
        """Test Stage 10: Best pair selection."""
        result = validate_stage(10, get_test_pdbs(10), JSON_DIR)
        assert result.failed == 0, f"Stage 10 failed for: {result.failed_pdbs[:5]}"
    
    def test_stage_11_step_params(self):
        """Test Stage 11: Step parameters."""
        result = validate_stage(11, get_test_pdbs(10), JSON_DIR)
        assert result.failed == 0, f"Stage 11 failed for: {result.failed_pdbs[:5]}"
    
    def test_stage_12_helical_params(self):
        """Test Stage 12: Helical parameters."""
        result = validate_stage(12, get_test_pdbs(10), JSON_DIR)
        assert result.failed == 0, f"Stage 12 failed for: {result.failed_pdbs[:5]}"


# --- CLI Entry Point ---

if __name__ == "__main__":
    cli_main()

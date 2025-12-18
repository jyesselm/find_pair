#!/usr/bin/env python3
"""
Unit tests for comparators with dummy data.

These tests verify that each comparator correctly:
1. Detects missing required fields
2. Detects mismatched values
3. Passes values within tolerance
4. Fails values outside tolerance
"""

import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from validation.comparators.atoms import compare_atoms
from validation.comparators.residues import compare_residues
from validation.comparators.frames import (
    compare_base_frame_calc,
    compare_ls_fitting,
    compare_frame_calc
)
from validation.comparators.pairs import (
    compare_pair_validation,
    compare_distance_checks,
    compare_base_pair,
    compare_best_pair_selection
)
from validation.comparators.hbonds import compare_hbonds
from validation.comparators.params import compare_step_params, compare_helical_params


class TestAtomComparator:
    """Test Stage 1: pdb_atoms comparator."""
    
    def test_matching_atoms_pass(self):
        """Identical atoms should pass."""
        legacy = [
            {"legacy_atom_idx": 1, "xyz": [1.0, 2.0, 3.0], 
             "atom_name": " N  ", "residue_name": "ADE", "chain_id": "A"}
        ]
        modern = [
            {"legacy_atom_idx": 1, "xyz": [1.0, 2.0, 3.0],
             "atom_name": " N  ", "residue_name": "ADE", "chain_id": "A"}
        ]
        passed, errors = compare_atoms(legacy, modern)
        assert passed, f"Should pass: {errors}"
    
    def test_xyz_mismatch_detected(self):
        """XYZ coordinate mismatch should be detected."""
        legacy = [{"legacy_atom_idx": 1, "xyz": [1.0, 2.0, 3.0],
                   "atom_name": " N  ", "residue_name": "ADE", "chain_id": "A"}]
        modern = [{"legacy_atom_idx": 1, "xyz": [1.1, 2.0, 3.0],  # 0.1 diff
                   "atom_name": " N  ", "residue_name": "ADE", "chain_id": "A"}]
        passed, errors = compare_atoms(legacy, modern)
        assert not passed, "Should fail on xyz mismatch"
        assert any("xyz" in e for e in errors), f"Should mention xyz: {errors}"
    
    def test_atom_name_mismatch_detected(self):
        """Atom name mismatch should be detected."""
        legacy = [{"legacy_atom_idx": 1, "xyz": [1.0, 2.0, 3.0],
                   "atom_name": " N  ", "residue_name": "ADE", "chain_id": "A"}]
        modern = [{"legacy_atom_idx": 1, "xyz": [1.0, 2.0, 3.0],
                   "atom_name": " C  ", "residue_name": "ADE", "chain_id": "A"}]
        passed, errors = compare_atoms(legacy, modern)
        assert not passed, "Should fail on atom_name mismatch"
        assert any("atom_name" in e for e in errors), f"Should mention atom_name: {errors}"
    
    def test_residue_name_mismatch_detected(self):
        """Residue name mismatch should be detected."""
        legacy = [{"legacy_atom_idx": 1, "xyz": [1.0, 2.0, 3.0],
                   "atom_name": " N  ", "residue_name": "ADE", "chain_id": "A"}]
        modern = [{"legacy_atom_idx": 1, "xyz": [1.0, 2.0, 3.0],
                   "atom_name": " N  ", "residue_name": "GUA", "chain_id": "A"}]
        passed, errors = compare_atoms(legacy, modern)
        assert not passed, "Should fail on residue_name mismatch"
        assert any("residue_name" in e for e in errors), f"Should mention residue_name: {errors}"
    
    def test_chain_id_mismatch_detected(self):
        """Chain ID mismatch should be detected."""
        legacy = [{"legacy_atom_idx": 1, "xyz": [1.0, 2.0, 3.0],
                   "atom_name": " N  ", "residue_name": "ADE", "chain_id": "A"}]
        modern = [{"legacy_atom_idx": 1, "xyz": [1.0, 2.0, 3.0],
                   "atom_name": " N  ", "residue_name": "ADE", "chain_id": "B"}]
        passed, errors = compare_atoms(legacy, modern)
        assert not passed, "Should fail on chain_id mismatch"
        assert any("chain_id" in e for e in errors), f"Should mention chain_id: {errors}"
    
    def test_missing_xyz_detected(self):
        """Missing xyz field should be detected."""
        legacy = [{"legacy_atom_idx": 1, "xyz": [1.0, 2.0, 3.0],
                   "atom_name": " N  ", "residue_name": "ADE", "chain_id": "A"}]
        modern = [{"legacy_atom_idx": 1,  # Missing xyz
                   "atom_name": " N  ", "residue_name": "ADE", "chain_id": "A"}]
        passed, errors = compare_atoms(legacy, modern)
        assert not passed, "Should fail on missing xyz"
        assert any("missing" in e.lower() and "xyz" in e for e in errors), f"Should mention missing xyz: {errors}"


class TestResidueComparator:
    """Test Stage 2: residue_indices comparator.
    
    Key: residue_idx
    Checks: start_atom, end_atom
    """
    
    def test_matching_residues_pass(self):
        """Identical residues should pass."""
        legacy = [{"residue_idx": 1, "start_atom": 1, "end_atom": 23}]
        modern = [{"residue_idx": 1, "start_atom": 1, "end_atom": 23}]
        passed, errors = compare_residues(legacy, modern)
        assert passed, f"Should pass: {errors}"
    
    def test_start_atom_mismatch_detected(self):
        """Start atom mismatch should be detected."""
        legacy = [{"residue_idx": 1, "start_atom": 1, "end_atom": 23}]
        modern = [{"residue_idx": 1, "start_atom": 5, "end_atom": 23}]
        passed, errors = compare_residues(legacy, modern)
        assert not passed, "Should fail on start_atom mismatch"
        assert any("start_atom" in e for e in errors), f"Should mention start_atom: {errors}"
    
    def test_end_atom_mismatch_detected(self):
        """End atom mismatch should be detected."""
        legacy = [{"residue_idx": 1, "start_atom": 1, "end_atom": 23}]
        modern = [{"residue_idx": 1, "start_atom": 1, "end_atom": 25}]
        passed, errors = compare_residues(legacy, modern)
        assert not passed, "Should fail on end_atom mismatch"
        assert any("end_atom" in e for e in errors), f"Should mention end_atom: {errors}"


class TestBaseFrameCalcComparator:
    """Test Stage 3: base_frame_calc comparator.
    
    Key: residue_idx (legacy) or legacy_residue_idx (modern)
    """
    
    def test_matching_base_frame_pass(self):
        """Identical base_frame_calc should pass."""
        legacy = [{"residue_idx": 1,
                   "base_type": "A", "rms_fit": 0.025, "num_matched_atoms": 10,
                   "matched_atoms": ["N1", "C2", "N3"], "standard_template": "Atomic_A.pdb"}]
        modern = [{"legacy_residue_idx": 1,
                   "base_type": "A", "rms_fit": 0.025, "num_matched_atoms": 10,
                   "matched_atoms": ["N1", "C2", "N3"], "standard_template": "data/templates/Atomic_A.pdb"}]
        passed, errors = compare_base_frame_calc(legacy, modern)
        assert passed, f"Should pass: {errors}"
    
    def test_base_type_mismatch_detected(self):
        """Base type mismatch should be detected."""
        legacy = [{"residue_idx": 1,
                   "base_type": "A", "rms_fit": 0.025, "num_matched_atoms": 10,
                   "matched_atoms": ["N1", "C2", "N3"]}]
        modern = [{"legacy_residue_idx": 1,
                   "base_type": "G", "rms_fit": 0.025, "num_matched_atoms": 10,
                   "matched_atoms": ["N1", "C2", "N3"]}]
        passed, errors = compare_base_frame_calc(legacy, modern)
        assert not passed, "Should fail on base_type mismatch"
        assert any("base_type" in e for e in errors), f"Should mention base_type: {errors}"
    
    def test_rms_fit_mismatch_detected(self):
        """RMS fit mismatch beyond tolerance should be detected."""
        legacy = [{"residue_idx": 1,
                   "base_type": "A", "rms_fit": 0.025, "num_matched_atoms": 10,
                   "matched_atoms": ["N1", "C2", "N3"]}]
        modern = [{"legacy_residue_idx": 1,
                   "base_type": "A", "rms_fit": 0.030, "num_matched_atoms": 10,  # 0.005 diff > 0.001
                   "matched_atoms": ["N1", "C2", "N3"]}]
        passed, errors = compare_base_frame_calc(legacy, modern)
        assert not passed, "Should fail on rms_fit mismatch"
        assert any("rms_fit" in e for e in errors), f"Should mention rms_fit: {errors}"
    
    def test_num_matched_atoms_mismatch_detected(self):
        """Num matched atoms mismatch should be detected."""
        legacy = [{"residue_idx": 1,
                   "base_type": "A", "rms_fit": 0.025, "num_matched_atoms": 10,
                   "matched_atoms": ["N1", "C2", "N3"]}]
        modern = [{"legacy_residue_idx": 1,
                   "base_type": "A", "rms_fit": 0.025, "num_matched_atoms": 9,
                   "matched_atoms": ["N1", "C2", "N3"]}]
        passed, errors = compare_base_frame_calc(legacy, modern)
        assert not passed, "Should fail on num_matched_atoms mismatch"
        assert any("num_matched_atoms" in e for e in errors), f"Should mention num_matched_atoms: {errors}"
    
    def test_matched_atoms_mismatch_detected(self):
        """Matched atoms set mismatch should be detected."""
        legacy = [{"residue_idx": 1,
                   "base_type": "A", "rms_fit": 0.025, "num_matched_atoms": 3,
                   "matched_atoms": ["N1", "C2", "N3"]}]
        modern = [{"legacy_residue_idx": 1,
                   "base_type": "A", "rms_fit": 0.025, "num_matched_atoms": 3,
                   "matched_atoms": ["N1", "C2", "C4"]}]  # C4 instead of N3
        passed, errors = compare_base_frame_calc(legacy, modern)
        assert not passed, "Should fail on matched_atoms mismatch"
        assert any("matched_atoms" in e for e in errors), f"Should mention matched_atoms: {errors}"
    
    def test_missing_base_type_detected(self):
        """Missing base_type should be detected."""
        legacy = [{"residue_idx": 1,
                   "base_type": "A", "rms_fit": 0.025, "num_matched_atoms": 10,
                   "matched_atoms": ["N1", "C2", "N3"]}]
        modern = [{"legacy_residue_idx": 1,
                   # Missing base_type
                   "rms_fit": 0.025, "num_matched_atoms": 10,
                   "matched_atoms": ["N1", "C2", "N3"]}]
        passed, errors = compare_base_frame_calc(legacy, modern)
        assert not passed, "Should fail on missing base_type"
        assert any("missing" in e.lower() and "base_type" in e for e in errors), f"Should mention missing base_type: {errors}"
    
    def test_standard_template_filename_only(self):
        """standard_template should compare filename only, not full path."""
        legacy = [{"residue_idx": 1, "base_type": "A", "rms_fit": 0.025, 
                   "num_matched_atoms": 10, "matched_atoms": ["N1", "C2", "N3"],
                   "standard_template": "/usr/local/x3dna/config/Atomic_A.pdb"}]
        modern = [{"legacy_residue_idx": 1, "base_type": "A", "rms_fit": 0.025,
                   "num_matched_atoms": 10, "matched_atoms": ["N1", "C2", "N3"],
                   "standard_template": "data/templates/Atomic_A.pdb"}]  # Different path, same filename
        passed, errors = compare_base_frame_calc(legacy, modern)
        assert passed, f"Should pass (same filename): {errors}"


class TestLSFittingComparator:
    """Test Stage 4: ls_fitting comparator.
    
    Key: residue_idx (legacy) or legacy_residue_idx (modern)
    Note: Legacy does NOT have base_type field.
    """
    
    def test_matching_ls_fitting_pass(self):
        """Identical ls_fitting should pass."""
        legacy = [{"residue_idx": 1, "rms_fit": 0.025, "num_points": 10,
                   "rotation_matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                   "translation": [1.0, 2.0, 3.0]}]
        modern = [{"legacy_residue_idx": 1, "rms_fit": 0.025, "num_points": 10,
                   "rotation_matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                   "translation": [1.0, 2.0, 3.0]}]
        passed, errors = compare_ls_fitting(legacy, modern)
        assert passed, f"Should pass: {errors}"
    
    def test_rms_fit_mismatch_detected(self):
        """RMS fit mismatch should be detected."""
        legacy = [{"residue_idx": 1, "rms_fit": 0.025, "num_points": 10,
                   "rotation_matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                   "translation": [1.0, 2.0, 3.0]}]
        modern = [{"legacy_residue_idx": 1, "rms_fit": 0.030, "num_points": 10,  # 0.005 diff > 0.001
                   "rotation_matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                   "translation": [1.0, 2.0, 3.0]}]
        passed, errors = compare_ls_fitting(legacy, modern)
        assert not passed, "Should fail on rms_fit mismatch"
        assert any("rms_fit" in e for e in errors), f"Should mention rms_fit: {errors}"
    
    def test_rotation_matrix_mismatch_detected(self):
        """Rotation matrix mismatch should be detected."""
        legacy = [{"residue_idx": 1, "rms_fit": 0.025, "num_points": 10,
                   "rotation_matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                   "translation": [1.0, 2.0, 3.0]}]
        modern = [{"legacy_residue_idx": 1, "rms_fit": 0.025, "num_points": 10,
                   "rotation_matrix": [[0.9, 0, 0], [0, 1, 0], [0, 0, 1]],  # 0.1 diff > 1e-4
                   "translation": [1.0, 2.0, 3.0]}]
        passed, errors = compare_ls_fitting(legacy, modern)
        assert not passed, "Should fail on rotation_matrix mismatch"
        assert any("rotation_matrix" in e for e in errors), f"Should mention rotation_matrix: {errors}"
    
    def test_translation_mismatch_detected(self):
        """Translation mismatch should be detected."""
        legacy = [{"residue_idx": 1, "rms_fit": 0.025, "num_points": 10,
                   "rotation_matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                   "translation": [1.0, 2.0, 3.0]}]
        modern = [{"legacy_residue_idx": 1, "rms_fit": 0.025, "num_points": 10,
                   "rotation_matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                   "translation": [1.1, 2.0, 3.0]}]  # 0.1 diff > 1e-6
        passed, errors = compare_ls_fitting(legacy, modern)
        assert not passed, "Should fail on translation mismatch"
        assert any("translation" in e for e in errors), f"Should mention translation: {errors}"


class TestFrameCalcComparator:
    """Test Stage 5: frame_calc comparator.
    
    Key: residue_idx (legacy) or legacy_residue_idx (modern)
    Checks: base_type, rms_fit, num_matched_atoms
    """
    
    def test_matching_frame_calc_pass(self):
        """Identical frame_calc should pass."""
        legacy = [{"residue_idx": 1, "base_type": "A", "rms_fit": 0.025,
                   "num_matched_atoms": 10}]
        modern = [{"legacy_residue_idx": 1, "base_type": "A", "rms_fit": 0.025,
                   "num_matched_atoms": 10}]
        passed, errors = compare_frame_calc(legacy, modern)
        assert passed, f"Should pass: {errors}"
    
    def test_base_type_mismatch_detected(self):
        """Base type mismatch should be detected."""
        legacy = [{"residue_idx": 1, "base_type": "g", "rms_fit": 0.025,
                   "num_matched_atoms": 10}]
        modern = [{"legacy_residue_idx": 1, "base_type": "?", "rms_fit": 0.025,
                   "num_matched_atoms": 10}]
        passed, errors = compare_frame_calc(legacy, modern)
        assert not passed, "Should fail on base_type mismatch"
        assert any("base_type" in e for e in errors), f"Should mention base_type: {errors}"
    
    def test_rms_fit_mismatch_detected(self):
        """RMS fit mismatch should be detected."""
        legacy = [{"residue_idx": 1, "base_type": "A", "rms_fit": 0.025,
                   "num_matched_atoms": 10}]
        modern = [{"legacy_residue_idx": 1, "base_type": "A", "rms_fit": 0.030,  # 0.005 diff > 0.001
                   "num_matched_atoms": 10}]
        passed, errors = compare_frame_calc(legacy, modern)
        assert not passed, "Should fail on rms_fit mismatch"
        assert any("rms_fit" in e for e in errors), f"Should mention rms_fit: {errors}"
    
    def test_num_matched_atoms_mismatch_detected(self):
        """num_matched_atoms mismatch should be detected."""
        legacy = [{"residue_idx": 1, "base_type": "A", "rms_fit": 0.025,
                   "num_matched_atoms": 10}]
        modern = [{"legacy_residue_idx": 1, "base_type": "A", "rms_fit": 0.025,
                   "num_matched_atoms": 9}]
        passed, errors = compare_frame_calc(legacy, modern)
        assert not passed, "Should fail on num_matched_atoms mismatch"
        assert any("num_matched_atoms" in e for e in errors), f"Should mention num_matched_atoms: {errors}"


class TestDistanceChecksComparator:
    """Test Stage 7: distance_checks comparator.
    
    Key: (base_i, base_j)
    Checks: values.{dorg, dNN, plane_angle, d_v, overlap_area}
    """
    
    def test_matching_distances_pass(self):
        """Identical distance checks should pass."""
        legacy = [{"base_i": 1, "base_j": 10,
                   "values": {"dorg": 5.5, "dNN": 8.9, "plane_angle": 15.0, "d_v": 0.5, "overlap_area": 0.0}}]
        modern = [{"base_i": 1, "base_j": 10,
                   "values": {"dorg": 5.5, "dNN": 8.9, "plane_angle": 15.0, "d_v": 0.5, "overlap_area": 0.0}}]
        passed, errors = compare_distance_checks(legacy, modern)
        assert passed, f"Should pass: {errors}"
    
    def test_dorg_mismatch_detected(self):
        """dorg mismatch should be detected."""
        legacy = [{"base_i": 1, "base_j": 10,
                   "values": {"dorg": 5.5, "dNN": 8.9, "plane_angle": 15.0, "d_v": 0.5, "overlap_area": 0.0}}]
        modern = [{"base_i": 1, "base_j": 10,
                   "values": {"dorg": 5.6, "dNN": 8.9, "plane_angle": 15.0, "d_v": 0.5, "overlap_area": 0.0}}]  # 0.1 diff
        passed, errors = compare_distance_checks(legacy, modern)
        assert not passed, "Should fail on dorg mismatch"
        assert any("dorg" in e for e in errors), f"Should mention dorg: {errors}"
    
    def test_dNN_mismatch_detected(self):
        """dNN mismatch should be detected."""
        legacy = [{"base_i": 1, "base_j": 10,
                   "values": {"dorg": 5.5, "dNN": 8.9, "plane_angle": 15.0, "d_v": 0.5, "overlap_area": 0.0}}]
        modern = [{"base_i": 1, "base_j": 10,
                   "values": {"dorg": 5.5, "dNN": 9.0, "plane_angle": 15.0, "d_v": 0.5, "overlap_area": 0.0}}]
        passed, errors = compare_distance_checks(legacy, modern)
        assert not passed, "Should fail on dNN mismatch"
        assert any("dNN" in e for e in errors), f"Should mention dNN: {errors}"
    
    def test_plane_angle_mismatch_detected(self):
        """plane_angle mismatch should be detected."""
        legacy = [{"base_i": 1, "base_j": 10,
                   "values": {"dorg": 5.5, "dNN": 8.9, "plane_angle": 15.0, "d_v": 0.5, "overlap_area": 0.0}}]
        modern = [{"base_i": 1, "base_j": 10,
                   "values": {"dorg": 5.5, "dNN": 8.9, "plane_angle": 16.0, "d_v": 0.5, "overlap_area": 0.0}}]
        passed, errors = compare_distance_checks(legacy, modern)
        assert not passed, "Should fail on plane_angle mismatch"
        assert any("plane_angle" in e for e in errors), f"Should mention plane_angle: {errors}"
    
    def test_d_v_mismatch_detected(self):
        """d_v mismatch should be detected."""
        legacy = [{"base_i": 1, "base_j": 10,
                   "values": {"dorg": 5.5, "dNN": 8.9, "plane_angle": 15.0, "d_v": 0.5, "overlap_area": 0.0}}]
        modern = [{"base_i": 1, "base_j": 10,
                   "values": {"dorg": 5.5, "dNN": 8.9, "plane_angle": 15.0, "d_v": 0.6, "overlap_area": 0.0}}]
        passed, errors = compare_distance_checks(legacy, modern)
        assert not passed, "Should fail on d_v mismatch"
        assert any("d_v" in e for e in errors), f"Should mention d_v: {errors}"


class TestHbondComparator:
    """Test Stage 8: hbond_list comparator.

    Key: (base_i, base_j) - residue indices
    Checks: hbonds[].{donor_atom, acceptor_atom, distance}
    """

    def test_matching_hbonds_pass(self):
        """Identical hbonds should pass."""
        legacy = [{"base_i": 1, "base_j": 10, "num_hbonds": 2,
                   "hbonds": [{"donor_atom": "N6", "acceptor_atom": "O4", "distance": 2.9},
                             {"donor_atom": "N1", "acceptor_atom": "N3", "distance": 2.8}]}]
        modern = [{"base_i": 1, "base_j": 10, "num_hbonds": 2,
                   "hbonds": [{"donor_atom": "N6", "acceptor_atom": "O4", "distance": 2.9},
                             {"donor_atom": "N1", "acceptor_atom": "N3", "distance": 2.8}]}]
        passed, errors = compare_hbonds(legacy, modern)
        assert passed, f"Should pass: {errors}"

    def test_num_hbonds_mismatch_detected(self):
        """hbond count mismatch should be detected (via actual hbonds list length)."""
        legacy = [{"base_i": 1, "base_j": 10, "num_hbonds": 2,
                   "hbonds": [{"donor_atom": "N6", "acceptor_atom": "O4", "distance": 2.9},
                             {"donor_atom": "N1", "acceptor_atom": "N3", "distance": 2.8}]}]
        modern = [{"base_i": 1, "base_j": 10, "num_hbonds": 3,
                   "hbonds": [{"donor_atom": "N6", "acceptor_atom": "O4", "distance": 2.9},
                             {"donor_atom": "N1", "acceptor_atom": "N3", "distance": 2.8},
                             {"donor_atom": "N2", "acceptor_atom": "O2", "distance": 3.0}]}]
        passed, errors = compare_hbonds(legacy, modern)
        assert not passed, "Should fail on hbond count mismatch"
        assert any("hbond" in e.lower() and "count" in e.lower() for e in errors), f"Should mention hbond count: {errors}"

    def test_donor_atom_mismatch_detected(self):
        """donor_atom mismatch should be detected."""
        legacy = [{"base_i": 1, "base_j": 10, "num_hbonds": 1,
                   "hbonds": [{"donor_atom": "N6", "acceptor_atom": "O4", "distance": 2.9}]}]
        modern = [{"base_i": 1, "base_j": 10, "num_hbonds": 1,
                   "hbonds": [{"donor_atom": "N1", "acceptor_atom": "O4", "distance": 2.9}]}]  # Different donor
        passed, errors = compare_hbonds(legacy, modern)
        assert not passed, "Should fail on donor_atom mismatch"
        assert any("atom" in e.lower() for e in errors), f"Should mention atom mismatch: {errors}"


class TestStepParamsComparator:
    """Test Stage 11: bpstep_params comparator.

    Key: (bp_idx1, bp_idx2) - step indices
    Legacy format: params dict with capitalized keys {"Shift": ..., "Slide": ...}
    Modern format: lowercase individual fields {"shift": ..., "slide": ...}
    """

    def test_matching_params_pass(self):
        """Identical step params should pass."""
        legacy = [{"bp_idx1": 1, "bp_idx2": 2,
                   "params": {"Shift": 0.1, "Slide": -1.5, "Rise": 3.3,
                              "Tilt": -2.0, "Roll": 8.5, "Twist": 31.0}}]
        modern = [{"bp_idx1": 1, "bp_idx2": 2,
                   "shift": 0.1, "slide": -1.5, "rise": 3.3,
                   "tilt": -2.0, "roll": 8.5, "twist": 31.0}]
        passed, errors = compare_step_params(legacy, modern)
        assert passed, f"Should pass: {errors}"

    def test_shift_mismatch_detected(self):
        """shift mismatch should be detected."""
        legacy = [{"bp_idx1": 1, "bp_idx2": 2,
                   "params": {"Shift": 0.1, "Slide": -1.5, "Rise": 3.3,
                              "Tilt": -2.0, "Roll": 8.5, "Twist": 31.0}}]
        modern = [{"bp_idx1": 1, "bp_idx2": 2,
                   "shift": 0.2, "slide": -1.5, "rise": 3.3,
                   "tilt": -2.0, "roll": 8.5, "twist": 31.0}]
        passed, errors = compare_step_params(legacy, modern)
        assert not passed, "Should fail on shift mismatch"
        assert any("shift" in e for e in errors), f"Should mention shift: {errors}"

    def test_slide_mismatch_detected(self):
        """slide mismatch should be detected."""
        legacy = [{"bp_idx1": 1, "bp_idx2": 2,
                   "params": {"Shift": 0.1, "Slide": -1.5, "Rise": 3.3,
                              "Tilt": -2.0, "Roll": 8.5, "Twist": 31.0}}]
        modern = [{"bp_idx1": 1, "bp_idx2": 2,
                   "shift": 0.1, "slide": -1.6, "rise": 3.3,
                   "tilt": -2.0, "roll": 8.5, "twist": 31.0}]
        passed, errors = compare_step_params(legacy, modern)
        assert not passed, "Should fail on slide mismatch"
        assert any("slide" in e for e in errors), f"Should mention slide: {errors}"

    def test_rise_mismatch_detected(self):
        """rise mismatch should be detected."""
        legacy = [{"bp_idx1": 1, "bp_idx2": 2,
                   "params": {"Shift": 0.1, "Slide": -1.5, "Rise": 3.3,
                              "Tilt": -2.0, "Roll": 8.5, "Twist": 31.0}}]
        modern = [{"bp_idx1": 1, "bp_idx2": 2,
                   "shift": 0.1, "slide": -1.5, "rise": 3.4,
                   "tilt": -2.0, "roll": 8.5, "twist": 31.0}]
        passed, errors = compare_step_params(legacy, modern)
        assert not passed, "Should fail on rise mismatch"
        assert any("rise" in e for e in errors), f"Should mention rise: {errors}"

    def test_tilt_mismatch_detected(self):
        """tilt mismatch should be detected."""
        legacy = [{"bp_idx1": 1, "bp_idx2": 2,
                   "params": {"Shift": 0.1, "Slide": -1.5, "Rise": 3.3,
                              "Tilt": -2.0, "Roll": 8.5, "Twist": 31.0}}]
        modern = [{"bp_idx1": 1, "bp_idx2": 2,
                   "shift": 0.1, "slide": -1.5, "rise": 3.3,
                   "tilt": -3.0, "roll": 8.5, "twist": 31.0}]
        passed, errors = compare_step_params(legacy, modern)
        assert not passed, "Should fail on tilt mismatch"
        assert any("tilt" in e for e in errors), f"Should mention tilt: {errors}"

    def test_roll_mismatch_detected(self):
        """roll mismatch should be detected."""
        legacy = [{"bp_idx1": 1, "bp_idx2": 2,
                   "params": {"Shift": 0.1, "Slide": -1.5, "Rise": 3.3,
                              "Tilt": -2.0, "Roll": 8.5, "Twist": 31.0}}]
        modern = [{"bp_idx1": 1, "bp_idx2": 2,
                   "shift": 0.1, "slide": -1.5, "rise": 3.3,
                   "tilt": -2.0, "roll": 9.5, "twist": 31.0}]
        passed, errors = compare_step_params(legacy, modern)
        assert not passed, "Should fail on roll mismatch"
        assert any("roll" in e for e in errors), f"Should mention roll: {errors}"

    def test_twist_mismatch_detected(self):
        """twist mismatch should be detected."""
        legacy = [{"bp_idx1": 1, "bp_idx2": 2,
                   "params": {"Shift": 0.1, "Slide": -1.5, "Rise": 3.3,
                              "Tilt": -2.0, "Roll": 8.5, "Twist": 31.0}}]
        modern = [{"bp_idx1": 1, "bp_idx2": 2,
                   "shift": 0.1, "slide": -1.5, "rise": 3.3,
                   "tilt": -2.0, "roll": 8.5, "twist": 32.0}]
        passed, errors = compare_step_params(legacy, modern)
        assert not passed, "Should fail on twist mismatch"
        assert any("twist" in e for e in errors), f"Should mention twist: {errors}"


class TestHelicalParamsComparator:
    """Test Stage 12: helical_params comparator.

    Key: (bp_idx1, bp_idx2) - step indices
    Legacy format: params array [x_disp, y_disp, rise, incl, tip, twist]
    Modern format: individual fields {"x_displacement": ..., "y_displacement": ...}
    """

    def test_matching_params_pass(self):
        """Identical helical params should pass."""
        legacy = [{"bp_idx1": 1, "bp_idx2": 2,
                   "params": [-0.5, 0.8, 2.8, 12.0, -4.5, 32.0]}]
        modern = [{"bp_idx1": 1, "bp_idx2": 2,
                   "x_displacement": -0.5, "y_displacement": 0.8,
                   "rise": 2.8, "inclination": 12.0, "tip": -4.5, "twist": 32.0}]
        passed, errors = compare_helical_params(legacy, modern)
        assert passed, f"Should pass: {errors}"

    def test_x_displacement_mismatch_detected(self):
        """x_displacement mismatch should be detected."""
        legacy = [{"bp_idx1": 1, "bp_idx2": 2,
                   "params": [-0.5, 0.8, 2.8, 12.0, -4.5, 32.0]}]
        modern = [{"bp_idx1": 1, "bp_idx2": 2,
                   "x_displacement": -0.6, "y_displacement": 0.8,
                   "rise": 2.8, "inclination": 12.0, "tip": -4.5, "twist": 32.0}]
        passed, errors = compare_helical_params(legacy, modern)
        assert not passed, "Should fail on x_displacement mismatch"
        assert any("x_displacement" in e for e in errors), f"Should mention x_displacement: {errors}"

    def test_y_displacement_mismatch_detected(self):
        """y_displacement mismatch should be detected."""
        legacy = [{"bp_idx1": 1, "bp_idx2": 2,
                   "params": [-0.5, 0.8, 2.8, 12.0, -4.5, 32.0]}]
        modern = [{"bp_idx1": 1, "bp_idx2": 2,
                   "x_displacement": -0.5, "y_displacement": 0.9,
                   "rise": 2.8, "inclination": 12.0, "tip": -4.5, "twist": 32.0}]
        passed, errors = compare_helical_params(legacy, modern)
        assert not passed, "Should fail on y_displacement mismatch"
        assert any("y_displacement" in e for e in errors), f"Should mention y_displacement: {errors}"

    def test_inclination_mismatch_detected(self):
        """inclination mismatch should be detected."""
        legacy = [{"bp_idx1": 1, "bp_idx2": 2,
                   "params": [-0.5, 0.8, 2.8, 12.0, -4.5, 32.0]}]
        modern = [{"bp_idx1": 1, "bp_idx2": 2,
                   "x_displacement": -0.5, "y_displacement": 0.8,
                   "rise": 2.8, "inclination": 13.0, "tip": -4.5, "twist": 32.0}]
        passed, errors = compare_helical_params(legacy, modern)
        assert not passed, "Should fail on inclination mismatch"
        assert any("inclination" in e for e in errors), f"Should mention inclination: {errors}"

    def test_tip_mismatch_detected(self):
        """tip mismatch should be detected."""
        legacy = [{"bp_idx1": 1, "bp_idx2": 2,
                   "params": [-0.5, 0.8, 2.8, 12.0, -4.5, 32.0]}]
        modern = [{"bp_idx1": 1, "bp_idx2": 2,
                   "x_displacement": -0.5, "y_displacement": 0.8,
                   "rise": 2.8, "inclination": 12.0, "tip": -5.5, "twist": 32.0}]
        passed, errors = compare_helical_params(legacy, modern)
        assert not passed, "Should fail on tip mismatch"
        assert any("tip" in e for e in errors), f"Should mention tip: {errors}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


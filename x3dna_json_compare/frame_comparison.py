"""
Frame calculation comparison utilities.

Provides functions to compare frame calculations between legacy and modern JSON outputs.
"""

import os
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
from .models import FrameComparison, FrameRecord, FrameMismatch, AtomLineInfo
from .pdb_utils import PdbFileReader


def trim_string(value: Any) -> Any:
    """Trim whitespace from string values, leave other types unchanged."""
    if isinstance(value, str):
        return value.strip()
    return value


def compare_frames(legacy_records: List[Dict], modern_records: List[Dict],
                   pdb_file: Path, pdb_reader: PdbFileReader = None,
                   legacy_atoms: List[Dict] = None, tolerance: float = 1e-6) -> FrameComparison:
    """
    Compare frame calculations between legacy and modern JSON.
    
    Uses legacy_residue_idx for direct matching when available, otherwise falls back
    to key-based matching by (chain_id, residue_seq, insertion).
    
    Args:
        legacy_records: List of legacy calculation records
        modern_records: List of modern calculation records
        pdb_file: Path to PDB file for line lookups
        pdb_reader: Optional pre-initialized PdbFileReader
        legacy_atoms: Optional list of legacy atoms for atom_idx lookup
        
    Returns:
        FrameComparison result
    """
    result = FrameComparison()
    
    # Build legacy atom lookup map by (chain_id, residue_seq, insertion, atom_name) -> atom_idx
    legacy_atom_idx_map = {}
    if legacy_atoms:
        for atom in legacy_atoms:
            key = (
                atom.get('chain_id', ''),
                atom.get('residue_seq', 0),
                atom.get('insertion', ' '),
                atom.get('atom_name', '')
            )
            legacy_atom_idx_map[key] = atom.get('atom_idx')
    
    # Use provided reader or create new one (skip if pdb_file is None)
    if pdb_reader is None and pdb_file is not None:
        try:
            pdb_reader = PdbFileReader(pdb_file)
        except Exception:
            pdb_reader = None
    
    # Deduplicate legacy records by residue_idx (legacy has duplicate record bug)
    # Keep first occurrence of each residue_idx
    legacy_deduped = []
    seen_residue_idx = set()
    for rec in legacy_records:
        residue_idx = rec.get('residue_idx')
        if residue_idx and residue_idx in seen_residue_idx:
            # Skip duplicate - legacy generates duplicate base_frame_calc and frame_calc records
            continue
        if residue_idx:
            seen_residue_idx.add(residue_idx)
        legacy_deduped.append(rec)
    
    # Separate legacy records by type and by residue_idx
    legacy_base_frame = {}  # base_frame_calc records
    legacy_ls_fitting = {}  # ls_fitting records
    legacy_frame_calc = {}  # frame_calc records
    legacy_by_residue_idx = {}  # Map residue_idx -> key for direct matching
    
    for rec in legacy_deduped:
        chain_id = rec.get('chain_id', '')
        residue_seq = rec.get('residue_seq', 0)
        insertion = rec.get('insertion', ' ')
        key = (chain_id, residue_seq, insertion)
        residue_idx = rec.get('residue_idx')
        
        if residue_idx:
            legacy_by_residue_idx[residue_idx] = key
        
        rec_type = rec.get('type')
        if rec_type == 'base_frame_calc':
            legacy_base_frame[key] = rec.copy() if rec else rec
        elif rec_type == 'ls_fitting':
            legacy_ls_fitting[key] = rec.copy() if rec else rec
        elif rec_type == 'frame_calc':
            legacy_frame_calc[key] = rec.copy() if rec else rec
    
    # Build modern maps by type and by legacy_residue_idx
    modern_base_frame = {}
    modern_ls_fitting = {}
    modern_frame_calc = {}
    modern_by_legacy_residue_idx_bf = {}  # For base_frame_calc
    modern_by_legacy_residue_idx_ls = {}  # For ls_fitting
    modern_by_legacy_residue_idx_fc = {}  # For frame_calc
    
    for rec in modern_records:
        chain_id = rec.get('chain_id', '')
        residue_seq = rec.get('residue_seq', 0)
        insertion = rec.get('insertion', ' ')
        key = (chain_id, residue_seq, insertion)
        
        rec_type = rec.get('type')
        if rec_type == 'base_frame_calc':
            if key not in modern_base_frame:
                modern_base_frame[key] = rec.copy() if rec else rec
            else:
                # Merge if multiple records - create a new dict to avoid mutation issues
                merged = modern_base_frame[key].copy()
                merged.update(rec)
                modern_base_frame[key] = merged
            
            # If modern record has legacy_residue_idx, use it for direct matching
            legacy_residue_idx = rec.get('legacy_residue_idx')
            if legacy_residue_idx:
                modern_by_legacy_residue_idx_bf[legacy_residue_idx] = modern_base_frame[key]
                
        elif rec_type == 'ls_fitting':
            if key not in modern_ls_fitting:
                modern_ls_fitting[key] = rec.copy() if rec else rec
            else:
                # Merge if multiple records - create a new dict to avoid mutation issues
                merged = modern_ls_fitting[key].copy()
                merged.update(rec)
                modern_ls_fitting[key] = merged
            
            # If modern record has legacy_residue_idx, use it for direct matching
            legacy_residue_idx = rec.get('legacy_residue_idx')
            if legacy_residue_idx:
                modern_by_legacy_residue_idx_ls[legacy_residue_idx] = modern_ls_fitting[key]
                
        elif rec_type == 'frame_calc':
            if key not in modern_frame_calc:
                modern_frame_calc[key] = rec.copy() if rec else rec
            else:
                # Merge if multiple records - create a new dict to avoid mutation issues
                merged = modern_frame_calc[key].copy()
                merged.update(rec)
                modern_frame_calc[key] = merged
            
            # If modern record has legacy_residue_idx, use it for direct matching
            legacy_residue_idx = rec.get('legacy_residue_idx')
            if legacy_residue_idx:
                modern_by_legacy_residue_idx_fc[legacy_residue_idx] = modern_frame_calc[key]
    
    # Count totals (use union of keys for all three types)
    all_legacy_keys = set(legacy_base_frame.keys()) | set(legacy_ls_fitting.keys()) | set(legacy_frame_calc.keys())
    all_modern_keys = set(modern_base_frame.keys()) | set(modern_ls_fitting.keys()) | set(modern_frame_calc.keys())
    result.total_legacy = len(all_legacy_keys)
    result.total_modern = len(all_modern_keys)
    
    # Helper function to find keys using legacy_residue_idx matching
    def get_matched_keys(legacy_map, modern_by_residue_idx_map):
        matched_keys = set()
        for legacy_residue_idx, mod_rec in modern_by_residue_idx_map.items():
            if legacy_residue_idx in legacy_by_residue_idx:
                key = legacy_by_residue_idx[legacy_residue_idx]
                if key in legacy_map:
                    matched_keys.add(key)
        return matched_keys
    
    # Match residues using legacy_residue_idx when available (PRIMARY matching method)
    matched_bf_keys = get_matched_keys(legacy_base_frame, modern_by_legacy_residue_idx_bf)
    matched_ls_keys = get_matched_keys(legacy_ls_fitting, modern_by_legacy_residue_idx_ls)
    matched_fc_keys = get_matched_keys(legacy_frame_calc, modern_by_legacy_residue_idx_fc)
    
    # Match by key (FALLBACK matching method - only use if residue_idx matching didn't work)
    # Exclude keys that were already matched by residue_idx to avoid duplicates
    common_bf_keys = (set(legacy_base_frame.keys()) & set(modern_base_frame.keys())) - matched_bf_keys
    common_ls_keys = (set(legacy_ls_fitting.keys()) & set(modern_ls_fitting.keys())) - matched_ls_keys
    common_fc_keys = (set(legacy_frame_calc.keys()) & set(modern_frame_calc.keys())) - matched_fc_keys
    
    # All keys to compare for each type (prioritize residue_idx matches)
    bf_keys_to_compare = matched_bf_keys | common_bf_keys
    ls_keys_to_compare = matched_ls_keys | common_ls_keys
    fc_keys_to_compare = matched_fc_keys | common_fc_keys
    
    # Find missing base_frame_calc residues (in legacy but not in modern)
    missing_bf_keys = set(legacy_base_frame.keys()) - bf_keys_to_compare
    for key in missing_bf_keys:
        bf_rec = legacy_base_frame[key]
        result.missing_residues.append(FrameRecord(
            residue_name=bf_rec.get('residue_name', ''),
            chain_id=bf_rec.get('chain_id', ''),
            residue_seq=bf_rec.get('residue_seq', 0),
            insertion=bf_rec.get('insertion', ' '),
            base_type=bf_rec.get('base_type', '?'),
            rms_fit=bf_rec.get('rms_fit', 0.0),
            num_matched_atoms=bf_rec.get('num_matched_atoms', 0),
            matched_atoms=bf_rec.get('matched_atoms', []),
            template_file=bf_rec.get('template_file', bf_rec.get('standard_template')),
            residue_idx=bf_rec.get('residue_idx'),
            legacy_residue_idx=bf_rec.get('residue_idx'),
        ))
    
    # Find missing ls_fitting residues (in legacy but not in modern)
    missing_ls_keys = set(legacy_ls_fitting.keys()) - ls_keys_to_compare
    for key in missing_ls_keys:
        ls_rec = legacy_ls_fitting[key]
        # Only add if not already added as missing base_frame_calc
        if key not in missing_bf_keys:
            result.missing_residues.append(FrameRecord(
                residue_name=ls_rec.get('residue_name', ''),
                chain_id=ls_rec.get('chain_id', ''),
                residue_seq=ls_rec.get('residue_seq', 0),
                insertion=ls_rec.get('insertion', ' '),
                base_type=ls_rec.get('base_type', '?'),
                rms_fit=ls_rec.get('rms_fit', 0.0),
                num_matched_atoms=ls_rec.get('num_points', 0),
                matched_atoms=[],
                template_file=ls_rec.get('template_file', ''),
                residue_idx=ls_rec.get('residue_idx'),
                legacy_residue_idx=ls_rec.get('residue_idx'),
            ))
    
    # Find missing frame_calc residues (in legacy but not in modern)
    missing_fc_keys = set(legacy_frame_calc.keys()) - fc_keys_to_compare
    for key in missing_fc_keys:
        fc_rec = legacy_frame_calc[key]
        # Only add if not already added as missing base_frame_calc or ls_fitting
        if key not in missing_bf_keys and key not in missing_ls_keys:
            result.missing_residues.append(FrameRecord(
                residue_name=fc_rec.get('residue_name', ''),
                chain_id=fc_rec.get('chain_id', ''),
                residue_seq=fc_rec.get('residue_seq', 0),
                insertion=fc_rec.get('insertion', ' '),
                base_type=fc_rec.get('base_type', '?'),
                rms_fit=fc_rec.get('rms_fit', 0.0),
                num_matched_atoms=fc_rec.get('num_matched_atoms', 0),
                matched_atoms=[],
                template_file=fc_rec.get('template_file', ''),
                residue_idx=fc_rec.get('residue_idx'),
                legacy_residue_idx=fc_rec.get('residue_idx'),
            ))
    
    # Compare base_frame_calc records
    for key in bf_keys_to_compare:
        leg_rec = legacy_base_frame.get(key)
        mod_rec = modern_base_frame.get(key)
        
        if not leg_rec or not mod_rec:
            continue
        
        mismatches = {}
        
        # CRITICAL: Compare by residue_idx - this is the primary matching key
        leg_residue_idx = leg_rec.get('residue_idx')
        mod_residue_idx = mod_rec.get('residue_idx')
        mod_legacy_residue_idx = mod_rec.get('legacy_residue_idx')
        
        # Primary check: legacy residue_idx should match modern legacy_residue_idx
        if leg_residue_idx is not None and mod_legacy_residue_idx is not None:
            if leg_residue_idx != mod_legacy_residue_idx:
                mismatches['residue_idx_mismatch'] = {
                    'legacy_residue_idx': leg_residue_idx,
                    'modern_legacy_residue_idx': mod_legacy_residue_idx,
                    'note': 'CRITICAL: residue_idx mismatch - records may not be for same residue'
                }
        
        # Secondary check: modern residue_idx vs legacy residue_idx (may differ if 0-based vs 1-based)
        if leg_residue_idx is not None and mod_residue_idx is not None and leg_residue_idx != mod_residue_idx:
            # Only warn if legacy_residue_idx also doesn't match (otherwise it's just 0-based vs 1-based)
            if mod_legacy_residue_idx is None or leg_residue_idx != mod_legacy_residue_idx:
                mismatches['residue_idx'] = {'legacy': leg_residue_idx, 'modern': mod_residue_idx}

        # Compare base_type (trim strings for consistent comparison)
        leg_base_type = trim_string(leg_rec.get('base_type', ''))
        mod_base_type = trim_string(mod_rec.get('base_type', ''))
        if leg_base_type != mod_base_type:
            mismatches['base_type'] = {'legacy': leg_base_type, 'modern': mod_base_type}

        # Compare residue_name (trim strings for consistent comparison)
        leg_resname = trim_string(leg_rec.get('residue_name', ''))
        mod_resname = trim_string(mod_rec.get('residue_name', ''))
        if leg_resname != mod_resname:
            mismatches['residue_name'] = {'legacy': leg_resname, 'modern': mod_resname}

        # Compare RMS
        leg_rms = leg_rec.get('rms_fit', 0.0)
        mod_rms = mod_rec.get('rms_fit', 0.0)
        if abs(leg_rms - mod_rms) > 0.001:
            mismatches['rms'] = {'legacy': leg_rms, 'modern': mod_rms, 'diff': abs(leg_rms - mod_rms)}

        # Compare num_matched_atoms
        leg_num = leg_rec.get('num_matched_atoms', 0)
        mod_num = mod_rec.get('num_matched_atoms', 0)
        if leg_num != mod_num:
            mismatches['num_matched_atoms'] = {'legacy': leg_num, 'modern': mod_num}

        # Compare matched_atoms
        leg_atoms = leg_rec.get('matched_atoms', [])
        mod_atoms = mod_rec.get('matched_atoms', [])
        if set(leg_atoms) != set(mod_atoms):
            mismatches['matched_atoms'] = {
                'legacy': leg_atoms,
                'modern': mod_atoms,
                'only_legacy': list(set(leg_atoms) - set(mod_atoms)),
                'only_modern': list(set(mod_atoms) - set(leg_atoms))
            }
        
        # Compare standard_template (compare only filename, not full path)
        leg_template = leg_rec.get('template_file', leg_rec.get('standard_template', ''))
        mod_template = mod_rec.get('template_file', mod_rec.get('standard_template', ''))
        # Normalize to just filename for comparison (paths differ between legacy/modern)
        leg_template_name = os.path.basename(leg_template) if leg_template else ''
        mod_template_name = os.path.basename(mod_template) if mod_template else ''
        if leg_template_name != mod_template_name:
            mismatches['template_file'] = {'legacy': leg_template, 'modern': mod_template}
        
        if mismatches:
            # Get PDB lines for matched atoms
            atom_line_info = []
            atoms_to_find = set(leg_atoms) | set(mod_atoms) if (leg_atoms or mod_atoms) else set()
            
            if atoms_to_find and pdb_reader:
                chain_id, residue_seq, insertion = key
                atom_lines = pdb_reader.get_atom_lines_by_names(
                    chain_id, residue_seq, insertion, list(atoms_to_find)
                )
                for atom_name, (line_num, line) in atom_lines.items():
                    atom_key = (chain_id, residue_seq, insertion, atom_name)
                    legacy_idx = legacy_atom_idx_map.get(atom_key)
                    atom_line_info.append(AtomLineInfo(
                        atom_name=atom_name,
                        line_number=line_num,
                        pdb_line=line,
                        legacy_atom_idx=legacy_idx
                    ))
            
            result.mismatched_calculations.append(FrameMismatch(
                residue_key=key,
                legacy_record=leg_rec,
                modern_record=mod_rec,
                mismatches={'base_frame_calc': mismatches},
                atom_pdb_lines=atom_line_info,
                legacy_matched_atoms=leg_atoms,
                modern_matched_atoms=mod_atoms
            ))
    
    # Compare ls_fitting records
    for key in ls_keys_to_compare:
        leg_rec = legacy_ls_fitting.get(key)
        mod_rec = modern_ls_fitting.get(key)
        
        if not leg_rec or not mod_rec:
            continue
        
        mismatches = {}
        
        # CRITICAL: Compare by residue_idx - this is the primary matching key
        leg_residue_idx = leg_rec.get('residue_idx')
        mod_residue_idx = mod_rec.get('residue_idx')
        mod_legacy_residue_idx = mod_rec.get('legacy_residue_idx')
        
        # Primary check: legacy residue_idx should match modern legacy_residue_idx
        if leg_residue_idx is not None and mod_legacy_residue_idx is not None:
            if leg_residue_idx != mod_legacy_residue_idx:
                mismatches['residue_idx_mismatch'] = {
                    'legacy_residue_idx': leg_residue_idx,
                    'modern_legacy_residue_idx': mod_legacy_residue_idx,
                    'note': 'CRITICAL: residue_idx mismatch - records may not be for same residue'
                }
        
        # Secondary check: modern residue_idx vs legacy residue_idx (may differ if 0-based vs 1-based)
        if leg_residue_idx is not None and mod_residue_idx is not None and leg_residue_idx != mod_residue_idx:
            # Only warn if legacy_residue_idx also doesn't match (otherwise it's just 0-based vs 1-based)
            if mod_legacy_residue_idx is None or leg_residue_idx != mod_legacy_residue_idx:
                mismatches['residue_idx'] = {'legacy': leg_residue_idx, 'modern': mod_residue_idx}
        
        # Compare residue_name (trim strings for consistent comparison)
        leg_resname = trim_string(leg_rec.get('residue_name', ''))
        mod_resname = trim_string(mod_rec.get('residue_name', ''))
        if leg_resname != mod_resname:
            mismatches['residue_name'] = {'legacy': leg_resname, 'modern': mod_resname}

        # Compare RMS
        leg_rms = leg_rec.get('rms_fit', 0.0)
        mod_rms = mod_rec.get('rms_fit', 0.0)
        if abs(leg_rms - mod_rms) > 0.001:
            mismatches['rms'] = {'legacy': leg_rms, 'modern': mod_rms, 'diff': abs(leg_rms - mod_rms)}

        # Compare num_points
        leg_num = leg_rec.get('num_points', 0)
        mod_num = mod_rec.get('num_points', 0)
        if leg_num != mod_num:
            mismatches['num_points'] = {'legacy': leg_num, 'modern': mod_num}
        
        # Compare rotation_matrix (3x3 matrix)
        leg_rot = leg_rec.get('rotation_matrix', [])
        mod_rot = mod_rec.get('rotation_matrix', [])
        if len(leg_rot) == 3 and len(mod_rot) == 3:
            rot_diff = False
            max_rot_diff = 0.0
            for i in range(3):
                if len(leg_rot[i]) == 3 and len(mod_rot[i]) == 3:
                    for j in range(3):
                        diff = abs(leg_rot[i][j] - mod_rot[i][j])
                        max_rot_diff = max(max_rot_diff, diff)
                        if diff > tolerance:
                            rot_diff = True
                else:
                    rot_diff = True
                    break
            if rot_diff:
                mismatches['rotation_matrix'] = {
                    'legacy': leg_rot,
                    'modern': mod_rot,
                    'max_diff': max_rot_diff
                }
        elif leg_rot != mod_rot:
            mismatches['rotation_matrix'] = {'legacy': leg_rot, 'modern': mod_rot}
        
        # Compare translation (3D vector)
        leg_trans = leg_rec.get('translation', [])
        mod_trans = mod_rec.get('translation', [])
        if len(leg_trans) == 3 and len(mod_trans) == 3:
            max_trans_diff = max(abs(leg_trans[i] - mod_trans[i]) for i in range(3))
            if max_trans_diff > tolerance:
                mismatches['translation'] = {
                    'legacy': leg_trans,
                    'modern': mod_trans,
                    'max_diff': max_trans_diff
                }
        elif leg_trans != mod_trans:
            mismatches['translation'] = {'legacy': leg_trans, 'modern': mod_trans}
        
        if mismatches:
            # Get matched atoms from base_frame_calc if available
            leg_atoms = []
            mod_atoms = []
            if key in legacy_base_frame:
                leg_atoms = legacy_base_frame[key].get('matched_atoms', [])
            if key in modern_base_frame:
                mod_atoms = modern_base_frame[key].get('matched_atoms', [])
            
            # Get PDB lines for matched atoms
            atom_line_info = []
            atoms_to_find = set(leg_atoms) | set(mod_atoms) if (leg_atoms or mod_atoms) else set()
            
            if atoms_to_find and pdb_reader:
                chain_id, residue_seq, insertion = key
                atom_lines = pdb_reader.get_atom_lines_by_names(
                    chain_id, residue_seq, insertion, list(atoms_to_find)
                )
                for atom_name, (line_num, line) in atom_lines.items():
                    atom_key = (chain_id, residue_seq, insertion, atom_name)
                    legacy_idx = legacy_atom_idx_map.get(atom_key)
                    atom_line_info.append(AtomLineInfo(
                        atom_name=atom_name,
                        line_number=line_num,
                        pdb_line=line,
                        legacy_atom_idx=legacy_idx
                    ))
            
            result.mismatched_calculations.append(FrameMismatch(
                residue_key=key,
                legacy_record=leg_rec,
                modern_record=mod_rec,
                mismatches={'ls_fitting': mismatches},
                atom_pdb_lines=atom_line_info,
                    legacy_matched_atoms=leg_atoms,
                    modern_matched_atoms=mod_atoms
            ))
    
    # Compare frame_calc records
    for key in fc_keys_to_compare:
        leg_rec = legacy_frame_calc.get(key)
        mod_rec = modern_frame_calc.get(key)
        
        if not leg_rec or not mod_rec:
            continue
        
        mismatches = {}
        tolerance = 1e-6  # Tolerance for coordinate comparisons
        
        # CRITICAL: Compare by residue_idx - this is the primary matching key
        leg_residue_idx = leg_rec.get('residue_idx')
        mod_residue_idx = mod_rec.get('residue_idx')
        mod_legacy_residue_idx = mod_rec.get('legacy_residue_idx')
        
        # Primary check: legacy residue_idx should match modern legacy_residue_idx
        if leg_residue_idx is not None and mod_legacy_residue_idx is not None:
            if leg_residue_idx != mod_legacy_residue_idx:
                mismatches['residue_idx_mismatch'] = {
                    'legacy_residue_idx': leg_residue_idx,
                    'modern_legacy_residue_idx': mod_legacy_residue_idx,
                    'note': 'CRITICAL: residue_idx mismatch - records may not be for same residue'
                }
        
        # Secondary check: modern residue_idx vs legacy residue_idx (may differ if 0-based vs 1-based)
        if leg_residue_idx is not None and mod_residue_idx is not None and leg_residue_idx != mod_residue_idx:
            # Only warn if legacy_residue_idx also doesn't match (otherwise it's just 0-based vs 1-based)
            if mod_legacy_residue_idx is None or leg_residue_idx != mod_legacy_residue_idx:
                mismatches['residue_idx'] = {'legacy': leg_residue_idx, 'modern': mod_residue_idx}
        
        # Compare base_type (trim strings for consistent comparison)
        leg_base_type = trim_string(leg_rec.get('base_type', ''))
        mod_base_type = trim_string(mod_rec.get('base_type', ''))
        if leg_base_type != mod_base_type:
            mismatches['base_type'] = {'legacy': leg_base_type, 'modern': mod_base_type}

        # Compare residue_name (trim strings for consistent comparison)
        leg_resname = trim_string(leg_rec.get('residue_name', ''))
        mod_resname = trim_string(mod_rec.get('residue_name', ''))
        if leg_resname != mod_resname:
            mismatches['residue_name'] = {'legacy': leg_resname, 'modern': mod_resname}

        # Compare RMS
        leg_rms = leg_rec.get('rms_fit', 0.0)
        mod_rms = mod_rec.get('rms_fit', 0.0)
        if abs(leg_rms - mod_rms) > 0.001:
            mismatches['rms'] = {'legacy': leg_rms, 'modern': mod_rms, 'diff': abs(leg_rms - mod_rms)}

        # Compare num_matched_atoms
        leg_num = leg_rec.get('num_matched_atoms', 0)
        mod_num = mod_rec.get('num_matched_atoms', 0)
        if leg_num != mod_num:
            mismatches['num_matched_atoms'] = {'legacy': leg_num, 'modern': mod_num}

        # Compare template_file (compare only filename, not full path)
        leg_template = leg_rec.get('template_file', '')
        mod_template = mod_rec.get('template_file', '')
        # Normalize to just filename for comparison (paths differ between legacy/modern)
        leg_template_name = os.path.basename(leg_template) if leg_template else ''
        mod_template_name = os.path.basename(mod_template) if mod_template else ''
        if leg_template_name != mod_template_name:
            mismatches['template_file'] = {'legacy': leg_template, 'modern': mod_template}
        
        # Compare matched_coordinates
        leg_coords = leg_rec.get('matched_coordinates', [])
        mod_coords = mod_rec.get('matched_coordinates', [])
        
        if len(leg_coords) != len(mod_coords):
            mismatches['matched_coordinates'] = {
                'legacy_count': len(leg_coords),
                'modern_count': len(mod_coords),
                'note': 'array length mismatch'
            }
        else:
            coord_diffs = []
            for i, (leg_coord, mod_coord) in enumerate(zip(leg_coords, mod_coords)):
                coord_mismatches = {}
                
                # Compare atom_idx
                leg_atom_idx = leg_coord.get('atom_idx')
                mod_atom_idx = mod_coord.get('atom_idx')
                if leg_atom_idx != mod_atom_idx:
                    coord_mismatches['atom_idx'] = {'legacy': leg_atom_idx, 'modern': mod_atom_idx}
                
                # Compare std_xyz
                leg_std_xyz = leg_coord.get('std_xyz', [])
                mod_std_xyz = mod_coord.get('std_xyz', [])
                if len(leg_std_xyz) == 3 and len(mod_std_xyz) == 3:
                    max_std_diff = max(abs(leg_std_xyz[j] - mod_std_xyz[j]) for j in range(3))
                    if max_std_diff > tolerance:
                        coord_mismatches['std_xyz'] = {
                            'legacy': leg_std_xyz,
                            'modern': mod_std_xyz,
                            'max_diff': max_std_diff
                        }
                elif leg_std_xyz != mod_std_xyz:
                    coord_mismatches['std_xyz'] = {'legacy': leg_std_xyz, 'modern': mod_std_xyz}
                
                # Compare exp_xyz
                leg_exp_xyz = leg_coord.get('exp_xyz', [])
                mod_exp_xyz = mod_coord.get('exp_xyz', [])
                if len(leg_exp_xyz) == 3 and len(mod_exp_xyz) == 3:
                    max_exp_diff = max(abs(leg_exp_xyz[j] - mod_exp_xyz[j]) for j in range(3))
                    if max_exp_diff > tolerance:
                        coord_mismatches['exp_xyz'] = {
                            'legacy': leg_exp_xyz,
                            'modern': mod_exp_xyz,
                            'max_diff': max_exp_diff
                        }
                elif leg_exp_xyz != mod_exp_xyz:
                    coord_mismatches['exp_xyz'] = {'legacy': leg_exp_xyz, 'modern': mod_exp_xyz}
                
                if coord_mismatches:
                    coord_diffs.append({
                        'index': i,
                        'atom_idx': leg_atom_idx or mod_atom_idx,
                        'mismatches': coord_mismatches
                    })
            
            if coord_diffs:
                mismatches['matched_coordinates'] = {
                    'total_pairs': len(leg_coords),
                    'mismatched_pairs': len(coord_diffs),
                    'differences': coord_diffs[:10]  # Limit to first 10 for reporting
                }
        
        if mismatches:
            # Get matched atoms from base_frame_calc if available for PDB line lookup
            leg_atoms = []
            mod_atoms = []
            if key in legacy_base_frame:
                leg_atoms = legacy_base_frame[key].get('matched_atoms', [])
            if key in modern_base_frame:
                mod_atoms = modern_base_frame[key].get('matched_atoms', [])
            
            # Get PDB lines for matched atoms
            atom_line_info = []
            atoms_to_find = set(leg_atoms) | set(mod_atoms) if (leg_atoms or mod_atoms) else set()
            
            if atoms_to_find and pdb_reader:
                chain_id, residue_seq, insertion = key
                atom_lines = pdb_reader.get_atom_lines_by_names(
                    chain_id, residue_seq, insertion, list(atoms_to_find)
                )
                for atom_name, (line_num, line) in atom_lines.items():
                    atom_key = (chain_id, residue_seq, insertion, atom_name)
                    legacy_idx = legacy_atom_idx_map.get(atom_key)
                    atom_line_info.append(AtomLineInfo(
                        atom_name=atom_name,
                        line_number=line_num,
                        pdb_line=line,
                        legacy_atom_idx=legacy_idx
                    ))
            
            result.mismatched_calculations.append(FrameMismatch(
                residue_key=key,
                legacy_record=leg_rec,
                modern_record=mod_rec,
                mismatches={'frame_calc': mismatches},
                atom_pdb_lines=atom_line_info,
                legacy_matched_atoms=leg_atoms,
                modern_matched_atoms=mod_atoms
            ))
    
    return result


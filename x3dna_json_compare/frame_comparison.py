"""
Frame calculation comparison utilities.

Provides functions to compare frame calculations between legacy and modern JSON outputs.
"""

from typing import Dict, List, Tuple, Optional
from pathlib import Path
from .models import FrameComparison, FrameRecord, FrameMismatch, AtomLineInfo
from .pdb_utils import PdbFileReader


def compare_frames(legacy_records: List[Dict], modern_records: List[Dict],
                   pdb_file: Path, pdb_reader: PdbFileReader = None,
                   legacy_atoms: List[Dict] = None) -> FrameComparison:
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
    
    # Separate legacy records by type and by residue_idx
    legacy_base_frame = {}  # For matched_atoms lookup
    legacy_calc = {}  # For calculations (frame_calc or ls_fitting)
    legacy_by_residue_idx = {}  # Map residue_idx -> key for direct matching
    
    for rec in legacy_records:
        chain_id = rec.get('chain_id', '')
        residue_seq = rec.get('residue_seq', 0)
        insertion = rec.get('insertion', ' ')
        key = (chain_id, residue_seq, insertion)
        residue_idx = rec.get('residue_idx')
        
        if residue_idx:
            legacy_by_residue_idx[residue_idx] = key
        
        if rec.get('type') == 'base_frame_calc':
            legacy_base_frame[key] = rec
        elif rec.get('type') in ['frame_calc', 'ls_fitting']:
            if key not in legacy_calc:
                legacy_calc[key] = rec
            else:
                # Merge multiple calculation types
                legacy_calc[key].update(rec)
    
    # Build modern map by key and by legacy_residue_idx
    modern_map = {}
    modern_by_legacy_residue_idx = {}
    
    for rec in modern_records:
        chain_id = rec.get('chain_id', '')
        residue_seq = rec.get('residue_seq', 0)
        insertion = rec.get('insertion', ' ')
        key = (chain_id, residue_seq, insertion)
        
        if key not in modern_map:
            modern_map[key] = rec
        else:
            modern_map[key].update(rec)
        
        # If modern record has legacy_residue_idx, use it for direct matching
        legacy_residue_idx = rec.get('legacy_residue_idx')
        if legacy_residue_idx:
            modern_by_legacy_residue_idx[legacy_residue_idx] = rec
    
    result.total_legacy = len(legacy_calc)
    result.total_modern = len(modern_map)
    
    # Match residues using legacy_residue_idx when available
    matched_by_residue_idx = set()
    matched_by_key = set()
    
    # First, match by legacy_residue_idx (direct matching)
    for legacy_residue_idx, mod_rec in modern_by_legacy_residue_idx.items():
        if legacy_residue_idx in legacy_by_residue_idx:
            key = legacy_by_residue_idx[legacy_residue_idx]
            if key in legacy_calc:
                matched_by_residue_idx.add(legacy_residue_idx)
                matched_by_key.add(key)
    
    # Then, match remaining residues by key
    legacy_keys = set(legacy_calc.keys())
    modern_keys = set(modern_map.keys())
    common_keys = legacy_keys & modern_keys
    
    # Find missing residues (in legacy but not matched in modern)
    for key in legacy_keys:
        if key not in matched_by_key and key not in common_keys:
            calc_rec = legacy_calc[key].copy()
            
            # Get matched_atoms and residue info from base_frame_calc if available
            if key in legacy_base_frame:
                bfc_rec = legacy_base_frame[key]
                calc_rec['matched_atoms'] = bfc_rec.get('matched_atoms', [])
                if 'residue_name' not in calc_rec:
                    calc_rec['residue_name'] = bfc_rec.get('residue_name', '')
                if 'chain_id' not in calc_rec:
                    calc_rec['chain_id'] = bfc_rec.get('chain_id', '')
                if 'residue_seq' not in calc_rec:
                    calc_rec['residue_seq'] = bfc_rec.get('residue_seq', 0)
                if 'insertion' not in calc_rec:
                    calc_rec['insertion'] = bfc_rec.get('insertion', ' ')
            
            result.missing_residues.append(FrameRecord(
                residue_name=calc_rec.get('residue_name', ''),
                chain_id=calc_rec.get('chain_id', ''),
                residue_seq=calc_rec.get('residue_seq', 0),
                insertion=calc_rec.get('insertion', ' '),
                base_type=calc_rec.get('base_type', '?'),
                rms_fit=calc_rec.get('rms_fit', 0.0),
                num_matched_atoms=calc_rec.get('num_matched_atoms', calc_rec.get('num_points', 0)),
                matched_atoms=calc_rec.get('matched_atoms', []),
                template_file=calc_rec.get('template_file', calc_rec.get('standard_template')),
                residue_idx=calc_rec.get('residue_idx'),
                legacy_residue_idx=calc_rec.get('residue_idx'),  # Legacy uses residue_idx directly
            ))
    
    # Compare matching residues (matched by legacy_residue_idx or by key)
    residues_to_compare = matched_by_key | common_keys
    for key in residues_to_compare:
        leg_rec = legacy_calc[key]
        mod_rec = modern_map[key]
        
        mismatches = {}
        
        # Compare RMS
        leg_rms = leg_rec.get('rms_fit', 0.0)
        mod_rms = mod_rec.get('rms_fit', 0.0)
        if abs(leg_rms - mod_rms) > 0.001:
            mismatches['rms'] = {'legacy': leg_rms, 'modern': mod_rms, 'diff': abs(leg_rms - mod_rms)}
        
        # Compare num_matched
        leg_num = leg_rec.get('num_matched_atoms', leg_rec.get('num_points', 0))
        mod_num = mod_rec.get('num_matched_atoms', mod_rec.get('num_points', 0))
        if leg_num != mod_num:
            mismatches['num_matched'] = {'legacy': leg_num, 'modern': mod_num}
        
        # Get matched atoms from base_frame_calc if available
        leg_atoms = []
        if key in legacy_base_frame:
            leg_atoms = legacy_base_frame[key].get('matched_atoms', [])
        mod_atoms = mod_rec.get('matched_atoms', [])
        
        if leg_atoms and mod_atoms and set(leg_atoms) != set(mod_atoms):
            mismatches['matched_atoms'] = {
                'legacy': leg_atoms,
                'modern': mod_atoms,
                'only_legacy': list(set(leg_atoms) - set(mod_atoms)),
                'only_modern': list(set(mod_atoms) - set(leg_atoms))
            }
        
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
                    # Look up legacy atom index
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
                mismatches=mismatches,
                atom_pdb_lines=atom_line_info,
                legacy_matched_atoms=leg_atoms,
                modern_matched_atoms=mod_atoms
            ))
    
    return result


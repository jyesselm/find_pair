"""
Frame calculation comparison utilities.

Provides functions to compare frame calculations between legacy and modern JSON outputs.
"""

from typing import Dict, List, Tuple, Optional
from pathlib import Path
from .models import FrameComparison, FrameRecord, FrameMismatch, AtomLineInfo
from .pdb_utils import PdbFileReader


def compare_frames(legacy_records: List[Dict], modern_records: List[Dict],
                   pdb_file: Path, pdb_reader: PdbFileReader = None) -> FrameComparison:
    """
    Compare frame calculations between legacy and modern JSON.
    
    Args:
        legacy_records: List of legacy calculation records
        modern_records: List of modern calculation records
        pdb_file: Path to PDB file for line lookups
        pdb_reader: Optional pre-initialized PdbFileReader
        
    Returns:
        FrameComparison result
    """
    result = FrameComparison()
    
    # Use provided reader or create new one
    if pdb_reader is None:
        try:
            pdb_reader = PdbFileReader(pdb_file)
        except Exception:
            pdb_reader = None
    
    # Separate legacy records by type
    legacy_base_frame = {}  # For matched_atoms lookup
    legacy_calc = {}  # For calculations (frame_calc or ls_fitting)
    
    for rec in legacy_records:
        chain_id = rec.get('chain_id', '')
        residue_seq = rec.get('residue_seq', 0)
        insertion = rec.get('insertion', ' ')
        key = (chain_id, residue_seq, insertion)
        
        if rec.get('type') == 'base_frame_calc':
            legacy_base_frame[key] = rec
        elif rec.get('type') in ['frame_calc', 'ls_fitting']:
            if key not in legacy_calc:
                legacy_calc[key] = rec
            else:
                # Merge multiple calculation types
                legacy_calc[key].update(rec)
    
    # Build modern map
    modern_map = {}
    for rec in modern_records:
        chain_id = rec.get('chain_id', '')
        residue_seq = rec.get('residue_seq', 0)
        insertion = rec.get('insertion', ' ')
        key = (chain_id, residue_seq, insertion)
        if key not in modern_map:
            modern_map[key] = rec
        else:
            modern_map[key].update(rec)
    
    result.total_legacy = len(legacy_calc)
    result.total_modern = len(modern_map)
    
    # Find missing residues
    for key in set(legacy_calc.keys()) - set(modern_map.keys()):
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
        ))
    
    # Compare matching residues
    for key in set(legacy_calc.keys()) & set(modern_map.keys()):
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
                    atom_line_info.append(AtomLineInfo(
                        atom_name=atom_name,
                        line_number=line_num,
                        pdb_line=line
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


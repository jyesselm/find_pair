"""
Atom-level comparison utilities.

Provides functions to compare atoms between legacy and modern JSON outputs.
"""

from typing import Dict, List, Tuple
from pathlib import Path
from .models import AtomInfo, AtomComparison, FieldMismatch
from .pdb_utils import PdbFileReader


def create_atom_key(atom: Dict) -> Tuple[str, int, str, str]:
    """
    Create a unique key for an atom.
    
    Args:
        atom: Atom dictionary from JSON
        
    Returns:
        Tuple of (chain_id, residue_seq, insertion, atom_name)
    """
    return (
        atom.get('chain_id', ''),
        atom.get('residue_seq', 0),
        atom.get('insertion', ' '),
        atom.get('atom_name', '')
    )


def compare_atoms(legacy_atoms: List[Dict], modern_atoms: List[Dict],
                 pdb_file: Path, pdb_reader: PdbFileReader = None) -> AtomComparison:
    """
    Compare atoms between legacy and modern JSON.
    
    Args:
        legacy_atoms: List of legacy atom dictionaries
        modern_atoms: List of modern atom dictionaries
        pdb_file: Path to PDB file for line lookups
        pdb_reader: Optional pre-initialized PdbFileReader
        
    Returns:
        AtomComparison result
    """
    result = AtomComparison(
        total_legacy=len(legacy_atoms),
        total_modern=len(modern_atoms),
        count_difference=len(legacy_atoms) != len(modern_atoms)
    )
    
    # Use provided reader or create new one
    if pdb_reader is None:
        try:
            pdb_reader = PdbFileReader(pdb_file)
        except Exception:
            pdb_reader = None
    
    # Create maps by atom key
    legacy_map = {create_atom_key(atom): atom for atom in legacy_atoms}
    modern_map = {create_atom_key(atom): atom for atom in modern_atoms}
    
    legacy_keys = set(legacy_map.keys())
    modern_keys = set(modern_map.keys())
    
    result.common_count = len(legacy_keys & modern_keys)
    
    # Find missing in modern
    for key in legacy_keys - modern_keys:
        leg_atom = legacy_map[key]
        line_number = leg_atom.get('line_number', 0)
        pdb_line = None
        
        if pdb_reader and line_number > 0:
            pdb_line = pdb_reader.get_line(line_number)
        elif pdb_reader:
            # Try to find by atom info
            found = pdb_reader.find_atom_line(
                key[0], key[1], key[2], key[3]
            )
            if found:
                line_number, pdb_line = found
        
        result.missing_in_modern.append(AtomInfo(
            atom_name=leg_atom.get('atom_name', ''),
            chain_id=leg_atom.get('chain_id', ''),
            residue_seq=leg_atom.get('residue_seq', 0),
            insertion=leg_atom.get('insertion', ' '),
            residue_name=leg_atom.get('residue_name', ''),
            line_number=line_number if line_number > 0 else None,
            pdb_line=pdb_line,
            atom_data=leg_atom
        ))
    
    # Find extra in modern
    for key in modern_keys - legacy_keys:
        mod_atom = modern_map[key]
        line_number = mod_atom.get('line_number', 0)
        pdb_line = None
        
        if pdb_reader and line_number > 0:
            pdb_line = pdb_reader.get_line(line_number)
        elif pdb_reader:
            # Try to find by atom info
            found = pdb_reader.find_atom_line(
                key[0], key[1], key[2], key[3]
            )
            if found:
                line_number, pdb_line = found
        
        result.extra_in_modern.append(AtomInfo(
            atom_name=mod_atom.get('atom_name', ''),
            chain_id=mod_atom.get('chain_id', ''),
            residue_seq=mod_atom.get('residue_seq', 0),
            insertion=mod_atom.get('insertion', ' '),
            residue_name=mod_atom.get('residue_name', ''),
            line_number=line_number if line_number > 0 else None,
            pdb_line=pdb_line,
            atom_data=mod_atom
        ))
    
    # Find mismatched fields in common atoms
    fields_to_check = ['xyz', 'residue_name', 'record_type', 'alt_loc']
    for key in legacy_keys & modern_keys:
        leg_atom = legacy_map[key]
        mod_atom = modern_map[key]
        
        mismatches = {}
        for field in fields_to_check:
            leg_val = leg_atom.get(field)
            mod_val = mod_atom.get(field)
            if leg_val != mod_val:
                mismatches[field] = {'legacy': leg_val, 'modern': mod_val}
        
        if mismatches:
            leg_line_num = leg_atom.get('line_number', 0)
            mod_line_num = mod_atom.get('line_number', 0)
            leg_pdb_line = None
            mod_pdb_line = None
            
            if pdb_reader:
                if leg_line_num > 0:
                    leg_pdb_line = pdb_reader.get_line(leg_line_num)
                if mod_line_num > 0:
                    mod_pdb_line = pdb_reader.get_line(mod_line_num)
            
            result.mismatched_fields.append(FieldMismatch(
                field_name=', '.join(mismatches.keys()),
                legacy_value=mismatches,
                modern_value=mismatches,
                atom_key=key,
                legacy_line_number=leg_line_num if leg_line_num > 0 else None,
                modern_line_number=mod_line_num if mod_line_num > 0 else None,
                legacy_pdb_line=leg_pdb_line,
                modern_pdb_line=mod_pdb_line
            ))
    
    return result


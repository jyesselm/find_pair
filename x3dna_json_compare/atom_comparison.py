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
    
    Uses legacy_atom_idx for direct matching when available, otherwise falls back
    to key-based matching by (chain_id, residue_seq, insertion, atom_name).
    
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
    
    # Build legacy map by atom_idx (for direct matching) and by key (for fallback)
    legacy_by_idx = {atom.get('atom_idx'): atom for atom in legacy_atoms if atom.get('atom_idx')}
    legacy_map = {create_atom_key(atom): atom for atom in legacy_atoms}
    
    # Build modern map by legacy_atom_idx (if available) and by key
    modern_by_legacy_idx = {}
    modern_map = {}
    for atom in modern_atoms:
        key = create_atom_key(atom)
        modern_map[key] = atom
        # If modern atom has legacy_atom_idx, use it for direct matching
        legacy_idx = atom.get('legacy_atom_idx')
        if legacy_idx:
            modern_by_legacy_idx[legacy_idx] = atom
    
    # Match atoms using legacy indices when available
    matched_by_idx = set()
    matched_by_key = set()
    
    # First, match by legacy_atom_idx (direct matching)
    for legacy_idx, mod_atom in modern_by_legacy_idx.items():
        if legacy_idx in legacy_by_idx:
            leg_atom = legacy_by_idx[legacy_idx]
            matched_by_idx.add(legacy_idx)
            matched_by_key.add(create_atom_key(leg_atom))
    
    # Then, match remaining atoms by key
    legacy_keys = set(legacy_map.keys())
    modern_keys = set(modern_map.keys())
    common_keys = legacy_keys & modern_keys
    
    # Count common atoms (matched by idx or by key)
    result.common_count = len(matched_by_idx) + len(common_keys - matched_by_key)
    
    # Find missing in modern (atoms in legacy but not matched in modern)
    for legacy_idx, leg_atom in legacy_by_idx.items():
        if legacy_idx not in matched_by_idx:
            # Check if it's also not matched by key
            key = create_atom_key(leg_atom)
            if key not in common_keys:
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
                    legacy_atom_idx=legacy_idx,
                    atom_data=leg_atom
                ))
    
    # Find extra in modern (atoms in modern but not matched in legacy)
    for mod_atom in modern_atoms:
        legacy_idx = mod_atom.get('legacy_atom_idx')
        key = create_atom_key(mod_atom)
        
        # Check if not matched by legacy_idx or by key
        if (not legacy_idx or legacy_idx not in legacy_by_idx) and key not in common_keys:
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
                legacy_atom_idx=legacy_idx or mod_atom.get('atom_idx'),  # Use legacy_idx if available
                atom_data=mod_atom
            ))
    
    # Find mismatched fields in common atoms
    # Check atoms matched by legacy_idx first, then by key
    # Compare all relevant fields
    fields_to_check = [
        'xyz',                    # Coordinates
        'residue_name',           # Residue name
        'record_type',            # PDB record type
        # 'alt_loc' intentionally excluded - modern doesn't store alternate location
        'atom_name',              # Atom name
        'chain_id',               # Chain identifier
        'residue_seq',            # Residue sequence number
        'line_number',            # PDB line number
        # 'pdb_line' intentionally excluded - modern doesn't store raw PDB lines
    ]
    
    # Special comparison: atom_idx (legacy) vs legacy_atom_idx (modern)
    # Also check that modern atom_idx matches legacy atom_idx when legacy_atom_idx is present
    
    # Compare atoms matched by legacy_idx
    for legacy_idx in matched_by_idx:
        leg_atom = legacy_by_idx[legacy_idx]
        mod_atom = modern_by_legacy_idx[legacy_idx]
        key = create_atom_key(leg_atom)
        
        mismatches = {}
        for field in fields_to_check:
            leg_val = leg_atom.get(field)
            mod_val = mod_atom.get(field)
            if leg_val != mod_val:
                mismatches[field] = {'legacy': leg_val, 'modern': mod_val}
        
        # Special check: legacy atom_idx should match modern legacy_atom_idx
        # This is the key comparison - legacy's atom_idx should equal modern's legacy_atom_idx
        leg_atom_idx = leg_atom.get('atom_idx')
        mod_legacy_atom_idx = mod_atom.get('legacy_atom_idx')
        if leg_atom_idx is not None and mod_legacy_atom_idx is not None:
            if leg_atom_idx != mod_legacy_atom_idx:
                mismatches['atom_idx_vs_legacy_atom_idx'] = {
                    'legacy_atom_idx': leg_atom_idx,
                    'modern_legacy_atom_idx': mod_legacy_atom_idx,
                    'note': 'Legacy atom_idx should match modern legacy_atom_idx'
                }
        
        # Note: modern atom_idx can be different from legacy_atom_idx - that's expected
        # The modern code uses its own sequential indexing, while legacy_atom_idx preserves
        # the original legacy index for direct comparison
        
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
    
    # Compare atoms matched by key (but not by legacy_idx)
    for key in common_keys - matched_by_key:
        leg_atom = legacy_map[key]
        mod_atom = modern_map[key]
        
        mismatches = {}
        for field in fields_to_check:
            leg_val = leg_atom.get(field)
            mod_val = mod_atom.get(field)
            if leg_val != mod_val:
                mismatches[field] = {'legacy': leg_val, 'modern': mod_val}
        
        # Special check: legacy atom_idx should match modern legacy_atom_idx (if present)
        # This is the key comparison - legacy's atom_idx should equal modern's legacy_atom_idx
        leg_atom_idx = leg_atom.get('atom_idx')
        mod_legacy_atom_idx = mod_atom.get('legacy_atom_idx')
        if leg_atom_idx is not None and mod_legacy_atom_idx is not None:
            if leg_atom_idx != mod_legacy_atom_idx:
                mismatches['atom_idx_vs_legacy_atom_idx'] = {
                    'legacy_atom_idx': leg_atom_idx,
                    'modern_legacy_atom_idx': mod_legacy_atom_idx,
                    'note': 'Legacy atom_idx should match modern legacy_atom_idx'
                }
        
        # Note: modern atom_idx can be different from legacy_atom_idx - that's expected
        # The modern code uses its own sequential indexing, while legacy_atom_idx preserves
        # the original legacy index for direct comparison
        
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


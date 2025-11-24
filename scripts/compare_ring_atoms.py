#!/usr/bin/env python3
"""
Compare ring atom matching between legacy and modern code.

This script compares which atoms are used for ring atom matching in base frame
calculations between the legacy (org) code and modern code. It uses legacy_atom_idx
to ensure we're comparing the same atoms.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

# Add parent directory to path to import pdb_utils
sys.path.insert(0, str(Path(__file__).parent.parent))
from x3dna_json_compare.pdb_utils import PdbFileReader


@dataclass
class RingAtomMatch:
    """Information about matched ring atoms for a residue."""
    residue_key: Tuple[str, int, str]  # (chain_id, residue_seq, insertion)
    residue_name: str
    base_type: str
    legacy_matched_atoms: List[str]
    modern_matched_atoms: List[str]
    legacy_atom_indices: List[int]  # legacy_atom_idx values
    modern_atom_indices: List[int]  # legacy_atom_idx values
    num_matched_legacy: int
    num_matched_modern: int
    is_match: bool
    differences: List[str]  # Description of differences
    legacy_atom_details: Optional[Dict[int, Dict]] = None  # legacy_atom_idx -> atom details with PDB line
    modern_atom_details: Optional[Dict[int, Dict]] = None  # legacy_atom_idx -> atom details with PDB line
    legacy_atom_details: Dict[int, Dict] = None  # legacy_atom_idx -> atom details with PDB line
    modern_atom_details: Dict[int, Dict] = None  # legacy_atom_idx -> atom details with PDB line


def load_json_file(json_file: Path) -> Optional[Dict]:
    """Load JSON file, handling potential extra data after the main JSON object/array."""
    if not json_file.exists():
        return None
    
    try:
        # Try fast path: just parse the whole file
        with open(json_file, 'r') as f:
            try:
                return json.load(f)
            except json.JSONDecodeError:
                # If that fails, read content and find first complete JSON
                f.seek(0)
                content = f.read()
        
        # Fast method: use JSONDecoder to parse incrementally
        decoder = json.JSONDecoder()
        try:
            obj, idx = decoder.raw_decode(content)
            return obj
        except (json.JSONDecodeError, ValueError):
            # Fallback to character-by-character (slower but more robust)
            brace_count = 0
            bracket_count = 0
            in_string = False
            escape_next = False
            end_pos = -1
            
            for i, char in enumerate(content):
                if escape_next:
                    escape_next = False
                    continue
                if char == '\\':
                    escape_next = True
                    continue
                if char == '"' and not escape_next:
                    in_string = not in_string
                    continue
                if not in_string:
                    if char == '{':
                        brace_count += 1
                    elif char == '}':
                        brace_count -= 1
                        if brace_count == 0 and bracket_count == 0:
                            end_pos = i + 1
                            break
                    elif char == '[':
                        bracket_count += 1
                    elif char == ']':
                        bracket_count -= 1
                        if bracket_count == 0 and brace_count == 0:
                            end_pos = i + 1
                            break
            
            if end_pos > 0:
                return json.loads(content[:end_pos])
            else:
                return json.loads(content)
    except Exception as e:
        print(f"Error loading JSON file {json_file}: {e}", file=sys.stderr)
        return None


def extract_atoms_map(json_data: Dict, json_file: Optional[Path] = None, is_legacy: bool = False) -> Dict[int, Dict]:
    """Extract atom records and create map by legacy_atom_idx."""
    atoms_map = {}
    atoms_by_name = {}  # Also create name->atom map for lookup
    
    # Check if split files are indicated
    is_split_files = False
    calculations = json_data.get('calculations', {})
    
    if isinstance(calculations, list):
        for calc in calculations:
            if isinstance(calc, dict) and calc.get('_split_files'):
                is_split_files = True
                break
    
    # Load from split file if available
    if is_split_files and json_file:
        pdb_name = json_data.get('pdb_name', json_file.stem)
        split_file = json_file.parent / f"{pdb_name}_pdb_atoms.json"
        if split_file.exists():
            split_data = load_json_file(split_file)
            if isinstance(split_data, list) and len(split_data) > 0:
                atoms = split_data[0].get('atoms', [])
                for atom in atoms:
                    # Legacy uses atom_idx, modern uses legacy_atom_idx
                    legacy_idx = atom.get('atom_idx') if is_legacy else atom.get('legacy_atom_idx')
                    if legacy_idx:
                        atoms_map[legacy_idx] = atom
                    # Also index by name for lookup
                    atom_name = atom.get('atom_name', '').strip()
                    if atom_name:
                        atoms_by_name[atom_name] = atom
    
    # Also try grouped format
    if not atoms_map and isinstance(calculations, dict):
        pdb_atoms_group = calculations.get('pdb_atoms', [])
        if pdb_atoms_group and len(pdb_atoms_group) > 0:
            atoms = pdb_atoms_group[0].get('atoms', [])
            for atom in atoms:
                legacy_idx = atom.get('atom_idx') if is_legacy else atom.get('legacy_atom_idx')
                if legacy_idx:
                    atoms_map[legacy_idx] = atom
                atom_name = atom.get('atom_name', '').strip()
                if atom_name:
                    atoms_by_name[atom_name] = atom
    
    # Also try array format
    if not atoms_map and isinstance(calculations, list):
        for calc in calculations:
            if calc.get('type') == 'pdb_atoms':
                atoms = calc.get('atoms', [])
                for atom in atoms:
                    legacy_idx = atom.get('atom_idx') if is_legacy else atom.get('legacy_atom_idx')
                    if legacy_idx:
                        atoms_map[legacy_idx] = atom
                    atom_name = atom.get('atom_name', '').strip()
                    if atom_name:
                        atoms_by_name[atom_name] = atom
                break
    
    # Store name map in atoms_map for lookup
    atoms_map['_by_name'] = atoms_by_name
    return atoms_map


def extract_frame_calculations(json_data: Dict, json_file: Optional[Path] = None) -> List[Dict]:
    """Extract base_frame_calc records from JSON."""
    calculations = json_data.get('calculations', {})
    frame_calcs = []
    seen_keys = set()  # Track (chain_id, residue_seq, insertion) to avoid duplicates
    
    # Check if split files
    is_split_files = False
    if isinstance(calculations, list):
        for calc in calculations:
            if isinstance(calc, dict) and calc.get('_split_files'):
                is_split_files = True
                break
    
    # Load from split file if available
    if is_split_files and json_file:
        pdb_name = json_data.get('pdb_name', json_file.stem)
        split_file = json_file.parent / f"{pdb_name}_base_frame_calc.json"
        if split_file.exists():
            split_data = load_json_file(split_file)
            if isinstance(split_data, list):
                for calc in split_data:
                    key = (calc.get('chain_id', ''), calc.get('residue_seq', 0), calc.get('insertion', ' '))
                    if key not in seen_keys:
                        frame_calcs.append(calc)
                        seen_keys.add(key)
            elif isinstance(split_data, dict) and 'base_frame_calc' in split_data:
                for calc in split_data['base_frame_calc']:
                    key = (calc.get('chain_id', ''), calc.get('residue_seq', 0), calc.get('insertion', ' '))
                    if key not in seen_keys:
                        frame_calcs.append(calc)
                        seen_keys.add(key)
    
    # Also try grouped format
    if not frame_calcs and isinstance(calculations, dict):
        for calc in calculations.get('base_frame_calc', []):
            key = (calc.get('chain_id', ''), calc.get('residue_seq', 0), calc.get('insertion', ' '))
            if key not in seen_keys:
                frame_calcs.append(calc)
                seen_keys.add(key)
    
    # Also try array format
    if not frame_calcs and isinstance(calculations, list):
        for calc in calculations:
            if calc.get('type') == 'base_frame_calc':
                key = (calc.get('chain_id', ''), calc.get('residue_seq', 0), calc.get('insertion', ' '))
                if key not in seen_keys:
                    frame_calcs.append(calc)
                    seen_keys.add(key)
    
    return frame_calcs


def get_atom_by_legacy_idx(atoms_map: Dict[int, Dict], legacy_idx: int) -> Optional[Dict]:
    """Get atom by legacy_atom_idx, using atom name matching as fallback."""
    if legacy_idx in atoms_map:
        return atoms_map[legacy_idx]
    return None


def find_atom_by_name_in_residue(atoms_map: Dict, atom_name: str, residue_key: Tuple[str, int, str]) -> Optional[Dict]:
    """Find atom by name within a specific residue."""
    atom_name_clean = atom_name.strip()
    chain_id, residue_seq, insertion = residue_key
    
    # Search all atoms in the map
    for key, atom in atoms_map.items():
        if key == '_by_name':
            continue
        if atom.get('atom_name', '').strip() == atom_name_clean:
            # Check if this atom belongs to the specified residue
            atom_chain = atom.get('chain_id', '')
            atom_res_seq = atom.get('residue_seq', 0)
            atom_ins = atom.get('insertion', ' ')
            if (atom_chain, atom_res_seq, atom_ins) == residue_key:
                return atom
    return None


def compare_ring_atoms(legacy_file: Path, modern_file: Path, pdb_file: Optional[Path] = None) -> List[RingAtomMatch]:
    """Compare ring atom matching between legacy and modern JSON."""
    legacy_data = load_json_file(legacy_file)
    modern_data = load_json_file(modern_file)
    
    if not legacy_data or not modern_data:
        return []
    
    # Load PDB file reader if available
    pdb_reader = None
    if pdb_file and pdb_file.exists():
        try:
            pdb_reader = PdbFileReader(pdb_file)
        except Exception as e:
            print(f"Warning: Could not load PDB file {pdb_file}: {e}", file=sys.stderr)
    
    # Extract atoms maps
    legacy_atoms = extract_atoms_map(legacy_data, legacy_file, is_legacy=True)
    modern_atoms = extract_atoms_map(modern_data, modern_file, is_legacy=False)
    
    # Extract frame calculations
    legacy_frames = extract_frame_calculations(legacy_data, legacy_file)
    modern_frames = extract_frame_calculations(modern_data, modern_file)
    
    # Build maps by residue key and by legacy_residue_idx
    legacy_by_residue = {}
    legacy_by_residue_idx = {}
    for frame in legacy_frames:
        chain_id = frame.get('chain_id', '')
        residue_seq = frame.get('residue_seq', 0)
        insertion = frame.get('insertion', ' ')
        key = (chain_id, residue_seq, insertion)
        legacy_by_residue[key] = frame
        # Also index by legacy_residue_idx if available (legacy uses residue_idx)
        residue_idx = frame.get('residue_idx')
        if residue_idx:
            legacy_by_residue_idx[residue_idx] = frame
    
    modern_by_residue = {}
    modern_by_residue_idx = {}
    for frame in modern_frames:
        chain_id = frame.get('chain_id', '')
        residue_seq = frame.get('residue_seq', 0)
        insertion = frame.get('insertion', ' ')
        key = (chain_id, residue_seq, insertion)
        modern_by_residue[key] = frame
        # Also index by legacy_residue_idx
        legacy_residue_idx = frame.get('legacy_residue_idx')
        if legacy_residue_idx:
            modern_by_residue_idx[legacy_residue_idx] = frame
    
    results = []
    
    # Compare residues - try to match by legacy_residue_idx first, then by key
    matched_residue_indices = set()
    processed_keys = set()
    
    # First, match by legacy_residue_idx
    for leg_residue_idx, leg_frame in legacy_by_residue_idx.items():
        if leg_residue_idx in modern_by_residue_idx:
            mod_frame = modern_by_residue_idx[leg_residue_idx]
            matched_residue_indices.add(leg_residue_idx)
            
            # Process this matched pair
            chain_id = leg_frame.get('chain_id', '')
            residue_seq = leg_frame.get('residue_seq', 0)
            insertion = leg_frame.get('insertion', ' ')
            key = (chain_id, residue_seq, insertion)
            processed_keys.add(key)
            
            leg_matched = leg_frame.get('matched_atoms', []) if leg_frame else []
            mod_matched = mod_frame.get('matched_atoms', []) if mod_frame else []
            
            # Get atom indices for matched atoms
            leg_indices = []
            mod_indices = []
            differences = []
            
            # For legacy: find atoms by name within this residue and get their atom_idx
            leg_residue_key = (leg_frame.get('chain_id', ''), leg_frame.get('residue_seq', 0), leg_frame.get('insertion', ' '))
            for atom_name in leg_matched:
                atom = find_atom_by_name_in_residue(legacy_atoms, atom_name, leg_residue_key)
                if atom:
                    leg_indices.append(atom.get('atom_idx'))  # Legacy uses atom_idx
                else:
                    differences.append(f"Legacy atom '{atom_name}' not found in residue {leg_residue_key}")
            
            # For modern: find atoms by name within this residue and get their legacy_atom_idx
            mod_residue_key = (mod_frame.get('chain_id', ''), mod_frame.get('residue_seq', 0), mod_frame.get('insertion', ' '))
            for atom_name in mod_matched:
                atom = find_atom_by_name_in_residue(modern_atoms, atom_name, mod_residue_key)
                if atom:
                    mod_indices.append(atom.get('legacy_atom_idx'))
                else:
                    differences.append(f"Modern atom '{atom_name}' not found in residue {mod_residue_key}")
            
            # Compare indices (AFTER collecting all atoms - only add each difference once)
            leg_set = set(leg_indices)
            mod_set = set(mod_indices)
            
            if leg_set != mod_set:
                only_legacy = sorted(leg_set - mod_set)
                only_modern = sorted(mod_set - leg_set)
                if only_legacy:
                    differences.append(f"Only in legacy: legacy_atom_idx {only_legacy}")
                if only_modern:
                    differences.append(f"Only in modern: legacy_atom_idx {only_modern}")
            
            # Compare order (only if sets match)
            if leg_indices != mod_indices and leg_set == mod_set:
                differences.append(f"Same atoms but different order")
            
            # Compare atom names (only add each difference once)
            if leg_matched != mod_matched:
                only_legacy_names = sorted(set(leg_matched) - set(mod_matched))
                only_modern_names = sorted(set(mod_matched) - set(leg_matched))
                if only_legacy_names:
                    differences.append(f"Only in legacy names: {only_legacy_names}")
                if only_modern_names:
                    differences.append(f"Only in modern names: {only_modern_names}")
            
            is_match = len(differences) == 0
            
            # Collect atom details with PDB lines for atoms that differ
            leg_atom_details = {}
            mod_atom_details = {}
            
            # Get details for legacy atoms
            for idx in leg_indices:
                if idx in legacy_atoms:
                    atom = legacy_atoms[idx]
                    line_num = atom.get('line_number', 0)
                    pdb_line = atom.get('pdb_line', '')
                    if not pdb_line and pdb_reader and line_num > 0:
                        pdb_line = pdb_reader.get_line(line_num) or ''
                    leg_atom_details[idx] = {
                        'atom_name': atom.get('atom_name', ''),
                        'line_number': line_num,
                        'pdb_line': pdb_line,
                        'xyz': atom.get('xyz', [])
                    }
            
            # Get details for modern atoms
            for idx in mod_indices:
                if idx in modern_atoms:
                    atom = modern_atoms[idx]
                    line_num = atom.get('line_number', 0)
                    pdb_line = atom.get('pdb_line', '')
                    if not pdb_line and pdb_reader and line_num > 0:
                        pdb_line = pdb_reader.get_line(line_num) or ''
                    mod_atom_details[idx] = {
                        'atom_name': atom.get('atom_name', ''),
                        'line_number': line_num,
                        'pdb_line': pdb_line,
                        'xyz': atom.get('xyz', [])
                    }
            
            residue_name = leg_frame.get('residue_name', '')
            base_type = leg_frame.get('base_type', '')
            
            results.append(RingAtomMatch(
                residue_key=key,
                residue_name=residue_name,
                base_type=base_type,
                legacy_matched_atoms=leg_matched,
                modern_matched_atoms=mod_matched,
                legacy_atom_indices=leg_indices,
                modern_atom_indices=mod_indices,
                num_matched_legacy=len(leg_matched),
                num_matched_modern=len(mod_matched),
                is_match=is_match,
                differences=differences,
                legacy_atom_details=leg_atom_details,
                modern_atom_details=mod_atom_details
            ))
    
    # Also handle residues matched by key (not by legacy_residue_idx)
    all_keys = set(legacy_by_residue.keys()) | set(modern_by_residue.keys())
    for key in sorted(all_keys):
        # Skip if already processed
        if key in processed_keys:
            continue
            
        leg_frame = legacy_by_residue.get(key)
        mod_frame = modern_by_residue.get(key)
        
        if not leg_frame or not mod_frame:
            continue
        
        # Process this pair
        chain_id, residue_seq, insertion = key
        leg_matched = leg_frame.get('matched_atoms', [])
        mod_matched = mod_frame.get('matched_atoms', [])
        
        # Get atom indices (same logic as above)
        leg_indices = []
        mod_indices = []
        differences = []
        
        leg_residue_key = (leg_frame.get('chain_id', ''), leg_frame.get('residue_seq', 0), leg_frame.get('insertion', ' '))
        for atom_name in leg_matched:
            atom = find_atom_by_name_in_residue(legacy_atoms, atom_name, leg_residue_key)
            if atom:
                leg_indices.append(atom.get('atom_idx'))
            else:
                differences.append(f"Legacy atom '{atom_name}' not found in residue {leg_residue_key}")
        
        mod_residue_key = (mod_frame.get('chain_id', ''), mod_frame.get('residue_seq', 0), mod_frame.get('insertion', ' '))
        for atom_name in mod_matched:
            atom = find_atom_by_name_in_residue(modern_atoms, atom_name, mod_residue_key)
            if atom:
                mod_indices.append(atom.get('legacy_atom_idx'))
            else:
                differences.append(f"Modern atom '{atom_name}' not found in residue {mod_residue_key}")
        
        # Compare (same as above) - only add each difference once
        leg_set = set(leg_indices)
        mod_set = set(mod_indices)
        
        if leg_set != mod_set:
            only_legacy = sorted(leg_set - mod_set)
            only_modern = sorted(mod_set - leg_set)
            if only_legacy:
                differences.append(f"Only in legacy: legacy_atom_idx {only_legacy}")
            if only_modern:
                differences.append(f"Only in modern: legacy_atom_idx {only_modern}")
        
        if leg_indices != mod_indices and leg_set == mod_set:
            differences.append(f"Same atoms but different order")
        
        if leg_matched != mod_matched:
            only_legacy_names = sorted(set(leg_matched) - set(mod_matched))
            only_modern_names = sorted(set(mod_matched) - set(leg_matched))
            if only_legacy_names:
                differences.append(f"Only in legacy names: {only_legacy_names}")
            if only_modern_names:
                differences.append(f"Only in modern names: {only_modern_names}")
        
        is_match = len(differences) == 0
        
        # Collect atom details with PDB lines for atoms that differ
        leg_atom_details = {}
        mod_atom_details = {}
        
        # Get details for legacy atoms
        for idx in leg_indices:
            if idx in legacy_atoms:
                atom = legacy_atoms[idx]
                line_num = atom.get('line_number', 0)
                pdb_line = atom.get('pdb_line', '')
                if not pdb_line and pdb_reader and line_num > 0:
                    pdb_line = pdb_reader.get_line(line_num) or ''
                leg_atom_details[idx] = {
                    'atom_name': atom.get('atom_name', ''),
                    'line_number': line_num,
                    'pdb_line': pdb_line,
                    'xyz': atom.get('xyz', [])
                }
        
        # Get details for modern atoms
        for idx in mod_indices:
            if idx in modern_atoms:
                atom = modern_atoms[idx]
                line_num = atom.get('line_number', 0)
                pdb_line = atom.get('pdb_line', '')
                if not pdb_line and pdb_reader and line_num > 0:
                    pdb_line = pdb_reader.get_line(line_num) or ''
                mod_atom_details[idx] = {
                    'atom_name': atom.get('atom_name', ''),
                    'line_number': line_num,
                    'pdb_line': pdb_line,
                    'xyz': atom.get('xyz', [])
                }
        
        residue_name = leg_frame.get('residue_name', '')
        base_type = leg_frame.get('base_type', '')
        
        results.append(RingAtomMatch(
            residue_key=key,
            residue_name=residue_name,
            base_type=base_type,
            legacy_matched_atoms=leg_matched,
            modern_matched_atoms=mod_matched,
            legacy_atom_indices=leg_indices,
            modern_atom_indices=mod_indices,
            num_matched_legacy=len(leg_matched),
            num_matched_modern=len(mod_matched),
            is_match=is_match,
            differences=differences,
            legacy_atom_details=leg_atom_details,
            modern_atom_details=mod_atom_details
        ))
    
    return results


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Compare ring atom matching between legacy and modern code')
    parser.add_argument('pdb_id', help='PDB ID to compare (e.g., 1H4S)')
    parser.add_argument('--legacy-dir', type=Path, default=Path('data/json_legacy'),
                       help='Directory containing legacy JSON files')
    parser.add_argument('--modern-dir', type=Path, default=Path('data/json'),
                       help='Directory containing modern JSON files')
    parser.add_argument('--diff-only', action='store_true',
                       help='Show only residues with differences')
    parser.add_argument('--pdb-dir', type=Path, default=Path('data/pdb'),
                       help='Directory containing PDB files')
    parser.add_argument('--show-pdb-lines', action='store_true',
                       help='Show PDB lines for atoms that differ')
    
    args = parser.parse_args()
    
    legacy_file = args.legacy_dir / f"{args.pdb_id}.json"
    modern_file = args.modern_dir / f"{args.pdb_id}.json"
    pdb_file = args.pdb_dir / f"{args.pdb_id}.pdb"
    
    if not legacy_file.exists():
        print(f"Error: Legacy file not found: {legacy_file}", file=sys.stderr)
        sys.exit(1)
    
    if not modern_file.exists():
        print(f"Error: Modern file not found: {modern_file}", file=sys.stderr)
        sys.exit(1)
    
    results = compare_ring_atoms(legacy_file, modern_file, pdb_file if args.show_pdb_lines else None)
    
    # Print report
    print("=" * 80)
    print("RING ATOM MATCHING COMPARISON")
    print("=" * 80)
    print(f"PDB ID: {args.pdb_id}")
    print()
    
    matches = sum(1 for r in results if r.is_match)
    mismatches = len(results) - matches
    
    print(f"SUMMARY")
    print(f"  Total residues: {len(results)}")
    print(f"  Matches: {matches}")
    print(f"  Mismatches: {mismatches}")
    print()
    
    if mismatches > 0 or not args.diff_only:
        print("=" * 80)
        print("DETAILED RESULTS")
        print("=" * 80)
        print()
        
        for result in results:
            if args.diff_only and result.is_match:
                continue
            
            chain_id, residue_seq, insertion = result.residue_key
            status = "✓" if result.is_match else "✗"
            
            print(f"{status} {result.residue_name} {chain_id}:{residue_seq}{insertion} ({result.base_type})")
            print(f"  Legacy: {result.num_matched_legacy} atoms - {result.legacy_matched_atoms}")
            print(f"  Modern: {result.num_matched_modern} atoms - {result.modern_matched_atoms}")
            
            if result.legacy_atom_indices or result.modern_atom_indices:
                print(f"  Legacy indices: {result.legacy_atom_indices}")
                print(f"  Modern indices: {result.modern_atom_indices}")
            
            if result.differences:
                print(f"  Differences:")
                for diff in result.differences:
                    print(f"    - {diff}")
            
            # Show PDB lines for atoms that differ
            if args.show_pdb_lines and (result.legacy_atom_details or result.modern_atom_details):
                leg_set = set(result.legacy_atom_indices)
                mod_set = set(result.modern_atom_indices)
                only_legacy = leg_set - mod_set
                only_modern = mod_set - leg_set
                
                if only_legacy or only_modern:
                    print(f"  PDB Lines for differing atoms:")
                    
                    if only_legacy:
                        print(f"    Only in Legacy:")
                        for idx in sorted(only_legacy):
                            if idx in result.legacy_atom_details:
                                details = result.legacy_atom_details[idx]
                                atom_name = details.get('atom_name', '?')
                                line_num = details.get('line_number', 0)
                                pdb_line = details.get('pdb_line', '')
                                if pdb_line:
                                    print(f"      legacy_atom_idx {idx} ({atom_name}) line {line_num}: {pdb_line[:80]}")
                                else:
                                    print(f"      legacy_atom_idx {idx} ({atom_name}) line {line_num}: (PDB line not available)")
                    
                    if only_modern:
                        print(f"    Only in Modern:")
                        for idx in sorted(only_modern):
                            if idx in result.modern_atom_details:
                                details = result.modern_atom_details[idx]
                                atom_name = details.get('atom_name', '?')
                                line_num = details.get('line_number', 0)
                                pdb_line = details.get('pdb_line', '')
                                if pdb_line:
                                    print(f"      legacy_atom_idx {idx} ({atom_name}) line {line_num}: {pdb_line[:80]}")
                                else:
                                    print(f"      legacy_atom_idx {idx} ({atom_name}) line {line_num}: (PDB line not available)")
            print()


if __name__ == '__main__':
    main()


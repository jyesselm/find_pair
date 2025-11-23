#!/usr/bin/env python3
"""
Comprehensive analysis of ring fitting differences between legacy and modern code.

This script analyzes JSON files and PDB files to understand why different atoms
are matched during ring fitting, including actual PDB lines for all relevant atoms.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional
from collections import defaultdict
import argparse


# Ring atom definitions (from RA_LIST)
RING_ATOMS_PURINE = [" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "]
RING_ATOMS_PYRIMIDINE = [" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "]


class PdbReader:
    """Simple PDB file reader with indexing."""
    
    def __init__(self, pdb_file: Path):
        self.pdb_file = pdb_file
        self.lines = []
        self.atom_index = {}  # (chain_id, residue_seq, insertion, atom_name) -> (line_num, line)
        
        if not pdb_file.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        self._load()
    
    def _load(self):
        """Load PDB file and build atom index."""
        with open(self.pdb_file, 'r') as f:
            self.lines = f.readlines()
        
        for line_num, line in enumerate(self.lines, 1):
            if line.startswith(('ATOM', 'HETATM')) and len(line) >= 16:
                chain_id = line[21] if len(line) > 21 else ' '
                try:
                    residue_seq = int(line[22:26].strip())
                except (ValueError, IndexError):
                    continue
                insertion = line[26] if len(line) > 26 else ' '
                atom_name = line[12:16] if len(line) >= 16 else ''
                
                if atom_name:
                    key = (chain_id, residue_seq, insertion, atom_name)
                    self.atom_index[key] = (line_num, line.rstrip('\n\r'))
    
    def find_atom(self, chain_id: str, residue_seq: int, insertion: str, atom_name: str) -> Optional[Tuple[int, str]]:
        """Find atom in PDB file. Returns (line_number, line_content) or None."""
        key = (chain_id, residue_seq, insertion, atom_name)
        return self.atom_index.get(key)
    
    def get_all_residue_atoms(self, chain_id: str, residue_seq: int, insertion: str = ' ') -> Dict[str, Tuple[int, str]]:
        """Get all atoms for a residue. Returns dict of atom_name -> (line_num, line)."""
        result = {}
        for (ch, seq, ins, atom), (line_num, line) in self.atom_index.items():
            if ch == chain_id and seq == residue_seq and ins == insertion:
                result[atom] = (line_num, line)
        return result


def load_json(file_path: Path) -> Optional[Dict]:
    """Load JSON file."""
    if not file_path.exists():
        return None
    with open(file_path, 'r') as f:
        return json.load(f)


def get_ring_atoms_for_base(base_type: str) -> List[str]:
    """Get list of ring atoms for a base type."""
    base_type_upper = base_type.upper()
    if base_type_upper in ['A', 'G']:
        return RING_ATOMS_PURINE
    else:
        return RING_ATOMS_PYRIMIDINE


def analyze_ring_fitting_differences(pdb_id: str, project_root: Path) -> Dict:
    """Analyze ring fitting differences for a single PDB."""
    pdb_file = project_root / 'data' / 'pdb' / f'{pdb_id}.pdb'
    legacy_json_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
    modern_json_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
    
    result = {
        'pdb_id': pdb_id,
        'pdb_file': str(pdb_file),
        'errors': [],
        'residue_analyses': []
    }
    
    # Check files exist
    if not pdb_file.exists():
        result['errors'].append(f"PDB file not found: {pdb_file}")
        return result
    
    if not legacy_json_file.exists():
        result['errors'].append(f"Legacy JSON not found: {legacy_json_file}")
        return result
    
    if not modern_json_file.exists():
        result['errors'].append(f"Modern JSON not found: {modern_json_file}")
        return result
    
    # Load files
    try:
        legacy_json = load_json(legacy_json_file)
        modern_json = load_json(modern_json_file)
        pdb_reader = PdbReader(pdb_file)
    except Exception as e:
        result['errors'].append(f"Error loading files: {e}")
        return result
    
    # Extract frame calculation records
    legacy_base_frame = {}  # For matched_atoms
    legacy_calc = {}  # For calculations
    
    for rec in legacy_json.get('calculations', []):
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
                legacy_calc[key].update(rec)
    
    modern_map = {}
    for rec in modern_json.get('calculations', []):
        if rec.get('type') in ['base_frame_calc', 'frame_calc', 'ls_fitting']:
            chain_id = rec.get('chain_id', '')
            residue_seq = rec.get('residue_seq', 0)
            insertion = rec.get('insertion', ' ')
            key = (chain_id, residue_seq, insertion)
            if key not in modern_map:
                modern_map[key] = rec
            else:
                modern_map[key].update(rec)
    
    # Analyze each residue
    all_keys = set(legacy_calc.keys()) | set(modern_map.keys())
    
    for key in all_keys:
        chain_id, residue_seq, insertion = key
        
        # Get legacy data
        leg_matched_atoms = []
        leg_base_type = '?'
        leg_rms = 0.0
        leg_num_matched = 0
        leg_residue_name = ''
        
        if key in legacy_base_frame:
            leg_rec = legacy_base_frame[key]
            leg_matched_atoms = leg_rec.get('matched_atoms', [])
            leg_base_type = leg_rec.get('base_type', '?')
            leg_residue_name = leg_rec.get('residue_name', '')
        
        if key in legacy_calc:
            leg_calc = legacy_calc[key]
            leg_rms = leg_calc.get('rms_fit', 0.0)
            leg_num_matched = leg_calc.get('num_matched_atoms', leg_calc.get('num_points', 0))
            if not leg_residue_name:
                leg_residue_name = leg_calc.get('residue_name', '')
            if leg_base_type == '?':
                leg_base_type = leg_calc.get('base_type', '?')
        
        # Get modern data
        mod_matched_atoms = []
        mod_base_type = '?'
        mod_rms = 0.0
        mod_num_matched = 0
        mod_residue_name = ''
        
        if key in modern_map:
            mod_rec = modern_map[key]
            mod_matched_atoms = mod_rec.get('matched_atoms', [])
            mod_base_type = mod_rec.get('base_type', '?')
            mod_rms = mod_rec.get('rms_fit', 0.0)
            mod_num_matched = mod_rec.get('num_matched_atoms', mod_rec.get('num_points', 0))
            mod_residue_name = mod_rec.get('residue_name', '')
        
        # Check if there are differences
        leg_atoms_set = set(leg_matched_atoms) if leg_matched_atoms else set()
        mod_atoms_set = set(mod_matched_atoms) if mod_matched_atoms else set()
        has_difference = (
            leg_atoms_set != mod_atoms_set or
            abs(leg_rms - mod_rms) > 0.001 or
            leg_num_matched != mod_num_matched or
            key not in legacy_calc or
            key not in modern_map
        )
        
        if not has_difference:
            continue  # Skip residues with no differences
        
        # Determine base type for ring atom list
        base_type = leg_base_type if leg_base_type != '?' else mod_base_type
        if base_type == '?':
            base_type = leg_residue_name.strip() if leg_residue_name else '?'
        
        expected_ring_atoms = get_ring_atoms_for_base(base_type)
        
        # Get all atoms in PDB for this residue
        pdb_atoms = pdb_reader.get_all_residue_atoms(chain_id, residue_seq, insertion)
        
        # Build atom analysis
        atom_analysis = []
        all_atoms_to_check = set(leg_matched_atoms) | set(mod_matched_atoms) | set(expected_ring_atoms)
        
        for atom_name in sorted(all_atoms_to_check):
            in_pdb = atom_name in pdb_atoms
            in_legacy = atom_name in leg_atoms_set
            in_modern = atom_name in mod_atoms_set
            is_expected = atom_name in expected_ring_atoms
            
            atom_info = {
                'atom_name': atom_name,
                'in_pdb': in_pdb,
                'in_legacy_matched': in_legacy,
                'in_modern_matched': in_modern,
                'is_expected_ring_atom': is_expected,
                'pdb_line': None,
                'pdb_line_number': None
            }
            
            if in_pdb:
                line_num, line = pdb_atoms[atom_name]
                atom_info['pdb_line'] = line
                atom_info['pdb_line_number'] = line_num
            
            atom_analysis.append(atom_info)
        
        # Build residue analysis
        residue_analysis = {
            'residue_key': {
                'chain_id': chain_id,
                'residue_seq': residue_seq,
                'insertion': insertion
            },
            'residue_name': leg_residue_name or mod_residue_name,
            'base_type': base_type,
            'legacy': {
                'matched_atoms': leg_matched_atoms,
                'num_matched': leg_num_matched,
                'rms_fit': leg_rms,
                'present': key in legacy_calc
            },
            'modern': {
                'matched_atoms': mod_matched_atoms,
                'num_matched': mod_num_matched,
                'rms_fit': mod_rms,
                'present': key in modern_map
            },
            'differences': {
                'only_in_legacy': list(leg_atoms_set - mod_atoms_set),
                'only_in_modern': list(mod_atoms_set - leg_atoms_set),
                'rms_diff': abs(leg_rms - mod_rms),
                'num_matched_diff': leg_num_matched != mod_num_matched
            },
            'expected_ring_atoms': expected_ring_atoms,
            'atom_analysis': atom_analysis,
            'all_pdb_atoms': {atom: {'line_num': num, 'line': line} 
                            for atom, (num, line) in pdb_atoms.items()}
        }
        
        result['residue_analyses'].append(residue_analysis)
    
    return result


def format_report(analysis: Dict) -> str:
    """Format analysis as a readable report."""
    lines = []
    lines.append("=" * 80)
    lines.append(f"RING FITTING DIFFERENCE ANALYSIS")
    lines.append("=" * 80)
    lines.append(f"PDB ID: {analysis['pdb_id']}")
    lines.append(f"PDB File: {analysis['pdb_file']}")
    lines.append("")
    
    if analysis['errors']:
        lines.append("ERRORS:")
        for error in analysis['errors']:
            lines.append(f"  - {error}")
        lines.append("")
        return "\n".join(lines)
    
    if not analysis['residue_analyses']:
        lines.append("No differences found.")
        lines.append("")
        return "\n".join(lines)
    
    lines.append(f"Found {len(analysis['residue_analyses'])} residues with differences\n")
    
    for i, residue in enumerate(analysis['residue_analyses'], 1):
        key = residue['residue_key']
        lines.append("-" * 80)
        lines.append(f"RESIDUE {i}: {residue['residue_name']} {key['chain_id']}:{key['residue_seq']}{key['insertion']}")
        lines.append("-" * 80)
        lines.append(f"Base Type: {residue['base_type']}")
        lines.append(f"Expected Ring Atoms: {', '.join(residue['expected_ring_atoms'])}")
        lines.append("")
        
        # Legacy vs Modern comparison
        lines.append("MATCHED ATOMS:")
        lines.append(f"  Legacy ({residue['legacy']['num_matched']} atoms): {residue['legacy']['matched_atoms']}")
        lines.append(f"  Modern ({residue['modern']['num_matched']} atoms): {residue['modern']['matched_atoms']}")
        lines.append("")
        
        if residue['differences']['only_in_legacy']:
            lines.append(f"  Only in Legacy: {residue['differences']['only_in_legacy']}")
        if residue['differences']['only_in_modern']:
            lines.append(f"  Only in Modern: {residue['differences']['only_in_modern']}")
        lines.append("")
        
        # RMS comparison
        if residue['differences']['rms_diff'] > 0.001:
            lines.append(f"RMS Fit: Legacy={residue['legacy']['rms_fit']:.6f}, "
                        f"Modern={residue['modern']['rms_fit']:.6f}, "
                        f"Diff={residue['differences']['rms_diff']:.6f}")
            lines.append("")
        
        # Atom analysis table
        lines.append("ATOM ANALYSIS:")
        lines.append("  Atom Name | In PDB | Expected | Legacy | Modern | PDB Line")
        lines.append("  " + "-" * 75)
        
        for atom_info in residue['atom_analysis']:
            atom_name = atom_info['atom_name']
            in_pdb = "✓" if atom_info['in_pdb'] else "✗"
            expected = "✓" if atom_info['is_expected_ring_atom'] else " "
            in_legacy = "✓" if atom_info['in_legacy_matched'] else "✗"
            in_modern = "✓" if atom_info['in_modern_matched'] else "✗"
            
            line_info = ""
            if atom_info['pdb_line']:
                line_info = f"Line {atom_info['pdb_line_number']}: {atom_info['pdb_line']}"
            
            lines.append(f"  {atom_name:10} |   {in_pdb}   |    {expected}    |   {in_legacy}   |   {in_modern}   | {line_info}")
        
        lines.append("")
        
        # All PDB atoms for this residue (for completeness)
        if residue['all_pdb_atoms']:
            lines.append(f"ALL ATOMS IN PDB FOR THIS RESIDUE ({len(residue['all_pdb_atoms'])} total):")
            for atom_name in sorted(residue['all_pdb_atoms'].keys()):
                atom_data = residue['all_pdb_atoms'][atom_name]
                lines.append(f"  Line {atom_data['line_num']}: {atom_data['line']}")
            lines.append("")
    
    lines.append("=" * 80)
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description='Analyze ring fitting differences between legacy and modern code')
    parser.add_argument('pdb_id', help='PDB identifier (e.g., 1H4S)')
    parser.add_argument('--project-root', type=Path, default=Path(__file__).parent.parent,
                       help='Project root directory (default: parent of script directory)')
    parser.add_argument('--output', type=Path, help='Output file path (default: print to stdout)')
    parser.add_argument('--json', action='store_true', help='Output as JSON instead of text')
    
    args = parser.parse_args()
    
    # Analyze
    analysis = analyze_ring_fitting_differences(args.pdb_id, args.project_root)
    
    # Format output
    if args.json:
        output = json.dumps(analysis, indent=2)
    else:
        output = format_report(analysis)
    
    # Write or print
    if args.output:
        with open(args.output, 'w') as f:
            f.write(output)
        print(f"Report written to: {args.output}")
    else:
        print(output)


if __name__ == '__main__':
    main()


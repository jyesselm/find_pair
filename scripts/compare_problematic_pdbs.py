#!/usr/bin/env python3
"""
Compare problematic PDBs between legacy and modern code, generating detailed reports
with actual PDB line content for differences.
"""

import os
import sys
import json
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from typing import Dict, List, Tuple, Optional
import time

from x3dna_json_compare import (
    JsonValidator,
    JsonComparator,
    PdbFileReader,
    ParallelExecutor,
)

def load_problematic_pdbs():
    """Load PDB IDs from problematic_pdbs.txt."""
    pdbs = []
    filepath = Path('docs/problematic_pdbs.txt')
    if not filepath.exists():
        return pdbs
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            pdb_id = line.split()[0] if line.split() else None
            if pdb_id:
                pdbs.append(pdb_id)
    return pdbs

def regenerate_legacy_json(pdb_id, executable_path, project_root):
    """Regenerate legacy JSON for a single PDB."""
    if isinstance(project_root, str):
        project_root = Path(project_root)
    if isinstance(executable_path, str):
        executable_path = Path(executable_path)
    
    pdb_file = project_root / 'data' / 'pdb' / f'{pdb_id}.pdb'
    json_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
    
    if not pdb_file.exists():
        return {'pdb_id': pdb_id, 'status': 'error', 'message': f'PDB file not found'}
    
    try:
        cmd = [str(executable_path), str(pdb_file)]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
            cwd=str(project_root / 'org')
        )
        
        if json_file.exists():
            # Validate the generated JSON
            is_valid, error_msg, was_removed = JsonValidator.validate_and_clean(json_file, remove_invalid=True)
            if not is_valid:
                return {
                    'pdb_id': pdb_id,
                    'status': 'error',
                    'message': f'Invalid JSON removed: {error_msg}',
                    'removed': was_removed
                }
            return {'pdb_id': pdb_id, 'status': 'success', 'json_file': str(json_file)}
        else:
            return {
                'pdb_id': pdb_id,
                'status': 'error',
                'message': 'JSON file not created',
                'stderr': result.stderr[-500:] if result.stderr else ''
            }
    except subprocess.TimeoutExpired:
        return {'pdb_id': pdb_id, 'status': 'error', 'message': 'Timeout'}
    except Exception as e:
        return {'pdb_id': pdb_id, 'status': 'error', 'message': str(e)}

def regenerate_modern_json_batch(pdb_ids, executable_path, project_root):
    """Regenerate modern JSON for multiple PDBs using the filtered test."""
    if isinstance(project_root, str):
        project_root = Path(project_root)
    if isinstance(executable_path, str):
        executable_path = Path(executable_path)
    
    try:
        # Create a temp file with PDB IDs
        temp_list_file = project_root / 'temp_pdb_list.txt'
        with open(temp_list_file, 'w') as f:
            for pdb_id in pdb_ids:
                f.write(f'{pdb_id}\n')
        
        cmd = [str(executable_path), str(temp_list_file)]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800,  # 30 minutes for batch
            cwd=str(project_root)
        )
        
        # Clean up temp file
        if temp_list_file.exists():
            temp_list_file.unlink()
        
        # Check which JSON files were created
        results = []
        for pdb_id in pdb_ids:
            json_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
            if json_file.exists():
                results.append({'pdb_id': pdb_id, 'status': 'success', 'json_file': str(json_file)})
            else:
                results.append({'pdb_id': pdb_id, 'status': 'error', 'message': 'JSON file not created'})
        
        return results
    except subprocess.TimeoutExpired:
        return [{'pdb_id': pdb_id, 'status': 'error', 'message': 'Timeout'} for pdb_id in pdb_ids]
    except Exception as e:
        return [{'pdb_id': pdb_id, 'status': 'error', 'message': str(e)} for pdb_id in pdb_ids]

# Use PdbFileReader from lib instead of custom function

def compare_atoms_with_pdb_lines(legacy_atoms: List[Dict], modern_atoms: List[Dict],
                                 pdb_file: Path) -> Dict:
    """Compare atoms and include actual PDB line content."""
    differences = {
        'missing_in_modern': [],
        'extra_in_modern': [],
        'mismatched_fields': [],
        'count_difference': len(legacy_atoms) != len(modern_atoms)
    }
    
    # Create maps by atom key (chain_id, residue_seq, atom_name)
    legacy_map = {}
    for atom in legacy_atoms:
        key = (
            atom.get('chain_id', ''),
            atom.get('residue_seq', 0),
            atom.get('atom_name', ''),
            atom.get('insertion', ' ')
        )
        legacy_map[key] = atom
    
    modern_map = {}
    for atom in modern_atoms:
        key = (
            atom.get('chain_id', ''),
            atom.get('residue_seq', 0),
            atom.get('atom_name', ''),
            atom.get('insertion', ' ')
        )
        modern_map[key] = atom
    
    # Find missing in modern
    for key, leg_atom in legacy_map.items():
        if key not in modern_map:
            leg_line_num = leg_atom.get('line_number', 0)
            leg_pdb_line = get_pdb_line(pdb_file, leg_line_num) if leg_line_num > 0 else None
            differences['missing_in_modern'].append({
                'atom': leg_atom,
                'line_number': leg_line_num,
                'pdb_line': leg_pdb_line
            })
    
    # Find extra in modern
    for key, mod_atom in modern_map.items():
        if key not in legacy_map:
            mod_line_num = mod_atom.get('line_number', 0)
            mod_pdb_line = get_pdb_line(pdb_file, mod_line_num) if mod_line_num > 0 else None
            differences['extra_in_modern'].append({
                'atom': mod_atom,
                'line_number': mod_line_num,
                'pdb_line': mod_pdb_line
            })
    
    # Find mismatched fields
    for key in set(legacy_map.keys()) & set(modern_map.keys()):
        leg_atom = legacy_map[key]
        mod_atom = modern_map[key]
        
        # Compare fields
        fields_to_check = ['xyz', 'residue_name', 'record_type', 'alt_loc']
        mismatches = {}
        for field in fields_to_check:
            leg_val = leg_atom.get(field)
            mod_val = mod_atom.get(field)
            if leg_val != mod_val:
                mismatches[field] = {'legacy': leg_val, 'modern': mod_val}
        
        if mismatches:
            leg_line_num = leg_atom.get('line_number', 0)
            mod_line_num = mod_atom.get('line_number', 0)
            leg_pdb_line = get_pdb_line(pdb_file, leg_line_num) if leg_line_num > 0 else None
            mod_pdb_line = get_pdb_line(pdb_file, mod_line_num) if mod_line_num > 0 else None
            differences['mismatched_fields'].append({
                'atom_key': key,
                'legacy_atom': leg_atom,
                'modern_atom': mod_atom,
                'mismatches': mismatches,
                'legacy_line_number': leg_line_num,
                'modern_line_number': mod_line_num,
                'legacy_pdb_line': leg_pdb_line,
                'modern_pdb_line': mod_pdb_line
            })
    
    return differences

def compare_frame_calculations_with_pdb_lines(legacy_records: List[Dict],
                                             modern_records: List[Dict],
                                             pdb_file: Path) -> Dict:
    """Compare frame calculation records and include PDB line info for matched atoms."""
    differences = {
        'missing_residues': [],
        'mismatched_calculations': []
    }
    
    # Separate by type - we need base_frame_calc for matched_atoms, and frame_calc/ls_fitting for calculations
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
    
    # Find missing residues (use legacy_calc for missing calculations)
    for key in set(legacy_calc.keys()) - set(modern_map.keys()):
        # Get matched_atoms from base_frame_calc if available
        calc_rec = legacy_calc[key].copy()
        # Also get residue info from base_frame_calc
        if key in legacy_base_frame:
            bfc_rec = legacy_base_frame[key]
            calc_rec['matched_atoms'] = bfc_rec.get('matched_atoms', [])
            # Ensure residue identification is included
            if 'residue_name' not in calc_rec:
                calc_rec['residue_name'] = bfc_rec.get('residue_name', '')
            if 'chain_id' not in calc_rec:
                calc_rec['chain_id'] = bfc_rec.get('chain_id', '')
            if 'residue_seq' not in calc_rec:
                calc_rec['residue_seq'] = bfc_rec.get('residue_seq', 0)
            if 'insertion' not in calc_rec:
                calc_rec['insertion'] = bfc_rec.get('insertion', ' ')
        differences['missing_residues'].append(calc_rec)
    
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
            if atoms_to_find:
                chain_id, residue_seq, insertion = key
                with open(pdb_file, 'r') as f:
                    lines = f.readlines()
                    for atom_name in atoms_to_find:
                        atom_name_str = str(atom_name).strip()
                        # Try to find this atom in the PDB
                        for line_num, line in enumerate(lines, 1):
                            if line.startswith(('ATOM', 'HETATM')):
                                line_chain = line[21] if len(line) > 21 else ' '
                                try:
                                    line_seq = int(line[22:26].strip())
                                except:
                                    continue
                                line_insertion = line[26] if len(line) > 26 else ' '
                                line_atom = line[12:16] if len(line) >= 16 else ''
                                
                                if (line_chain == chain_id and line_seq == residue_seq and
                                    line_insertion == insertion and line_atom == atom_name_str):
                                    atom_line_info.append({
                                        'atom_name': atom_name_str,
                                        'line_number': line_num,
                                        'pdb_line': line.rstrip('\n')
                                    })
                                    break
            
            differences['mismatched_calculations'].append({
                'residue_key': key,
                'legacy_record': leg_rec,
                'modern_record': mod_rec,
                'mismatches': mismatches,
                'atom_pdb_lines': atom_line_info,
                'legacy_matched_atoms': leg_atoms,
                'modern_matched_atoms': mod_atoms
            })
    
    return differences

def analyze_single_pdb(pdb_id: str, project_root: Path, legacy_exe: Path, modern_exe: Path) -> Dict:
    """Analyze a single PDB and generate detailed comparison."""
    print(f"Processing {pdb_id}...")
    
    pdb_file = project_root / 'data' / 'pdb' / f'{pdb_id}.pdb'
    legacy_json_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
    modern_json_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
    
    result = {
        'pdb_id': pdb_id,
        'status': 'unknown',
        'atom_differences': None,
        'frame_differences': None,
        'errors': []
    }
    
    # Check files exist
    if not pdb_file.exists():
        result['status'] = 'error'
        result['errors'].append('PDB file not found')
        return result
    
    if not legacy_json_file.exists():
        result['status'] = 'error'
        result['errors'].append('Legacy JSON not found')
        return result
    
    if not modern_json_file.exists():
        result['status'] = 'error'
        result['errors'].append('Modern JSON not found')
        return result
    
    try:
        # Load JSON files
        with open(legacy_json_file) as f:
            legacy_json = json.load(f)
        with open(modern_json_file) as f:
            modern_json = json.load(f)
        
        # Extract pdb_atoms records
        legacy_atoms_rec = None
        modern_atoms_rec = None
        
        for rec in legacy_json.get('calculations', []):
            if rec.get('type') == 'pdb_atoms':
                legacy_atoms_rec = rec
                break
        
        for rec in modern_json.get('calculations', []):
            if rec.get('type') == 'pdb_atoms':
                modern_atoms_rec = rec
                break
        
        if legacy_atoms_rec and modern_atoms_rec:
            legacy_atoms = legacy_atoms_rec.get('atoms', [])
            modern_atoms = modern_atoms_rec.get('atoms', [])
            result['atom_differences'] = compare_atoms_with_pdb_lines(
                legacy_atoms, modern_atoms, pdb_file
            )
        
        # Extract frame calculation records (include all types for complete comparison)
        legacy_frame_recs = [r for r in legacy_json.get('calculations', [])
                            if r.get('type') in ['base_frame_calc', 'ls_fitting', 'frame_calc']]
        modern_frame_recs = [r for r in modern_json.get('calculations', [])
                            if r.get('type') in ['base_frame_calc', 'ls_fitting', 'frame_calc']]
        
        if legacy_frame_recs or modern_frame_recs:
            result['frame_differences'] = compare_frame_calculations_with_pdb_lines(
                legacy_frame_recs, modern_frame_recs, pdb_file
            )
        
        # Store pdb_file path for report generation
        result['pdb_file_path'] = str(pdb_file)
        
        # Determine overall status
        has_diff = False
        if result.get('atom_differences'):
            atom_diff = result['atom_differences']
            if (atom_diff.get('count_difference') or
                atom_diff.get('missing_in_modern') or
                atom_diff.get('extra_in_modern') or
                atom_diff.get('mismatched_fields')):
                has_diff = True
        
        if result.get('frame_differences'):
            frame_diff = result['frame_differences']
            if (frame_diff.get('missing_residues') or
                frame_diff.get('mismatched_calculations')):
                has_diff = True
        
        result['status'] = 'diff' if has_diff else 'match'
        
    except Exception as e:
        result['status'] = 'error'
        result['errors'].append(str(e))
    
    return result

def generate_report(all_results: List[Dict], output_file: Path):
    """Generate detailed HTML/text report with PDB line content."""
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("DETAILED PDB COMPARISON REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total PDBs analyzed: {len(all_results)}\n\n")
        
        # Summary
        status_counts = defaultdict(int)
        for r in all_results:
            status_counts[r['status']] += 1
        
        f.write("SUMMARY\n")
        f.write("-" * 80 + "\n")
        for status, count in sorted(status_counts.items()):
            f.write(f"{status}: {count}\n")
        f.write("\n")
        
        # Detailed results
        f.write("\n" + "=" * 80 + "\n")
        f.write("DETAILED RESULTS\n")
        f.write("=" * 80 + "\n\n")
        
        for result in all_results:
            if result['status'] == 'error':
                continue  # Skip errors for now
            
            f.write(f"\n{'=' * 80}\n")
            f.write(f"PDB: {result['pdb_id']}\n")
            f.write(f"{'=' * 80}\n\n")
            
            # Atom differences
            if result['atom_differences']:
                diff = result['atom_differences']
                f.write("ATOM DIFFERENCES:\n")
                f.write("-" * 80 + "\n")
                
                if diff['count_difference']:
                    f.write(f"  Count difference: Legacy={len(diff.get('legacy_count', 0))}, "
                           f"Modern={len(diff.get('modern_count', 0))}\n")
                
                if diff['missing_in_modern']:
                    f.write(f"\n  Missing in modern ({len(diff['missing_in_modern'])}):\n")
                    for item in diff['missing_in_modern'][:10]:  # Limit to 10
                        atom = item['atom']
                        f.write(f"    {atom.get('atom_name')} "
                               f"({atom.get('chain_id')}:{atom.get('residue_seq')})\n")
                        if item['pdb_line']:
                            f.write(f"      Line {item['line_number']}: {item['pdb_line']}\n")
                
                if diff['extra_in_modern']:
                    f.write(f"\n  Extra in modern ({len(diff['extra_in_modern'])}):\n")
                    for item in diff['extra_in_modern'][:10]:  # Limit to 10
                        atom = item['atom']
                        f.write(f"    {atom.get('atom_name')} "
                               f"({atom.get('chain_id')}:{atom.get('residue_seq')})\n")
                        if item['pdb_line']:
                            f.write(f"      Line {item['line_number']}: {item['pdb_line']}\n")
                
                if diff['mismatched_fields']:
                    f.write(f"\n  Mismatched fields ({len(diff['mismatched_fields'])}):\n")
                    for item in diff['mismatched_fields'][:5]:  # Limit to 5
                        key = item['atom_key']
                        f.write(f"    Atom: {key[2]} ({key[0]}:{key[1]})\n")
                        for field, mismatch in item['mismatches'].items():
                            f.write(f"      {field}: Legacy={mismatch['legacy']}, "
                                   f"Modern={mismatch['modern']}\n")
                        if item['legacy_pdb_line']:
                            f.write(f"      Legacy line {item['legacy_line_number']}: "
                                   f"{item['legacy_pdb_line']}\n")
                        if item['modern_pdb_line']:
                            f.write(f"      Modern line {item['modern_line_number']}: "
                                   f"{item['modern_pdb_line']}\n")
            
            # Frame calculation differences
            if result['frame_differences']:
                diff = result['frame_differences']
                f.write("\nFRAME CALCULATION DIFFERENCES:\n")
                f.write("-" * 80 + "\n")
                
                if diff.get('missing_residues'):
                    f.write(f"\n  Missing in modern ({len(diff['missing_residues'])} residues):\n")
                    for rec in diff['missing_residues'][:10]:  # Limit to 10
                        chain_id = rec.get('chain_id', '?')
                        residue_seq = rec.get('residue_seq', 0)
                        residue_name = rec.get('residue_name', '?')
                        insertion = rec.get('insertion', ' ')
                        base_type = rec.get('base_type', '?')
                        num_matched = rec.get('num_matched_atoms', 0)
                        rms = rec.get('rms_fit', 0.0)
                        
                        f.write(f"    {residue_name} {chain_id}:{residue_seq}"
                               f"{insertion if insertion != ' ' else ''} (base: {base_type})\n")
                        f.write(f"      Matched atoms: {num_matched}, RMS: {rms:.6f}\n")
                        
                        # Get PDB lines for matched atoms if available
                        matched_atoms = rec.get('matched_atoms', [])
                        pdb_file_path = Path(result.get('pdb_file_path', ''))
                        if matched_atoms and pdb_file_path.exists():
                            pdb_file = pdb_file_path
                            f.write("      Matched atoms PDB lines:\n")
                            with open(pdb_file, 'r') as pdb_f:
                                lines = pdb_f.readlines()
                                for atom_name in matched_atoms[:6]:  # Limit to 6 atoms
                                    # Search for this atom in the PDB
                                    for line_num, line in enumerate(lines, 1):
                                        if line.startswith(('ATOM', 'HETATM')) and len(line) >= 16:
                                            line_chain = line[21] if len(line) > 21 else ' '
                                            try:
                                                line_seq = int(line[22:26].strip())
                                            except:
                                                continue
                                            line_insertion = line[26] if len(line) > 26 else ' '
                                            line_atom = line[12:16] if len(line) >= 16 else ''
                                            
                                            if (line_chain == chain_id and line_seq == residue_seq and
                                                line_insertion == insertion and line_atom == atom_name):
                                                f.write(f"        {atom_name} (line {line_num}): "
                                                       f"{line.rstrip()[:80]}\n")
                                                break
                
                if diff.get('mismatched_calculations'):
                    f.write(f"\n  Mismatched calculations ({len(diff['mismatched_calculations'])}):\n")
                    for item in diff['mismatched_calculations'][:5]:  # Limit to 5
                        key = item['residue_key']
                        f.write(f"    Residue: {key[0]}:{key[1]}{key[2] if key[2] != ' ' else ''}\n")
                        for mismatch_type, mismatch_data in item['mismatches'].items():
                            f.write(f"      {mismatch_type}: {mismatch_data}\n")
                        
                        # Show matched atoms from legacy
                        if item.get('legacy_matched_atoms'):
                            leg_atoms = item['legacy_matched_atoms']
                            f.write(f"      Legacy matched atoms: {leg_atoms[:10]}\n")  # Limit to 10
                        
                        # Show matched atoms from modern
                        if item.get('modern_matched_atoms'):
                            mod_atoms = item['modern_matched_atoms']
                            f.write(f"      Modern matched atoms: {mod_atoms[:10]}\n")  # Limit to 10
                        
                        if item.get('atom_pdb_lines'):
                            f.write("      Matched atoms PDB lines:\n")
                            for atom_info in item['atom_pdb_lines'][:10]:  # Limit to 10
                                f.write(f"        {atom_info['atom_name']} "
                                       f"(line {atom_info['line_number']}): "
                                       f"{atom_info['pdb_line'][:80]}\n")
            
            f.write("\n")

def main():
    project_root = Path(__file__).parent.parent.absolute()
    os.chdir(project_root)
    
    # Find executables
    legacy_exe = project_root / 'org' / 'build' / 'bin' / 'find_pair_analyze'
    if not legacy_exe.exists():
        legacy_exe = project_root / 'org' / 'find_pair_analyze'
    
    modern_exe = project_root / 'build' / 'tests' / 'integration' / 'test_json_generation_filtered'
    
    if not legacy_exe.exists():
        print(f"Error: Legacy executable not found: {legacy_exe}")
        print("Build it with: cd org && make")
        return 1
    
    if not modern_exe.exists():
        print(f"Error: Modern executable not found: {modern_exe}")
        print("Build it with: cd build && make test_json_generation_filtered")
        return 1
    
    print(f"Using legacy executable: {legacy_exe}")
    print(f"Using modern executable: {modern_exe}")
    
    # Load problematic PDBs
    print("Loading problematic PDBs...")
    problematic_pdbs = load_problematic_pdbs()
    print(f"Found {len(problematic_pdbs)} problematic PDBs")
    
    # Limit to first 20 for testing
    if '--all' not in sys.argv:
        problematic_pdbs = problematic_pdbs[:20]
        print(f"Limiting to first 20 PDBs (use --all for all)")
    
    # Step 1: Regenerate JSON files
    print("\nStep 1: Regenerating JSON files...")
    
    # Regenerate legacy JSON
    print("Regenerating legacy JSON...")
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = []
        for pdb_id in problematic_pdbs:
            future = executor.submit(regenerate_legacy_json, pdb_id, legacy_exe, project_root)
            futures.append(future)
        
        legacy_results = []
        for future in as_completed(futures):
            result = future.result()
            legacy_results.append(result)
            if result['status'] == 'success':
                print(f"  ✓ {result['pdb_id']}")
            else:
                print(f"  ✗ {result['pdb_id']}: {result.get('message', 'error')}")
    
    # Regenerate modern JSON (batch processing via gtest)
    print("\nRegenerating modern JSON (batch via test)...")
    print("Note: test_json_generation_filtered reads from docs/problematic_pdbs.txt")
    print("Running test...")
    
    # The test reads from docs/problematic_pdbs.txt, so we need to ensure that file has our PDBs
    # For now, let's just run it - the file should already have them
    try:
        result = subprocess.run(
            [str(modern_exe), '--gtest_filter=JsonGenerationFilteredTest.GenerateProblematicPdbs'],
            capture_output=True,
            text=True,
            timeout=1800,
            cwd=str(project_root)
        )
        print("Test output:")
        print(result.stdout[-1000:] if len(result.stdout) > 1000 else result.stdout)
        if result.stderr:
            print("Test errors:")
            print(result.stderr[-500:] if len(result.stderr) > 500 else result.stderr)
        
        # Check which JSON files were created
        modern_results = []
        for pdb_id in problematic_pdbs:
            json_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
            if json_file.exists():
                modern_results.append({'pdb_id': pdb_id, 'status': 'success', 'json_file': str(json_file)})
                print(f"  ✓ {pdb_id}")
            else:
                modern_results.append({'pdb_id': pdb_id, 'status': 'error', 'message': 'JSON file not created'})
                print(f"  ✗ {pdb_id}: JSON not created")
    except Exception as e:
        print(f"Error running modern JSON generation: {e}")
        modern_results = [{'pdb_id': pdb_id, 'status': 'error', 'message': str(e)} for pdb_id in problematic_pdbs]
    
    # Step 2: Compare and analyze
    print("\nStep 2: Comparing JSON files and extracting PDB lines...")
    all_results = []
    
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = []
        for pdb_id in problematic_pdbs:
            future = executor.submit(analyze_single_pdb, pdb_id, project_root, legacy_exe, modern_exe)
            futures.append(future)
        
        for future in as_completed(futures):
            result = future.result()
            all_results.append(result)
            status_symbol = '✓' if result['status'] == 'match' else '✗' if result['status'] == 'diff' else '!'
            print(f"  {status_symbol} {result['pdb_id']}: {result['status']}")
    
    # Step 3: Generate report
    print("\nStep 3: Generating detailed report...")
    output_file = project_root / 'docs' / 'detailed_pdb_comparison_report.txt'
    generate_report(all_results, output_file)
    print(f"Report written to: {output_file}")
    
    # Also save JSON version
    json_output = project_root / 'docs' / 'detailed_pdb_comparison_report.json'
    with open(json_output, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"JSON report written to: {json_output}")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())


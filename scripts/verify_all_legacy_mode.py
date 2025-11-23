#!/usr/bin/env python3
"""
Comprehensive legacy mode verification for all PDB files.
Regenerates JSON with --legacy flag and compares with legacy JSON.
"""

import json
import sys
import subprocess
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import time
from typing import Dict, List, Tuple

def get_all_pdb_ids(project_root: Path) -> List[str]:
    """Get all PDB IDs from data/pdb directory."""
    pdb_dir = project_root / 'data' / 'pdb'
    if not pdb_dir.exists():
        return []
    
    pdb_ids = []
    for pdb_file in pdb_dir.glob('*.pdb'):
        pdb_id = pdb_file.stem
        pdb_ids.append(pdb_id)
    
    return sorted(pdb_ids)

def regenerate_modern_json_legacy(pdb_id: str, executable_path: Path, project_root: Path) -> Dict:
    """Regenerate modern JSON with --legacy flag for a single PDB."""
    pdb_file = project_root / 'data' / 'pdb' / f'{pdb_id}.pdb'
    json_file = project_root / 'data' / 'json' / f'{pdb_id}_legacy.json'
    legacy_json_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
    
    if not pdb_file.exists():
        return {
            'pdb_id': pdb_id,
            'status': 'error',
            'message': f'PDB file not found',
            'skipped': True
        }
    
    if not legacy_json_file.exists():
        return {
            'pdb_id': pdb_id,
            'status': 'skipped',
            'message': 'No legacy JSON to compare',
            'skipped': True
        }
    
    try:
        # Run generate_modern_json with --legacy flag
        cmd = [
            str(executable_path),
            str(pdb_file),
            str(json_file),
            '--legacy'
        ]
        
        result = subprocess.run(
            cmd,
            cwd=str(project_root),
            capture_output=True,
            text=True,
            timeout=120
        )
        
        if result.returncode != 0:
            return {
                'pdb_id': pdb_id,
                'status': 'error',
                'message': f'Generation failed: {result.stderr[:200]}',
                'skipped': False
            }
        
        if not json_file.exists():
            return {
                'pdb_id': pdb_id,
                'status': 'error',
                'message': 'JSON file not created',
                'skipped': False
            }
        
        return {
            'pdb_id': pdb_id,
            'status': 'success',
            'message': 'Generated successfully',
            'skipped': False
        }
        
    except subprocess.TimeoutExpired:
        return {
            'pdb_id': pdb_id,
            'status': 'error',
            'message': 'Timeout',
            'skipped': False
        }
    except Exception as e:
        return {
            'pdb_id': pdb_id,
            'status': 'error',
            'message': f'Exception: {str(e)[:200]}',
            'skipped': False
        }

def compare_json_files(pdb_id: str, project_root: Path) -> Dict:
    """Compare legacy and modern JSON files for a single PDB."""
    legacy_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
    modern_file = project_root / 'data' / 'json' / f'{pdb_id}_legacy.json'
    
    if not legacy_file.exists() or not modern_file.exists():
        return {
            'pdb_id': pdb_id,
            'status': 'skipped',
            'total': 0,
            'exact_matches': 0,
            'set_matches': 0,
            'real_diffs': 0,
            'examples': []
        }
    
    try:
        with open(legacy_file) as f:
            legacy_data = json.load(f)
        with open(modern_file) as f:
            modern_data = json.load(f)
    except Exception as e:
        return {
            'pdb_id': pdb_id,
            'status': 'error',
            'message': f'JSON parse error: {str(e)[:200]}',
            'total': 0,
            'exact_matches': 0,
            'set_matches': 0,
            'real_diffs': 0
        }
    
    # Extract base_frame_calc entries
    leg_calcs = {
        (c.get('chain_id', ''), c.get('residue_seq', 0)): c 
        for c in legacy_data.get('calculations', []) 
        if c.get('type') == 'base_frame_calc'
    }
    mod_calcs = {
        (c.get('chain_id', ''), c.get('residue_seq', 0)): c 
        for c in modern_data.get('calculations', []) 
        if c.get('type') == 'base_frame_calc'
    }
    
    total = 0
    exact_matches = 0
    set_matches = 0
    real_diffs = 0
    examples = []
    
    for key, leg_calc in leg_calcs.items():
        if key not in mod_calcs:
            continue
        
        mod_calc = mod_calcs[key]
        total += 1
        
        leg_atoms = leg_calc.get('matched_atoms', [])
        mod_atoms = mod_calc.get('matched_atoms', [])
        
        if leg_atoms == mod_atoms:
            exact_matches += 1
            set_matches += 1
        elif set(leg_atoms) == set(mod_atoms):
            set_matches += 1
            if len(examples) < 3:
                examples.append({
                    'residue': f"{leg_calc.get('base_type', '?')} {key[0]}:{key[1]}",
                    'issue': 'order_only',
                    'legacy': leg_atoms,
                    'modern': mod_atoms
                })
        else:
            real_diffs += 1
            if len(examples) < 5:
                leg_only = sorted(set(leg_atoms) - set(mod_atoms))
                mod_only = sorted(set(mod_atoms) - set(leg_atoms))
                examples.append({
                    'residue': f"{leg_calc.get('base_type', '?')} {key[0]}:{key[1]}",
                    'issue': 'set_difference',
                    'legacy': leg_atoms,
                    'modern': mod_atoms,
                    'legacy_only': leg_only,
                    'modern_only': mod_only
                })
    
    return {
        'pdb_id': pdb_id,
        'status': 'success',
        'total': total,
        'exact_matches': exact_matches,
        'set_matches': set_matches,
        'real_diffs': real_diffs,
        'examples': examples
    }

def process_pdb(pdb_id: str, executable_path: Path, project_root: Path) -> Dict:
    """Process a single PDB: regenerate and compare."""
    # First regenerate
    regen_result = regenerate_modern_json_legacy(pdb_id, executable_path, project_root)
    
    if regen_result.get('skipped', False) or regen_result['status'] != 'success':
        return regen_result
    
    # Then compare
    compare_result = compare_json_files(pdb_id, project_root)
    compare_result['regen_status'] = regen_result['status']
    compare_result['regen_message'] = regen_result['message']
    
    return compare_result

def main():
    project_root = Path(__file__).parent.parent.absolute()
    
    # Find executable
    executable = project_root / 'build' / 'generate_modern_json'
    if not executable.exists():
        print(f"ERROR: Executable not found: {executable}")
        print("Please build the project first: cd build && make -j8")
        sys.exit(1)
    
    print("=" * 80)
    print("COMPREHENSIVE LEGACY MODE VERIFICATION")
    print("=" * 80)
    print(f"Project root: {project_root}")
    print(f"Executable: {executable}")
    print()
    
    # Get all PDB IDs
    all_pdbs = get_all_pdb_ids(project_root)
    print(f"Found {len(all_pdbs)} PDB files")
    print()
    
    # Filter to only PDBs that have legacy JSON
    legacy_dir = project_root / 'data' / 'json_legacy'
    pdbs_with_legacy = []
    for pdb_id in all_pdbs:
        if (legacy_dir / f'{pdb_id}.json').exists():
            pdbs_with_legacy.append(pdb_id)
    
    print(f"Found {len(pdbs_with_legacy)} PDBs with legacy JSON to verify")
    print()
    
    # Process all PDBs
    print("=" * 80)
    print("PROCESSING PDBs")
    print("=" * 80)
    print()
    
    num_workers = min(len(pdbs_with_legacy), 8)  # Limit concurrent processes
    start_time = time.time()
    results = []
    
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        future_to_pdb = {
            executor.submit(process_pdb, pdb_id, executable, project_root): pdb_id
            for pdb_id in pdbs_with_legacy
        }
        
        completed = 0
        for future in as_completed(future_to_pdb):
            pdb_id = future_to_pdb[future]
            try:
                result = future.result()
                results.append(result)
                completed += 1
                
                if result.get('skipped', False):
                    status = "⊘"
                elif result.get('status') == 'success':
                    if result.get('real_diffs', 0) == 0:
                        status = "✓"
                    else:
                        status = "✗"
                else:
                    status = "✗"
                
                total = result.get('total', 0)
                exact = result.get('exact_matches', 0)
                diffs = result.get('real_diffs', 0)
                
                print(f"{status} [{completed}/{len(pdbs_with_legacy)}] {pdb_id}: "
                      f"{total} residues, {exact} exact, {diffs} diffs")
                
            except Exception as e:
                results.append({
                    'pdb_id': pdb_id,
                    'status': 'error',
                    'message': f'Exception: {str(e)}',
                    'total': 0,
                    'exact_matches': 0,
                    'real_diffs': 0
                })
                completed += 1
                print(f"✗ [{completed}/{len(pdbs_with_legacy)}] {pdb_id}: Exception: {str(e)}")
    
    elapsed = time.time() - start_time
    
    # Generate summary
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    
    successful = [r for r in results if r.get('status') == 'success' and not r.get('skipped', False)]
    skipped = [r for r in results if r.get('skipped', False)]
    errors = [r for r in results if r.get('status') == 'error' and not r.get('skipped', False)]
    
    total_residues = sum(r.get('total', 0) for r in successful)
    total_exact = sum(r.get('exact_matches', 0) for r in successful)
    total_set_match = sum(r.get('set_matches', 0) for r in successful)
    total_diffs = sum(r.get('real_diffs', 0) for r in successful)
    
    perfect_pdbs = [r for r in successful if r.get('real_diffs', 0) == 0 and r.get('total', 0) > 0]
    diff_pdbs = [r for r in successful if r.get('real_diffs', 0) > 0]
    
    print(f"Processing time: {elapsed:.1f} seconds")
    print(f"PDBs processed: {len(successful)}")
    print(f"PDBs skipped: {len(skipped)}")
    print(f"PDBs with errors: {len(errors)}")
    print()
    print(f"Total residues compared: {total_residues}")
    print(f"Exact matches: {total_exact}/{total_residues} ({total_exact*100/max(total_residues,1):.1f}%)")
    print(f"Set matches: {total_set_match}/{total_residues} ({total_set_match*100/max(total_residues,1):.1f}%)")
    print(f"Real differences: {total_diffs}/{total_residues} ({total_diffs*100/max(total_residues,1):.1f}%)")
    print()
    print(f"Perfect PDBs (0 differences): {len(perfect_pdbs)}/{len(successful)}")
    print(f"PDBs with differences: {len(diff_pdbs)}/{len(successful)}")
    print()
    
    # Show examples of differences
    all_examples = []
    for r in diff_pdbs[:10]:  # Top 10 PDBs with differences
        for ex in r.get('examples', []):
            if ex.get('issue') == 'set_difference':
                all_examples.append({
                    'pdb_id': r['pdb_id'],
                    **ex
                })
    
    if all_examples:
        print("=" * 80)
        print("EXAMPLE DIFFERENCES")
        print("=" * 80)
        print()
        for i, ex in enumerate(all_examples[:10], 1):
            print(f"{i}. {ex['pdb_id']} - {ex['residue']}:")
            print(f"   Legacy: {ex['legacy']}")
            print(f"   Modern: {ex['modern']}")
            if ex.get('legacy_only'):
                print(f"   Only in legacy: {ex['legacy_only']}")
            if ex.get('modern_only'):
                print(f"   Only in modern: {ex['modern_only']}")
            print()
    
    # Generate report file
    report_file = project_root / 'docs' / 'LEGACY_MODE_VERIFICATION_REPORT.md'
    with open(report_file, 'w') as f:
        f.write("# Legacy Mode Verification Report\n\n")
        f.write(f"**Generated**: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("## Summary\n\n")
        f.write(f"- **Total PDBs processed**: {len(successful)}\n")
        f.write(f"- **PDBs skipped**: {len(skipped)}\n")
        f.write(f"- **PDBs with errors**: {len(errors)}\n")
        f.write(f"- **Total residues compared**: {total_residues}\n")
        f.write(f"- **Exact matches**: {total_exact}/{total_residues} ({total_exact*100/max(total_residues,1):.1f}%)\n")
        f.write(f"- **Set matches**: {total_set_match}/{total_residues} ({total_set_match*100/max(total_residues,1):.1f}%)\n")
        f.write(f"- **Real differences**: {total_diffs}/{total_residues} ({total_diffs*100/max(total_residues,1):.1f}%)\n")
        f.write(f"- **Perfect PDBs**: {len(perfect_pdbs)}/{len(successful)}\n")
        f.write(f"- **PDBs with differences**: {len(diff_pdbs)}/{len(successful)}\n\n")
        
        if perfect_pdbs:
            f.write("## Perfect Matches (0 differences)\n\n")
            f.write(f"Total: {len(perfect_pdbs)} PDBs\n\n")
            for r in perfect_pdbs[:20]:
                f.write(f"- {r['pdb_id']}: {r.get('total', 0)} residues, {r.get('exact_matches', 0)} exact matches\n")
            if len(perfect_pdbs) > 20:
                f.write(f"\n... and {len(perfect_pdbs) - 20} more\n")
            f.write("\n")
        
        if diff_pdbs:
            f.write("## PDBs with Differences\n\n")
            for r in diff_pdbs:
                f.write(f"### {r['pdb_id']}\n\n")
                f.write(f"- Total residues: {r.get('total', 0)}\n")
                f.write(f"- Exact matches: {r.get('exact_matches', 0)}\n")
                f.write(f"- Set matches: {r.get('set_matches', 0)}\n")
                f.write(f"- Real differences: {r.get('real_diffs', 0)}\n\n")
                if r.get('examples'):
                    f.write("Examples:\n")
                    for ex in r['examples']:
                        f.write(f"- {ex['residue']}: {ex.get('issue', 'unknown')}\n")
                f.write("\n")
    
    print(f"Report saved to: {report_file}")
    print()
    
    # Final verdict
    print("=" * 80)
    print("VERDICT")
    print("=" * 80)
    
    if total_diffs == 0:
        print("✅ PERFECT MATCH! All residues match exactly!")
        print(f"   {len(perfect_pdbs)} PDBs verified with 0 differences")
        sys.exit(0)
    else:
        print(f"⚠️  {total_diffs} residues have differences across {len(diff_pdbs)} PDBs")
        print(f"   {len(perfect_pdbs)} PDBs match perfectly")
        if total_diffs * 100 / max(total_residues, 1) < 1.0:
            print("   (Less than 1% difference rate - very good!)")
        sys.exit(0 if total_diffs * 100 / max(total_residues, 1) < 5.0 else 1)

if __name__ == "__main__":
    main()


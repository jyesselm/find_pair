#!/usr/bin/env python3
"""
Regenerate legacy JSON files for problematic PDBs using the original C code.
This ensures the legacy JSON has the removed atom tracking information.
"""

import os
import sys
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import time

def load_problematic_pdbs():
    """Load problematic PDB list."""
    problem_file = Path('docs/problematic_pdbs.txt')
    if not problem_file.exists():
        return []
    
    pdbs = []
    with open(problem_file, 'r') as f:
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
    # Convert strings back to Path objects if needed (for ProcessPoolExecutor)
    if isinstance(project_root, str):
        project_root = Path(project_root)
    if isinstance(executable_path, str):
        executable_path = Path(executable_path)
    
    pdb_file = project_root / 'data' / 'pdb' / f'{pdb_id}.pdb'
    json_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
    
    if not pdb_file.exists():
        return {
            'pdb_id': pdb_id,
            'status': 'error',
            'message': f'PDB file not found: {pdb_file}'
        }
    
    try:
        # Run find_pair_analyze on the PDB
        # This will generate JSON with removed atom tracking
        # Use absolute path for PDB file
        cmd = [str(executable_path), str(pdb_file)]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout per PDB
            cwd=str(project_root / 'org')
        )
        
        # Check if JSON file was created
        if json_file.exists():
            return {
                'pdb_id': pdb_id,
                'status': 'success',
                'message': 'JSON regenerated successfully'
            }
        else:
            return {
                'pdb_id': pdb_id,
                'status': 'error',
                'message': f'JSON file not created: {json_file}',
                'stderr': result.stderr[-500:] if result.stderr else ''
            }
    
    except subprocess.TimeoutExpired:
        return {
            'pdb_id': pdb_id,
            'status': 'error',
            'message': 'Timeout after 5 minutes'
        }
    except Exception as e:
        return {
            'pdb_id': pdb_id,
            'status': 'error',
            'message': str(e)
        }

def main():
    # Find executable - use absolute path from project root
    project_root = Path(__file__).parent.parent.absolute()
    executable = project_root / 'org' / 'build' / 'bin' / 'find_pair_analyze'
    
    if not executable.exists():
        # Try alternative path
        executable = Path('/Users/jyesselman2/Dropbox/2_code/cpp/find_pair_2/org/build/bin/find_pair_analyze')
        if not executable.exists():
            print(f"Error: find_pair_analyze executable not found")
            print(f"  Tried: {project_root / 'org' / 'build' / 'bin' / 'find_pair_analyze'}")
            print(f"  Tried: {executable}")
            sys.exit(1)
    
    print(f"Using executable: {executable}")
    print(f"Project root: {project_root}")
    print()
    
    # Load problematic PDBs
    problematic_pdbs = load_problematic_pdbs()
    
    if not problematic_pdbs:
        print("Error: No problematic PDBs found in docs/problematic_pdbs.txt")
        sys.exit(1)
    
    print(f"Found {len(problematic_pdbs)} problematic PDBs to regenerate")
    print()
    
    # Filter to only PDBs that exist
    existing_pdbs = []
    for pdb_id in problematic_pdbs:
        pdb_file = project_root / 'data' / 'pdb' / f'{pdb_id}.pdb'
        if pdb_file.exists():
            existing_pdbs.append(pdb_id)
        else:
            print(f"Warning: PDB file not found: {pdb_file}")
    
    print(f"Found {len(existing_pdbs)} PDB files to process")
    print()
    
    # Process in parallel using processes (better for subprocess calls)
    num_workers = min(len(existing_pdbs), os.cpu_count() or 4)
    print(f"Processing with {num_workers} processes...")
    print()
    
    start_time = time.time()
    results = []
    
    # Use ProcessPoolExecutor for better parallelism with subprocess calls
    # Convert Path objects to strings for pickle compatibility
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit all tasks
        # Convert Path objects to strings for pickle compatibility with ProcessPoolExecutor
        future_to_pdb = {
            executor.submit(regenerate_legacy_json, pdb_id, str(executable), str(project_root)): pdb_id
            for pdb_id in existing_pdbs
        }
        
        completed = 0
        for future in as_completed(future_to_pdb):
            pdb_id = future_to_pdb[future]
            try:
                result = future.result()
                results.append(result)
                completed += 1
                
                if result['status'] == 'success':
                    status = "✓"
                else:
                    status = "✗"
                
                print(f"{status} [{completed}/{len(existing_pdbs)}] {pdb_id}: {result['message']}")
                
            except Exception as e:
                results.append({
                    'pdb_id': pdb_id,
                    'status': 'error',
                    'message': f'Exception: {str(e)}'
                })
                completed += 1
                print(f"✗ [{completed}/{len(existing_pdbs)}] {pdb_id}: Exception: {str(e)}")
    
    elapsed = time.time() - start_time
    
    # Summary
    print()
    print("=" * 80)
    print("REGENERATION SUMMARY")
    print("=" * 80)
    print(f"Total processed: {len(results)}")
    print(f"Success: {sum(1 for r in results if r['status'] == 'success')}")
    print(f"Errors: {sum(1 for r in results if r['status'] == 'error')}")
    print(f"Time elapsed: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
    print()
    
    # Show errors
    errors = [r for r in results if r['status'] == 'error']
    if errors:
        print("Errors:")
        for r in errors[:10]:
            print(f"  {r['pdb_id']}: {r['message']}")
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more errors")
    
    print()
    print("Legacy JSON files regenerated with removed atom tracking!")
    print("You can now compare with modern JSON files.")

if __name__ == '__main__':
    main()


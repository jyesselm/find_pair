#!/usr/bin/env python3
"""
Refresh all JSON files (both legacy and modern) with validation.
Removes corrupted or empty JSON files instead of saving them.
"""

import os
import sys
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import time
from typing import Dict, Optional

# Import reusable validation from lib
from x3dna_json_compare import JsonValidator

def regenerate_legacy_json(pdb_id: str, executable_path: str, project_root: str) -> Dict:
    """Regenerate legacy JSON for a single PDB with validation."""
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
        cmd = [str(executable_path), str(pdb_file)]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout per PDB
            cwd=str(project_root / 'org')
        )
        
        # Check if JSON file was created
        if not json_file.exists():
            return {
                'pdb_id': pdb_id,
                'status': 'error',
                'message': f'JSON file not created: {json_file}',
                'stderr': result.stderr[-500:] if result.stderr else ''
            }
        
        # Validate the generated JSON
        is_valid, error_msg, was_removed = JsonValidator.validate_and_clean(json_file, remove_invalid=True)
        
        if not is_valid:
            return {
                'pdb_id': pdb_id,
                'status': 'error',
                'message': f'Invalid JSON removed: {error_msg}',
                'removed': was_removed
            }
        
        return {
            'pdb_id': pdb_id,
            'status': 'success',
            'message': 'JSON regenerated and validated successfully'
        }
    
    except subprocess.TimeoutExpired:
        # Remove any partial JSON file
        if json_file.exists():
            JsonValidator.remove_invalid_file(json_file, "Timeout during generation")
        return {
            'pdb_id': pdb_id,
            'status': 'error',
            'message': 'Timeout after 5 minutes',
            'removed': True
        }
    except Exception as e:
        # Remove any partial JSON file
        if json_file.exists():
            JsonValidator.remove_invalid_file(json_file, f"Exception: {str(e)}")
        return {
            'pdb_id': pdb_id,
            'status': 'error',
            'message': str(e),
            'removed': True
        }

def regenerate_modern_json(pdb_id: str, executable_path: Path, project_root: Path) -> Dict:
    """Regenerate modern JSON for a single PDB using test_json_generation."""
    pdb_file = project_root / 'data' / 'pdb' / f'{pdb_id}.pdb'
    json_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
    
    if not pdb_file.exists():
        return {
            'pdb_id': pdb_id,
            'status': 'error',
            'message': f'PDB file not found: {pdb_file}'
        }
    
    try:
        # Use test_json_generation to generate JSON for a specific PDB
        # We'll use test_single_pdb as a helper, or create a filtered test
        # For now, let's use the test_json_generation executable with filtering
        
        # Create a temporary input file with just this PDB
        # Actually, test_json_generation processes all PDBs, so we need a different approach
        # Let's use test_single_pdb if available, or we need to modify the approach
        
        # Check if test_single_pdb exists
        test_single_pdb = project_root / 'build' / 'tests' / 'integration' / 'test_single_pdb'
        if test_single_pdb.exists():
            # Use test_single_pdb to generate JSON for this PDB
            cmd = [str(test_single_pdb), pdb_id]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
                cwd=str(project_root / 'build')
            )
            
            if result.returncode != 0:
                # Remove any partial JSON
                if json_file.exists():
                    remove_invalid_json(json_file, f"test_single_pdb failed: {result.stderr[-200:]}")
                return {
                    'pdb_id': pdb_id,
                    'status': 'error',
                    'message': f'test_single_pdb failed: {result.stderr[-200:] if result.stderr else "Unknown error"}',
                    'removed': True
                }
        else:
            # Fallback: use test_json_generation (will process all, but we'll validate this one)
            # This is less efficient but works
            test_json_gen = project_root / 'build' / 'tests' / 'integration' / 'test_json_generation'
            if not test_json_gen.exists():
                return {
                    'pdb_id': pdb_id,
                    'status': 'error',
                    'message': 'test_json_generation executable not found'
                }
            
            # Run test_json_generation (it will process all, but we only care about this one)
            cmd = [str(test_json_gen), '--gtest_filter=JsonGenerationTest.GenerateAllJsonFiles']
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,  # Longer timeout since it processes all
                cwd=str(project_root / 'build')
            )
        
        # Validate the generated JSON
        if json_file.exists():
            is_valid, error_msg, was_removed = JsonValidator.validate_and_clean(json_file, remove_invalid=True)
            
            if not is_valid:
                return {
                    'pdb_id': pdb_id,
                    'status': 'error',
                    'message': f'Invalid JSON removed: {error_msg}',
                    'removed': was_removed
                }
            
            return {
                'pdb_id': pdb_id,
                'status': 'success',
                'message': 'JSON regenerated and validated successfully'
            }
        else:
            return {
                'pdb_id': pdb_id,
                'status': 'error',
                'message': f'JSON file not created: {json_file}'
            }
    
    except subprocess.TimeoutExpired:
        if json_file.exists():
            JsonValidator.remove_invalid_file(json_file, "Timeout during generation")
        return {
            'pdb_id': pdb_id,
            'status': 'error',
            'message': 'Timeout during generation',
            'removed': True
        }
    except Exception as e:
        if json_file.exists():
            JsonValidator.remove_invalid_file(json_file, f"Exception: {str(e)}")
        return {
            'pdb_id': pdb_id,
            'status': 'error',
            'message': str(e),
            'removed': True
        }

def validate_existing_json(json_file: Path) -> Dict:
    """Validate an existing JSON file and remove if invalid."""
    pdb_id = json_file.stem
    
    is_valid, error_msg, was_removed = JsonValidator.validate_and_clean(json_file, remove_invalid=True)
    
    if not is_valid:
        return {
            'pdb_id': pdb_id,
            'status': 'removed' if was_removed else 'error',
            'message': f'Invalid JSON removed: {error_msg}',
            'removed': was_removed
        }
    
    return {
        'pdb_id': pdb_id,
        'status': 'valid',
        'message': 'JSON is valid'
    }

def get_all_pdb_ids(project_root: Path) -> list:
    """Get all PDB IDs from the data/pdb directory."""
    pdb_dir = project_root / 'data' / 'pdb'
    if not pdb_dir.exists():
        return []
    
    pdb_ids = []
    for pdb_file in pdb_dir.glob('*.pdb'):
        pdb_ids.append(pdb_file.stem)
    
    return sorted(pdb_ids)

def main():
    project_root = Path(__file__).parent.parent.absolute()
    
    # Find executables
    legacy_exe = project_root / 'org' / 'build' / 'bin' / 'find_pair_analyze'
    if not legacy_exe.exists():
        legacy_exe = Path('/Users/jyesselman2/Dropbox/2_code/cpp/find_pair_2/org/build/bin/find_pair_analyze')
        if not legacy_exe.exists():
            print(f"Warning: Legacy executable not found: {legacy_exe}")
            print("  Legacy JSON regeneration will be skipped")
            legacy_exe = None
    
    modern_exe = project_root / 'build' / 'tests' / 'integration' / 'test_single_pdb'
    if not modern_exe.exists():
        modern_exe = project_root / 'build' / 'tests' / 'integration' / 'test_json_generation'
        if not modern_exe.exists():
            print(f"Warning: Modern executable not found")
            print("  Modern JSON regeneration will be skipped")
            modern_exe = None
    
    print("=" * 80)
    print("JSON FILE REFRESH")
    print("=" * 80)
    print(f"Project root: {project_root}")
    print(f"Legacy executable: {legacy_exe if legacy_exe else 'NOT FOUND'}")
    print(f"Modern executable: {modern_exe if modern_exe else 'NOT FOUND'}")
    print()
    
    # Get all PDB IDs
    all_pdbs = get_all_pdb_ids(project_root)
    print(f"Found {len(all_pdbs)} PDB files")
    print()
    
    # Parse command line arguments
    validate_only = '--validate-only' in sys.argv
    legacy_only = '--legacy-only' in sys.argv
    modern_only = '--modern-only' in sys.argv
    
    if validate_only:
        print("VALIDATION MODE: Validating existing JSON files only")
        print()
        
        # Validate legacy JSON files
        legacy_dir = project_root / 'data' / 'json_legacy'
        if legacy_dir.exists():
            print("Validating legacy JSON files...")
            legacy_json_files = list(legacy_dir.glob('*.json'))
            print(f"Found {len(legacy_json_files)} legacy JSON files")
            
            num_workers = min(len(legacy_json_files), os.cpu_count() or 4)
            with ThreadPoolExecutor(max_workers=num_workers) as executor:
                future_to_file = {
                    executor.submit(validate_existing_json, json_file): json_file
                    for json_file in legacy_json_files
                }
                
                removed_count = 0
                valid_count = 0
                for future in as_completed(future_to_file):
                    result = future.result()
                    if result['status'] == 'removed':
                        removed_count += 1
                        print(f"  ✗ Removed: {result['pdb_id']} - {result['message']}")
                    elif result['status'] == 'valid':
                        valid_count += 1
                
                print(f"\nLegacy JSON: {valid_count} valid, {removed_count} removed")
        
        # Validate modern JSON files
        modern_dir = project_root / 'data' / 'json'
        if modern_dir.exists():
            print("\nValidating modern JSON files...")
            modern_json_files = list(modern_dir.glob('*.json'))
            print(f"Found {len(modern_json_files)} modern JSON files")
            
            num_workers = min(len(modern_json_files), os.cpu_count() or 4)
            with ThreadPoolExecutor(max_workers=num_workers) as executor:
                future_to_file = {
                    executor.submit(validate_existing_json, json_file): json_file
                    for json_file in modern_json_files
                }
                
                removed_count = 0
                valid_count = 0
                for future in as_completed(future_to_file):
                    result = future.result()
                    if result['status'] == 'removed':
                        removed_count += 1
                        print(f"  ✗ Removed: {result['pdb_id']} - {result['message']}")
                    elif result['status'] == 'valid':
                        valid_count += 1
                
                print(f"\nModern JSON: {valid_count} valid, {removed_count} removed")
        
        return
    
    # Regeneration mode
    print("REGENERATION MODE: Regenerating all JSON files with validation")
    print()
    
    # Regenerate legacy JSON
    if not modern_only and legacy_exe:
        print("=" * 80)
        print("REGENERATING LEGACY JSON FILES")
        print("=" * 80)
        print(f"Processing {len(all_pdbs)} PDBs...")
        print()
        
        num_workers = min(len(all_pdbs), os.cpu_count() or 4)
        start_time = time.time()
        results = []
        
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            future_to_pdb = {
                executor.submit(regenerate_legacy_json, pdb_id, str(legacy_exe), str(project_root)): pdb_id
                for pdb_id in all_pdbs
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
                    
                    print(f"{status} [{completed}/{len(all_pdbs)}] {pdb_id}: {result['message']}")
                    
                except Exception as e:
                    results.append({
                        'pdb_id': pdb_id,
                        'status': 'error',
                        'message': f'Exception: {str(e)}'
                    })
                    completed += 1
                    print(f"✗ [{completed}/{len(all_pdbs)}] {pdb_id}: Exception: {str(e)}")
        
        elapsed = time.time() - start_time
        success_count = sum(1 for r in results if r['status'] == 'success')
        error_count = sum(1 for r in results if r['status'] == 'error')
        removed_count = sum(1 for r in results if r.get('removed', False))
        
        print()
        print("=" * 80)
        print("LEGACY JSON REGENERATION SUMMARY")
        print("=" * 80)
        print(f"Total processed: {len(results)}")
        print(f"Success: {success_count}")
        print(f"Errors: {error_count}")
        print(f"Removed (invalid): {removed_count}")
        print(f"Time elapsed: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
        print()
    
    # Regenerate modern JSON
    if not legacy_only and modern_exe:
        print("=" * 80)
        print("REGENERATING MODERN JSON FILES")
        print("=" * 80)
        print(f"Processing {len(all_pdbs)} PDBs...")
        print()
        print("Note: Modern JSON generation processes all PDBs at once")
        print("      Individual validation will be performed after generation")
        print()
        
        # For modern JSON, we'll use test_json_generation which processes all at once
        # Then validate each file individually
        test_json_gen = project_root / 'build' / 'tests' / 'integration' / 'test_json_generation'
        if test_json_gen.exists():
            print("Running test_json_generation...")
            start_time = time.time()
            
            try:
                result = subprocess.run(
                    [str(test_json_gen), '--gtest_filter=JsonGenerationTest.GenerateAllJsonFiles'],
                    capture_output=True,
                    text=True,
                    timeout=3600,  # 1 hour timeout
                    cwd=str(project_root / 'build')
                )
                
                if result.returncode != 0:
                    print(f"Error running test_json_generation:")
                    print(result.stderr[-1000:])
                else:
                    print("✓ test_json_generation completed")
            except subprocess.TimeoutExpired:
                print("✗ test_json_generation timed out after 1 hour")
            except Exception as e:
                print(f"✗ Error running test_json_generation: {e}")
            
            elapsed = time.time() - start_time
            print(f"Generation time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
            print()
            
            # Now validate all generated files
            print("Validating generated modern JSON files...")
            modern_dir = project_root / 'data' / 'json'
            if modern_dir.exists():
                modern_json_files = list(modern_dir.glob('*.json'))
                print(f"Found {len(modern_json_files)} modern JSON files to validate")
                
                num_workers = min(len(modern_json_files), os.cpu_count() or 4)
                validation_results = []
                
                with ThreadPoolExecutor(max_workers=num_workers) as executor:
                    future_to_file = {
                        executor.submit(validate_existing_json, json_file): json_file
                        for json_file in modern_json_files
                    }
                    
                    completed = 0
                    for future in as_completed(future_to_file):
                        result = future.result()
                        validation_results.append(result)
                        completed += 1
                        
                        if result['status'] == 'removed':
                            print(f"  ✗ [{completed}/{len(modern_json_files)}] Removed: {result['pdb_id']} - {result['message']}")
                        elif result['status'] == 'valid':
                            if completed % 100 == 0:
                                print(f"  ✓ [{completed}/{len(modern_json_files)}] Validated...")
                
                valid_count = sum(1 for r in validation_results if r['status'] == 'valid')
                removed_count = sum(1 for r in validation_results if r['status'] == 'removed')
                
                print()
                print("=" * 80)
                print("MODERN JSON REGENERATION SUMMARY")
                print("=" * 80)
                print(f"Total files: {len(validation_results)}")
                print(f"Valid: {valid_count}")
                print(f"Removed (invalid): {removed_count}")
                print()
    
    print("=" * 80)
    print("JSON REFRESH COMPLETE")
    print("=" * 80)

if __name__ == '__main__':
    main()


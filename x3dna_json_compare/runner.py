"""
Unified validation runner.

Single implementation for all validation types - replaces the many validate_*.py scripts.

Features:
- Checkpoint/resume: Save progress to JSON, resume from where you left off
- Clean on match: Delete modern JSON files that match legacy (save disk space)
- Stage-specific: Only validate/generate specific stages
"""

import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Dict, Any, Set
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from .json_comparison import JsonComparator
from .output import OutputFormatter, ValidationSummary
from .pdb_list import get_pdb_list


# Linear stages in legacy execution order (12 stages)
# Each stage maps to a comparison flag and the JSON directories it validates
STAGE_CONFIG = {
    # Stage 1: Atom parsing
    '1': {'name': 'pdb_atoms', 'compare_flag': 'compare_atoms', 'json_dirs': ['pdb_atoms']},
    'pdb_atoms': {'name': 'pdb_atoms', 'compare_flag': 'compare_atoms', 'json_dirs': ['pdb_atoms']},
    
    # Stage 2: Residue indices
    '2': {'name': 'residue_indices', 'compare_flag': 'compare_residue_indices', 'json_dirs': ['residue_indices']},
    'residue_indices': {'name': 'residue_indices', 'compare_flag': 'compare_residue_indices', 'json_dirs': ['residue_indices']},
    
    # Stage 3: Base frame calculation
    '3': {'name': 'base_frame_calc', 'compare_flag': 'compare_frames', 'json_dirs': ['base_frame_calc']},
    'base_frame_calc': {'name': 'base_frame_calc', 'compare_flag': 'compare_frames', 'json_dirs': ['base_frame_calc']},
    
    # Stage 4: LS fitting
    '4': {'name': 'ls_fitting', 'compare_flag': 'compare_frames', 'json_dirs': ['ls_fitting']},
    'ls_fitting': {'name': 'ls_fitting', 'compare_flag': 'compare_frames', 'json_dirs': ['ls_fitting']},
    
    # Stage 5: Frame calculation  
    '5': {'name': 'frame_calc', 'compare_flag': 'compare_frames', 'json_dirs': ['frame_calc']},
    'frame_calc': {'name': 'frame_calc', 'compare_flag': 'compare_frames', 'json_dirs': ['frame_calc']},
    
    # Stage 6: Pair validation
    '6': {'name': 'pair_validation', 'compare_flag': 'compare_pairs', 'json_dirs': ['pair_validation']},
    'pair_validation': {'name': 'pair_validation', 'compare_flag': 'compare_pairs', 'json_dirs': ['pair_validation']},
    
    # Stage 7: Distance checks
    '7': {'name': 'distance_checks', 'compare_flag': 'compare_pairs', 'json_dirs': ['distance_checks']},
    'distance_checks': {'name': 'distance_checks', 'compare_flag': 'compare_pairs', 'json_dirs': ['distance_checks']},
    
    # Stage 8: H-bond list
    '8': {'name': 'hbond_list', 'compare_flag': 'compare_hbond_list', 'json_dirs': ['hbond_list']},
    'hbond_list': {'name': 'hbond_list', 'compare_flag': 'compare_hbond_list', 'json_dirs': ['hbond_list']},
    
    # Stage 9: Base pair
    '9': {'name': 'base_pair', 'compare_flag': 'compare_pairs', 'json_dirs': ['base_pair']},
    'base_pair': {'name': 'base_pair', 'compare_flag': 'compare_pairs', 'json_dirs': ['base_pair']},
    
    # Stage 10: Find bestpair selection
    '10': {'name': 'find_bestpair_selection', 'compare_flag': 'compare_pairs', 'json_dirs': ['find_bestpair_selection']},
    'find_bestpair_selection': {'name': 'find_bestpair_selection', 'compare_flag': 'compare_pairs', 'json_dirs': ['find_bestpair_selection']},
    
    # Stage 11: Step parameters
    '11': {'name': 'bpstep_params', 'compare_flag': 'compare_steps', 'json_dirs': ['bpstep_params']},
    'bpstep_params': {'name': 'bpstep_params', 'compare_flag': 'compare_steps', 'json_dirs': ['bpstep_params']},
    
    # Stage 12: Helical parameters  
    '12': {'name': 'helical_params', 'compare_flag': 'compare_steps', 'json_dirs': ['helical_params']},
    'helical_params': {'name': 'helical_params', 'compare_flag': 'compare_steps', 'json_dirs': ['helical_params']},
}

# Grouped stages (for backward compatibility)
STAGE_GROUPS = {
    'atoms': ['1'],
    'frames': ['3', '4', '5'],
    'hbonds': ['8'],
    'pairs': ['6', '7', '9', '10'],
    'steps': ['11', '12'],
}

# Stage to JSON subdirectory mapping (for backward compatibility)
STAGE_JSON_DIRS = {
    'atoms': ['pdb_atoms'],
    'frames': ['base_frame_calc', 'frame_calc', 'ls_fitting'],
    'hbonds': ['hbond_list'],
    'pairs': ['base_pair', 'pair_validation', 'find_bestpair_selection', 'distance_checks'],
    'steps': ['bpstep_params', 'helical_params'],
    'residue_indices': ['residue_indices'],
}


def _validate_single_pdb(args) -> Dict[str, Any]:
    """Worker function for parallel validation (must be top-level for pickling)."""
    pdb_id, project_root_str, comparator_kwargs = args
    project_root = Path(project_root_str)
    
    # Create comparator for this process
    comparator = JsonComparator(**comparator_kwargs)
    
    legacy_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
    modern_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
    pdb_file = project_root / 'data' / 'pdb' / f'{pdb_id}.pdb'
    
    result = comparator.compare_files(legacy_file, modern_file, pdb_file, pdb_id)
    
    # Convert to dict for pickling
    return {
        'pdb_id': pdb_id,
        'status': result.status,
        'has_differences': result.has_differences(),
        'errors': result.errors,
        'atom_comparison': _summarize_atom_comparison(result.atom_comparison),
        'frame_comparison': _summarize_frame_comparison(result.frame_comparison),
        'hbond_comparison': _summarize_hbond_comparison(result.hbond_list_comparison),
    }


def _summarize_atom_comparison(ac) -> Optional[Dict]:
    """Summarize atom comparison for serialization."""
    if ac is None:
        return None
    return {
        'total_legacy': ac.total_legacy,
        'total_modern': ac.total_modern,
        'missing_in_modern': len(ac.missing_in_modern),
        'extra_in_modern': len(ac.extra_in_modern),
        'mismatched_fields': len(ac.mismatched_fields),
    }


def _summarize_frame_comparison(fc) -> Optional[Dict]:
    """Summarize frame comparison for serialization."""
    if fc is None:
        return None
    # Compute matched as total - missing - mismatched
    matched = fc.total_legacy - len(fc.missing_residues) - len(fc.mismatched_calculations)
    return {
        'total_legacy': fc.total_legacy,
        'matched': matched,
        'missing_residues': len(fc.missing_residues),
        'mismatched_calculations': len(fc.mismatched_calculations),
    }


def _summarize_hbond_comparison(hc) -> Optional[Dict]:
    """Summarize hbond comparison for serialization."""
    if hc is None:
        return None
    return {
        'total_legacy': hc.total_legacy,
        'total_modern': hc.total_modern,
        'missing_in_modern': len(hc.missing_in_modern),
        'extra_in_modern': len(hc.extra_in_modern),
    }


class ValidationRunner:
    """Core validation runner - single implementation for all validation types.
    
    Features:
    - Checkpoint/resume: Save progress to JSON file, skip already-passed PDBs
    - Clean on match: Delete modern JSON files that match legacy
    - Stage-specific validation
    """
    
    def __init__(self,
                 project_root: Path = None,
                 stages: List[str] = None,
                 workers: int = None,
                 quiet: bool = False,
                 verbose: bool = False,
                 stop_on_first: bool = False,
                 document_differences: bool = False,
                 diff_file: Path = None,
                 checkpoint_file: Path = None,
                 resume: bool = False,
                 clean_on_match: bool = False,
                 **comparator_kwargs):
        """Initialize validation runner.
        
        Args:
            project_root: Project root directory (default: current working directory)
            stages: List of stages to validate (atoms, frames, hbonds, pairs, steps)
            workers: Number of parallel workers (default: CPU count)
            quiet: Suppress all output except exit code
            verbose: Show per-PDB results
            stop_on_first: Stop at first failure
            document_differences: Save differences to file
            diff_file: Custom path for differences file
            checkpoint_file: Path to checkpoint file for resume support
            resume: If True, skip PDBs that already passed in checkpoint
            clean_on_match: If True, delete modern JSON files that match legacy
            **comparator_kwargs: Additional arguments for JsonComparator
        """
        self.project_root = Path(project_root) if project_root else Path.cwd()
        self.stages = stages or ['all']
        self.workers = workers or min(multiprocessing.cpu_count(), 10)
        self.quiet = quiet
        self.verbose = verbose
        self.stop_on_first = stop_on_first
        self.document_differences = document_differences
        self.diff_file = diff_file
        self.checkpoint_file = checkpoint_file
        self.resume = resume
        self.clean_on_match = clean_on_match
        
        # Build comparator kwargs from stages
        self.comparator_kwargs = self._build_comparator_kwargs(stages, comparator_kwargs)
        self.output = OutputFormatter(quiet=quiet, verbose=verbose)
        
        # Load checkpoint if resuming
        self.checkpoint_data = self._load_checkpoint() if resume and checkpoint_file else None
    
    def _build_comparator_kwargs(self, stages: List[str], extra_kwargs: Dict) -> Dict:
        """Build JsonComparator kwargs from stage list."""
        # Default: all off, cache disabled (cache doesn't track which stages were compared)
        kwargs = {
            'compare_atoms': False,
            'compare_frames': False,
            'compare_hbond_list': False,
            'compare_pairs': False,
            'compare_steps': False,
            'compare_residue_indices': False,
            'enable_cache': False,  # Disable cache - it doesn't account for comparison flags
        }
        
        if not stages or 'all' in stages:
            # Enable all comparison flags (but keep enable_cache=False)
            for k in kwargs:
                if k != 'enable_cache':
                    kwargs[k] = True
        else:
            # Expand stage groups and individual stages
            expanded_stages = set()
            for stage in stages:
                # Check if it's a group name
                if stage in STAGE_GROUPS:
                    expanded_stages.update(STAGE_GROUPS[stage])
                # Check if it's an individual stage
                elif stage in STAGE_CONFIG:
                    expanded_stages.add(stage)
                # Legacy support for old stage names
                elif stage in ['atoms', 'frames', 'hbonds', 'pairs', 'steps', 'residue_indices']:
                    if stage in STAGE_GROUPS:
                        expanded_stages.update(STAGE_GROUPS[stage])
                    else:
                        # Map old stage name to new config
                        stage_map = {
                            'atoms': '1',
                            'residue_indices': '2',
                        }
                        if stage in stage_map:
                            expanded_stages.add(stage_map[stage])
            
            # Enable comparison flags for all expanded stages
            for stage in expanded_stages:
                if stage in STAGE_CONFIG:
                    compare_flag = STAGE_CONFIG[stage]['compare_flag']
                    kwargs[compare_flag] = True
        
        # Override with any explicit kwargs
        kwargs.update(extra_kwargs)
        return kwargs
    
    def _load_checkpoint(self) -> Optional[Dict]:
        """Load checkpoint file if it exists."""
        if not self.checkpoint_file:
            return None
        checkpoint_path = Path(self.checkpoint_file)
        if not checkpoint_path.exists():
            return None
        try:
            with open(checkpoint_path, 'r') as f:
                return json.load(f)
        except Exception:
            return None
    
    def _save_checkpoint(self, results: Dict[str, Dict]):
        """Save checkpoint to file."""
        if not self.checkpoint_file:
            return
        
        checkpoint_path = Path(self.checkpoint_file)
        checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
        
        checkpoint = {
            'timestamp': datetime.now().isoformat(),
            'stages': self.stages,
            'results': results
        }
        
        with open(checkpoint_path, 'w') as f:
            json.dump(checkpoint, f, indent=2, default=str)
    
    def _get_passed_pdbs(self) -> Set[str]:
        """Get set of PDBs that already passed (from checkpoint)."""
        if not self.checkpoint_data:
            return set()
        
        results = self.checkpoint_data.get('results', {})
        return {pdb_id for pdb_id, result in results.items() 
                if result.get('status') == 'match'}
    
    def _clean_modern_json(self, pdb_id: str):
        """Delete modern JSON files for a PDB that matched."""
        if not self.clean_on_match:
            return
        
        # Determine which directories to clean based on stages
        dirs_to_clean = set()
        if 'all' in self.stages:
            for stage_info in STAGE_CONFIG.values():
                dirs_to_clean.update(stage_info['json_dirs'])
        else:
            for stage in self.stages:
                # Check individual stage config
                if stage in STAGE_CONFIG:
                    dirs_to_clean.update(STAGE_CONFIG[stage]['json_dirs'])
                # Check stage groups
                elif stage in STAGE_GROUPS:
                    for sub_stage in STAGE_GROUPS[stage]:
                        if sub_stage in STAGE_CONFIG:
                            dirs_to_clean.update(STAGE_CONFIG[sub_stage]['json_dirs'])
                # Legacy support
                elif stage in STAGE_JSON_DIRS:
                    dirs_to_clean.update(STAGE_JSON_DIRS[stage])
        
        # Delete files in each directory
        for subdir in dirs_to_clean:
            json_file = self.project_root / 'data' / 'json' / subdir / f'{pdb_id}.json'
            if json_file.exists():
                json_file.unlink()
        
        # Also try main json file
        main_json = self.project_root / 'data' / 'json' / f'{pdb_id}.json'
        if main_json.exists():
            main_json.unlink()
    
    def get_pdb_list(self,
                     specific: List[str] = None,
                     max_count: int = None,
                     test_set: int = None) -> List[str]:
        """Get list of PDBs to validate.
        
        Uses data/valid_pdbs_fast.json by default.
        If resume is enabled, skips PDBs that already passed.
        """
        pdb_ids = get_pdb_list(
            self.project_root,
            specific=specific,
            max_count=max_count,
            test_set=test_set
        )
        
        # Filter out already-passed PDBs if resuming
        if self.resume and self.checkpoint_data:
            passed_pdbs = self._get_passed_pdbs()
            original_count = len(pdb_ids)
            pdb_ids = [pdb_id for pdb_id in pdb_ids if pdb_id not in passed_pdbs]
            if not self.quiet and passed_pdbs:
                skipped = original_count - len(pdb_ids)
                print(f"Resuming: skipping {skipped} already-passed PDBs")
        
        return pdb_ids
    
    def validate(self, pdb_ids: List[str]) -> ValidationSummary:
        """Validate all PDBs. Returns summary.
        
        Args:
            pdb_ids: List of PDB IDs to validate
            
        Returns:
            ValidationSummary with results
        """
        passed = 0
        failed = 0
        skipped = 0
        differences = []
        first_failure_pdb = None
        all_results = {}  # For checkpoint
        
        # Load existing results from checkpoint if resuming
        if self.checkpoint_data:
            all_results = self.checkpoint_data.get('results', {}).copy()
        
        self.output.start(len(pdb_ids), self.stages)
        
        if self.stop_on_first or len(pdb_ids) == 1:
            # Sequential processing for stop-on-first or single PDB
            for i, pdb_id in enumerate(pdb_ids):
                result = self._validate_single_sequential(pdb_id)
                all_results[pdb_id] = result
                
                if result['status'] == 'match':
                    passed += 1
                    self._clean_modern_json(pdb_id)  # Clean matched files
                elif result['status'] == 'skip' or result['status'] == 'error':
                    if 'not found' in str(result.get('errors', [])):
                        skipped += 1
                    else:
                        failed += 1
                        if self.document_differences:
                            differences.append(result)
                        if self.stop_on_first and first_failure_pdb is None:
                            first_failure_pdb = pdb_id
                            self.output.first_failure(pdb_id, result)
                            # Save checkpoint before stopping
                            self._save_checkpoint(all_results)
                            break
                else:  # diff
                    failed += 1
                    if self.document_differences:
                        differences.append(result)
                    if self.stop_on_first and first_failure_pdb is None:
                        first_failure_pdb = pdb_id
                        self.output.first_failure(pdb_id, result)
                        # Save checkpoint before stopping
                        self._save_checkpoint(all_results)
                        break
                
                self.output.progress(i + 1, len(pdb_ids), pdb_id, result)
                
                # Save checkpoint periodically (every 100 PDBs)
                if self.checkpoint_file and (i + 1) % 100 == 0:
                    self._save_checkpoint(all_results)
        else:
            # Parallel processing
            results = self._validate_parallel(pdb_ids)
            
            for i, (pdb_id, result) in enumerate(results):
                all_results[pdb_id] = result
                
                if result['status'] == 'match':
                    passed += 1
                    self._clean_modern_json(pdb_id)  # Clean matched files
                elif result['status'] == 'skip' or result['status'] == 'error':
                    if 'not found' in str(result.get('errors', [])):
                        skipped += 1
                    else:
                        failed += 1
                        if self.document_differences:
                            differences.append(result)
                else:  # diff
                    failed += 1
                    if self.document_differences:
                        differences.append(result)
                
                self.output.progress(i + 1, len(pdb_ids), pdb_id, result)
            
            # Save checkpoint after parallel run
            if self.checkpoint_file:
                self._save_checkpoint(all_results)
        
        # Final checkpoint save
        if self.checkpoint_file:
            self._save_checkpoint(all_results)
        
        # Write differences file if requested
        diff_file_path = None
        if differences and self.document_differences:
            diff_file_path = self._write_differences(differences)
        
        summary = ValidationSummary(
            total=len(pdb_ids),
            passed=passed,
            failed=failed,
            skipped=skipped,
            stages_tested=self.stages,
            differences_file=diff_file_path,
            first_failure_pdb=first_failure_pdb,
            checkpoint_file=self.checkpoint_file
        )
        
        self.output.summary(summary)
        return summary
    
    def _validate_single_sequential(self, pdb_id: str) -> Dict[str, Any]:
        """Validate single PDB (sequential mode)."""
        comparator = JsonComparator(**self.comparator_kwargs)
        
        legacy_file = self.project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
        modern_file = self.project_root / 'data' / 'json' / f'{pdb_id}.json'
        pdb_file = self.project_root / 'data' / 'pdb' / f'{pdb_id}.pdb'
        
        result = comparator.compare_files(legacy_file, modern_file, pdb_file, pdb_id)
        
        return {
            'pdb_id': pdb_id,
            'status': result.status,
            'has_differences': result.has_differences(),
            'errors': result.errors,
            'atom_comparison': _summarize_atom_comparison(result.atom_comparison),
            'frame_comparison': _summarize_frame_comparison(result.frame_comparison),
            'hbond_comparison': _summarize_hbond_comparison(result.hbond_list_comparison),
        }
    
    def _validate_parallel(self, pdb_ids: List[str]) -> List[tuple]:
        """Validate PDBs in parallel."""
        results = []
        
        # Prepare arguments for worker function
        args_list = [
            (pdb_id, str(self.project_root), self.comparator_kwargs)
            for pdb_id in pdb_ids
        ]
        
        with ProcessPoolExecutor(max_workers=self.workers) as executor:
            future_to_pdb = {
                executor.submit(_validate_single_pdb, args): args[0]
                for args in args_list
            }
            
            for future in as_completed(future_to_pdb):
                pdb_id = future_to_pdb[future]
                try:
                    result = future.result()
                    results.append((pdb_id, result))
                except Exception as e:
                    results.append((pdb_id, {
                        'pdb_id': pdb_id,
                        'status': 'error',
                        'has_differences': True,
                        'errors': [f'Exception: {str(e)}'],
                    }))
        
        # Sort by original order
        pdb_order = {pdb_id: i for i, pdb_id in enumerate(pdb_ids)}
        results.sort(key=lambda x: pdb_order.get(x[0], len(pdb_ids)))
        
        return results
    
    def _write_differences(self, differences: List[Dict]) -> Path:
        """Write differences to JSON file."""
        diff_file = self.diff_file
        if diff_file is None:
            diff_file = self.project_root / 'data' / 'validation_results' / 'differences.json'
        
        diff_file = Path(diff_file)
        diff_file.parent.mkdir(parents=True, exist_ok=True)
        
        output = {
            'total_differences': len(differences),
            'stages_tested': self.stages,
            'differences': differences
        }
        
        with open(diff_file, 'w') as f:
            json.dump(output, f, indent=2, default=str)
        
        return diff_file


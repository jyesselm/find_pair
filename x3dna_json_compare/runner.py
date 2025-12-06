"""
Unified validation runner.

Single implementation for all validation types - replaces the many validate_*.py scripts.
"""

import json
from pathlib import Path
from typing import List, Optional, Dict, Any
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from .json_comparison import JsonComparator
from .output import OutputFormatter, ValidationSummary
from .pdb_list import get_pdb_list


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
    return {
        'total_legacy': fc.total_legacy,
        'matched': fc.matched,
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
    """Core validation runner - single implementation for all validation types."""
    
    def __init__(self,
                 project_root: Path = None,
                 stages: List[str] = None,
                 workers: int = None,
                 quiet: bool = False,
                 verbose: bool = False,
                 stop_on_first: bool = False,
                 document_differences: bool = False,
                 diff_file: Path = None,
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
        
        # Build comparator kwargs from stages
        self.comparator_kwargs = self._build_comparator_kwargs(stages, comparator_kwargs)
        self.output = OutputFormatter(quiet=quiet, verbose=verbose)
    
    def _build_comparator_kwargs(self, stages: List[str], extra_kwargs: Dict) -> Dict:
        """Build JsonComparator kwargs from stage list."""
        # Default: all off
        kwargs = {
            'compare_atoms': False,
            'compare_frames': False,
            'compare_hbond_list': False,
            'compare_pairs': False,
            'compare_steps': False,
            'compare_residue_indices': False,
        }
        
        if not stages or 'all' in stages:
            # Enable all
            kwargs = {k: True for k in kwargs}
        else:
            # Enable only requested stages
            stage_map = {
                'atoms': 'compare_atoms',
                'frames': 'compare_frames',
                'hbonds': 'compare_hbond_list',
                'pairs': 'compare_pairs',
                'steps': 'compare_steps',
                'residue_indices': 'compare_residue_indices',
            }
            for stage in stages:
                if stage in stage_map:
                    kwargs[stage_map[stage]] = True
        
        # Override with any explicit kwargs
        kwargs.update(extra_kwargs)
        return kwargs
    
    def get_pdb_list(self,
                     specific: List[str] = None,
                     max_count: int = None,
                     test_set: int = None) -> List[str]:
        """Get list of PDBs to validate.
        
        Uses data/valid_pdbs_fast.json by default.
        """
        return get_pdb_list(
            self.project_root,
            specific=specific,
            max_count=max_count,
            test_set=test_set
        )
    
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
        
        self.output.start(len(pdb_ids), self.stages)
        
        if self.stop_on_first or len(pdb_ids) == 1:
            # Sequential processing for stop-on-first or single PDB
            for i, pdb_id in enumerate(pdb_ids):
                result = self._validate_single_sequential(pdb_id)
                
                if result['status'] == 'match':
                    passed += 1
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
                            break
                else:  # diff
                    failed += 1
                    if self.document_differences:
                        differences.append(result)
                    if self.stop_on_first and first_failure_pdb is None:
                        first_failure_pdb = pdb_id
                        self.output.first_failure(pdb_id, result)
                        break
                
                self.output.progress(i + 1, len(pdb_ids), pdb_id, result)
        else:
            # Parallel processing
            results = self._validate_parallel(pdb_ids)
            
            for i, (pdb_id, result) in enumerate(results):
                if result['status'] == 'match':
                    passed += 1
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
            first_failure_pdb=first_failure_pdb
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


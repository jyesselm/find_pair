"""
Main JSON comparison engine.

Provides the core comparison functionality with caching and parallel processing.
"""

import json
from pathlib import Path
from typing import Dict, List, Optional
import time

from .models import ComparisonResult, AtomComparison, FrameComparison
from .atom_comparison import compare_atoms
from .frame_comparison import compare_frames
from .pdb_utils import PdbFileReader
from .result_cache import ComparisonCache
from .parallel_executor import ParallelExecutor, print_progress


class JsonComparator:
    """Main comparison engine with caching support."""
    
    def __init__(self, cache_dir: Path = Path(".comparison_cache"),
                 tolerance: float = 1e-6, enable_cache: bool = True,
                 force_recompute: bool = False):
        """
        Initialize comparator with caching support.
        
        Args:
            cache_dir: Directory to store cache files
            tolerance: Tolerance for floating point comparisons
            enable_cache: Enable result caching
            force_recompute: Force recomputation even if cache exists
        """
        self.cache_dir = Path(cache_dir)
        self.tolerance = tolerance
        self.enable_cache = enable_cache
        self.force_recompute = force_recompute
        self.cache = ComparisonCache(self.cache_dir) if enable_cache else None
    
    def _load_json(self, json_file: Path) -> Optional[Dict]:
        """Load JSON file."""
        if not json_file.exists():
            return None
        
        try:
            with open(json_file, 'r') as f:
                return json.load(f)
        except Exception:
            return None
    
    def _extract_atoms(self, json_data: Dict) -> List[Dict]:
        """Extract atom records from JSON."""
        atoms = []
        for calc in json_data.get('calculations', []):
            if calc.get('type') == 'pdb_atoms':
                atoms = calc.get('atoms', [])
                break
        return atoms
    
    def _extract_frame_records(self, json_data: Dict) -> List[Dict]:
        """Extract frame calculation records from JSON."""
        records = []
        for calc in json_data.get('calculations', []):
            if calc.get('type') in ['base_frame_calc', 'ls_fitting', 'frame_calc']:
                records.append(calc)
        return records
    
    def compare_files(self, legacy_file: Path, modern_file: Path,
                     pdb_file: Path, pdb_id: str) -> ComparisonResult:
        """
        Compare two JSON files, using cache if available.
        
        Args:
            legacy_file: Path to legacy JSON file
            modern_file: Path to modern JSON file
            pdb_file: Path to PDB file
            pdb_id: PDB identifier
            
        Returns:
            ComparisonResult
        """
        # Check cache first
        if self.cache and not self.force_recompute:
            cached_result = self.cache.get_cached_result(pdb_id)
            if cached_result:
                return cached_result
        
        # Load JSON files
        legacy_json = self._load_json(legacy_file)
        modern_json = self._load_json(modern_file)
        
        result = ComparisonResult(
            pdb_id=pdb_id,
            status='error',
            timestamp=time.time(),
            pdb_file_path=str(pdb_file)
        )
        
        if legacy_json is None:
            result.errors.append(f"Legacy JSON not found: {legacy_file}")
            return result
        
        if modern_json is None:
            result.errors.append(f"Modern JSON not found: {modern_file}")
            return result
        
        # Generate cache key
        result.cache_key = self.cache.get_cache_key(
            legacy_file, modern_file, pdb_file
        ) if self.cache else None
        
        # Initialize PDB reader
        try:
            pdb_reader = PdbFileReader(pdb_file)
        except Exception as e:
            pdb_reader = None
            result.errors.append(f"Error loading PDB file: {e}")
        
        # Compare atoms
        try:
            legacy_atoms = self._extract_atoms(legacy_json)
            modern_atoms = self._extract_atoms(modern_json)
            
            atom_comparison = compare_atoms(
                legacy_atoms, modern_atoms, pdb_file, pdb_reader
            )
            result.atom_comparison = atom_comparison
        except Exception as e:
            result.errors.append(f"Error comparing atoms: {e}")
        
        # Compare frames
        try:
            legacy_frames = self._extract_frame_records(legacy_json)
            modern_frames = self._extract_frame_records(modern_json)
            
            frame_comparison = compare_frames(
                legacy_frames, modern_frames, pdb_file, pdb_reader
            )
            result.frame_comparison = frame_comparison
        except Exception as e:
            result.errors.append(f"Error comparing frames: {e}")
        
        # Determine status
        if result.errors:
            result.status = 'error'
        elif result.has_differences():
            result.status = 'diff'
        else:
            result.status = 'match'
        
        # Cache result
        if self.cache and result.status != 'error':
            self.cache.cache_result(result, legacy_file, modern_file, pdb_file)
        
        return result
    
    def batch_compare(self, pdb_list: List[str], project_root: Path = Path("."),
                     max_workers: Optional[int] = None,
                     progress_callback: Optional[callable] = None) -> Dict[str, ComparisonResult]:
        """
        Compare multiple PDBs in parallel.
        
        Args:
            pdb_list: List of PDB IDs to compare
            project_root: Project root directory
            max_workers: Maximum number of parallel workers
            progress_callback: Optional callback(count, total) for progress
            
        Returns:
            Dictionary mapping PDB ID -> ComparisonResult
        """
        project_root = Path(project_root)
        
        # Prepare comparison tasks
        def compare_pdb(pdb_id: str) -> tuple[str, ComparisonResult]:
            legacy_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
            modern_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
            pdb_file = project_root / 'data' / 'pdb' / f'{pdb_id}.pdb'
            
            result = self.compare_files(legacy_file, modern_file, pdb_file, pdb_id)
            return (pdb_id, result)
        
        # Execute in parallel
        executor = ParallelExecutor(max_workers=max_workers, use_processes=True)
        results_list = executor.map(
            compare_pdb,
            pdb_list,
            progress_callback=progress_callback or print_progress
        )
        
        return dict(results_list)


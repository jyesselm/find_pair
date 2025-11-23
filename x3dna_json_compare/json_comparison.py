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
                 force_recompute: bool = False,
                 compare_atoms: bool = True,
                 compare_frames: bool = True):
        """
        Initialize comparator with caching support.
        
        Args:
            cache_dir: Directory to store cache files
            tolerance: Tolerance for floating point comparisons
            enable_cache: Enable result caching
            force_recompute: Force recomputation even if cache exists
            compare_atoms: If True, compare atom records (pdb_atoms)
            compare_frames: If True, compare frame calculations (base_frame_calc, frame_calc, ls_fitting)
        """
        self.cache_dir = Path(cache_dir)
        self.tolerance = tolerance
        self.enable_cache = enable_cache
        self.force_recompute = force_recompute
        self.compare_atoms = compare_atoms
        self.compare_frames = compare_frames
        self.cache = ComparisonCache(self.cache_dir) if enable_cache else None
    
    def _load_json(self, json_file: Path) -> Optional[Dict]:
        """Load JSON file, handling both objects and arrays, and files with extra data."""
        if not json_file.exists():
            return None
        
        try:
            with open(json_file, 'r') as f:
                content = f.read().strip()
                
                # Try to parse as complete JSON first
                try:
                    return json.loads(content)
                except json.JSONDecodeError:
                    pass
                
                # If that fails, try to extract the first complete JSON structure
                # (object or array)
                brace_count = 0
                bracket_count = 0
                in_string = False
                escape_next = False
                end_pos = 0
                is_array = False
                
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
                            if i == 0 or (i > 0 and content[i-1].isspace()):
                                is_array = True
                        elif char == ']':
                            bracket_count -= 1
                            if bracket_count == 0 and brace_count == 0:
                                end_pos = i + 1
                                break
                
                if end_pos > 0:
                    # Parse only the first complete JSON structure
                    return json.loads(content[:end_pos])
                else:
                    # Last resort: try to parse what we can
                    raise json.JSONDecodeError("Could not find complete JSON structure", content, 0)
        except Exception as e:
            # Log the error for debugging
            import sys
            print(f"Error loading JSON file {json_file}: {e}", file=sys.stderr)
            return None
    
    def _extract_atoms(self, json_data: Dict, json_file: Optional[Path] = None) -> List[Dict]:
        """Extract atom records from JSON (supports both array, grouped, and split file formats)."""
        atoms = []
        calculations = json_data.get('calculations', {})
        
        # Check if split files are indicated
        is_split_files = False
        if isinstance(calculations, list):
            for calc in calculations:
                if isinstance(calc, dict) and calc.get('_split_files'):
                    is_split_files = True
                    break
        
        # If split files, load from split file
        if is_split_files and json_file:
            pdb_name = json_data.get('pdb_name', json_file.stem)
            split_file = json_file.parent / f"{pdb_name}_pdb_atoms.json"
            if split_file.exists():
                try:
                    split_data = self._load_json(split_file)
                    if isinstance(split_data, list) and len(split_data) > 0:
                        # Split file contains array of entries, first entry has atoms
                        atoms = split_data[0].get('atoms', [])
                except Exception:
                    pass
        
        # Handle grouped format: calculations is a dict with type keys
        if not atoms and isinstance(calculations, dict):
            pdb_atoms_group = calculations.get('pdb_atoms', [])
            if pdb_atoms_group and len(pdb_atoms_group) > 0:
                # First entry in the group should have the atoms array
                atoms = pdb_atoms_group[0].get('atoms', [])
        # Handle array format: calculations is a list
        elif not atoms and isinstance(calculations, list):
            for calc in calculations:
                if isinstance(calc, dict) and calc.get('type') == 'pdb_atoms':
                    atoms = calc.get('atoms', [])
                    break
        
        # If still no atoms found and we have a file path, try split file format as fallback
        if not atoms and json_file:
            split_file = json_file.parent / f"{json_file.stem}_pdb_atoms.json"
            if split_file.exists():
                try:
                    split_data = self._load_json(split_file)
                    if isinstance(split_data, list) and len(split_data) > 0:
                        # Split file contains array of entries, first entry has atoms
                        atoms = split_data[0].get('atoms', [])
                except Exception:
                    pass
        
        return atoms
    
    def _extract_frame_records(self, json_data: Dict, json_file: Optional[Path] = None) -> List[Dict]:
        """Extract frame calculation records from JSON (supports array, grouped, and split file formats)."""
        records = []
        calculations = json_data.get('calculations', {})
        
        # Handle grouped format: calculations is a dict with type keys
        if isinstance(calculations, dict):
            for calc_type in ['base_frame_calc', 'ls_fitting', 'frame_calc']:
                calc_group = calculations.get(calc_type, [])
                if isinstance(calc_group, list):
                    records.extend(calc_group)
        # Handle array format: calculations is a list
        elif isinstance(calculations, list):
            for calc in calculations:
                if calc.get('type') in ['base_frame_calc', 'ls_fitting', 'frame_calc']:
                    records.append(calc)
        
        # If no records found and we have a file path, try split file format
        if not records and json_file:
            for calc_type in ['base_frame_calc', 'ls_fitting', 'frame_calc']:
                split_file = json_file.parent / f"{json_file.stem}_{calc_type}.json"
                if split_file.exists():
                    try:
                        split_data = self._load_json(split_file)
                        if isinstance(split_data, list):
                            # Split file contains array of entries (already in correct format)
                            records.extend(split_data)
                    except Exception:
                        pass
        
        return records
    
    def compare_files(self, legacy_file: Path, modern_file: Path,
                     pdb_file: Path, pdb_id: str) -> ComparisonResult:
        """
        Compare two JSON files, using cache if available.
        
        The types of comparisons performed are controlled by the compare_atoms
        and compare_frames options set during initialization.
        
        Args:
            legacy_file: Path to legacy JSON file
            modern_file: Path to modern JSON file
            pdb_file: Path to PDB file
            pdb_id: PDB identifier
            
        Returns:
            ComparisonResult with atom_comparison and/or frame_comparison
            based on the enabled comparison options
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
            if legacy_file.exists():
                result.errors.append(f"Legacy JSON file exists but could not be parsed: {legacy_file} (may be malformed)")
            else:
                result.errors.append(f"Legacy JSON not found: {legacy_file}")
            return result
        
        if modern_json is None:
            if modern_file and modern_file.exists():
                result.errors.append(f"Modern JSON file exists but could not be parsed: {modern_file} (may be malformed)")
            else:
                result.errors.append(f"Modern JSON not found: {modern_file}")
            return result
        
        # Generate cache key
        result.cache_key = self.cache.get_cache_key(
            legacy_file, modern_file, pdb_file
        ) if self.cache else None
        
        # Initialize PDB reader (lazy - only if file exists)
        pdb_reader = None
        if pdb_file and pdb_file.exists():
            try:
                pdb_reader = PdbFileReader(pdb_file)
            except Exception as e:
                # Don't fail comparison if PDB can't be loaded
                pdb_reader = None
        
        # Extract legacy atoms (needed for both atom and frame comparison)
        legacy_atoms = []
        if self.compare_atoms or self.compare_frames:
            legacy_atoms = self._extract_atoms(legacy_json, legacy_file)
        
        # Compare atoms
        if self.compare_atoms:
            try:
                modern_atoms = self._extract_atoms(modern_json, modern_file)
                
                if legacy_atoms or modern_atoms:
                    atom_comparison = compare_atoms(
                        legacy_atoms, modern_atoms, pdb_file, pdb_reader
                    )
                    result.atom_comparison = atom_comparison
            except Exception as e:
                # Don't fail on atom comparison errors
                pass
        
        # Compare frames
        if self.compare_frames:
            try:
                legacy_frames = self._extract_frame_records(legacy_json, legacy_file)
                modern_frames = self._extract_frame_records(modern_json, modern_file)
                
                frame_comparison = compare_frames(
                    legacy_frames, modern_frames, pdb_file, pdb_reader, legacy_atoms
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


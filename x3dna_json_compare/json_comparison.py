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
from .step_comparison import compare_step_parameters
from .pair_comparison import compare_pair_validations, compare_distance_checks
from .base_pair_comparison import compare_base_pairs
from .find_bestpair_comparison import compare_find_bestpair_selections
from .hbond_comparison import compare_hbond_lists
from .residue_indices_comparison import compare_residue_indices
from .pdb_utils import PdbFileReader
from .parallel_executor import ParallelExecutor, print_progress
from .json_file_finder import find_json_file


class JsonComparator:
    """Main comparison engine."""
    
    def __init__(self,
                 tolerance: float = 2e-5,  # 2e-5 handles normal FP variations
                 compare_atoms: bool = True,
                 compare_frames: bool = True,
                 compare_steps: bool = True,
                 compare_pairs: bool = True,
                 compare_hbond_list: bool = True,
                 compare_residue_indices: bool = True,
                 hbond_ignore_count_mismatch: bool = True):
        """
        Initialize comparator.
        
        Args:
            tolerance: Tolerance for floating point comparisons
            compare_atoms: If True, compare atom records (pdb_atoms)
            compare_frames: If True, compare frame calculations (base_frame_calc, frame_calc, ls_fitting)
            compare_steps: If True, compare step parameters (bpstep_params, helical_params)
            compare_pairs: If True, compare pair validation and distance checks
            compare_hbond_list: If True, compare hbond_list records
            compare_residue_indices: If True, compare residue_indices records
            hbond_ignore_count_mismatch: If True, don't flag H-bond count differences as mismatches
        """
        self.tolerance = tolerance
        self.compare_atoms = compare_atoms
        self.compare_frames = compare_frames
        self.compare_steps = compare_steps
        self.compare_pairs = compare_pairs
        self.compare_hbond_list = compare_hbond_list
        self.compare_residue_indices = compare_residue_indices
        self.hbond_ignore_count_mismatch = hbond_ignore_count_mismatch
    
    def _load_json(self, json_file: Path) -> Optional[Dict]:
        """Load JSON file, handling both objects and arrays, and files with extra data."""
        if not json_file.exists():
            return None
        
        # Skip if path is a directory (segmented JSON structure)
        if json_file.is_dir():
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
            base_dir = json_file.parent
            
            # Try new structure first: <record_type>/<PDB_ID>.json
            split_file = find_json_file(base_dir, pdb_name, "pdb_atoms")
            if not split_file or not split_file.exists():
                # Fall back to old structure: <PDB_ID>_<record_type>.json
                split_file = base_dir / f"{pdb_name}_pdb_atoms.json"
            
            if split_file and split_file.exists():
                try:
                    split_data = self._load_json(split_file)
                    if isinstance(split_data, dict):
                        # Modern format: dict with 'atoms' key directly
                        atoms = split_data.get('atoms', [])
                    elif isinstance(split_data, list) and len(split_data) > 0:
                        # Legacy format: array where first entry has atoms
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
            pdb_id = json_file.stem
            base_dir = json_file.parent
            
            # Try new structure first: <record_type>/<PDB_ID>.json
            split_file = find_json_file(base_dir, pdb_id, "pdb_atoms")
            if not split_file or not split_file.exists():
                # Fall back to old structure: <PDB_ID>_<record_type>.json
                split_file = base_dir / f"{pdb_id}_pdb_atoms.json"
            
            if split_file and split_file.exists():
                try:
                    split_data = self._load_json(split_file)
                    if isinstance(split_data, dict):
                        # Modern format: dict with 'atoms' key directly
                        atoms = split_data.get('atoms', [])
                    elif isinstance(split_data, list) and len(split_data) > 0:
                        # Legacy format: array where first entry has atoms
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
        # Try new structure first, then fall back to old structure
        if not records and json_file:
            pdb_id = json_file.stem
            base_dir = json_file.parent
            
            for calc_type in ['base_frame_calc', 'ls_fitting', 'frame_calc']:
                # Each calc_type has its own directory
                dir_name = calc_type
                
                # Try new structure: <record_type>/<PDB_ID>.json
                split_file = find_json_file(base_dir, pdb_id, dir_name)
                if split_file and split_file.exists():
                    try:
                        split_data = self._load_json(split_file)
                        if isinstance(split_data, list):
                            # Add type field to each record if missing
                            for rec in split_data:
                                if 'type' not in rec:
                                    rec['type'] = calc_type
                            records.extend(split_data)
                    except Exception:
                        pass
        
        return records
    
    def _extract_step_records(self, json_data: Dict, json_file: Optional[Path] = None) -> List[Dict]:
        """Extract step parameter records from JSON (supports array, grouped, and split file formats)."""
        records = []
        calculations = json_data.get('calculations', {})
        
        # Handle grouped format: calculations is a dict with type keys
        if isinstance(calculations, dict):
            for calc_type in ['bpstep_params', 'helical_params']:
                calc_group = calculations.get(calc_type, [])
                if isinstance(calc_group, list):
                    records.extend(calc_group)
        # Handle array format: calculations is a list
        elif isinstance(calculations, list):
            for calc in calculations:
                if calc.get('type') in ['bpstep_params', 'helical_params']:
                    records.append(calc)
        
        # If no records found and we have a file path, try split file format
        # Try new structure first, then fall back to old structure
        if not records and json_file:
            pdb_id = json_file.stem
            base_dir = json_file.parent
            
            for calc_type in ['bpstep_params', 'helical_params']:
                # Try new structure: <record_type>/<PDB_ID>.json
                split_file = find_json_file(base_dir, pdb_id, calc_type)
                if split_file and split_file.exists():
                    try:
                        split_data = self._load_json(split_file)
                        if isinstance(split_data, list):
                            # Split file contains array of entries (already in correct format)
                            records.extend(split_data)
                    except Exception:
                        pass
                else:
                    # Fall back to old structure: <PDB_ID>_<record_type>.json
                    split_file = base_dir / f"{pdb_id}_{calc_type}.json"
                    if split_file.exists():
                        try:
                            split_data = self._load_json(split_file)
                            if isinstance(split_data, list):
                                records.extend(split_data)
                        except Exception:
                            pass
        
        return records
    
    def _extract_find_bestpair_selection_records(self, json_data: Dict, json_file: Optional[Path] = None) -> List[Dict]:
        """Extract find_bestpair_selection records from JSON (actual selected pairs from find_bestpair)."""
        records = []
        
        # If json_file is a split file (ends with _find_bestpair_selection.json), load it directly
        if json_file and json_file.name.endswith('_find_bestpair_selection.json'):
            try:
                split_data = self._load_json(json_file)
                if isinstance(split_data, list):
                    records.extend(split_data)
                    return records
            except Exception:
                pass
        
        # Handle split file format (array of records) - this is the primary format
        # Try new structure first, then fall back to old structure
        if json_file:
            # Extract PDB ID from file
            if json_file.name.endswith('_find_bestpair_selection.json'):
                pdb_id = json_file.name.replace('_find_bestpair_selection.json', '')
            else:
                pdb_id = json_file.stem
            base_dir = json_file.parent
            
            # Try new structure: find_bestpair_selection/<PDB_ID>.json
            split_file = find_json_file(base_dir, pdb_id, "find_bestpair_selection")
            if split_file and split_file.exists() and split_file != json_file:
                try:
                    split_data = self._load_json(split_file)
                    if isinstance(split_data, list):
                        records.extend(split_data)
                        return records
                except Exception:
                    pass
            
            # Fall back to old structure: <PDB_ID>_find_bestpair_selection.json
            old_split_file = base_dir / f"{pdb_id}_find_bestpair_selection.json"
            if old_split_file.exists() and old_split_file != json_file:
                try:
                    split_data = self._load_json(old_split_file)
                    if isinstance(split_data, list):
                        records.extend(split_data)
                        return records
                except Exception:
                    pass
        
        # Handle main JSON file format
        if isinstance(json_data, dict):
            calculations = json_data.get('calculations', {})
            
            # Handle grouped format: calculations is a dict with type keys
            if isinstance(calculations, dict):
                calc_group = calculations.get('find_bestpair_selection', [])
                if isinstance(calc_group, list):
                    records.extend(calc_group)
            # Handle array format: calculations is a list
            elif isinstance(calculations, list):
                for calc in calculations:
                    if isinstance(calc, dict) and calc.get('type') == 'find_bestpair_selection':
                        records.append(calc)
        elif isinstance(json_data, list):
            # Direct array format (split file loaded directly)
            records.extend(json_data)
        
        return records
    
    def _extract_base_pair_records(self, json_data: Dict, json_file: Optional[Path] = None) -> List[Dict]:
        """Extract base_pair records from JSON (original identification from find_pair phase).
        
        Legacy writes base_pair records to the main JSON file (not split files).
        Modern writes base_pair records to split files.
        """
        records = []
        calculations = json_data.get('calculations', {})
        
        # Handle grouped format: calculations is a dict with type keys
        if isinstance(calculations, dict):
            calc_group = calculations.get('base_pair', [])
            if isinstance(calc_group, list):
                records.extend(calc_group)
        # Handle array format: calculations is a list
        elif isinstance(calculations, list):
            for calc in calculations:
                if isinstance(calc, dict) and calc.get('type') == 'base_pair':
                    records.append(calc)
        
        # If no records found in main file and we have a file path, try split file format
        # (Modern code writes to split files, legacy writes to main file)
        # Try new structure first, then fall back to old structure
        if not records and json_file:
            pdb_id = json_file.stem
            base_dir = json_file.parent
            
            # Try new structure: base_pair/<PDB_ID>.json
            split_file = find_json_file(base_dir, pdb_id, "base_pair")
            if split_file and split_file.exists():
                try:
                    split_data = self._load_json(split_file)
                    if isinstance(split_data, list):
                        records.extend(split_data)
                except Exception:
                    pass
        
        return records
    
    def _extract_pair_validation_records(self, json_data: Dict, json_file: Optional[Path] = None) -> List[Dict]:
        """Extract pair_validation records from JSON (supports array, grouped, and split file formats)."""
        records = []
        calculations = json_data.get('calculations', {})
        
        # Handle grouped format: calculations is a dict with type keys
        if isinstance(calculations, dict):
            calc_group = calculations.get('pair_validation', [])
            if isinstance(calc_group, list):
                records.extend(calc_group)
        # Handle array format: calculations is a list
        elif isinstance(calculations, list):
            for calc in calculations:
                if calc.get('type') == 'pair_validation':
                    records.append(calc)
        
        # If no records found and we have a file path, try split file format
        # Try new structure first, then fall back to old structure
        if not records and json_file:
            pdb_id = json_file.stem
            base_dir = json_file.parent
            
            # Try new structure: pair_validation/<PDB_ID>.json
            split_file = find_json_file(base_dir, pdb_id, "pair_validation")
            if split_file and split_file.exists():
                try:
                    split_data = self._load_json(split_file)
                    if isinstance(split_data, list):
                        records.extend(split_data)
                except Exception:
                    pass
        
        return records
    
    def _extract_distance_checks_records(self, json_data: Dict, json_file: Optional[Path] = None) -> List[Dict]:
        """Extract distance_checks records from JSON (supports array, grouped, and split file formats)."""
        records = []
        calculations = json_data.get('calculations', {})
        
        # Handle grouped format: calculations is a dict with type keys
        if isinstance(calculations, dict):
            calc_group = calculations.get('distance_checks', [])
            if isinstance(calc_group, list):
                records.extend(calc_group)
        # Handle array format: calculations is a list
        elif isinstance(calculations, list):
            for calc in calculations:
                if calc.get('type') == 'distance_checks':
                    records.append(calc)
        
        # If no records found and we have a file path, try split file format
        # Try new structure first, then fall back to old structure
        if not records and json_file:
            pdb_id = json_file.stem
            base_dir = json_file.parent
            
            # Try new structure: distance_checks/<PDB_ID>.json
            split_file = find_json_file(base_dir, pdb_id, "distance_checks")
            if split_file and split_file.exists():
                try:
                    split_data = self._load_json(split_file)
                    if isinstance(split_data, list):
                        # Split files don't have type field - add it for consistency
                        for record in split_data:
                            if isinstance(record, dict):
                                record_with_type = record.copy()
                                record_with_type['type'] = 'distance_checks'
                                records.append(record_with_type)
                except Exception:
                    pass
        
        return records
    
    def _extract_hbond_list_records(self, json_data: Dict, json_file: Optional[Path] = None) -> List[Dict]:
        """Extract hbond_list records from JSON (supports array, grouped, and split file formats)."""
        records = []
        calculations = json_data.get('calculations', {})
        
        # Handle grouped format: calculations is a dict with type keys
        if isinstance(calculations, dict):
            calc_group = calculations.get('hbond_list', [])
            if isinstance(calc_group, list):
                records.extend(calc_group)
        # Handle array format: calculations is a list
        elif isinstance(calculations, list):
            for calc in calculations:
                if isinstance(calc, dict) and calc.get('type') == 'hbond_list':
                    records.append(calc)
        
        # If no records found and we have a file path, try split file format
        # Try new structure first, then fall back to old structure
        if not records and json_file:
            pdb_id = json_file.stem
            base_dir = json_file.parent
            
            # Try new structure: hbond_list/<PDB_ID>.json
            split_file = find_json_file(base_dir, pdb_id, "hbond_list")
            if split_file and split_file.exists():
                try:
                    split_data = self._load_json(split_file)
                    if isinstance(split_data, list):
                        records.extend(split_data)
                except Exception:
                    pass
        
        # If still no records and file has split files indicator, try loading from main file directly
        # (legacy may have records embedded in main file even if split_files is true)
        if not records and json_file and json_file.exists():
            # Check if calculations indicates split files but we didn't find records
            split_indicator = False
            if isinstance(calculations, list):
                split_indicator = any(isinstance(c, dict) and c.get('_split_files') for c in calculations)
            
            # If split file doesn't exist, parse all JSON objects from main file
            if split_indicator:
                split_file = json_file.parent / f"{json_file.stem}_hbond_list.json"
                if not split_file.exists():
                    # Parse all JSON objects from main file
                    try:
                        with open(json_file, 'r') as f:
                            content = f.read()
                        decoder = json.JSONDecoder()
                        idx = 0
                        while idx < len(content):
                            # Skip whitespace
                            while idx < len(content) and content[idx].isspace():
                                idx += 1
                            if idx >= len(content):
                                break
                            try:
                                obj, new_idx = decoder.raw_decode(content, idx)
                                if isinstance(obj, dict) and obj.get('type') == 'hbond_list':
                                    records.append(obj)
                                idx = new_idx
                            except json.JSONDecodeError:
                                # Try to find next object start
                                next_brace = content.find('{', idx + 1)
                                if next_brace == -1:
                                    break
                                idx = next_brace
                    except Exception:
                        pass
        
        return records
    
    def _extract_residue_indices_records(self, json_data: Optional[Dict], json_file: Optional[Path] = None) -> List[Dict]:
        """Extract residue_indices records from JSON (supports array, grouped, and split file formats)."""
        records = []
        
        # Try to extract from main JSON data if available
        if json_data:
            calculations = json_data.get('calculations', {})
            
            # Handle grouped format: calculations is a dict with type keys
            if isinstance(calculations, dict):
                calc_group = calculations.get('residue_indices', [])
                if isinstance(calc_group, list):
                    records.extend(calc_group)
            # Handle array format: calculations is a list
            elif isinstance(calculations, list):
                for calc in calculations:
                    if isinstance(calc, dict) and calc.get('type') == 'residue_indices':
                        records.append(calc)
        
        # If no records found and we have a file path, try split file format
        # Try new structure first, then fall back to old structure
        if not records and json_file:
            pdb_id = json_file.stem
            base_dir = json_file.parent
            
            # Try new structure: residue_indices/<PDB_ID>.json
            split_file = find_json_file(base_dir, pdb_id, "residue_indices")
            if split_file and split_file.exists():
                try:
                    split_data = self._load_json(split_file)
                    if isinstance(split_data, list):
                        records.extend(split_data)
                    elif isinstance(split_data, dict) and split_data.get('type') == 'residue_indices':
                        records.append(split_data)
                except Exception:
                    pass
            
            # Fall back to old structure: <PDB_ID>_residue_indices.json
            if not records:
                old_split_file = base_dir / f"{pdb_id}_residue_indices.json"
                if old_split_file.exists():
                    try:
                        split_data = self._load_json(old_split_file)
                        if isinstance(split_data, list):
                            records.extend(split_data)
                        elif isinstance(split_data, dict) and split_data.get('type') == 'residue_indices':
                            records.append(split_data)
                    except Exception:
                        pass
        
        # Handle main JSON file directly if it's an array or dict
        if not records and json_file and json_file.exists():
            try:
                with open(json_file, 'r') as f:
                    content = f.read()
                decoder = json.JSONDecoder()
                idx = 0
                while idx < len(content):
                    # Skip whitespace
                    while idx < len(content) and content[idx].isspace():
                        idx += 1
                    if idx >= len(content):
                        break
                    try:
                        obj, new_idx = decoder.raw_decode(content, idx)
                        if isinstance(obj, dict) and obj.get('type') == 'residue_indices':
                            records.append(obj)
                        idx = new_idx
                    except json.JSONDecodeError:
                        # Try to find next object start
                        next_brace = content.find('{', idx + 1)
                        if next_brace == -1:
                            break
                        idx = next_brace
            except Exception:
                pass
        
        return records
    
    def compare_files(self, legacy_file: Path, modern_file: Path,
                     pdb_file: Path, pdb_id: str) -> ComparisonResult:
        """
        Compare two JSON files.
        
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
        
        # Load JSON files
        # If main legacy file doesn't exist, check for split files and create minimal wrapper
        legacy_json = self._load_json(legacy_file)
        if legacy_json is None and not legacy_file.exists():
            # Check if split files exist in new structure: <record_type>/<PDB_ID>.json
            # or old structure: <PDB_ID>_<record_type>.json
            pdb_id = legacy_file.stem
            base_dir = legacy_file.parent
            
            # Try new structure first - check essential record types
            split_file_found = None
            for record_type in ["find_bestpair_selection", "base_pair", "hbond_list", "base_frame_calc", "pdb_atoms", "pair_validation", "distance_checks"]:
                split_file = find_json_file(base_dir, pdb_id, record_type)
                if split_file and split_file.exists():
                    split_file_found = split_file
                    break
            
            # Fall back to old structure
            if not split_file_found:
                for record_type in ["find_bestpair_selection", "base_pair", "hbond_list", "base_frame_calc", "pdb_atoms", "pair_validation", "distance_checks"]:
                    old_file = base_dir / f"{pdb_id}_{record_type}.json"
                    if old_file.exists():
                        split_file_found = old_file
                        break
            
            if split_file_found:
                # Create minimal JSON structure that indicates split files exist
                legacy_json = {
                    'pdb_name': pdb_id,
                    'calculations': [{'_split_files': True}]
                }
            else:
                legacy_json = None
        
        modern_json = self._load_json(modern_file)
        if modern_json is None and modern_file and not modern_file.exists():
            # Check if split files exist in new structure: <record_type>/<PDB_ID>.json
            # or old structure: <PDB_ID>_<record_type>.json
            pdb_id_from_file = modern_file.stem
            base_dir = modern_file.parent
            
            # Try new structure first - check essential record types
            split_file_found = None
            for record_type in ["find_bestpair_selection", "base_pair", "hbond_list", "base_frame_calc", "pdb_atoms", "pair_validation", "distance_checks"]:
                split_file = find_json_file(base_dir, pdb_id_from_file, record_type)
                if split_file and split_file.exists():
                    split_file_found = split_file
                    break
            
            # Fall back to old structure
            if not split_file_found:
                for record_type in ["find_bestpair_selection", "base_pair", "hbond_list", "base_frame_calc", "pdb_atoms", "pair_validation", "distance_checks"]:
                    old_file = base_dir / f"{pdb_id_from_file}_{record_type}.json"
                    if old_file.exists():
                        split_file_found = old_file
                        break
            
            if split_file_found:
                # Create minimal JSON structure that indicates split files exist
                modern_json = {
                    'pdb_name': pdb_id_from_file,
                    'calculations': [{'_split_files': True}]
                }
            else:
                modern_json = None
        
        result = ComparisonResult(
            pdb_id=pdb_id,
            status='pending',  # Use 'pending' initially so has_differences() works correctly
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
        
        # No cache - removed for reliability
        
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
        
        # Compare frames FIRST (step parameters depend on frames, so verify frames are correct)
        # ALWAYS check frames if comparing steps, even if compare_frames is False
        frame_comparison = None
        frames_match = True
        need_frame_check = self.compare_frames or self.compare_steps
        
        if need_frame_check:
            try:
                legacy_frames = self._extract_frame_records(legacy_json, legacy_file)
                modern_frames = self._extract_frame_records(modern_json, modern_file)
                
                frame_comparison = compare_frames(
                    legacy_frames, modern_frames, pdb_file, pdb_reader, legacy_atoms
                )
                
                # Only store in result if compare_frames was explicitly requested
                if self.compare_frames:
                    result.frame_comparison = frame_comparison
                
                # Check if frames match (no mismatches)
                if frame_comparison.mismatched_calculations or frame_comparison.missing_residues:
                    frames_match = False
                    if self.compare_steps:
                        result.errors.append(
                            "WARNING: Frame mismatches detected. Step parameter comparison "
                            "may be unreliable since step parameters are calculated from frames. "
                            "Fix frame differences first before comparing step parameters."
                        )
            except Exception as e:
                result.errors.append(f"Error comparing frames: {e}")
                frames_match = False
        
        # Compare step parameters (warn if frames don't match)
        if self.compare_steps:
            try:
                legacy_steps = self._extract_step_records(legacy_json, legacy_file)
                modern_steps = self._extract_step_records(modern_json, modern_file)
                
                if legacy_steps or modern_steps:
                    # If frames don't match, warn but still do step comparison
                    if not frames_match and frame_comparison:
                        result.errors.append(
                            f"Step parameter comparison proceeding despite frame mismatches "
                            f"({len(frame_comparison.mismatched_calculations)} mismatched frames, "
                            f"{len(frame_comparison.missing_residues)} missing residues). "
                            f"Results may be unreliable."
                        )
                    
                    # Compare bpstep_params
                    legacy_bpstep = [r for r in legacy_steps if r.get('type') == 'bpstep_params']
                    modern_bpstep = [r for r in modern_steps if r.get('type') == 'bpstep_params']
                    
                    bpstep_comparison = None
                    if legacy_bpstep or modern_bpstep:
                        bpstep_comparison = compare_step_parameters(
                            legacy_bpstep, modern_bpstep, 
                            parameter_type="bpstep_params",
                            frame_comparison=frame_comparison  # Pass frame comparison to verify frames first
                        )
                    
                    # Compare helical_params
                    legacy_helical = [r for r in legacy_steps if r.get('type') == 'helical_params']
                    modern_helical = [r for r in modern_steps if r.get('type') == 'helical_params']
                    
                    helical_comparison = None
                    if legacy_helical or modern_helical:
                        helical_comparison = compare_step_parameters(
                            legacy_helical, modern_helical, 
                            parameter_type="helical_params",
                            frame_comparison=frame_comparison  # Pass frame comparison to verify frames first
                        )
                    
                    # Store both comparison types separately (unified approach)
                    if bpstep_comparison:
                        result.step_comparison = bpstep_comparison
                    if helical_comparison:
                        result.helical_comparison = helical_comparison
            except Exception as e:
                # Don't fail comparison on step parameter errors (they may not be generated yet)
                result.errors.append(f"Error comparing step parameters: {e}")
                pass
        
        # Compare find_bestpair_selection records (PRIMARY - actual selected pairs from find_bestpair)
        if self.compare_pairs:
            try:
                legacy_selection = self._extract_find_bestpair_selection_records(legacy_json, legacy_file)
                modern_selection = self._extract_find_bestpair_selection_records(modern_json, modern_file)
                
                if legacy_selection or modern_selection:
                    find_bestpair_comparison = compare_find_bestpair_selections(
                        legacy_selection, modern_selection
                    )
                    result.find_bestpair_comparison = find_bestpair_comparison
            except Exception as e:
                result.errors.append(f"Error comparing find_bestpair_selection records: {e}")

        # Compare base_pair records (from calculate_more_bppars - pairs that pass all checks)
        if self.compare_pairs:
            try:
                legacy_base_pairs = self._extract_base_pair_records(legacy_json, legacy_file)
                modern_base_pairs = self._extract_base_pair_records(modern_json, modern_file)
                
                if legacy_base_pairs or modern_base_pairs:
                    base_pair_comparison = compare_base_pairs(
                        legacy_base_pairs, modern_base_pairs, self.tolerance
                    )
                    result.base_pair_comparison = base_pair_comparison
            except Exception as e:
                result.errors.append(f"Error comparing base_pair records: {e}")

        # Compare pair validations and distance checks (SECONDARY - for debugging)
        if self.compare_pairs:
            try:
                legacy_pair_validations = self._extract_pair_validation_records(legacy_json, legacy_file)
                modern_pair_validations = self._extract_pair_validation_records(modern_json, modern_file)
                
                if legacy_pair_validations or modern_pair_validations:
                    pair_validation_comparison = compare_pair_validations(
                        legacy_pair_validations, modern_pair_validations, self.tolerance
                    )
                    result.pair_validation_comparison = pair_validation_comparison
                
                legacy_distance_checks = self._extract_distance_checks_records(legacy_json, legacy_file)
                modern_distance_checks = self._extract_distance_checks_records(modern_json, modern_file)
                
                if legacy_distance_checks or modern_distance_checks:
                    distance_checks_comparison = compare_distance_checks(
                        legacy_distance_checks, modern_distance_checks, self.tolerance
                    )
                    result.distance_checks_comparison = distance_checks_comparison
            except Exception as e:
                # Don't fail comparison on pair validation errors
                result.errors.append(f"Error comparing pair validations: {e}")
        
        # Compare hbond_list records
        if self.compare_hbond_list:
            try:
                legacy_hbond_lists = self._extract_hbond_list_records(legacy_json, legacy_file)
                modern_hbond_lists = self._extract_hbond_list_records(modern_json, modern_file)
                
                if legacy_hbond_lists or modern_hbond_lists:
                    hbond_list_comparison = compare_hbond_lists(
                        legacy_hbond_lists, modern_hbond_lists, self.tolerance,
                        ignore_count_mismatch=self.hbond_ignore_count_mismatch
                    )
                    result.hbond_list_comparison = hbond_list_comparison
            except Exception as e:
                result.errors.append(f"Error comparing hbond_list records: {e}")
        
        # Compare residue_indices records
        if self.compare_residue_indices:
            try:
                legacy_residue_indices = self._extract_residue_indices_records(legacy_json, legacy_file)
                modern_residue_indices = self._extract_residue_indices_records(modern_json, modern_file)
                
                if legacy_residue_indices or modern_residue_indices:
                    residue_indices_comparison = compare_residue_indices(
                        legacy_residue_indices, modern_residue_indices, self.tolerance
                    )
                    result.residue_indices_comparison = residue_indices_comparison
            except Exception as e:
                result.errors.append(f"Error comparing residue_indices records: {e}")
        
        # Determine status
        if result.errors:
            result.status = 'error'
        elif result.has_differences():
            result.status = 'diff'
        else:
            result.status = 'match'
        
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
            # Use dummy file paths - JsonComparator will find split files in new structure
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


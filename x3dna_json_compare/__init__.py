"""
X3DNA JSON Comparison Library

A Python package for comparing legacy and modern JSON outputs from 3DNA/X3DNA.
Provides utilities for comparing atom records, frame calculations, and PDB files.

Usage:
    from x3dna_json_compare import JsonComparator, ComparisonResult
    
    comparator = JsonComparator()
    result = comparator.compare_files(legacy_file, modern_file, pdb_file, pdb_id)
"""

from .models import (
    AtomInfo,
    FieldMismatch,
    AtomComparison,
    FrameRecord,
    AtomLineInfo,
    FrameMismatch,
    FrameComparison,
    ComparisonResult,
)
from .pdb_utils import PdbFileReader, get_pdb_line
from .atom_comparison import compare_atoms
from .frame_comparison import compare_frames
from .result_cache import ComparisonCache
from .json_comparison import JsonComparator
from .json_validator import JsonValidator
from .parallel_executor import ParallelExecutor, print_progress
from .hbond_comparison import compare_hbond_lists, HBondListComparison
from .config import load_config, get_comparison_flags

__all__ = [
    # Models
    'AtomInfo',
    'FieldMismatch',
    'AtomComparison',
    'FrameRecord',
    'AtomLineInfo',
    'FrameMismatch',
    'FrameComparison',
    'ComparisonResult',
    # PDB utilities
    'PdbFileReader',
    'get_pdb_line',
    # Comparison functions
    'compare_atoms',
    'compare_frames',
    'compare_hbond_lists',
    # Comparison result models
    'HBondListComparison',
    # Main classes
    'ComparisonCache',
    'JsonComparator',
    'JsonValidator',
    'ParallelExecutor',
    'print_progress',
    # Configuration
    'load_config',
    'get_comparison_flags',
]

__version__ = '1.0.0'

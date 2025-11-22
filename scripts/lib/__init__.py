"""
X3DNA Comparison Library

Reusable modules for comparing legacy and modern JSON outputs.
Provides caching, parallel processing, and unified comparison interfaces.
"""

from .models import (
    AtomInfo,
    AtomComparison,
    FrameMismatch,
    FrameComparison,
    ComparisonResult,
)
from .pdb_utils import PdbFileReader
from .result_cache import ComparisonCache
from .json_comparison import JsonComparator
from .parallel_executor import ParallelExecutor

__all__ = [
    'AtomInfo',
    'AtomComparison',
    'FrameMismatch',
    'FrameComparison',
    'ComparisonResult',
    'PdbFileReader',
    'ComparisonCache',
    'JsonComparator',
    'ParallelExecutor',
]

__version__ = '1.0.0'


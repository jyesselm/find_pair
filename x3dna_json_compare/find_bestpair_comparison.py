"""
Find bestpair selection comparison utilities.

Compares find_bestpair_selection records (the actual pairs selected by find_bestpair)
between legacy and modern JSON outputs.
"""

from typing import Dict, List
from dataclasses import dataclass, field


@dataclass
class FindBestpairComparison:
    """Result of find_bestpair_selection comparison."""
    missing_in_modern: List[tuple] = field(default_factory=list)
    extra_in_modern: List[tuple] = field(default_factory=list)
    total_legacy: int = 0
    total_modern: int = 0
    common_count: int = 0


def compare_find_bestpair_selections(
    legacy_records: List[Dict],
    modern_records: List[Dict]
) -> FindBestpairComparison:
    """
    Compare find_bestpair_selection records between legacy and modern JSON.
    
    These records represent the actual pairs selected by find_bestpair
    (mutual best matches), which is the original base pair identification.
    
    Args:
        legacy_records: List of legacy find_bestpair_selection records
        modern_records: List of modern find_bestpair_selection records
        
    Returns:
        FindBestpairComparison result
    """
    result = FindBestpairComparison()
    
    # Extract pairs from records
    legacy_pairs = set()
    modern_pairs = set()
    
    for rec in legacy_records:
        if rec.get('type') != 'find_bestpair_selection':
            continue
        pairs = rec.get('pairs', [])
        for pair in pairs:
            if isinstance(pair, list) and len(pair) >= 2:
                # Normalize pair order (always i < j) for comparison
                i, j = pair[0], pair[1]
                normalized_pair = (min(i, j), max(i, j))
                legacy_pairs.add(normalized_pair)
    
    for rec in modern_records:
        if rec.get('type') != 'find_bestpair_selection':
            continue
        pairs = rec.get('pairs', [])
        for pair in pairs:
            if isinstance(pair, list) and len(pair) >= 2:
                # Normalize pair order (always i < j) for comparison
                i, j = pair[0], pair[1]
                normalized_pair = (min(i, j), max(i, j))
                modern_pairs.add(normalized_pair)
    
    result.total_legacy = len(legacy_pairs)
    result.total_modern = len(modern_pairs)
    
    common_pairs = legacy_pairs & modern_pairs
    result.common_count = len(common_pairs)
    
    # Find missing in modern
    result.missing_in_modern = sorted(legacy_pairs - modern_pairs)
    
    # Find extra in modern
    result.extra_in_modern = sorted(modern_pairs - legacy_pairs)
    
    return result


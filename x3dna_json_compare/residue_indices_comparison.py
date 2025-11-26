"""
Residue indices comparison utilities.

Compares residue_indices records between legacy and modern JSON outputs.
Each residue_indices record contains the mapping of residues to atom ranges (seidx).
"""

from typing import Dict, List, Optional
from dataclasses import dataclass, field


@dataclass
class ResidueIndexEntry:
    """Single residue index entry from seidx array."""
    residue_idx: int  # 1-based residue index
    start_atom: int   # First atom index (1-based)
    end_atom: int     # Last atom index (1-based)


@dataclass
class ResidueIndicesComparison:
    """Result of residue_indices comparison."""

    missing_in_modern: bool = False
    extra_in_modern: bool = False
    num_residue_match: bool = True
    legacy_num_residue: int = 0
    modern_num_residue: int = 0
    mismatched_entries: List[Dict] = field(default_factory=list)
    legacy_seidx: List[ResidueIndexEntry] = field(default_factory=list)
    modern_seidx: List[ResidueIndexEntry] = field(default_factory=list)


def parse_seidx_array(seidx_array: List[Dict]) -> List[ResidueIndexEntry]:
    """Parse seidx array from JSON record into list of ResidueIndexEntry objects."""
    entries = []
    for entry_dict in seidx_array:
        residue_idx = entry_dict.get("residue_idx", 0)
        start_atom = entry_dict.get("start_atom", 0)
        end_atom = entry_dict.get("end_atom", 0)
        entries.append(ResidueIndexEntry(
            residue_idx=residue_idx,
            start_atom=start_atom,
            end_atom=end_atom
        ))
    return entries


def compare_residue_indices(
    legacy_records: List[Dict], modern_records: List[Dict], tolerance: float = 1e-6
) -> ResidueIndicesComparison:
    """
    Compare residue_indices records between legacy and modern.
    
    Args:
        legacy_records: List of legacy residue_indices records (should be single record)
        modern_records: List of modern residue_indices records (should be single record)
        tolerance: Tolerance for comparisons (not used for residue indices, kept for API consistency)
    
    Returns:
        ResidueIndicesComparison result
    """
    result = ResidueIndicesComparison()
    
    # residue_indices should be a single record, not a list of multiple records
    legacy_record = legacy_records[0] if legacy_records else None
    modern_record = modern_records[0] if modern_records else None
    
    if not legacy_record:
        result.extra_in_modern = True
        if modern_record:
            result.modern_num_residue = modern_record.get("num_residue", 0)
            seidx_array = modern_record.get("seidx", [])
            result.modern_seidx = parse_seidx_array(seidx_array)
        return result
    
    if not modern_record:
        result.missing_in_modern = True
        result.legacy_num_residue = legacy_record.get("num_residue", 0)
        seidx_array = legacy_record.get("seidx", [])
        result.legacy_seidx = parse_seidx_array(seidx_array)
        return result
    
    # Both records exist - compare them
    result.legacy_num_residue = legacy_record.get("num_residue", 0)
    result.modern_num_residue = modern_record.get("num_residue", 0)
    result.num_residue_match = result.legacy_num_residue == result.modern_num_residue
    
    legacy_seidx_array = legacy_record.get("seidx", [])
    modern_seidx_array = modern_record.get("seidx", [])
    
    result.legacy_seidx = parse_seidx_array(legacy_seidx_array)
    result.modern_seidx = parse_seidx_array(modern_seidx_array)
    
    # Compare seidx arrays entry by entry
    max_len = max(len(result.legacy_seidx), len(result.modern_seidx))
    
    for i in range(max_len):
        legacy_entry = result.legacy_seidx[i] if i < len(result.legacy_seidx) else None
        modern_entry = result.modern_seidx[i] if i < len(result.modern_seidx) else None
        
        if not legacy_entry:
            result.mismatched_entries.append({
                "index": i,
                "issue": "missing_in_legacy",
                "modern": {
                    "residue_idx": modern_entry.residue_idx if modern_entry else None,
                    "start_atom": modern_entry.start_atom if modern_entry else None,
                    "end_atom": modern_entry.end_atom if modern_entry else None,
                }
            })
        elif not modern_entry:
            result.mismatched_entries.append({
                "index": i,
                "issue": "missing_in_modern",
                "legacy": {
                    "residue_idx": legacy_entry.residue_idx,
                    "start_atom": legacy_entry.start_atom,
                    "end_atom": legacy_entry.end_atom,
                }
            })
        else:
            # Both entries exist - compare values
            if (legacy_entry.residue_idx != modern_entry.residue_idx or
                legacy_entry.start_atom != modern_entry.start_atom or
                legacy_entry.end_atom != modern_entry.end_atom):
                result.mismatched_entries.append({
                    "index": i,
                    "issue": "mismatch",
                    "legacy": {
                        "residue_idx": legacy_entry.residue_idx,
                        "start_atom": legacy_entry.start_atom,
                        "end_atom": legacy_entry.end_atom,
                    },
                    "modern": {
                        "residue_idx": modern_entry.residue_idx,
                        "start_atom": modern_entry.start_atom,
                        "end_atom": modern_entry.end_atom,
                    }
                })
    
    return result


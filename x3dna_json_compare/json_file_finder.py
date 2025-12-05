"""
Utility functions for finding JSON files in the new directory structure.

New structure:
- data/json/<record_type>/<PDB_ID>.json
- data/json_legacy/<record_type>/<PDB_ID>.json

Old structure (for backward compatibility):
- data/json/<PDB_ID>_<record_type>.json
- data/json/<PDB_ID>.json
"""

from pathlib import Path
from typing import Optional, Dict, List


def find_json_file(base_dir: Path, pdb_id: str, record_type: str) -> Optional[Path]:
    """
    Find a JSON file for a specific record type.
    
    Tries new structure first, then falls back to old structure for compatibility.
    
    Args:
        base_dir: Base directory (data/json or data/json_legacy)
        pdb_id: PDB identifier
        record_type: Record type (e.g., 'pdb_atoms', 'base_pair')
        
    Returns:
        Path to JSON file if found, None otherwise
    """
    # Try new structure: <record_type>/<PDB_ID>.json
    new_path = base_dir / record_type / f"{pdb_id}.json"
    if new_path.exists():
        return new_path
    
    # Fall back to old structure: <PDB_ID>_<record_type>.json
    old_path = base_dir / f"{pdb_id}_{record_type}.json"
    if old_path.exists():
        return old_path
    
    return None


def find_all_record_files(base_dir: Path, pdb_id: str, record_types: List[str]) -> Dict[str, Optional[Path]]:
    """
    Find all JSON files for a PDB across multiple record types.
    
    Args:
        base_dir: Base directory (data/json or data/json_legacy)
        pdb_id: PDB identifier
        record_types: List of record types to find
        
    Returns:
        Dictionary mapping record_type -> Path (or None if not found)
    """
    result = {}
    for record_type in record_types:
        result[record_type] = find_json_file(base_dir, pdb_id, record_type)
    return result


def get_record_types() -> List[str]:
    """Get list of all record types."""
    return [
        "pdb_atoms",
        "base_frame_calc",
        "frame_calc",
        "base_pair",
        "pair_validation",
        "distance_checks",
        "hbond_list",
        "find_bestpair_selection",
    ]


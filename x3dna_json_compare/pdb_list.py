"""
PDB list management utilities.

Handles loading PDB lists from various sources (valid_pdbs_fast.json, test sets, etc.)
"""

import json
from pathlib import Path
from typing import List, Optional


def load_valid_pdbs_fast(project_root: Path) -> List[str]:
    """Load valid_pdbs_fast.json and return list of PDB IDs.
    
    Args:
        project_root: Project root directory
        
    Returns:
        List of PDB IDs
        
    Raises:
        FileNotFoundError: If valid_pdbs_fast.json doesn't exist
        ValueError: If file is invalid or empty
    """
    fast_pdbs_file = project_root / "data" / "valid_pdbs_fast.json"
    
    if not fast_pdbs_file.exists():
        raise FileNotFoundError(
            f"{fast_pdbs_file} not found! Run create_fast_pdbs_json.py first."
        )
    
    with open(fast_pdbs_file) as f:
        pdbs_data = json.load(f)
        pdb_list = pdbs_data.get("valid_pdbs_with_atoms_and_frames", [])
    
    if not pdb_list:
        raise ValueError(f"No PDBs found in {fast_pdbs_file.name}")
    
    return pdb_list


def load_test_set(project_root: Path, size: int) -> Optional[List[str]]:
    """Load a test set of the specified size.
    
    Args:
        project_root: Project root directory
        size: Test set size (10, 50, 100, 500, 1000)
        
    Returns:
        List of PDB IDs, or None if test set doesn't exist
    """
    # Check resources/test_sets first (new location), fall back to data/test_sets
    resources_dir = project_root / "resources" / "test_sets"
    data_dir = project_root / "data" / "test_sets"
    
    test_sets_dir = resources_dir if resources_dir.exists() else data_dir
    test_set_file = test_sets_dir / f"test_set_{size}.json"
    
    if not test_set_file.exists():
        return None
    
    try:
        with open(test_set_file, 'r') as f:
            data = json.load(f)
            return data.get('pdb_ids', [])
    except Exception:
        return None


def get_pdb_list(
    project_root: Path,
    specific: Optional[List[str]] = None,
    max_count: Optional[int] = None,
    test_set: Optional[int] = None
) -> List[str]:
    """Get list of PDBs to process.
    
    Priority:
    1. Specific PDBs if provided
    2. Test set if specified
    3. All valid fast PDBs
    
    Args:
        project_root: Project root directory
        specific: Specific PDB IDs to use
        max_count: Maximum number of PDBs
        test_set: Test set size to load
        
    Returns:
        List of PDB IDs
    """
    if specific:
        return list(specific)
    
    if test_set:
        pdb_ids = load_test_set(project_root, test_set)
        if pdb_ids is None:
            raise FileNotFoundError(
                f"Test set of size {test_set} not found. "
                f"Generate it with: fp2-validate generate-test-sets"
            )
        return pdb_ids
    
    # Default: load valid_pdbs_fast.json
    pdb_ids = load_valid_pdbs_fast(project_root)
    
    if max_count:
        pdb_ids = pdb_ids[:max_count]
    
    return pdb_ids


"""
Base comparison utilities shared across all comparators.
"""

from typing import Any, Dict, List, Optional, Tuple

from ..config import Tolerance

# Type aliases
CompareResult = Tuple[bool, List[str]]
RecordKey = Tuple[str, int, str]  # (chain_id, residue_seq, insertion)


def make_residue_key(rec: Dict[str, Any]) -> Optional[RecordKey]:
    """
    Create a residue key from a record.
    
    Args:
        rec: Record with chain_id, residue_seq, and optional insertion
        
    Returns:
        Tuple key or None if required fields missing
    """
    chain = rec.get("chain_id")
    seq = rec.get("residue_seq")
    if chain is None or seq is None:
        return None
    return (chain, seq, rec.get("insertion", " "))


def build_residue_lookup(
    records: List[Dict[str, Any]],
    errors: List[str],
    label: str
) -> Dict[RecordKey, Dict[str, Any]]:
    """
    Build a lookup dict keyed by residue identity.
    
    Args:
        records: List of records to index
        errors: List to append error messages to
        label: Label for error messages ("legacy" or "modern")
        
    Returns:
        Dict mapping residue keys to records
    """
    lookup = {}
    for rec in records:
        key = make_residue_key(rec)
        if key is None:
            errors.append(f"{label} record missing chain_id or residue_seq")
            continue
        lookup[key] = rec
    return lookup


def check_required_field(
    legacy_rec: Dict[str, Any],
    modern_rec: Dict[str, Any],
    field: str,
    key: Any,
    errors: List[str]
) -> Tuple[Optional[Any], Optional[Any]]:
    """
    Check that a required field exists in both records.
    
    Args:
        legacy_rec: Legacy record
        modern_rec: Modern record
        field: Field name to check
        key: Key for error messages
        errors: List to append errors to
        
    Returns:
        Tuple of (legacy_value, modern_value) - either may be None if missing
    """
    leg_val = legacy_rec.get(field)
    mod_val = modern_rec.get(field)
    
    if leg_val is None:
        errors.append(f"Key {key} legacy missing {field}")
    if mod_val is None:
        errors.append(f"Key {key} modern missing {field}")
    
    return leg_val, mod_val


def compare_scalar(
    leg_val: Any,
    mod_val: Any,
    field: str,
    key: Any,
    errors: List[str],
    tolerance: float = 0.0
) -> None:
    """
    Compare scalar values (exact or with tolerance).
    
    Args:
        leg_val: Legacy value
        mod_val: Modern value
        field: Field name for error message
        key: Key for error message
        errors: List to append errors to
        tolerance: If > 0, allows numeric difference within tolerance
    """
    if leg_val is None or mod_val is None:
        return
    
    if tolerance > 0 and isinstance(leg_val, (int, float)):
        diff = abs(leg_val - mod_val)
        if diff > tolerance:
            errors.append(f"Key {key} {field}: {leg_val} vs {mod_val} (diff={diff})")
    elif leg_val != mod_val:
        errors.append(f"Key {key} {field}: {leg_val} vs {mod_val}")


def compare_vector(
    leg_vec: Optional[List[float]],
    mod_vec: Optional[List[float]],
    field: str,
    key: Any,
    errors: List[str],
    tolerance: float = Tolerance.COORDINATE,
    expected_len: int = 3
) -> None:
    """
    Compare two vectors element-wise with tolerance.
    
    Args:
        leg_vec: Legacy vector
        mod_vec: Modern vector
        field: Field name for error message
        key: Key for error message
        errors: List to append errors to
        tolerance: Maximum allowed difference per element
        expected_len: Expected vector length
    """
    if leg_vec is None or mod_vec is None:
        return
    
    if len(leg_vec) != expected_len or len(mod_vec) != expected_len:
        errors.append(f"Key {key} {field} wrong length: {len(leg_vec)} vs {len(mod_vec)}")
        return
    
    max_diff = max(abs(leg_vec[i] - mod_vec[i]) for i in range(expected_len))
    if max_diff > tolerance:
        errors.append(f"Key {key} {field} max_diff: {max_diff}")


def compare_matrix(
    leg_mat: Optional[List[List[float]]],
    mod_mat: Optional[List[List[float]]],
    field: str,
    key: Any,
    errors: List[str],
    tolerance: float = Tolerance.MATRIX
) -> None:
    """
    Compare two 3x3 matrices element-wise with tolerance.
    
    Args:
        leg_mat: Legacy matrix (3x3)
        mod_mat: Modern matrix (3x3)
        field: Field name for error message
        key: Key for error message
        errors: List to append errors to
        tolerance: Maximum allowed difference per element
    """
    if leg_mat is None or mod_mat is None:
        return
    
    if len(leg_mat) != 3 or len(mod_mat) != 3:
        errors.append(f"Key {key} {field} not 3x3")
        return
    
    for i in range(3):
        if len(leg_mat[i]) != 3 or len(mod_mat[i]) != 3:
            errors.append(f"Key {key} {field} row {i} not length 3")
            return
    
    max_diff = 0.0
    for i in range(3):
        for j in range(3):
            diff = abs(leg_mat[i][j] - mod_mat[i][j])
            max_diff = max(max_diff, diff)
    
    if max_diff > tolerance:
        errors.append(f"Key {key} {field} max_diff: {max_diff}")


"""
Residue ID utilities for JSON comparison.

Provides functions to construct and extract res_id strings for matching records
between legacy and modern JSON outputs.

res_id format: chain_id-residue_name-residue_seq[insertion]
Examples:
    A-G-5       (chain A, guanine, position 5, no insertion)
    A-C-10A     (chain A, cytosine, position 10, insertion code A)
    B-PSU-25    (chain B, pseudouridine, position 25)
"""

from typing import Dict, Any, Optional


def get_res_id(record: Dict[str, Any]) -> str:
    """
    Extract or construct res_id from a JSON record.

    If the record has a 'res_id' field, returns it directly.
    Otherwise, constructs res_id from (chain_id, residue_name, residue_seq, insertion).

    Args:
        record: JSON record dictionary

    Returns:
        res_id string in format "chain_id-residue_name-residue_seq[insertion]"
    """
    # If res_id is already present, use it directly
    if 'res_id' in record:
        return record['res_id']

    # Construct res_id from fields
    chain = record.get('chain_id', '')
    name = record.get('residue_name', '')
    seq = record.get('residue_seq', 0)
    ins = record.get('insertion', '')

    # Trim whitespace from name (PDB format uses padded names)
    if isinstance(name, str):
        name = name.strip()

    # Trim whitespace from insertion
    if isinstance(ins, str):
        ins = ins.strip()

    # Construct res_id
    if ins:
        return f"{chain}-{name}-{seq}{ins}"
    else:
        return f"{chain}-{name}-{seq}"


def get_res_id_i(record: Dict[str, Any]) -> Optional[str]:
    """
    Extract or construct res_id_i from a pair-related JSON record.

    For records with two residues (base_pair, pair_validation, etc.),
    this extracts the res_id for the first residue (i).

    Args:
        record: JSON record dictionary

    Returns:
        res_id_i string or None if not available
    """
    # If res_id_i is already present, use it directly
    if 'res_id_i' in record:
        return record['res_id_i']

    # For step params, try res_id_1i (pair 1, residue i)
    if 'res_id_1i' in record:
        return record['res_id_1i']

    # Cannot construct - base_i/base_j are indices, not residue info
    return None


def get_res_id_j(record: Dict[str, Any]) -> Optional[str]:
    """
    Extract or construct res_id_j from a pair-related JSON record.

    For records with two residues (base_pair, pair_validation, etc.),
    this extracts the res_id for the second residue (j).

    Args:
        record: JSON record dictionary

    Returns:
        res_id_j string or None if not available
    """
    # If res_id_j is already present, use it directly
    if 'res_id_j' in record:
        return record['res_id_j']

    # For step params, try res_id_1j (pair 1, residue j)
    if 'res_id_1j' in record:
        return record['res_id_1j']

    # Cannot construct - base_i/base_j are indices, not residue info
    return None


def get_pair_res_ids(record: Dict[str, Any]) -> tuple:
    """
    Extract both res_ids from a pair record.

    Args:
        record: JSON record dictionary

    Returns:
        Tuple of (res_id_i, res_id_j) or (None, None) if not available
    """
    return (get_res_id_i(record), get_res_id_j(record))


def get_step_res_ids(record: Dict[str, Any]) -> tuple:
    """
    Extract all 4 res_ids from a step record (2 pairs x 2 residues each).

    Args:
        record: JSON record dictionary

    Returns:
        Tuple of (res_id_1i, res_id_1j, res_id_2i, res_id_2j)
        or (None, None, None, None) if not available
    """
    return (
        record.get('res_id_1i'),
        record.get('res_id_1j'),
        record.get('res_id_2i'),
        record.get('res_id_2j'),
    )


def make_pair_key(res_id_i: Optional[str], res_id_j: Optional[str]) -> Optional[tuple]:
    """
    Create a normalized pair key from two res_ids.

    The key is always (smaller_res_id, larger_res_id) for consistent comparison.

    Args:
        res_id_i: First residue id
        res_id_j: Second residue id

    Returns:
        Normalized tuple (res_id_1, res_id_2) or None if either is None
    """
    if res_id_i is None or res_id_j is None:
        return None

    # Normalize order
    if res_id_i <= res_id_j:
        return (res_id_i, res_id_j)
    else:
        return (res_id_j, res_id_i)

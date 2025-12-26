"""Residue identifier conversion utilities.

Handles conversion between different ID formats:
- DSSR format: "A.G1", "0.C530", "D.DG10"
- Slot/res_id format: "A-G-1", "0-C-530", "D-DG-10"
- Legacy format: chain + residue index
"""

import re
from typing import Optional, Tuple


def dssr_to_res_id(dssr_nt: str) -> Optional[str]:
    """Convert DSSR nucleotide ID to res_id format.

    Examples:
        "A.G1" -> "A-G-1"
        "0.C530" -> "0-C-530"
        "D.DG10" -> "D-DG-10"
        "A.G1^A" -> "A-G-1A" (alternate conformer)
        "A.G190^K" -> "A-G-190K" (alternate conformer)

    Args:
        dssr_nt: DSSR format nucleotide ID

    Returns:
        res_id format string, or None if parsing fails
    """
    if not dssr_nt or "." not in dssr_nt:
        return None

    # Pattern: chain.base+number[insertion][^altconf]
    # Chain can be letter or number
    # Base can be multi-letter (DG, DC, PSU, etc.)
    # Number can have insertion code suffix
    # ^X indicates alternate conformer (should be included in res_id)
    match = re.match(
        r"^([A-Za-z0-9]+)\.([A-Za-z]+)(\d+)([A-Za-z]?)(?:\^([A-Za-z]))?$", dssr_nt
    )

    if not match:
        return None

    chain, base, num, ins_code, alt_conf = match.groups()

    # Build res_id - include both insertion code and alternate conformer
    seq_num = num
    if ins_code:
        seq_num += ins_code
    if alt_conf:
        seq_num += alt_conf
    return f"{chain}-{base}-{seq_num}"


def res_id_to_dssr(res_id: str) -> Optional[str]:
    """Convert res_id format to DSSR nucleotide ID.

    Examples:
        "A-G-1" -> "A.G1"
        "0-C-530" -> "0.C530"
        "D-DG-10" -> "D.DG10"

    Args:
        res_id: res_id format string

    Returns:
        DSSR format string, or None if parsing fails
    """
    if not res_id:
        return None

    parts = res_id.split("-")
    if len(parts) < 3:
        return None

    chain = parts[0]
    base = parts[1]
    seq_num = "-".join(parts[2:])  # Handle negative residue numbers

    return f"{chain}.{base}{seq_num}"


def parse_res_id(res_id: str) -> Optional[Tuple[str, str, str]]:
    """Parse res_id into components.

    Args:
        res_id: Format "chain-base-seqnum" (e.g., "A-G-1")

    Returns:
        Tuple of (chain, base_type, seq_num) or None if invalid
    """
    if not res_id:
        return None

    parts = res_id.split("-")
    if len(parts) < 3:
        return None

    chain = parts[0]
    base = parts[1]
    seq_num = "-".join(parts[2:])

    return chain, base, seq_num


def parse_dssr_id(dssr_nt: str) -> Optional[Tuple[str, str, str]]:
    """Parse DSSR ID into components.

    Handles DSSR format "chain.basename" where basename may have
    insertion codes or alternate location indicators.

    Args:
        dssr_nt: DSSR format like "A.G1", "0.C530", "D.DG10"

    Returns:
        Tuple of (chain, base_type, seq_num) or None if invalid

    Examples:
        "A.G1" -> ("A", "G", "1")
        "0.C530" -> ("0", "C", "530")
        "D.DG10" -> ("D", "DG", "10")
        "A.G1^A" -> ("A", "G", "1")  # Ignores alt location
    """
    if not dssr_nt or "." not in dssr_nt:
        return None

    parts = dssr_nt.split(".")
    if len(parts) != 2:
        return None

    chain = parts[0]
    nt_part = parts[1]

    # Strip alternate location indicator (^A)
    if "^" in nt_part:
        nt_part = nt_part.split("^")[0]

    # Find where digits start
    i = 0
    while i < len(nt_part) and not nt_part[i].isdigit() and nt_part[i] != "-":
        i += 1

    if i == len(nt_part):
        return None

    base_type = nt_part[:i]
    seq_num = nt_part[i:]

    return chain, base_type, seq_num


def make_res_id(chain: str, base_type: str, seq_num: str) -> str:
    """Create res_id from components.

    Args:
        chain: Chain identifier
        base_type: Base type (single letter)
        seq_num: Sequence number (may include insertion code)

    Returns:
        res_id format string
    """
    return f"{chain}-{base_type}-{seq_num}"


def normalize_chain(chain: str) -> str:
    """Normalize chain identifier.

    Handles empty chains, whitespace, etc.
    """
    chain = chain.strip()
    return chain if chain else "A"


def extract_sequence(res_id1: str, res_id2: str) -> Optional[str]:
    """Extract two-letter sequence from pair of res_ids.

    Args:
        res_id1: First residue ID
        res_id2: Second residue ID

    Returns:
        Two-letter sequence (e.g., "GC") or None
    """
    parsed1 = parse_res_id(res_id1)
    parsed2 = parse_res_id(res_id2)

    if parsed1 is None or parsed2 is None:
        return None

    return parsed1[1] + parsed2[1]


def is_standard_wc_sequence(sequence: str) -> bool:
    """Check if sequence is a standard Watson-Crick pair.

    Standard WC: GC, CG, AU, UA (and AT, TA for DNA)
    """
    return sequence.upper() in {"GC", "CG", "AU", "UA", "AT", "TA"}


def is_canonical_pair(sequence: str) -> bool:
    """Check if sequence could form a canonical base pair.

    Includes WC and wobble pairs.
    """
    canonical = {"GC", "CG", "AU", "UA", "AT", "TA", "GU", "UG", "GT", "TG"}
    return sequence.upper() in canonical

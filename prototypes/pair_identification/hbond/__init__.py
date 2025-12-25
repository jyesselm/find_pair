"""H-bond detection and analysis for base pair identification.

This package provides:
- Consolidated H-bond pattern definitions for all LW classes
- H-bond detection with donor/acceptor validation
- Slot-based saturation tracking for NH2/O donors and acceptors
- Scoring and quality assessment
"""

from .patterns import (
    # Type aliases
    HBondPattern,

    # Pattern dictionaries
    CWW_PATTERNS,
    HBOND_PATTERNS,

    # Pattern lookup functions
    get_expected_hbonds,
    get_cww_expected,
    is_hbond_match,

    # Atom sets
    get_all_donor_atoms,
    get_all_acceptor_atoms,
    get_base_atoms,

    # Constants
    LW_CLASSES,
    STANDARD_WC_SEQUENCES,
    CANONICAL_SEQUENCES,
)

__all__ = [
    # Type aliases
    "HBondPattern",

    # Pattern dictionaries
    "CWW_PATTERNS",
    "HBOND_PATTERNS",

    # Pattern lookup functions
    "get_expected_hbonds",
    "get_cww_expected",
    "is_hbond_match",

    # Atom sets
    "get_all_donor_atoms",
    "get_all_acceptor_atoms",
    "get_base_atoms",

    # Constants
    "LW_CLASSES",
    "STANDARD_WC_SEQUENCES",
    "CANONICAL_SEQUENCES",
]

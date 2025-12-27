"""H-bond detection and analysis for base pair identification.

This package provides:
- Consolidated H-bond pattern definitions for all LW classes
- H-bond detection with donor/acceptor validation
- Slot-based saturation tracking for NH2/O donors and acceptors
- Scoring and quality assessment
"""

from .finder import HBond, HBondCandidate, HBondFinder
from .geometry import (
    ACCEPTOR_CAPACITY,
    BASE_CONNECTIVITY,
    DONOR_CAPACITY,
    HSlot,
    LPSlot,
    angle_between,
    compute_base_normal,
    normalize,
    predict_h_slots,
    predict_lp_slots,
    rotate_vector,
)
from .patterns import (
    HBOND_PATTERNS,
    HBondPattern,
    get_all_lw_classes,
    get_canonical_sequences,
    get_cww_expected,
    get_expected_hbonds,
    get_standard_wc_sequences,
    is_hbond_match,
)

__all__ = [
    # Finder classes
    "HBond",
    "HBondCandidate",
    "HBondFinder",
    # Geometry classes
    "HSlot",
    "LPSlot",
    # Geometry constants
    "DONOR_CAPACITY",
    "ACCEPTOR_CAPACITY",
    "BASE_CONNECTIVITY",
    # Geometry functions
    "normalize",
    "angle_between",
    "rotate_vector",
    "compute_base_normal",
    "predict_h_slots",
    "predict_lp_slots",
    # Pattern type aliases
    "HBondPattern",
    # Pattern dictionaries
    "HBOND_PATTERNS",
    # Pattern lookup functions
    "get_expected_hbonds",
    "get_cww_expected",
    "is_hbond_match",
    "get_all_lw_classes",
    "get_standard_wc_sequences",
    "get_canonical_sequences",
]

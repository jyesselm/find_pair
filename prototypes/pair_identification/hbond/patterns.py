"""Consolidated H-bond pattern definitions for all LW classes.

This module is the single source of truth for expected H-bond patterns
in base pair classification.
"""

from typing import Dict, List, Tuple, Set

# Type alias for H-bond pattern: (donor_atom, acceptor_atom)
HBondPattern = Tuple[str, str]


# =============================================================================
# CANONICAL WATSON-CRICK PATTERNS (cWW)
# =============================================================================

CWW_PATTERNS: Dict[str, List[HBondPattern]] = {
    # Standard Watson-Crick
    "GC": [("N1", "N3"), ("N2", "O2"), ("O6", "N4")],  # 3 H-bonds
    "CG": [("N4", "O6"), ("N3", "N1"), ("O2", "N2")],  # Reverse of GC
    "AU": [("N1", "N3"), ("N6", "O4")],                # 2 H-bonds
    "UA": [("N3", "N1"), ("O4", "N6")],                # Reverse of AU
    "AT": [("N1", "N3"), ("N6", "O4")],                # DNA
    "TA": [("N3", "N1"), ("O4", "N6")],                # DNA reverse
    # Wobble pairs
    "GU": [("N1", "O2"), ("O6", "N3")],                # G-U wobble
    "UG": [("N3", "O6"), ("O2", "N1")],                # U-G wobble
    "GT": [("N1", "O2"), ("O6", "N3")],                # G-T wobble (DNA)
    "TG": [("N3", "O6"), ("O2", "N1")],                # T-G wobble (DNA)
}


# =============================================================================
# ALL LEONTIS-WESTHOF CLASS PATTERNS
# =============================================================================

# Full patterns for all 12 LW classes
# Format: LW_CLASS -> SEQUENCE -> [(donor, acceptor), ...]

HBOND_PATTERNS: Dict[str, Dict[str, List[HBondPattern]]] = {
    "cWW": CWW_PATTERNS,

    "tWW": {
        "GC": [("N2", "O2")],
        "CG": [("O2", "N2")],
        "AU": [("N6", "O4")],
        "UA": [("O4", "N6")],
        "GU": [("N1", "O4"), ("N2", "O2")],
        "UG": [("O4", "N1"), ("O2", "N2")],
    },

    "cWH": {
        "GA": [("N2", "N7")],
        "AG": [("N7", "N2")],
        "GG": [("N2", "N7"), ("N1", "O6")],
        "AA": [("N6", "N7")],
        "CA": [("N4", "N7")],
        "AC": [("N7", "N4")],
    },

    "tWH": {
        "GA": [("O6", "N6"), ("N1", "N7")],
        "AG": [("N6", "O6"), ("N7", "N1")],
        "UA": [("O4", "N6")],
        "AU": [("N6", "O4")],
    },

    "cWS": {
        "GC": [("N2", "O2")],  # Only 1 H-bond expected
        "CG": [("O2", "N2")],
        "GA": [("N2", "N3")],
        "AG": [("N3", "N2")],
    },

    "tWS": {
        "GU": [("N2", "O2")],
        "UG": [("O2", "N2")],
        "AU": [("N6", "O2")],
        "UA": [("O2", "N6")],
    },

    "cHH": {
        "AA": [("N6", "N7")],
        "GG": [("N2", "N7")],
    },

    "tHH": {
        "AA": [("N6", "N7")],
        "GG": [("O6", "N7")],
    },

    "cHS": {
        "GA": [("N7", "N2")],
        "AG": [("N2", "N7")],
        "GG": [("N7", "N2")],
    },

    "tHS": {
        "AU": [("N7", "N6")],
        "UA": [("N6", "N7")],
        "GU": [("N7", "N1"), ("O6", "N3")],
        "UG": [("N1", "N7"), ("N3", "O6")],
    },

    "cSS": {
        "GG": [("N2", "N3"), ("N3", "N2")],
        "AA": [("N3", "N1")],
    },

    "tSS": {
        "AA": [("N1", "N1")],
        "GG": [("N2", "N2")],
    },
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_expected_hbonds(lw_class: str, sequence: str) -> List[HBondPattern]:
    """Get expected H-bond patterns for a LW class and sequence.

    Args:
        lw_class: Leontis-Westhof class (e.g., "cWW", "tWH")
        sequence: Two-letter sequence (e.g., "GC", "AU")

    Returns:
        List of (donor_atom, acceptor_atom) tuples
    """
    class_patterns = HBOND_PATTERNS.get(lw_class, {})
    return class_patterns.get(sequence.upper(), [])


def get_cww_expected(sequence: str) -> List[HBondPattern]:
    """Get expected H-bonds for canonical Watson-Crick pair.

    Convenience function for the common case.
    """
    return CWW_PATTERNS.get(sequence.upper(), [])


def is_hbond_match(
    found: HBondPattern,
    expected: List[HBondPattern],
    bidirectional: bool = True,
) -> bool:
    """Check if a found H-bond matches any expected pattern.

    Args:
        found: (donor, acceptor) tuple from detection
        expected: List of expected (donor, acceptor) patterns
        bidirectional: If True, also check reversed direction

    Returns:
        True if found matches any expected pattern
    """
    donor, acceptor = found

    for exp_donor, exp_acceptor in expected:
        if donor == exp_donor and acceptor == exp_acceptor:
            return True
        if bidirectional and donor == exp_acceptor and acceptor == exp_donor:
            return True

    return False


def get_all_donor_atoms() -> Set[str]:
    """Get set of all atoms that can be H-bond donors."""
    return {"N1", "N2", "N3", "N4", "N6", "N7", "O2'"}  # O2' for ribose


def get_all_acceptor_atoms() -> Set[str]:
    """Get set of all atoms that can be H-bond acceptors."""
    return {"N1", "N3", "N7", "O2", "O4", "O6", "O2'", "O3'", "O4'", "O5'"}


def get_base_atoms() -> Set[str]:
    """Get set of base atoms (not sugar/phosphate)."""
    return {"N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6",
            "C2", "C4", "C5", "C6", "C8"}


# List of all 12 LW classes
LW_CLASSES = ["cWW", "tWW", "cWH", "tWH", "cWS", "tWS",
              "cHH", "tHH", "cHS", "tHS", "cSS", "tSS"]

# Standard Watson-Crick sequences
STANDARD_WC_SEQUENCES = {"GC", "CG", "AU", "UA", "AT", "TA"}

# Canonical pair sequences (WC + wobble)
CANONICAL_SEQUENCES = STANDARD_WC_SEQUENCES | {"GU", "UG", "GT", "TG"}

"""Hydrogen bond patterns for Leontis-Westhof base pair classification.

This module is the single authoritative source for expected H-bond donor-acceptor
pairs across all 12 Leontis-Westhof (LW) base pair classes. Each pattern defines
the canonical hydrogen bonding interactions that characterize a specific LW class
and sequence combination.

Pattern Format:
    Each pattern is a list of (donor_atom, acceptor_atom) tuples representing
    expected hydrogen bonds. Atom names follow PDB nomenclature (N1, N2, O6, etc.).

LW Classification System:
    - Edge notation: W (Watson-Crick), H (Hoogsteen), S (Sugar)
    - Orientation: c (cis), t (trans)
    - Classes: cWW, tWW, cWH, tWH, cWS, tWS, cHH, tHH, cHS, tHS, cSS, tSS

Usage:
    from hbond.patterns import HBOND_PATTERNS, get_expected_hbonds

    # Get expected H-bonds for G-C Watson-Crick pair
    expected = get_expected_hbonds("cWW", "GC")
    # Returns: [("N1", "N3"), ("N2", "O2"), ("O6", "N4")]
"""

from typing import Dict, List, Tuple, Set

# Type alias for clarity
HBondPattern = Tuple[str, str]


# =============================================================================
# WATSON-CRICK CIS (cWW) - Canonical Base Pairing
# =============================================================================
# The most common base pairing mode in double helices. Watson-Crick edges
# of both bases face each other in cis orientation.
#
# Key patterns:
# - GC: 3 H-bonds (strongest canonical pair)
# - AU/AT: 2 H-bonds (standard RNA/DNA pairing)
# - GU: Wobble pair (2 H-bonds, less stable than GC/AU)
# - Non-canonical: GA sheared, AA imino, etc.

HBOND_PATTERNS: Dict[str, Dict[str, List[HBondPattern]]] = {
    "cWW": {
        # Standard Watson-Crick pairs (RNA)
        "GC": [("N1", "N3"), ("N2", "O2"), ("O6", "N4")],  # 3 H-bonds
        "CG": [("N4", "O6"), ("N3", "N1"), ("O2", "N2")],
        "AU": [("N6", "O4"), ("N1", "N3")],  # 2 H-bonds
        "UA": [("O4", "N6"), ("N3", "N1")],

        # DNA Watson-Crick
        "AT": [("N6", "O4"), ("N1", "N3")],
        "TA": [("O4", "N6"), ("N3", "N1")],

        # Wobble pairs (G-U/T)
        "GU": [("N1", "O2")],  # O6-N3 optional
        "UG": [("O2", "N1")],
        "GT": [("N1", "O2")],
        "TG": [("O2", "N1")],

        # Non-canonical cWW (sheared, imino pairs)
        "GA": [("O6", "N6"), ("N1", "N1")],  # Sheared G-A
        "AG": [("N6", "O6"), ("N1", "N1")],
        "AA": [("N6", "N1"), ("N1", "N6")],  # A-A imino
        "GG": [("O6", "N1"), ("N1", "O6")],  # G-G pairs
        "UC": [("N3", "N3")],  # U-C pairs
        "CU": [("N3", "N3")],
        "UU": [("N3", "O4"), ("O4", "N3")],  # U-U pairs
        "CC": [("N4", "N3"), ("N3", "N4")],  # C-C pairs
    },

    # =============================================================================
    # WATSON-CRICK TRANS (tWW)
    # =============================================================================
    # Watson-Crick edges face each other in trans orientation. Common in
    # parallel duplexes and certain tertiary motifs.

    "tWW": {
        "GC": [("N2", "O2"), ("N1", "N3")],
        "CG": [("O2", "N2"), ("N3", "N1")],
        "AU": [("N6", "O4")],
        "UA": [("O4", "N6")],
    },

    # =============================================================================
    # WATSON-HOOGSTEEN CIS (cWH)
    # =============================================================================
    # Watson-Crick edge of first base pairs with Hoogsteen edge of second base
    # in cis orientation. Common in loop regions and tertiary interactions.

    "cWH": {
        "AU": [("N6", "N7")],
        "UA": [("O4", "N6")],
        "GC": [("N2", "N7")],
    },

    # =============================================================================
    # WATSON-HOOGSTEEN TRANS (tWH)
    # =============================================================================
    # Watson-Crick edge pairs with Hoogsteen edge in trans orientation.
    # Frequently seen in base triples and higher-order structures.

    "tWH": {
        "AU": [("N6", "N7"), ("N1", "N6")],
        "UA": [("N3", "N7")],
        "GU": [("N1", "O6"), ("O6", "N7")],
        "UG": [("N3", "N7")],
        "GG": [("N1", "N7"), ("O6", "N2")],
    },

    # =============================================================================
    # WATSON-SUGAR CIS (cWS)
    # =============================================================================
    # Watson-Crick edge pairs with sugar edge in cis orientation.
    # Sugar edge uses 2'-OH and other ribose groups for interaction.

    "cWS": {
        "GC": [("N2", "O2")],
        "CG": [("O2", "N2")],
        "AU": [("N6", "O2")],
        "UA": [("O2", "N6")],
    },

    # =============================================================================
    # WATSON-SUGAR TRANS (tWS)
    # =============================================================================
    # Watson-Crick edge pairs with sugar edge in trans orientation.

    "tWS": {
        "GU": [("N2", "O2")],
        "UG": [("O2", "N2")],
        "GA": [("N2", "N3")],
        "AG": [("N3", "N2")],
    },

    # =============================================================================
    # HOOGSTEEN-HOOGSTEEN CIS (cHH)
    # =============================================================================
    # Both bases use Hoogsteen edges in cis orientation.
    # Hoogsteen edge uses N7 (purines) for interactions.

    "cHH": {
        "AA": [("N6", "N7")],
        "GG": [("N2", "N7"), ("O6", "N7")],
    },

    # =============================================================================
    # HOOGSTEEN-HOOGSTEEN TRANS (tHH)
    # =============================================================================
    # Both bases use Hoogsteen edges in trans orientation.

    "tHH": {
        "AA": [("N6", "N7")],
        "GG": [("N2", "N7")],
    },

    # =============================================================================
    # HOOGSTEEN-SUGAR CIS (cHS)
    # =============================================================================
    # Hoogsteen edge pairs with sugar edge in cis orientation.

    "cHS": {
        "GA": [("N2", "N3")],
        "AG": [("N7", "N6")],
        "GG": [("N2", "N3")],
    },

    # =============================================================================
    # HOOGSTEEN-SUGAR TRANS (tHS) - Sheared Pairs
    # =============================================================================
    # Hoogsteen edge pairs with sugar edge in trans orientation.
    # Includes sheared G-A pairs common in GNRA tetraloops.

    "tHS": {
        "GA": [("N2", "N7"), ("N3", "N6")],  # Sheared G-A
        "AG": [("N7", "N2"), ("N6", "N3")],
        "AA": [("N6", "N1")],
        "GG": [("N2", "N7")],
    },

    # =============================================================================
    # SUGAR-SUGAR CIS (cSS)
    # =============================================================================
    # Both bases use sugar edges in cis orientation.

    "cSS": {
        "GC": [("N2", "O2")],
        "AU": [("O2", "O2")],  # Sugar-sugar interaction
    },

    # =============================================================================
    # SUGAR-SUGAR TRANS (tSS)
    # =============================================================================
    # Both bases use sugar edges in trans orientation.

    "tSS": {
        "GA": [("N3", "N3")],
        "AG": [("N3", "N3")],
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
    """Get expected H-bonds for Watson-Crick cis pair.

    Convenience function for the most common pairing mode.

    Args:
        sequence: Two-letter sequence (e.g., "GC", "AU")

    Returns:
        List of (donor_atom, acceptor_atom) tuples for cWW pattern
    """
    return get_expected_hbonds("cWW", sequence)


def is_hbond_match(
    found: HBondPattern,
    expected: List[HBondPattern],
    bidirectional: bool = True,
) -> bool:
    """Check if found H-bond matches any expected pattern.

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


def get_all_lw_classes() -> List[str]:
    """Get list of all 12 Leontis-Westhof classes."""
    return ["cWW", "tWW", "cWH", "tWH", "cWS", "tWS",
            "cHH", "tHH", "cHS", "tHS", "cSS", "tSS"]


def get_standard_wc_sequences() -> Set[str]:
    """Get standard Watson-Crick sequences (no wobbles)."""
    return {"GC", "CG", "AU", "UA", "AT", "TA"}


def get_canonical_sequences() -> Set[str]:
    """Get canonical sequences including wobbles."""
    return get_standard_wc_sequences() | {"GU", "UG", "GT", "TG"}

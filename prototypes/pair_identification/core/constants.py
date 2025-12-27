"""Constants for nucleic acid base pair identification.

Centralizes definitions for ring atoms, base atoms, glycosidic nitrogen
positions, and other structural constants used throughout the package.
"""
from typing import Dict, FrozenSet

# =============================================================================
# Ring Atom Definitions
# =============================================================================

PURINE_RING_ATOMS: FrozenSet[str] = frozenset({
    "N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"
})
"""Ring atoms for purine bases (A, G, I) - 9-membered fused ring system."""

PYRIMIDINE_RING_ATOMS: FrozenSet[str] = frozenset({
    "N1", "C2", "N3", "C4", "C5", "C6"
})
"""Ring atoms for pyrimidine bases (C, U, T) - 6-membered ring."""

# =============================================================================
# Glycosidic Nitrogen (N1 for pyrimidines, N9 for purines)
# =============================================================================

GLYCOSIDIC_N: Dict[str, str] = {
    "A": "N9",
    "G": "N9",
    "I": "N9",
    "C": "N1",
    "U": "N1",
    "T": "N1",
    "P": "C5",  # Pseudouridine has C-glycosidic bond at C5
}
"""Glycosidic attachment point for each base type."""

# =============================================================================
# Base Atom Sets (exclude sugar/phosphate backbone)
# =============================================================================

BASE_ATOMS: Dict[str, FrozenSet[str]] = {
    "A": frozenset({"N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"}),
    "G": frozenset({"N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"}),
    "C": frozenset({"N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"}),
    "U": frozenset({"N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"}),
    "T": frozenset({"N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6", "C7"}),
    "I": frozenset({"N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N3", "C4"}),
    "P": frozenset({"N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"}),
}
"""Base atoms for standard nucleotides (excludes sugar/phosphate)."""

# =============================================================================
# Base Classification
# =============================================================================

PURINES: FrozenSet[str] = frozenset({"A", "G", "I"})
"""Purine bases: Adenine, Guanine, Inosine."""

PYRIMIDINES: FrozenSet[str] = frozenset({"C", "U", "T"})
"""Pyrimidine bases: Cytosine, Uracil, Thymine."""

DNA_BASES: FrozenSet[str] = frozenset({"A", "G", "C", "T"})
"""Standard DNA bases."""

RNA_BASES: FrozenSet[str] = frozenset({"A", "G", "C", "U"})
"""Standard RNA bases."""

MODIFIED_BASES: FrozenSet[str] = frozenset({"I", "P"})
"""Common modified bases: Inosine, Pseudouridine."""

# =============================================================================
# Watson-Crick Sequences
# =============================================================================

CANONICAL_WC_SEQUENCES: FrozenSet[str] = frozenset({
    "GC", "CG",  # G-C Watson-Crick
    "AU", "UA",  # A-U Watson-Crick (RNA)
    "AT", "TA",  # A-T Watson-Crick (DNA)
})
"""Canonical Watson-Crick base pair sequences."""

WOBBLE_SEQUENCES: FrozenSet[str] = frozenset({
    "GU", "UG",  # G-U wobble (RNA)
    "GT", "TG",  # G-T wobble (DNA)
})
"""Wobble base pair sequences (G-U/G-T)."""

CANONICAL_SEQUENCES: FrozenSet[str] = CANONICAL_WC_SEQUENCES | WOBBLE_SEQUENCES
"""All canonical sequences (WC + wobble)."""

# =============================================================================
# Geometric Validation Constants
# =============================================================================

MAX_DORG: float = 15.0
"""Maximum origin distance between paired bases (Angstroms)."""

MAX_D_V: float = 2.5
"""Maximum vertical distance along helix axis (Angstroms)."""

MAX_PLANE_ANGLE: float = 65.0
"""Maximum angle between base plane normals (degrees)."""

MIN_DNN: float = 4.5
"""Minimum N1/N9 glycosidic nitrogen distance (Angstroms)."""

OVERLAP_THRESHOLD: float = 0.01
"""Threshold for base overlap detection."""

D_V_WEIGHT: float = 1.5
"""Weight for vertical distance in quality score calculation."""

PLANE_ANGLE_DIVISOR: float = 180.0
"""Divisor for plane angle in quality score calculation."""

# =============================================================================
# Modified Base Normalization Map
# =============================================================================

MODIFIED_BASE_MAP: Dict[str, str] = {
    # Guanine modifications
    "2MG": "G", "7MG": "G", "M2G": "G", "OMG": "G", "YG": "G",
    # Uracil modifications
    "PSU": "U", "H2U": "U", "5MU": "U", "4SU": "U",
    # Cytosine modifications
    "5MC": "C", "OMC": "C",
    # Adenine modifications
    "1MA": "A", "MIA": "A", "6MA": "A",
    # Inosine (purine)
    "INO": "I",
}
"""Mapping from modified base names to parent base types."""

# =============================================================================
# DNA Prefix Map (DA, DG, DC, DT -> A, G, C, T)
# =============================================================================

DNA_PREFIX_MAP: Dict[str, str] = {
    "DA": "A",
    "DG": "G",
    "DC": "C",
    "DT": "T",
}
"""Mapping from DNA residue names to single-letter codes."""

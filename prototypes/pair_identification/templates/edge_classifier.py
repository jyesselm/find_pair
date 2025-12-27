"""Edge-based classification for Leontis-Westhof base pairs.

This module groups LW classes by their edge pairs (WW, WH, WS, HH, HS, SS)
since cis/trans variants within an edge pair are geometrically similar and
cannot be reliably distinguished by ring atom RMSD alone.

The cis/trans distinction is determined separately by glycosidic bond
orientation (C1' positions relative to the base pair plane).
"""

from typing import Dict, List, Optional, Tuple

import numpy as np

# Map edge pairs to their cis/trans variants
EDGE_PAIR_CLASSES: Dict[str, List[str]] = {
    "WW": ["cWW", "tWW"],
    "WH": ["cWH", "tWH"],
    "WS": ["cWS", "tWS"],
    "HW": ["cHW", "tHW"],
    "HH": ["cHH", "tHH"],
    "HS": ["cHS", "tHS"],
    "SW": ["cSW", "tSW"],
    "SH": ["cSH", "tSH"],
    "SS": ["cSS", "tSS"],
}

# Reverse mapping: LW class -> edge pair
LW_TO_EDGE_PAIR: Dict[str, str] = {}
for edge_pair, lw_classes in EDGE_PAIR_CLASSES.items():
    for lw_class in lw_classes:
        LW_TO_EDGE_PAIR[lw_class] = edge_pair


def get_edge_pair(lw_class: str) -> str:
    """Get the edge pair for an LW class.

    Args:
        lw_class: Leontis-Westhof class (e.g., "cWW", "tWW").

    Returns:
        Edge pair (e.g., "WW").

    Example:
        >>> get_edge_pair("cWW")
        'WW'
        >>> get_edge_pair("tWW")
        'WW'
    """
    if lw_class in LW_TO_EDGE_PAIR:
        return LW_TO_EDGE_PAIR[lw_class]

    # Try extracting from string (e.g., "cWW" -> "WW")
    if len(lw_class) >= 3:
        return lw_class[1:3]

    return lw_class


def get_lw_classes_for_edge(edge_pair: str) -> List[str]:
    """Get all LW classes for an edge pair.

    Args:
        edge_pair: Edge pair (e.g., "WW").

    Returns:
        List of LW classes (e.g., ["cWW", "tWW"]).
    """
    return EDGE_PAIR_CLASSES.get(edge_pair, [])


def is_same_edge_pair(lw1: str, lw2: str) -> bool:
    """Check if two LW classes share the same edge pair.

    Args:
        lw1: First LW class.
        lw2: Second LW class.

    Returns:
        True if same edge pair.

    Example:
        >>> is_same_edge_pair("cWW", "tWW")
        True
        >>> is_same_edge_pair("cWW", "cWS")
        False
    """
    return get_edge_pair(lw1) == get_edge_pair(lw2)


def determine_cis_trans(
    c1_prime_1: np.ndarray,
    c1_prime_2: np.ndarray,
    glycosidic_n1: np.ndarray,
    glycosidic_n2: np.ndarray,
    ring_center_1: np.ndarray,
    ring_center_2: np.ndarray,
) -> str:
    """Determine cis/trans orientation from C1' positions.

    Cis: Both C1' atoms on same side of the base pair plane
    Trans: C1' atoms on opposite sides of the base pair plane

    Args:
        c1_prime_1: C1' position of residue 1.
        c1_prime_2: C1' position of residue 2.
        glycosidic_n1: N1/N9 position of residue 1.
        glycosidic_n2: N1/N9 position of residue 2.
        ring_center_1: Ring center of residue 1.
        ring_center_2: Ring center of residue 2.

    Returns:
        "cis" or "trans".
    """
    # Compute base pair plane normal from ring centers and glycosidic Ns
    pair_center = (ring_center_1 + ring_center_2) / 2
    v1 = ring_center_1 - pair_center
    v2 = glycosidic_n1 - pair_center

    # Cross product gives plane normal
    normal = np.cross(v1, v2)
    norm = np.linalg.norm(normal)
    if norm < 1e-6:
        # Fallback: use vector between ring centers
        v2 = ring_center_2 - ring_center_1
        normal = np.cross(v1, v2)
        norm = np.linalg.norm(normal)
        if norm < 1e-6:
            return "cis"  # Default if geometry is degenerate

    normal = normal / norm

    # Project C1' positions onto normal
    c1_vec_1 = c1_prime_1 - pair_center
    c1_vec_2 = c1_prime_2 - pair_center

    proj1 = np.dot(c1_vec_1, normal)
    proj2 = np.dot(c1_vec_2, normal)

    # Same side = cis, opposite sides = trans
    if proj1 * proj2 > 0:
        return "cis"
    else:
        return "trans"


def determine_cis_trans_simple(
    c1_prime_1: np.ndarray,
    c1_prime_2: np.ndarray,
    plane_normal: np.ndarray,
    plane_center: np.ndarray,
) -> str:
    """Simplified cis/trans determination with pre-computed plane.

    Args:
        c1_prime_1: C1' position of residue 1.
        c1_prime_2: C1' position of residue 2.
        plane_normal: Normal vector of base pair plane.
        plane_center: Center point of base pair plane.

    Returns:
        "cis" or "trans".
    """
    c1_vec_1 = c1_prime_1 - plane_center
    c1_vec_2 = c1_prime_2 - plane_center

    proj1 = np.dot(c1_vec_1, plane_normal)
    proj2 = np.dot(c1_vec_2, plane_normal)

    return "cis" if proj1 * proj2 > 0 else "trans"


class EdgeClassifier:
    """Classifies base pairs by edge pair with combined RMSD scoring.

    This class wraps a TemplateAligner and provides edge-pair-based
    classification where cWW/tWW (and other cis/trans pairs) are treated
    as equivalent for RMSD comparison.

    Example:
        classifier = EdgeClassifier(aligner)
        result = classifier.classify_edge_pair(res1, res2, "WW")
        print(f"Best RMSD for WW: {result.best_rmsd}")
        print(f"Orientation: {result.cis_trans}")
    """

    def __init__(self, template_aligner):
        """Initialize with a TemplateAligner.

        Args:
            template_aligner: TemplateAligner instance for RMSD computation.
        """
        self.aligner = template_aligner

    def get_best_edge_rmsd(
        self, res1, res2, edge_pair: str
    ) -> Tuple[float, str]:
        """Get best RMSD across all LW classes for an edge pair.

        Args:
            res1: First residue.
            res2: Second residue.
            edge_pair: Edge pair to check (e.g., "WW").

        Returns:
            Tuple of (best_rmsd, best_lw_class).
        """
        lw_classes = get_lw_classes_for_edge(edge_pair)
        if not lw_classes:
            return float("inf"), ""

        best_rmsd = float("inf")
        best_lw = ""

        for lw_class in lw_classes:
            try:
                result = self.aligner.align_to_class(res1, res2, lw_class)
                if result and result.rmsd < best_rmsd:
                    best_rmsd = result.rmsd
                    best_lw = lw_class
            except Exception:
                continue

        return best_rmsd, best_lw

    def get_best_overall(
        self, res1, res2, edge_pairs: Optional[List[str]] = None
    ) -> Tuple[float, str, str]:
        """Get best RMSD across multiple edge pairs.

        Args:
            res1: First residue.
            res2: Second residue.
            edge_pairs: Edge pairs to check (default: all).

        Returns:
            Tuple of (best_rmsd, best_lw_class, best_edge_pair).
        """
        if edge_pairs is None:
            edge_pairs = list(EDGE_PAIR_CLASSES.keys())

        best_rmsd = float("inf")
        best_lw = ""
        best_edge = ""

        for edge_pair in edge_pairs:
            rmsd, lw_class = self.get_best_edge_rmsd(res1, res2, edge_pair)
            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_lw = lw_class
                best_edge = edge_pair

        return best_rmsd, best_lw, best_edge

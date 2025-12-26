"""Base pair quality scoring for cWW classification.

Computes a composite score (0-1) based on:
- RMSD to expected template
- H-bond coverage (found vs expected)
- H-bond quality (distance and alignment)

For pairs with good geometry but sparse H-bonds, recalculates H-bonds
with extended distance (5Å) and weights alignment more heavily.
"""

from typing import List, Dict, Any, Tuple, Optional
import numpy as np


# Expected H-bond counts per sequence type
EXPECTED_HBOND_COUNT = {
    "GC": 3, "CG": 3,
    "AU": 2, "UA": 2,
    "GU": 2, "UG": 2,
}

# Ideal H-bond distance range
IDEAL_DISTANCE_MIN = 2.7
IDEAL_DISTANCE_MAX = 3.2


def compute_bp_score(
    sequence: str,
    rmsd: float,
    found_hbonds: List[Dict[str, Any]],
    base_base_only: bool = True,
    res1_atoms: Optional[Dict[str, np.ndarray]] = None,
    res2_atoms: Optional[Dict[str, np.ndarray]] = None,
    interbase_angle: float = 0.0,
) -> Tuple[float, Dict[str, float]]:
    """Compute composite base pair quality score.

    Args:
        sequence: Base pair sequence (e.g., "GC", "AU")
        rmsd: RMSD to cWW template in Angstroms
        found_hbonds: List of detected H-bonds with distance, alignment, context
        base_base_only: Only count base-base H-bonds (default True)
        res1_atoms: Optional atom coordinates for extended H-bond search
        res2_atoms: Optional atom coordinates for extended H-bond search
        interbase_angle: Angle between base planes (for planarity check)

    Returns:
        Tuple of (total_score, component_scores_dict)
        total_score: 0.0 (worst) to 1.0 (best)
        components: Dict with 'rmsd', 'coverage', 'quality' scores
    """
    # Filter to base-base H-bonds if requested
    if base_base_only:
        hbonds = [hb for hb in found_hbonds if hb.get("context") == "base_base"]
    else:
        hbonds = found_hbonds

    # If geometry is good but H-bonds are sparse, try extended search
    extended_search_used = False
    extended_hbond_count = 0  # Count of H-bonds found only via extended search
    if res1_atoms is not None and res2_atoms is not None:
        expected_count = EXPECTED_HBOND_COUNT.get(sequence, 2)
        original_count = len(hbonds)
        if rmsd < 1.0 and interbase_angle < 30.0 and original_count < expected_count:
            try:
                from .extended_hbond_finder import recalculate_hbonds_if_needed
                extended_hbonds, was_recalculated = recalculate_hbonds_if_needed(
                    res1_atoms, res2_atoms, sequence, rmsd, hbonds, interbase_angle
                )
                if was_recalculated:
                    hbonds = [hb for hb in extended_hbonds if hb.get("context") == "base_base"]
                    extended_search_used = True
                    extended_hbond_count = len(hbonds) - original_count
            except ImportError:
                pass  # Extended finder not available

    expected_count = EXPECTED_HBOND_COUNT.get(sequence, 2)
    found_count = len(hbonds)

    # Component 1: RMSD score (0-1, where 1 is perfect)
    # RMSD < 0.3A = 1.0, RMSD > 1.0A = 0.0, linear in between
    if rmsd <= 0.3:
        rmsd_score = 1.0
    elif rmsd >= 1.0:
        rmsd_score = 0.0
    else:
        rmsd_score = 1.0 - (rmsd - 0.3) / 0.7

    # Component 2: Coverage score (fraction of expected H-bonds found)
    # Extended H-bonds count at 85% - they're real but stretched
    if expected_count > 0:
        normal_count = found_count - extended_hbond_count
        effective_count = normal_count + (extended_hbond_count * 0.85)
        coverage_score = min(effective_count / expected_count, 1.0)
    else:
        coverage_score = 0.0

    # Component 3: H-bond quality score (distance and alignment)
    # Key insight: if geometry is good (RMSD < 0.8A), be very lenient on H-bond distances
    # Long H-bonds with good geometry = stretched but real cWW pair
    #
    # Geometry leniency: RMSD < 0.5A = full leniency (1.0)
    #                    RMSD 0.5-0.8A = partial leniency (linear)
    #                    RMSD > 0.8A = no leniency (0.0)
    if rmsd <= 0.5:
        geometry_leniency = 1.0
    elif rmsd >= 0.8:
        geometry_leniency = 0.0
    else:
        geometry_leniency = 1.0 - (rmsd - 0.5) / 0.3

    if hbonds:
        quality_scores = []
        for hb in hbonds:
            dist = hb.get("distance", 3.0)
            alignment = hb.get("alignment", 1.0)

            # Distance score with geometry-adjusted leniency
            if IDEAL_DISTANCE_MIN <= dist <= IDEAL_DISTANCE_MAX:
                dist_score = 1.0
            elif dist < IDEAL_DISTANCE_MIN:
                # Short distances are okay
                dist_score = max(0.5, 1.0 - (IDEAL_DISTANCE_MIN - dist) / 0.5)
            else:
                # Long distances: with good geometry, accept up to 4.2A without penalty
                # With geometry_leniency=1.0: H-bonds up to 4.2A score 1.0
                # With geometry_leniency=0.0: H-bonds penalized starting at 3.2A
                lenient_max = IDEAL_DISTANCE_MAX + (1.0 * geometry_leniency)  # 3.2A to 4.2A
                if dist <= lenient_max:
                    dist_score = 1.0
                else:
                    # Gradual penalty beyond lenient max
                    dist_score = max(0.0, 1.0 - (dist - lenient_max) / 0.5)

            # Alignment score: 0-2 scale, lower is better
            if alignment <= 1.0:
                align_score = 1.0
            elif alignment >= 2.0:
                align_score = 0.0
            else:
                align_score = 1.0 - (alignment - 1.0)

            # Combine distance and alignment (weighted)
            hb_quality = 0.7 * dist_score + 0.3 * align_score
            quality_scores.append(hb_quality)

        quality_score = sum(quality_scores) / len(quality_scores)
    else:
        quality_score = 0.0

    # Combine components (weighted average)
    # RMSD: 30%, Coverage: 40%, Quality: 30%
    total_score = 0.3 * rmsd_score + 0.4 * coverage_score + 0.3 * quality_score

    components = {
        "rmsd": round(rmsd_score, 3),
        "coverage": round(coverage_score, 3),
        "quality": round(quality_score, 3),
        "extended_search": extended_search_used,
    }

    return round(total_score, 3), components


def score_to_grade(score: float) -> str:
    """Convert numeric score to letter grade.

    Args:
        score: 0.0 to 1.0 score

    Returns:
        Letter grade: A, B, C, D, or F
    """
    if score >= 0.9:
        return "A"
    elif score >= 0.8:
        return "B"
    elif score >= 0.7:
        return "C"
    elif score >= 0.6:
        return "D"
    else:
        return "F"


def format_score_report(
    sequence: str,
    rmsd: float,
    found_hbonds: List[Dict[str, Any]],
) -> str:
    """Generate a human-readable score report.

    Args:
        sequence: Base pair sequence
        rmsd: RMSD to template
        found_hbonds: Detected H-bonds

    Returns:
        Multi-line string with score breakdown
    """
    total, components = compute_bp_score(sequence, rmsd, found_hbonds)
    grade = score_to_grade(total)

    expected = EXPECTED_HBOND_COUNT.get(sequence, 2)
    found = len([hb for hb in found_hbonds if hb.get("context") == "base_base"])

    lines = [
        f"BP Score: {total:.2f} ({grade})",
        f"  RMSD:     {components['rmsd']:.2f} (RMSD={rmsd:.2f}Å)",
        f"  Coverage: {components['coverage']:.2f} ({found}/{expected} H-bonds)",
        f"  Quality:  {components['quality']:.2f}",
    ]

    return "\n".join(lines)

"""Extended H-bond finder for pairs with good geometry but poor refinement.

When a base pair has good cWW geometry (RMSD < 1.0Å) but few/no H-bonds detected
at normal thresholds, this module recalculates H-bonds with extended distance
(up to 5Å) and adjusts scoring to weight alignment more heavily as distance increases.
"""

import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np

# Add parent paths for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from hbond_optimizer.geometry import (
    normalize, compute_base_normal, predict_h_slots, predict_lp_slots,
    score_hbond_alignment, DONOR_CAPACITY, ACCEPTOR_CAPACITY
)
from hbond_optimizer.optimizer import Residue, HBondCandidate, HBond


# Base atoms that can form Watson-Crick H-bonds
WC_DONOR_ATOMS = {'N1', 'N2', 'N3', 'N4', 'N6'}
WC_ACCEPTOR_ATOMS = {'N1', 'N3', 'O2', 'O4', 'O6'}

# Map single-letter codes to full base type for DNA
# DNA uses DA, DG, DC, DT in geometry definitions
DNA_BASE_MAP = {
    'A': 'A', 'G': 'G', 'C': 'C', 'U': 'U', 'T': 'T',
    'DA': 'DA', 'DG': 'DG', 'DC': 'DC', 'DT': 'DT',
}


def normalize_base_type(base_type: str, res_id: str = "") -> str:
    """Normalize base type, detecting DNA from residue ID if needed.

    Args:
        base_type: Single letter base type (A, G, C, U, T)
        res_id: Residue ID that may contain DNA indicator (e.g., 'A-DG-5')

    Returns:
        Normalized base type (e.g., 'DG' for DNA guanine)
    """
    # If already a DNA type, return as-is
    if base_type.upper() in ('DA', 'DG', 'DC', 'DT'):
        return base_type.upper()

    # Check if residue ID indicates DNA (has 'D' prefix in residue name)
    if res_id and '-D' in res_id.upper():
        # Extract the base letter after 'D'
        parts = res_id.split('-')
        for part in parts:
            if part.upper().startswith('D') and len(part) >= 2:
                return 'D' + part[1].upper()

    return base_type.upper()


def find_extended_hbonds(
    res1_atoms: Dict[str, np.ndarray],
    res2_atoms: Dict[str, np.ndarray],
    base1_type: str,
    base2_type: str,
    max_distance: float = 5.0,
    min_alignment: float = 0.5,
    res1_id: str = "",
    res2_id: str = "",
) -> List[Dict]:
    """Find H-bonds with extended distance threshold.

    Args:
        res1_atoms: Dict mapping atom name to coordinates for residue 1
        res2_atoms: Dict mapping atom name to coordinates for residue 2
        base1_type: Base type for residue 1 (A, G, C, U, DA, DG, DC, DT)
        base2_type: Base type for residue 2 (A, G, C, U, DA, DG, DC, DT)
        max_distance: Maximum donor-acceptor distance (default 5.0Å)
        min_alignment: Minimum alignment score to accept (0-2 scale, default 0.5)
        res1_id: Residue ID for residue 1 (for DNA detection)
        res2_id: Residue ID for residue 2 (for DNA detection)

    Returns:
        List of H-bond dicts with keys: donor_atom, acceptor_atom, distance,
        alignment, context, quality_score
    """
    # Normalize base types (detect DNA if needed)
    base1_type = normalize_base_type(base1_type, res1_id)
    base2_type = normalize_base_type(base2_type, res2_id)

    hbonds = []

    # Compute base normals
    normal1 = compute_base_normal(res1_atoms)
    normal2 = compute_base_normal(res2_atoms)

    # Find candidates: res1 donors -> res2 acceptors
    for donor_atom, donor_pos in res1_atoms.items():
        donor_atom = donor_atom.strip()
        if donor_atom not in WC_DONOR_ATOMS:
            continue
        donor_key = (base1_type.upper(), donor_atom)
        if donor_key not in DONOR_CAPACITY:
            continue

        for acceptor_atom, acceptor_pos in res2_atoms.items():
            acceptor_atom = acceptor_atom.strip()
            if acceptor_atom not in WC_ACCEPTOR_ATOMS:
                continue
            acceptor_key = (base2_type.upper(), acceptor_atom)
            if acceptor_key not in ACCEPTOR_CAPACITY:
                continue

            dist = np.linalg.norm(np.array(acceptor_pos) - np.array(donor_pos))
            if dist > max_distance:
                continue

            # Compute alignment score
            h_slots = predict_h_slots(base1_type, donor_atom, res1_atoms, normal1)
            lp_slots = predict_lp_slots(base2_type, acceptor_atom, res2_atoms, normal2)

            if h_slots and lp_slots:
                _, _, alignment = score_hbond_alignment(
                    np.array(donor_pos), np.array(acceptor_pos), h_slots, lp_slots
                )
            else:
                alignment = 0.0

            # Skip if alignment is too poor
            if alignment < min_alignment:
                continue

            # Compute quality score that weights alignment more as distance increases
            quality = compute_extended_quality(dist, alignment)

            hbonds.append({
                'donor_atom': donor_atom,
                'acceptor_atom': acceptor_atom,
                'distance': round(dist, 3),
                'alignment': round(alignment, 3),
                'context': 'base_base',
                'quality_score': round(quality, 3),
            })

    # Find candidates: res2 donors -> res1 acceptors
    for donor_atom, donor_pos in res2_atoms.items():
        donor_atom = donor_atom.strip()
        if donor_atom not in WC_DONOR_ATOMS:
            continue
        donor_key = (base2_type.upper(), donor_atom)
        if donor_key not in DONOR_CAPACITY:
            continue

        for acceptor_atom, acceptor_pos in res1_atoms.items():
            acceptor_atom = acceptor_atom.strip()
            if acceptor_atom not in WC_ACCEPTOR_ATOMS:
                continue
            acceptor_key = (base1_type.upper(), acceptor_atom)
            if acceptor_key not in ACCEPTOR_CAPACITY:
                continue

            dist = np.linalg.norm(np.array(acceptor_pos) - np.array(donor_pos))
            if dist > max_distance:
                continue

            # Compute alignment score
            h_slots = predict_h_slots(base2_type, donor_atom, res2_atoms, normal2)
            lp_slots = predict_lp_slots(base1_type, acceptor_atom, res1_atoms, normal1)

            if h_slots and lp_slots:
                _, _, alignment = score_hbond_alignment(
                    np.array(donor_pos), np.array(acceptor_pos), h_slots, lp_slots
                )
            else:
                alignment = 0.0

            # Skip if alignment is too poor
            if alignment < min_alignment:
                continue

            # Compute quality score
            quality = compute_extended_quality(dist, alignment)

            hbonds.append({
                'donor_atom': donor_atom,
                'acceptor_atom': acceptor_atom,
                'distance': round(dist, 3),
                'alignment': round(alignment, 3),
                'context': 'base_base',
                'quality_score': round(quality, 3),
            })

    return hbonds


def compute_extended_quality(distance: float, alignment: float) -> float:
    """Compute quality score that weights alignment more as distance increases.

    For short distances (< 3.2Å): distance dominates (70% distance, 30% alignment)
    For medium distances (3.2-4.0Å): balanced (50% distance, 50% alignment)
    For long distances (> 4.0Å): alignment dominates (30% distance, 70% alignment)

    Extended H-bonds get mild penalties - if geometry is good, long H-bonds
    are likely real but stretched due to poor refinement.

    Args:
        distance: H-bond distance in Angstroms
        alignment: Alignment score (0-2, where 2 is perfect)

    Returns:
        Quality score 0-1 (mildly penalized for long distances)
    """
    # Distance score: 1.0 if <= 2.8Å, 0.0 if >= 5.0Å
    if distance <= 2.8:
        dist_score = 1.0
    elif distance >= 5.0:
        dist_score = 0.0
    else:
        dist_score = 1.0 - (distance - 2.8) / 2.2

    # Alignment score: normalize from 0-2 to 0-1
    align_score = alignment / 2.0

    # Adjust weights based on distance
    if distance < 3.2:
        # Short: distance-dominated (these are clearly real H-bonds)
        dist_weight = 0.7
        align_weight = 0.3
        length_penalty = 1.0  # No penalty
    elif distance < 4.0:
        # Medium: balanced, mild penalty
        dist_weight = 0.5
        align_weight = 0.5
        length_penalty = 0.9  # 10% penalty
    else:
        # Long: alignment-dominated, moderate penalty
        dist_weight = 0.3
        align_weight = 0.7
        length_penalty = 0.8  # 20% penalty

    base_quality = dist_weight * dist_score + align_weight * align_score
    return base_quality * length_penalty


def recalculate_hbonds_if_needed(
    res1_atoms: Dict[str, np.ndarray],
    res2_atoms: Dict[str, np.ndarray],
    sequence: str,
    rmsd_cww: float,
    existing_hbonds: List[Dict],
    interbase_angle: float = 0.0,
    res1_id: str = "",
    res2_id: str = "",
) -> Tuple[List[Dict], bool]:
    """Recalculate H-bonds with extended threshold if geometry is good but H-bonds are sparse.

    Args:
        res1_atoms: Atom coordinates for first residue
        res2_atoms: Atom coordinates for second residue
        sequence: Two-letter sequence (e.g., "GC", "AU")
        rmsd_cww: RMSD to cWW template in Angstroms
        existing_hbonds: Currently detected H-bonds
        interbase_angle: Angle between base planes in degrees
        res1_id: Residue ID for first residue (for DNA detection)
        res2_id: Residue ID for second residue (for DNA detection)

    Returns:
        Tuple of (hbonds_list, was_recalculated)
    """
    # Only recalculate if:
    # 1. Geometry is good (RMSD < 1.0Å)
    # 2. Planarity is acceptable (angle < 30°)
    # 3. Few or no H-bonds found at normal threshold
    if rmsd_cww >= 1.0:
        return existing_hbonds, False

    if interbase_angle >= 30.0:
        return existing_hbonds, False

    # Count base-base H-bonds
    base_base_count = len([hb for hb in existing_hbonds if hb.get('context') == 'base_base'])

    # Expected H-bonds
    expected = 3 if sequence in ('GC', 'CG') else 2

    # If we have all expected H-bonds, no need to recalculate
    if base_base_count >= expected:
        return existing_hbonds, False

    # Extract base types
    base1 = sequence[0].upper()
    base2 = sequence[1].upper()

    # Find extended H-bonds
    extended_hbonds = find_extended_hbonds(
        res1_atoms, res2_atoms, base1, base2,
        max_distance=5.0,
        min_alignment=0.1,  # Very lenient - accept poor alignment if distance is OK
        res1_id=res1_id,
        res2_id=res2_id,
    )

    if not extended_hbonds:
        return existing_hbonds, False

    # Merge with existing (keep existing short ones, add new long ones)
    merged = list(existing_hbonds)
    existing_pairs = {(hb['donor_atom'], hb['acceptor_atom']) for hb in existing_hbonds}

    for hb in extended_hbonds:
        pair = (hb['donor_atom'], hb['acceptor_atom'])
        if pair not in existing_pairs:
            merged.append(hb)

    return merged, True

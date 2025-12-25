"""Kabsch alignment and RMSD computation utilities."""

from typing import Dict, List, Optional, Tuple
import numpy as np


def kabsch_align(
    P: np.ndarray,
    Q: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute optimal rotation and translation to align P onto Q.

    Uses the Kabsch algorithm to find the rotation matrix R and
    translation vector t that minimizes RMSD between P and Q.

    Args:
        P: Nx3 array of source points to be aligned
        Q: Nx3 array of target points

    Returns:
        Tuple of (R, t, centroid_Q) where:
        - R: 3x3 rotation matrix
        - t: Translation to apply after rotation (centroid of Q)
        - centroid_P: Centroid of P (for transforming other points)
    """
    assert P.shape == Q.shape, "Point sets must have same shape"
    assert P.shape[1] == 3, "Points must be 3D"

    # Center the points
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    # Compute covariance matrix
    H = P_centered.T @ Q_centered

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Compute rotation
    R = Vt.T @ U.T

    # Handle reflection case
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    return R, centroid_Q, centroid_P


def compute_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Compute RMSD between two point sets (must be pre-aligned).

    Args:
        P: Nx3 array of points
        Q: Nx3 array of points

    Returns:
        RMSD value
    """
    assert P.shape == Q.shape
    diff = P - Q
    return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))


def align_and_compute_rmsd(
    source: np.ndarray,
    target: np.ndarray,
) -> Tuple[np.ndarray, float]:
    """Align source to target and compute RMSD.

    Args:
        source: Nx3 array of points to align
        target: Nx3 array of target points

    Returns:
        Tuple of (aligned_source, rmsd)
    """
    R, centroid_target, centroid_source = kabsch_align(source, target)

    # Transform source
    aligned = (source - centroid_source) @ R + centroid_target

    rmsd = compute_rmsd(aligned, target)
    return aligned, rmsd


def align_atom_dicts(
    source_atoms: Dict[str, np.ndarray],
    target_atoms: Dict[str, np.ndarray],
    alignment_atoms: Optional[List[str]] = None,
) -> Tuple[Dict[str, np.ndarray], float, int]:
    """Align source atoms to target atoms using common atoms.

    Args:
        source_atoms: Dict mapping atom name to coordinates
        target_atoms: Dict mapping atom name to coordinates
        alignment_atoms: Optional list of atom names to use for alignment.
                        If None, uses all common atoms.

    Returns:
        Tuple of (transformed_source, rmsd, num_aligned_atoms)
    """
    # Find common atoms
    if alignment_atoms is None:
        common = set(source_atoms.keys()) & set(target_atoms.keys())
    else:
        common = set(alignment_atoms) & set(source_atoms.keys()) & set(target_atoms.keys())

    if len(common) < 3:
        return source_atoms.copy(), float('inf'), 0

    common_list = sorted(common)

    # Build coordinate arrays
    source_coords = np.array([source_atoms[name] for name in common_list])
    target_coords = np.array([target_atoms[name] for name in common_list])

    # Compute alignment
    R, centroid_target, centroid_source = kabsch_align(source_coords, target_coords)

    # Transform ALL source atoms
    transformed = {}
    for name, coords in source_atoms.items():
        transformed[name] = (coords - centroid_source) @ R + centroid_target

    # Compute RMSD on aligned atoms
    aligned_coords = np.array([transformed[name] for name in common_list])
    rmsd = compute_rmsd(aligned_coords, target_coords)

    return transformed, rmsd, len(common_list)


def transform_points(
    points: np.ndarray,
    R: np.ndarray,
    centroid_source: np.ndarray,
    centroid_target: np.ndarray,
) -> np.ndarray:
    """Apply a transformation to a set of points.

    Args:
        points: Nx3 array of points to transform
        R: 3x3 rotation matrix
        centroid_source: Centroid of source (to center before rotation)
        centroid_target: Centroid of target (to translate after rotation)

    Returns:
        Nx3 array of transformed points
    """
    return (points - centroid_source) @ R + centroid_target

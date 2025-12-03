#!/usr/bin/env python3
"""
Manually calculate RMSD for 2YR to understand the difference.
"""

import numpy as np
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize

# Standard ring geometry (from legacy xyz_ring and modern STANDARD_RING_GEOMETRY)
# Order: C4, N3, C2, N1, C6, C5, N7, C8, N9
STANDARD_COORDS = np.array([
    [-1.265, 3.177, 0.000],  # C4
    [-2.342, 2.364, 0.001],  # N3
    [-1.999, 1.087, 0.000],  # C2
    [-0.700, 0.641, 0.000],  # N1
    [ 0.424, 1.460, 0.000],  # C6
    [ 0.071, 2.833, 0.000],  # C5
    [ 0.870, 3.969, 0.000],  # N7
    [ 0.023, 4.962, 0.000],  # C8
    [-1.289, 4.551, 0.000],  # N9
])

# 2YR atoms from 9CJI (only pyrimidine ring - indices 0-5)
EXPERIMENTAL_2YR = np.array([
    [144.086, 142.749, 180.080],  # C4
    [143.559, 141.527, 180.167],  # N3
    [144.135, 140.500, 179.503],  # C2
    [145.268, 140.736, 178.736],  # N1
    [145.803, 141.983, 178.646],  # C6
    [145.248, 143.012, 179.296],  # C5
])

# Corresponding standard coords (only pyrimidine)
STANDARD_2YR = STANDARD_COORDS[0:6]

def kabsch(P, Q):
    """
    Kabsch algorithm to find optimal rotation matrix.
    P: experimental coords
    Q: standard coords
    """
    # Center both sets
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q
    
    # Compute covariance matrix
    H = P_centered.T @ Q_centered
    
    # SVD
    U, S, Vt = np.linalg.svd(H)
    
    # Rotation matrix
    R = Vt.T @ U.T
    
    # Ensure right-handed coordinate system
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    # Apply rotation and translation
    P_fitted = (R @ P_centered.T).T + centroid_Q
    
    # Calculate RMSD
    diff = P_fitted - Q_centered
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    return rmsd, R, centroid_P, centroid_Q

# Calculate RMSD
rmsd, R, cent_P, cent_Q = kabsch(EXPERIMENTAL_2YR, STANDARD_2YR)

print("=" * 70)
print("Manual RMSD Calculation for 2YR (9CJI chain C seq 7)")
print("=" * 70)
print(f"\nNumber of atoms matched: {len(EXPERIMENTAL_2YR)}")
print(f"Atoms: C4, N3, C2, N1, C6, C5 (pyrimidine ring only)")
print(f"\nCalculated RMSD: {rmsd:.6f}")
print(f"\nComparison:")
print(f"  Expected (legacy): ≤ 0.2618")
print(f"  Modern reported:   0.593301")
print(f"  Manual calculation: {rmsd:.6f}")

if rmsd > 0.2618:
    print(f"\n⚠️  RMSD > 0.2618 - Would be REJECTED")
else:
    print(f"\n✓  RMSD ≤ 0.2618 - Would be ACCEPTED")

# Show the fitted coordinates to understand the transformation
print(f"\nExperimental centroid: {cent_P}")
print(f"Standard centroid: {cent_Q}")


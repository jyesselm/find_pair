/**
 * @file validation_constants.hpp
 * @brief Named constants for base pair validation
 *
 * These constants match the legacy X3DNA code behavior.
 */

#pragma once

namespace x3dna {
namespace algorithms {
namespace validation_constants {

// Quality score formula weights
constexpr double D_V_WEIGHT = 2.0;           // Weight for vertical distance in quality score
constexpr double PLANE_ANGLE_DIVISOR = 20.0; // Divisor for plane angle in quality score

// Distance thresholds
constexpr double DNN_FALLBACK = 1e10;     // Large fallback value when N1/N9 not found
constexpr double BOND_DISTANCE = 2.0;     // Max covalent bond distance (Angstroms)
constexpr double MIN_ATOM_DISTANCE = 0.1; // Min distance to exclude same-atom matches

// Overlap calculation constants (from legacy x3dna.h)
constexpr double XBIG = 1.0e+18; // Very large number for bounds
constexpr double GAMUT = 5.0e8;  // Range for numeric stability

// Hydrogen bond distance thresholds (from legacy cmn_fncs.c)
constexpr double HB_DEFAULT_DIST2 = 4.5;   // Default hb_dist2 for find_hydrogen_bonds
constexpr double HB_GOOD_MIN = 2.5;        // Minimum distance for "good" H-bond
constexpr double HB_GOOD_MAX = 3.5;        // Maximum distance for "good" H-bond
constexpr double HB_FILTER_MAX = 3.6;      // Max distance before filtering
constexpr double HB_NONSTANDARD_MIN = 2.6; // Min distance for non-standard H-bond filtering
constexpr double HB_NONSTANDARD_MAX = 3.2; // Max distance for non-standard H-bond filtering
constexpr int HB_LINKAGE_CONFLICT = 18;    // Linkage type indicating conflict pair

// Nucleotide RMSD threshold (from legacy residue_ident)
constexpr double NT_RMSD_CUTOFF = 0.2618;  // RMSD threshold for nucleotide identification

// Watson-Crick bonus (from legacy find_bestpair)
constexpr double WC_QUALITY_BONUS = 2.0;   // Bonus subtracted from quality score for WC pairs

} // namespace validation_constants
} // namespace algorithms
} // namespace x3dna

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
constexpr double DNN_FALLBACK = 1e10;        // Large fallback value when N1/N9 not found
constexpr double BOND_DISTANCE = 2.0;        // Max covalent bond distance (Angstroms)
constexpr double MIN_ATOM_DISTANCE = 0.1;    // Min distance to exclude same-atom matches

// Overlap calculation constants (from legacy x3dna.h)
constexpr double XBIG = 1.0e+18;             // Very large number for bounds
constexpr double GAMUT = 5.0e8;              // Range for numeric stability

} // namespace validation_constants
} // namespace algorithms
} // namespace x3dna

/**
 * @file validation_constants.hpp
 * @brief Named constants for base pair validation
 *
 * This file provides backward-compatible access to validation constants.
 * All values are sourced from config/parameters.hpp (single source of truth).
 *
 * For new code, prefer using config::Parameters::instance() or config::defaults::
 */

#pragma once

#include <x3dna/config/parameters.hpp>

namespace x3dna {
namespace algorithms {
namespace validation_constants {

// Import all defaults from the centralized parameters file
using namespace config::defaults;

// Additional constants not in parameters.hpp (algorithm-internal)
constexpr double BOND_DISTANCE = 2.0;     // Max covalent bond distance (Angstroms)
constexpr double MIN_ATOM_DISTANCE = 0.1; // Min distance to exclude same-atom matches

} // namespace validation_constants
} // namespace algorithms
} // namespace x3dna

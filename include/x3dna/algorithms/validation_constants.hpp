/**
 * @file validation_constants.hpp
 * @brief Named constants for base pair validation
 *
 * This file provides backward-compatible access to validation constants.
 * All values are sourced from config/parameters_generated.hpp (auto-generated
 * from resources/config/parameters.json).
 *
 * For new code, prefer using config::params:: directly.
 */

#pragma once

#include <x3dna/config/parameters_generated.hpp>

namespace x3dna {
namespace algorithms {
namespace validation_constants {

// Import all constants from the generated parameters
using namespace config::params;

} // namespace validation_constants
} // namespace algorithms
} // namespace x3dna

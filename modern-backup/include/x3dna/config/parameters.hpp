/**
 * @file parameters.hpp
 * @brief Centralized parameter configuration
 *
 * This header includes the auto-generated compile-time constants from parameters.json.
 *
 * Usage:
 *   #include <x3dna/config/parameters.hpp>
 *   using namespace x3dna::config::params;
 *
 *   if (distance > MAX_DORG) { ... }
 *   if (rmsd > NT_RMSD_CUTOFF) { ... }
 *
 * To modify constants:
 *   1. Edit resources/config/parameters.json
 *   2. Rebuild (CMake automatically regenerates parameters_generated.hpp)
 *
 * The same parameters.json is also readable from Python:
 *   from x3dna_json_compare.parameters import params
 *   print(params.validation.distance.max_dorg)
 */

#pragma once

#include <x3dna/config/parameters_generated.hpp>

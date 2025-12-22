/**
 * @file common_types.hpp
 * @brief Common types and enumerations for X3DNA
 *
 * Note: Core types (BasePairType) are defined in their
 * respective headers (base_pair.hpp). This file only
 * contains types that aren't tied to a specific class.
 */

#pragma once

#include <cstddef>
#include <string>

namespace x3dna {
namespace algorithms {

/**
 * @enum PairFindingStrategy
 * @brief Strategy for finding base pairs
 */
enum class PairFindingStrategy {
    BEST_PAIR,     // Greedy mutual best match
    ALL_PAIRS,     // Exhaustive search
    DISTANCE_BASED // Simple distance-based
};

} // namespace algorithms

} // namespace x3dna

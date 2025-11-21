/**
 * @file common_types.hpp
 * @brief Common types and enumerations for X3DNA
 */

#pragma once

#include <cstddef>
#include <string>

namespace x3dna {
namespace core {

/**
 * @enum ResidueType
 * @brief Type of residue
 */
enum class ResidueType { ADENINE, CYTOSINE, GUANINE, THYMINE, URACIL, AMINO_ACID, OTHER };

/**
 * @enum PurinePyrimidine
 * @brief Purine/Pyrimidine classification
 */
enum class PurinePyrimidine { PURINE = 1, PYRIMIDINE = 0, NOT_BASE = -1 };

/**
 * @enum BasePairType
 * @brief Type of base pair
 */
enum class BasePairType { WATSON_CRICK, WOBBLE, HOOGSTEEN, REVERSE_HOOGSTEEN, UNKNOWN };

} // namespace core

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

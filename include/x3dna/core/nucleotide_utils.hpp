/**
 * @file nucleotide_utils.hpp
 * @brief Free functions for nucleotide-specific operations on Residue objects
 */

#pragma once

#include <vector>
#include <x3dna/core/atom.hpp>
#include <x3dna/core/typing/type_registry.hpp>

namespace x3dna {
namespace core {

// Forward declaration to avoid circular include
class Residue;

/**
 * @brief Extract ring atoms from a residue (base ring atoms for nucleotides)
 * @param residue The residue to extract ring atoms from
 * @return Vector of ring atoms
 */
[[nodiscard]] std::vector<Atom> ring_atoms(const Residue& residue);

/**
 * @brief Get one-letter code for a residue
 * @param residue The residue to get the code for
 * @return One-letter code (A, C, G, T, U for nucleotides, lowercase for modified)
 *
 * Recognizes standard nucleotide names:
 * - Single letter: A, C, G, T, U
 * - Three letter: ADE, CYT, GUA, THY, URA
 * - DNA format: DA, DC, DG, DT
 * - Modified nucleotides: PSU, 5MC, 5MU, H2U, A2M, etc.
 */
[[nodiscard]] char one_letter_code(const Residue& residue);

/**
 * @brief Get RY classification for a residue
 * @param residue The residue to classify
 * @return 1 for purines, 0 for pyrimidines, -1 for non-nucleotides
 */
[[nodiscard]] int ry_classification(const Residue& residue);

} // namespace core
} // namespace x3dna

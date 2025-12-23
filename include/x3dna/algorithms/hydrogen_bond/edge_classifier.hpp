/**
 * @file edge_classifier.hpp
 * @brief Leontis-Westhof edge classification for nucleotide base atoms
 *
 * Classifies which "edge" of a nucleotide base an H-bonding atom is on:
 * - Watson edge (W): The canonical Watson-Crick face
 * - Hoogsteen edge (H): The N7/major groove face (mainly purines)
 * - Sugar edge (S): The face toward the ribose sugar
 *
 * Reference: Leontis & Westhof (2001) RNA 7:499-512
 */

#pragma once

#include <string>
#include <vector>
#include <x3dna/core/hbond_types.hpp>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

/**
 * @brief Classifies nucleotide atoms by Leontis-Westhof edge
 *
 * Edge definitions:
 *
 * PURINES (A, G):
 *   Watson edge (W): N1, C2, N6 (A) / N1, C2, O6 (G)
 *   Hoogsteen edge (H): N7, C8, N6 (A) / N7, C8, O6 (G)
 *   Sugar edge (S): N3, C4, N2 (G only), O2'
 *
 * PYRIMIDINES (C, U, T):
 *   Watson edge (W): N3, C4, N4 (C) / N3, C4, O4 (U/T)
 *   Hoogsteen edge (H): C5, C6 (less common for H-bonding)
 *   Sugar edge (S): O2, N1, O2'
 *
 * Note: Some atoms like N6 (A) and O6 (G) participate in both W and H edges.
 * In such cases, we classify by the primary interaction pattern.
 */
class EdgeClassifier {
public:
    /**
     * @brief Classify which edge an atom is on
     * @param atom_name Atom name (trimmed, e.g., "N7", "O6")
     * @param base_type One-letter base code ('A', 'C', 'G', 'U', 'T')
     * @return BaseEdge classification
     *
     * Returns UNKNOWN for non-base atoms (backbone, sugar except O2')
     * or atoms that cannot participate in H-bonding.
     */
    [[nodiscard]] static core::BaseEdge classify(const std::string& atom_name, char base_type);

    /**
     * @brief Classify edge from residue name (handles modified bases)
     * @param atom_name Atom name (trimmed)
     * @param residue_name Residue name (e.g., "A", "PSU", "2MG")
     * @return BaseEdge classification
     *
     * Uses parent base type for modified nucleotides.
     */
    [[nodiscard]] static core::BaseEdge classify_from_residue(const std::string& atom_name,
                                                               const std::string& residue_name);

    /**
     * @brief Get all H-bonding atoms on a specific edge
     * @param base_type One-letter base code
     * @param edge The edge to query
     * @return Vector of atom names on that edge
     */
    [[nodiscard]] static std::vector<std::string> atoms_on_edge(char base_type, core::BaseEdge edge);

    /**
     * @brief Check if an atom is on the base (vs backbone/sugar)
     * @param atom_name Atom name
     * @return True if this is a base atom
     */
    [[nodiscard]] static bool is_base_atom(const std::string& atom_name);

private:
    /**
     * @brief Get base type from residue name
     * @param residue_name Residue name
     * @return One-letter code ('A', 'C', 'G', 'U', 'T') or '?' if unknown
     */
    [[nodiscard]] static char get_base_type(const std::string& residue_name);
};

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna

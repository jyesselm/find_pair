#pragma once

#include <string>

namespace x3dna::core::atom_classification {

/**
 * @brief Determines if an atom is part of the nucleobase.
 *
 * A nucleobase atom is one that is not part of the sugar or backbone.
 * Excludes atoms with prime notation (', *) and phosphate/hydrogen atoms.
 *
 * @param atom_name PDB-format atom name (4 characters, space-padded)
 * @return true if atom is part of the nucleobase
 */
[[nodiscard]] bool is_nucleobase_atom(const std::string& atom_name);

/**
 * @brief Determines if an atom is part of the backbone.
 *
 * Backbone atoms: P, OP1, OP2, O1P, O2P, O5', O3'
 *
 * @param atom_name PDB-format atom name (4 characters, space-padded)
 * @return true if atom is a backbone atom
 */
[[nodiscard]] bool is_backbone_atom(const std::string& atom_name);

/**
 * @brief Determines if an atom is part of the sugar.
 *
 * Sugar atoms: C1', C2', C3', C4', C5', O4', O2'
 *
 * @param atom_name PDB-format atom name (4 characters, space-padded)
 * @return true if atom is a sugar atom
 */
[[nodiscard]] bool is_sugar_atom(const std::string& atom_name);

/**
 * @brief Determines if an atom is a base atom for H-bond counting.
 *
 * Legacy pattern: " C5M" or " XD " where X is not H or P.
 * Used specifically in hydrogen bond counting algorithms.
 *
 * @param atom_name PDB-format atom name (4 characters, space-padded)
 * @return true if atom matches legacy H-bond counting pattern
 */
[[nodiscard]] bool is_base_atom_for_hbond(const std::string& atom_name);

/**
 * @brief Checks if two atoms can form a hydrogen bond.
 *
 * Verifies that both atoms have elements from the allowed set.
 * Default allowed elements: ".O.N." (oxygen and nitrogen).
 *
 * @param atom1 First atom name (PDB format)
 * @param atom2 Second atom name (PDB format)
 * @param allowed_elements String with period-delimited elements (e.g., ".O.N.")
 * @return true if both atoms can participate in H-bonding
 */
[[nodiscard]] bool can_form_hbond(
    const std::string& atom1,
    const std::string& atom2,
    const std::string& allowed_elements = ".O.N.");

/**
 * @brief Gets the legacy-compatible element index for an atom.
 *
 * Element indices (legacy asym_idx compatible):
 * - 0: Unknown
 * - 1: Carbon (C)
 * - 2: Oxygen (O)
 * - 3: Hydrogen (H)
 * - 4: Nitrogen (N)
 * - 5: Sulfur (S)
 * - 6: Phosphorus (P)
 *
 * @param atom_name PDB-format atom name (4 characters, space-padded)
 * @return Element index (0-6)
 */
[[nodiscard]] int get_element_index(const std::string& atom_name);

} // namespace x3dna::core::atom_classification

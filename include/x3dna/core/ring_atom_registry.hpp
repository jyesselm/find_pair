/**
 * @file ring_atom_registry.hpp
 * @brief Single source of truth for ring atom definitions
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>
#include <x3dna/core/residue_type.hpp>

namespace x3dna {
namespace core {

/**
 * @brief Single source of truth for ring atom definitions
 *
 * Ring atoms are the base ring atoms used for least-squares fitting
 * to calculate reference frames. Purines have 9 ring atoms (fused
 * 6+5 ring system), pyrimidines have 6 ring atoms (single 6-membered ring).
 *
 * This registry centralizes all ring atom logic that was previously
 * scattered across Atom::is_ring_atom(), RingAtomMatcher::get_ring_atom_names(),
 * and Residue::ring_atoms().
 */
class RingAtomRegistry {
public:
    /**
     * @brief Check if an atom name is a ring atom
     * @param atom_name Atom name (leading/trailing spaces are trimmed)
     * @return true if this is a base ring atom
     */
    [[nodiscard]] static bool is_ring_atom(std::string_view atom_name);

    /**
     * @brief Get ring atom names for purine bases (A, G, I)
     * @return 9 atom names: N1, C2, N3, C4, C5, C6, N7, C8, N9
     */
    [[nodiscard]] static const std::vector<std::string>& purine_atoms();

    /**
     * @brief Get ring atom names for pyrimidine bases (C, U, T, P)
     * @return 6 atom names: N1, C2, N3, C4, C5, C6
     */
    [[nodiscard]] static const std::vector<std::string>& pyrimidine_atoms();

    /**
     * @brief Get ring atom names for a residue type
     * @param type ResidueType enum value
     * @return Ring atom names for that type (purine or pyrimidine set)
     * @note Returns pyrimidine atoms for UNKNOWN type
     */
    [[nodiscard]] static const std::vector<std::string>& atoms_for_type(ResidueType type);

    /**
     * @brief Check if a residue type is a purine
     * @param type ResidueType enum value
     * @return true for ADENINE, GUANINE, INOSINE; false otherwise
     */
    [[nodiscard]] static bool is_purine(ResidueType type);

private:
    static const std::vector<std::string> PURINE_RING_ATOMS;
    static const std::vector<std::string> PYRIMIDINE_RING_ATOMS;
};

} // namespace core
} // namespace x3dna

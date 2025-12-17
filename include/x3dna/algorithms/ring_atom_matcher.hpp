/**
 * @file ring_atom_matcher.hpp
 * @brief Matches experimental ring atoms to standard template atoms
 */

#pragma once

#include <vector>
#include <string>
#include <optional>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>

namespace x3dna {
namespace algorithms {

/**
 * @struct MatchedAtoms
 * @brief Result of ring atom matching
 */
struct MatchedAtoms {
    std::vector<core::Atom> experimental; // Experimental atoms (from residue)
    std::vector<core::Atom> standard;     // Standard template atoms
    std::vector<std::string> atom_names;  // Names of matched atoms
    size_t num_matched = 0;               // Number of matched atom pairs

    [[nodiscard]] bool is_valid() const {
        // Check that we have matching experimental/standard vectors with at least 3 matched atoms
        // for fitting
        return experimental.size() == standard.size() && experimental.size() == atom_names.size() &&
               experimental.size() == num_matched && num_matched >= 3; // Minimum 3 atoms required for fitting
    }
};

/**
 * @class RingAtomMatcher
 * @brief Matches experimental residue atoms to standard base template atoms
 *
 * Ring atoms are used for least-squares fitting to calculate reference frames.
 * For purines (A, G): 9 ring atoms (C4, N3, C2, N1, C6, C5, N7, C8, N9)
 * For pyrimidines (C, T, U): 6 ring atoms (C4, N3, C2, N1, C6, C5)
 * Note: C1' is a sugar atom, not a ring atom, so it is never included.
 */
class RingAtomMatcher {
public:
    /**
     * @brief Match ring atoms between experimental residue and standard template
     * @param residue Experimental residue to match
     * @param standard_template Standard template structure
     * @param residue_type Optional residue type (if provided, overrides residue.residue_type())
     * @return MatchedAtoms structure with matched atom pairs
     */
    [[nodiscard]] static MatchedAtoms match(const core::Residue& residue, const core::Structure& standard_template,
                              std::optional<core::ResidueType> residue_type = std::nullopt);

    /**
     * @brief Get list of ring atom names for a residue type
     * @param residue_type Residue type (ADENINE, CYTOSINE, etc.)
     * @return Vector of atom names
     */
    [[nodiscard]] static std::vector<std::string> get_ring_atom_names(core::ResidueType residue_type);

private:
    /**
     * @brief Find first atom with given name in a residue
     * @param residue Residue to search
     * @param atom_name Atom name to find (must match exactly, including spaces)
     * @return Optional atom if found, nullopt otherwise
     */
    [[nodiscard]] static std::optional<core::Atom> find_atom_by_name(const core::Residue& residue, const std::string& atom_name);

    /**
     * @brief Find first atom with given name in a structure (checks all chains/residues)
     * @param structure Structure to search
     * @param atom_name Atom name to find
     * @return Optional atom if found, nullopt otherwise
     */
    [[nodiscard]] static std::optional<core::Atom> find_atom_by_name(const core::Structure& structure, const std::string& atom_name);

    /**
     * @brief Check if residue type is a purine
     * @param type Residue type
     * @return true if purine (A or G), false if pyrimidine
     */
    [[nodiscard]] static bool is_purine(core::ResidueType type);
};

} // namespace algorithms
} // namespace x3dna

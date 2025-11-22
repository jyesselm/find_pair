/**
 * @file ring_atom_matcher.hpp
 * @brief Matches experimental ring atoms to standard template atoms
 */

#pragma once

#include <vector>
#include <string>
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
    std::vector<core::Atom> experimental;  // Experimental atoms (from residue)
    std::vector<core::Atom> standard;      // Standard template atoms
    std::vector<std::string> atom_names;   // Names of matched atoms
    size_t num_matched = 0;                // Number of matched atom pairs
    
    bool is_valid() const {
        return experimental.size() == standard.size() &&
               experimental.size() == atom_names.size() &&
               num_matched >= 3; // Minimum 3 atoms required for fitting
    }
};

/**
 * @class RingAtomMatcher
 * @brief Matches experimental residue atoms to standard base template atoms
 *
 * Ring atoms are used for least-squares fitting to calculate reference frames.
 * For purines (A, G): 9 ring atoms (or 10 for RNA including C1')
 * For pyrimidines (C, T, U): 6 ring atoms (or 7 for RNA including C1')
 */
class RingAtomMatcher {
public:
    /**
     * @brief Match ring atoms between experimental residue and standard template
     * @param residue Experimental residue to match
     * @param standard_template Standard template structure
     * @param is_rna Whether this is RNA (includes C1' in matching)
     * @return MatchedAtoms structure with matched atom pairs
     */
    static MatchedAtoms match(const core::Residue& residue,
                              const core::Structure& standard_template,
                              bool is_rna = false);

    /**
     * @brief Get list of ring atom names for a residue type
     * @param residue_type Residue type (ADENINE, CYTOSINE, etc.)
     * @param is_rna Whether this is RNA (includes C1' in list)
     * @return Vector of atom names
     */
    static std::vector<std::string> get_ring_atom_names(core::ResidueType residue_type,
                                                         bool is_rna = false);

private:
    /**
     * @brief Find first atom with given name in a residue
     * @param residue Residue to search
     * @param atom_name Atom name to find (must match exactly, including spaces)
     * @return Optional atom if found, nullopt otherwise
     */
    static std::optional<core::Atom> find_atom_by_name(const core::Residue& residue,
                                                        const std::string& atom_name);

    /**
     * @brief Find first atom with given name in a structure (checks all chains/residues)
     * @param structure Structure to search
     * @param atom_name Atom name to find
     * @return Optional atom if found, nullopt otherwise
     */
    static std::optional<core::Atom> find_atom_by_name(const core::Structure& structure,
                                                        const std::string& atom_name);

    /**
     * @brief Check if residue type is a purine
     * @param type Residue type
     * @return true if purine (A or G), false if pyrimidine
     */
    static bool is_purine(core::ResidueType type);
};

} // namespace algorithms
} // namespace x3dna


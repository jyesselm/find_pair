/**
 * @file role_classifier.hpp
 * @brief Validates H-bonds and determines donor/acceptor roles
 */

#pragma once

#include <string>
#include <vector>
#include <x3dna/core/hbond_types.hpp>
#include <x3dna/core/hbond.hpp>

namespace x3dna {
namespace algorithms {

/**
 * @brief Validates H-bonds and determines donor/acceptor roles
 *
 * Extracted from legacy donor_acceptor() function logic.
 * Maps atom roles (donor/acceptor/either) for each base type.
 */
class HBondRoleClassifier {
public:
    /**
     * @brief Determine H-bond classification based on atom roles
     * @param base1 One-letter code of first base (e.g., 'A', 'G')
     * @param base2 One-letter code of second base
     * @param atom1 First atom name (4-char PDB format)
     * @param atom2 Second atom name (4-char PDB format)
     * @return Classification (STANDARD, NON_STANDARD, or INVALID)
     */
    [[nodiscard]] static core::HBondClassification classify_bond(char base1, char base2, const std::string& atom1,
                                                                 const std::string& atom2);

    /**
     * @brief Get the H-bond role of an atom in a specific base
     * @param base One-letter code (A, C, G, T, U, I)
     * @param atom_name Atom name (4-char PDB format)
     * @return Role (DONOR, ACCEPTOR, EITHER, or UNKNOWN)
     */
    [[nodiscard]] static core::HBondAtomRole get_atom_role(char base, const std::string& atom_name);

    /**
     * @brief Check if H-bond distance is in "good" range
     * @param distance H-bond distance in Angstroms
     * @param min_dist Minimum distance (default 2.5 A)
     * @param max_dist Maximum distance (default 3.5 A)
     * @return true if distance is in range
     */
    [[nodiscard]] static bool is_good_hbond_distance(double distance, double min_dist = 2.5, double max_dist = 3.5);

    /**
     * @brief Count good H-bonds (STANDARD with distance in good range)
     * @param bonds Vector of H-bonds to evaluate
     * @param min_dist Minimum distance (default 2.5 A)
     * @param max_dist Maximum distance (default 3.5 A)
     * @return Number of good H-bonds
     */
    [[nodiscard]] static int count_good_hbonds(const std::vector<core::HBond>& bonds, double min_dist = 2.5,
                                               double max_dist = 3.5);
};

} // namespace algorithms
} // namespace x3dna

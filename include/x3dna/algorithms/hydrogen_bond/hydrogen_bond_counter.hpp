/**
 * @file hydrogen_bond_counter.hpp
 * @brief Simple H-bond counting for pair validation (matches legacy check_pair counting)
 *
 * This class counts H-bonds WITHOUT validation - matches legacy's simple counting
 * in check_pair before validation filtering.
 */

#pragma once

#include <x3dna/core/residue.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <string>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

/**
 * @class HydrogenBondCounter
 * @brief Simple H-bond counting for pair validation (matches legacy check_pair counting)
 *
 * This class counts H-bonds WITHOUT validation - matches legacy's simple counting
 * in check_pair before validation filtering.
 */
class HydrogenBondCounter {
public:
    /**
     * @brief Count H-bonds simply (no validation) - matches legacy check_pair
     * @param res1 First residue
     * @param res2 Second residue
     * @param hb_lower Lower distance limit
     * @param hb_dist1 Upper distance limit
     * @param hb_atoms H-bond atom list (default ".O.N")
     * @param num_base_hb Output: count of base-base H-bonds
     * @param num_o2_hb Output: count of O2' H-bonds
     *
     * Matches legacy check_pair H-bond counting (lines 4605-4614 in cmn_fncs.c)
     * Counts H-bonds BEFORE validation - this is the key difference from validated counting
     */
    static void count_simple(const core::Residue& res1, const core::Residue& res2, double hb_lower, double hb_dist1,
                             const std::string& hb_atoms, int& num_base_hb, int& num_o2_hb);

private:
    /**
     * @brief Check if distance is within limits
     */
    static bool within_limits(const geometry::Vector3D& pos1, const geometry::Vector3D& pos2, double lower,
                              double upper);

    /**
     * @brief Check if atom is a base atom (matches legacy is_baseatom)
     */
    static bool is_base_atom(const std::string& atom_name);

    /**
     * @brief Check if two atoms can form a hydrogen bond (matches legacy good_hbatoms)
     */
    static bool good_hb_atoms(const std::string& atom1, const std::string& atom2, const std::string& hb_atoms);
};

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna

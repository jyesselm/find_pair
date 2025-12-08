/**
 * @file hydrogen_bond_counter.cpp
 * @brief Implementation of simple H-bond counting for pair validation
 */

#include <x3dna/algorithms/hydrogen_bond/hydrogen_bond_counter.hpp>
#include <x3dna/algorithms/hydrogen_bond/hydrogen_bond_utils.hpp>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

bool HydrogenBondCounter::within_limits(const geometry::Vector3D& pos1, const geometry::Vector3D& pos2, double lower,
                                        double upper) {
    geometry::Vector3D diff = pos1 - pos2;
    double dist = diff.length();
    return (dist >= lower && dist <= upper);
}

bool HydrogenBondCounter::is_base_atom(const std::string& atom_name) {
    return hydrogen_bond::is_base_atom(atom_name);
}

bool HydrogenBondCounter::good_hb_atoms(const std::string& atom1, const std::string& atom2,
                                        const std::string& hb_atoms) {
    return hydrogen_bond::good_hb_atoms(atom1, atom2, hb_atoms);
}

void HydrogenBondCounter::count_simple(const core::Residue& res1, const core::Residue& res2, double hb_lower,
                                       double hb_dist1, const std::string& hb_atoms, int& num_base_hb, int& num_o2_hb) {
    // Matches legacy check_pair H-bond counting (lines 4605-4614 in cmn_fncs.c)
    // Counts H-bonds BEFORE validation - this is the key difference

    num_base_hb = 0;
    num_o2_hb = 0;

    // Loop through all atom pairs (matches legacy nested loops)
    for (const auto& atom1 : res1.atoms()) {
        for (const auto& atom2 : res2.atoms()) {
            // Check distance using within_limits equivalent
            geometry::Vector3D diff = atom1.position() - atom2.position();
            double dist = diff.length();

            // Check if distance is in range [hb_lower, hb_dist1]
            if (dist < hb_lower || dist > hb_dist1) {
                continue;
            }

            // Check if both are base atoms and can form H-bond
            // Legacy: O2' is counted separately as num_o2_hb, not num_base_hb
            // So exclude O2' from base H-bond counting
            bool atom1_is_base = is_base_atom(atom1.name());
            bool atom2_is_base = is_base_atom(atom2.name());
            bool both_base = atom1_is_base && atom2_is_base;
            bool not_o2prime = (atom1.name() != " O2'" && atom2.name() != " O2'");

            if (both_base && not_o2prime) {
                if (good_hb_atoms(atom1.name(), atom2.name(), hb_atoms)) {
                    num_base_hb++;
                }
            }

            // Check if either atom is O2'
            if (atom1.name() == " O2'" || atom2.name() == " O2'") {
                num_o2_hb++;
            }
        }
    }
}

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna

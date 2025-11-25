/**
 * @file hydrogen_bond_finder.cpp
 * @brief Implementation of hydrogen bond finder - matches legacy get_hbond_ij
 */

#include <x3dna/algorithms/hydrogen_bond_finder.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cctype>

namespace x3dna {
namespace algorithms {

using namespace x3dna::core;
using namespace x3dna::geometry;

std::vector<HydrogenBondFinder::HydrogenBondResult>
HydrogenBondFinder::find_hydrogen_bonds(const Residue& res1, const Residue& res2, double hb_lower,
                                        double hb_dist1) {
    auto detailed = find_hydrogen_bonds_detailed(res1, res2, hb_lower, hb_dist1);
    return detailed.final_hbonds;
}

DetailedHBondResult HydrogenBondFinder::find_hydrogen_bonds_detailed(const Residue& res1,
                                                                      const Residue& res2,
                                                                      double hb_lower,
                                                                      double hb_dist1) {
    DetailedHBondResult result;

    // Step 1: Find all potential H-bonds (matches legacy get_hbond_ij initial loop)
    // Legacy: for (m = seidx[i][1]; m <= seidx[i][2]; m++)
    //         for (n = seidx[j][1]; n <= seidx[j][2]; n++)
    //             if (good_hbatoms(...) && within_limits(...))
    for (const auto& atom1 : res1.atoms()) {
        for (const auto& atom2 : res2.atoms()) {
            // Check distance using within_limits equivalent
            Vector3D diff = atom1.position() - atom2.position();
            double dist = diff.length();

            if (dist < hb_lower || dist > hb_dist1) {
                continue;
            }

            // Check if atoms can form H-bond (matches legacy good_hbatoms)
            if (good_hb_atoms(atom1.name(), atom2.name())) {
                HydrogenBondResult hbond;
                hbond.donor_atom = atom1.name();
                hbond.acceptor_atom = atom2.name();
                hbond.distance = dist;
                hbond.type = '-';       // Will be validated later
                hbond.linkage_type = 0; // Not yet implemented
                result.initial_hbonds.push_back(hbond);
            }
        }
    }

    if (result.initial_hbonds.empty()) {
        return result;
    }

    // Step 2: Resolve conflicts (matches legacy hb_atompair)
    // Legacy hb_atompair uses iterative algorithm with distance negation
    // We implement a simplified version that should match the result
    result.after_conflict_resolution = result.initial_hbonds;
    resolve_conflicts(result.after_conflict_resolution);

    // Step 3: Validate H-bonds (matches legacy validate_hbonds)
    result.after_validation = result.after_conflict_resolution;
    char base1 = res1.one_letter_code();
    char base2 = res2.one_letter_code();
    validate_hbonds(result.after_validation, base1, base2);

    // Step 4: Filter to only H-bonds with type != ' ' (matches legacy)
    for (const auto& hbond : result.after_validation) {
        if (hbond.type != ' ') {
            result.final_hbonds.push_back(hbond);

            // Count good H-bonds (type='-' and distance in [2.5, 3.5])
            if (hbond.type == '-' && hbond.distance >= 2.5 && hbond.distance <= 3.5) {
                result.num_good_hb++;
            }
        }
    }

    return result;
}

bool HydrogenBondFinder::good_hb_atoms(const std::string& atom1, const std::string& atom2) {
    // Reuse BasePairValidator's logic by creating a temporary validator
    // Note: good_hb_atoms is a member function, so we need to create an instance
    ValidationParameters params;
    params.hb_atoms = ".O.N"; // Default
    BasePairValidator validator(params);
    // Access via public interface - good_hb_atoms is public
    return validator.good_hb_atoms(atom1, atom2);
}

void HydrogenBondFinder::resolve_conflicts(std::vector<HydrogenBondResult>& hbonds) {
    // Matches legacy hb_atompair logic
    // Legacy uses iterative algorithm that marks conflicts by negating distances
    // Simplified version: if same atom has multiple H-bonds, keep shortest

    std::vector<bool> keep(hbonds.size(), true);

    for (size_t i = 0; i < hbonds.size(); ++i) {
        if (!keep[i]) {
            continue;
        }

        const auto& hbond_i = hbonds[i];

        // Check if donor atom has a shorter H-bond
        for (size_t j = 0; j < hbonds.size(); ++j) {
            if (i == j || !keep[j]) {
                continue;
            }
            const auto& hbond_j = hbonds[j];
            if (hbond_j.donor_atom == hbond_i.donor_atom && hbond_j.distance < hbond_i.distance) {
                keep[i] = false;
                break;
            }
        }

        if (!keep[i]) {
            continue;
        }

        // Check if acceptor atom has a shorter H-bond
        for (size_t j = 0; j < hbonds.size(); ++j) {
            if (i == j || !keep[j]) {
                continue;
            }
            const auto& hbond_j = hbonds[j];
            if (hbond_j.acceptor_atom == hbond_i.acceptor_atom &&
                hbond_j.distance < hbond_i.distance) {
                keep[i] = false;
                break;
            }
        }
    }

    // Remove conflicting H-bonds
    std::vector<HydrogenBondResult> resolved;
    for (size_t i = 0; i < hbonds.size(); ++i) {
        if (keep[i]) {
            resolved.push_back(hbonds[i]);
        }
    }
    hbonds = resolved;
}

void HydrogenBondFinder::validate_hbonds(std::vector<HydrogenBondResult>& hbonds, char base1,
                                         char base2) {
    // Matches legacy validate_hbonds logic
    // 1. Determine H-bond type using donor_acceptor
    // 2. Count good H-bonds (type='-' and distance in [2.5, 3.5])
    // 3. If there are good H-bonds, filter:
    //    - Remove H-bonds with distance > 3.6
    //    - Remove non-standard H-bonds (type='*') with distance outside [2.6, 3.2]
    //    - Note: Linkage type checking (lkg_type != 18) not yet implemented
    // 4. Only keep H-bonds with type != ' '

    int num_good_hb = 0;

    // First pass: determine types and count good H-bonds
    for (auto& hbond : hbonds) {
        hbond.type = donor_acceptor(base1, base2, hbond.donor_atom, hbond.acceptor_atom);
        if (hbond.type == '-' && hbond.distance >= 2.5 && hbond.distance <= 3.5) {
            num_good_hb++;
        }
    }

    // Second pass: apply filtering if there are good H-bonds
    if (num_good_hb > 0) {
        for (auto& hbond : hbonds) {
            if (hbond.type == ' ') {
                continue; // Already invalid
            }

            // Filter out H-bonds with distance > 3.6
            if (hbond.distance > 3.6) {
                hbond.type = ' ';
                continue;
            }

            // Filter out non-standard H-bonds with distance outside [2.6, 3.2]
            // Note: Linkage type checking (lkg_type != 18) not yet implemented
            if (hbond.type == '*' && (hbond.distance < 2.6 || hbond.distance > 3.2)) {
                hbond.type = ' ';
                continue;
            }
        }
    }
}

char HydrogenBondFinder::donor_acceptor(char base1, char base2, const std::string& atom1,
                                        const std::string& atom2) {
    // Use BasePairValidator's donor_acceptor function
    return BasePairValidator::donor_acceptor(base1, base2, atom1, atom2);
}

} // namespace algorithms
} // namespace x3dna

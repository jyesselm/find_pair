/**
 * @file hydrogen_bond_finder.cpp
 * @brief Implementation of hydrogen bond finder - matches legacy get_hbond_ij
 */

#include <x3dna/algorithms/hydrogen_bond_finder.hpp>
#include <x3dna/algorithms/hydrogen_bond/hydrogen_bond_utils.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cctype>
#include <iostream>

namespace x3dna {
namespace algorithms {

using namespace x3dna::core;
using namespace x3dna::geometry;

std::vector<HydrogenBondResult> HydrogenBondFinder::find_hydrogen_bonds(const Residue& res1, const Residue& res2,
                                                                        double hb_lower, double hb_dist1) {
    // Default hb_dist2 = 4.5 (typical value, can be overridden)
    auto detailed = find_hydrogen_bonds_detailed(res1, res2, hb_lower, hb_dist1, 4.5);
    return detailed.final_hbonds;
}

DetailedHBondResult HydrogenBondFinder::find_hydrogen_bonds_detailed(const Residue& res1, const Residue& res2,
                                                                     double hb_lower, double hb_dist1,
                                                                     double hb_dist2) {
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
            if (hydrogen_bond::good_hb_atoms(atom1.name(), atom2.name(), ".O.N")) {
                HydrogenBondResult hbond;
                hbond.donor_atom = atom1.name();
                hbond.acceptor_atom = atom2.name();
                hbond.distance = dist;
                hbond.type = '-';       // Will be validated later
                hbond.linkage_type = 0; // Will be calculated in resolve_conflicts
                result.initial_hbonds.push_back(hbond);
            }
        }
    }

    if (result.initial_hbonds.empty()) {
        return result;
    }

    // Step 2: Resolve conflicts (matches legacy hb_atompair)
    // Legacy hb_atompair uses iterative algorithm with distance negation and linkage type
    // calculation
    result.after_conflict_resolution = result.initial_hbonds;
    resolve_conflicts(result.after_conflict_resolution, hb_lower, hb_dist2);

    // Step 3: Validate H-bonds (matches legacy validate_hbonds)
    // Only processes H-bonds with positive distance (conflicts marked by negative distance)
    result.after_validation = result.after_conflict_resolution;

    // Get base types for H-bond validation
    // For modified nucleotides, one_letter_code() returns '?' but we need the actual base type
    // Legacy uses lowercase letters (a, c, g, t, u) for modified nucleotides, which get
    // converted to uppercase in donor_acceptor via toupper()
    // Use get_base_type_for_hbond() to handle modified nucleotides correctly
    char base1 = get_base_type_for_hbond(res1);
    char base2 = get_base_type_for_hbond(res2);

    validate_hbonds(result.after_validation, base1, base2);

    // Step 4: Filter to only H-bonds with type != ' ' for final_hbonds (used for counting)
    // NOTE: Legacy records ALL H-bonds (including type=' ') to JSON in get_hbond_ij
    // after_validation contains ALL H-bonds (including type=' ') for JSON recording
    // final_hbonds contains only type != ' ' for quality adjustment counting
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

void HydrogenBondFinder::resolve_conflicts(std::vector<HydrogenBondResult>& hbonds, double hb_lower, double hb_dist2) {
    // Matches legacy hb_atompair logic exactly
    // Uses iterative algorithm that marks conflicts by negating distances
    // Also calculates linkage types (lkg_type)

    if (hbonds.empty()) {
        return;
    }

    const size_t num_hbonds = hbonds.size();
    std::vector<bool> matched_idx(num_hbonds, false);
    std::vector<std::vector<int>> idx2(num_hbonds, std::vector<int>(2, 0));

    // Phase 1: Iterative conflict resolution (matches legacy lines 3932-3963)
    size_t num_iter = 0;
    size_t m = 0;

    while (true) {
        // Find next unmatched H-bond
        while (num_iter < num_hbonds && matched_idx[num_iter]) {
            num_iter++;
        }

        if (num_iter >= num_hbonds) {
            break;
        }

        // Find shortest H-bonds for donor and acceptor atoms of current H-bond
        // Use absolute distance for comparison (conflicts may be negated)
        double dtmp[3] = {0.0, 0.0, 0.0};
        size_t ddidx[3] = {num_hbonds, num_hbonds, num_hbonds};

        // Initialize with current H-bond (use absolute distance)
        double current_dist = std::abs(hbonds[num_iter].distance);
        dtmp[1] = current_dist;
        ddidx[1] = num_iter;
        dtmp[2] = current_dist;
        ddidx[2] = num_iter;

        // Find shorter H-bonds for same donor or acceptor
        // Note: Need to use absolute distance for comparison (conflicts may be negated)
        for (size_t n = 0; n < num_hbonds; ++n) {
            if (n == num_iter || matched_idx[n]) {
                continue;
            }

            double dist_n = std::abs(hbonds[n].distance);

            // Check if same donor atom and shorter distance
            if (hbonds[n].donor_atom == hbonds[num_iter].donor_atom && dist_n < dtmp[1]) {
                dtmp[1] = dist_n;
                ddidx[1] = n;
            }

            // Check if same acceptor atom and shorter distance
            if (hbonds[n].acceptor_atom == hbonds[num_iter].acceptor_atom && dist_n < dtmp[2]) {
                dtmp[2] = dist_n;
                ddidx[2] = n;
            }
        }

        // If donor and acceptor both point to same H-bond, mark it as conflict
        if (ddidx[1] == ddidx[2] && ddidx[1] < num_hbonds) {
            size_t k = ddidx[1];
            // Mark conflict by negating distance
            hbonds[k].distance = -hbonds[k].distance;

            // Mark all H-bonds sharing atoms with this conflict as matched
            num_iter = 0;
            for (size_t n = 0; n < num_hbonds; ++n) {
                if (matched_idx[n]) {
                    continue;
                }
                if (hbonds[n].donor_atom == hbonds[k].donor_atom ||
                    hbonds[n].acceptor_atom == hbonds[k].acceptor_atom) {
                    matched_idx[n] = true;
                    m++;
                }
            }

            if (m >= num_hbonds) {
                break;
            }
        } else {
            num_iter++;
        }
    }

    // Phase 2: Calculate linkage types (matches legacy lines 3964-3978)
    // For each conflicted H-bond (negative distance), mark which atoms it conflicts with
    for (size_t k = 0; k < num_hbonds; ++k) {
        if (hbonds[k].distance > 0.0) {
            continue; // Not a conflict
        }

        // Mark as conflict (9, 9)
        idx2[k][0] = 9;
        idx2[k][1] = 9;

        // Find all non-conflicted H-bonds that share atoms with this conflict
        for (size_t m_idx = 0; m_idx < num_hbonds; ++m_idx) {
            if (m_idx == k || hbonds[m_idx].distance < 0.0) {
                continue;
            }

            if (hbonds[m_idx].donor_atom == hbonds[k].donor_atom) {
                idx2[m_idx][0] = 1;
            }
            if (hbonds[m_idx].acceptor_atom == hbonds[k].acceptor_atom) {
                idx2[m_idx][1] = 1;
            }
        }
    }

    // Phase 3: Set linkage types and mark additional conflicts (matches legacy lines 3979-3984)
    for (size_t k = 0; k < num_hbonds; ++k) {
        int linkage_sum = idx2[k][0] + idx2[k][1];
        hbonds[k].linkage_type = linkage_sum;

        // Mark additional conflicts: if linkage_type != 18 and distance in range, negate
        if (linkage_sum != 18 && hbonds[k].distance > 0.0) {
            if (hbonds[k].distance >= hb_lower && hbonds[k].distance <= hb_dist2) {
                hbonds[k].distance = -hbonds[k].distance;
            }
        }
    }
}

void HydrogenBondFinder::validate_hbonds(std::vector<HydrogenBondResult>& hbonds, char base1, char base2) {
    // Matches legacy validate_hbonds logic exactly (lines 3989-4019 in cmn_fncs.c)
    // KEY: Legacy only processes H-bonds with NEGATIVE distance (conflicts marked by negative)
    //      Positive distances are skipped (type remains ' ')
    // 1. Initialize all types as ' '
    // 2. Only process H-bonds with negative distance (conflicts)
    // 3. Determine H-bond type using donor_acceptor
    // 4. Restore absolute distance
    // 5. Count good H-bonds (type='-' and distance in [2.5, 3.5])
    // 6. If there are good H-bonds, filter:
    //    - Remove H-bonds with distance > 3.6
    //    - Remove non-standard H-bonds (type='*') with lkg_type != 18 and distance outside
    //    [2.6, 3.2]

    int num_good_hb = 0;

    // First pass: determine types and count good H-bonds (ONLY for negative distances)
    // Legacy: if (hb_dist[k] > 0.0) continue;  // Skip positive distances!
    for (auto& hbond : hbonds) {
        hbond.type = ' '; // Initialize as invalid (matches legacy line 3994)

        // Only process H-bonds with NEGATIVE distance (conflicts marked by negative)
        // Legacy: if (hb_dist[k] > 0.0) continue;  // Positive = skip!
        if (hbond.distance > 0.0) {
            continue; // Skip positive distances, they remain type=' '
        }

        // Process negative distances (conflicts)
        hbond.type = BasePairValidator::donor_acceptor(base1, base2, hbond.donor_atom, hbond.acceptor_atom);

        // Restore absolute distance (matches legacy line 3998: hb_dist[k] = fabs(hb_dist[k]))
        hbond.distance = std::abs(hbond.distance);

        if (hbond.type == '-' && hbond.distance >= 2.5 && hbond.distance <= 3.5) {
            num_good_hb++;
        }
    }

    // Second pass: apply filtering if there are good H-bonds (matches legacy lines 4002-4011)
    if (num_good_hb > 0) {
        for (auto& hbond : hbonds) {
            if (hbond.type == ' ') {
                continue; // Already invalid
            }

            // Filter out H-bonds with distance > 3.6 (matches legacy line 4006)
            if (hbond.distance > 3.6) {
                hbond.type = ' ';
                continue;
            }

            // Filter out non-standard H-bonds with lkg_type != 18 and distance outside [2.6, 3.2]
            // (matches legacy lines 4007-4008)
            if (hbond.type == '*' && hbond.linkage_type != 18 && (hbond.distance < 2.6 || hbond.distance > 3.2)) {
                hbond.type = ' ';
                continue;
            }
        }
    }
}

char HydrogenBondFinder::get_base_type_for_hbond(const core::Residue& residue) {
    // Get base type for H-bond detection
    // Matches legacy behavior: uses one_letter_code() if available, otherwise uses residue_type()
    // Legacy assigns lowercase letters (a, c, g, t, u) to modified nucleotides, which get
    // converted to uppercase in donor_acceptor via toupper()

    char code = residue.one_letter_code();
    if (code != '?') {
        return code; // Standard nucleotide or already determined
    }

    // For modified nucleotides, use residue_type() to determine base type
    // This matches legacy's get_seq() behavior for modified nucleotides
    core::ResidueType res_type = residue.residue_type();
    switch (res_type) {
        case core::ResidueType::ADENINE:
            return 'A';
        case core::ResidueType::CYTOSINE:
            return 'C';
        case core::ResidueType::GUANINE:
            return 'G';
        case core::ResidueType::THYMINE:
            return 'T';
        case core::ResidueType::URACIL:
            return 'U';
        default:
            // For unknown/modified nucleotides, try to determine from atoms
            // Legacy uses lowercase letters for modified nucleotides
            // Check if it has purine or pyrimidine ring atoms
            bool has_purine = false, has_pyrimidine = false;
            for (const auto& atom : residue.atoms()) {
                std::string name = atom.name();
                if (name == " N9 ") {
                    has_purine = true;
                } else if (name == " N1 ") {
                    // Could be purine or pyrimidine, check for other atoms
                    bool has_c6 = false;
                    for (const auto& a2 : residue.atoms()) {
                        if (a2.name() == " C6 ") {
                            has_c6 = true;
                            break;
                        }
                    }
                    if (has_c6) {
                        has_purine = true;
                    } else {
                        has_pyrimidine = true;
                    }
                }
            }

            if (has_purine) {
                // Check for G vs A
                bool has_o6 = false, has_n6 = false;
                for (const auto& atom : residue.atoms()) {
                    if (atom.name() == " O6 ")
                        has_o6 = true;
                    if (atom.name() == " N6 ")
                        has_n6 = true;
                }
                return (has_o6 || (!has_n6)) ? 'G' : 'A';
            } else if (has_pyrimidine) {
                // Check for C vs T vs U
                bool has_n4 = false, has_c5m = false;
                for (const auto& atom : residue.atoms()) {
                    if (atom.name() == " N4 ")
                        has_n4 = true;
                    if (atom.name() == " C5M" || atom.name() == " C7 ")
                        has_c5m = true;
                }
                if (has_n4)
                    return 'C';
                if (has_c5m)
                    return 'T';
                return 'U'; // Default for pyrimidines
            }

            return '?'; // Cannot determine
    }
}

} // namespace algorithms
} // namespace x3dna

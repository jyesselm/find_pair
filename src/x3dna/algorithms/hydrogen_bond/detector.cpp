/**
 * @file detector.cpp
 * @brief Implementation of general-purpose H-bond detector
 */

#include <x3dna/algorithms/hydrogen_bond/detector.hpp>
#include <x3dna/algorithms/hydrogen_bond/geometry.hpp>
#include <x3dna/algorithms/hydrogen_bond/role_classifier.hpp>
#include <x3dna/algorithms/hydrogen_bond/hydrogen_bond_utils.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/core/typing/atom_classification.hpp>
#include <x3dna/core/typing/nucleotide_type.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <algorithm>
#include <cmath>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

using namespace x3dna::core;
using namespace x3dna::core::typing;
using namespace x3dna::geometry;

namespace {

/**
 * @brief Find the index of the H-bond with shortest distance sharing a specific atom
 */
[[nodiscard]] size_t find_shortest_hbond_sharing_atom(const std::vector<HBond>& bonds, size_t current_idx,
                                                      const std::vector<bool>& matched_idx, bool check_donor) {
    const auto& current = bonds[current_idx];
    const std::string& atom_to_match = check_donor ? current.donor_atom_name : current.acceptor_atom_name;
    double shortest_dist = current.distance;
    size_t shortest_idx = current_idx;

    for (size_t n = 0; n < bonds.size(); ++n) {
        if (n == current_idx || matched_idx[n]) {
            continue;
        }

        const auto& candidate = bonds[n];
        const std::string& candidate_atom = check_donor ? candidate.donor_atom_name : candidate.acceptor_atom_name;

        if (candidate_atom != atom_to_match) {
            continue;
        }

        if (candidate.distance < shortest_dist) {
            shortest_dist = candidate.distance;
            shortest_idx = n;
        }
    }

    return shortest_idx;
}

/**
 * @brief Mark H-bonds sharing atoms with a conflict as matched
 */
size_t mark_sharing_hbonds_as_matched(const std::vector<HBond>& bonds, size_t conflict_idx,
                                      std::vector<bool>& matched_idx) {
    const auto& conflict = bonds[conflict_idx];
    size_t count = 0;

    for (size_t n = 0; n < bonds.size(); ++n) {
        if (matched_idx[n]) {
            continue;
        }

        const bool shares_donor = bonds[n].donor_atom_name == conflict.donor_atom_name;
        const bool shares_acceptor = bonds[n].acceptor_atom_name == conflict.acceptor_atom_name;

        if (shares_donor || shares_acceptor) {
            matched_idx[n] = true;
            ++count;
        }
    }

    return count;
}

/**
 * @brief Helper to check if residue has a specific atom (trimmed name)
 */
bool has_atom(const Residue& residue, const std::string& name) {
    for (const auto& atom : residue.atoms()) {
        if (atom.name() == name) {
            return true;
        }
    }
    return false;
}

/**
 * @brief Determine purine type (A or G) from atoms
 */
char determine_purine_type(const Residue& residue) {
    const bool has_o6 = has_atom(residue, "O6");
    const bool has_n6 = has_atom(residue, "N6");
    return (has_o6 || !has_n6) ? 'G' : 'A';
}

/**
 * @brief Determine pyrimidine type (C, T, or U) from atoms
 */
char determine_pyrimidine_type(const Residue& residue) {
    if (has_atom(residue, "N4")) {
        return 'C';
    }
    if (has_atom(residue, "C5M") || has_atom(residue, "C7")) {
        return 'T';
    }
    return 'U';
}

/**
 * @brief Determine base type from atoms for unknown residues
 */
char determine_base_type_from_atoms(const Residue& residue) {
    const bool has_n9 = has_atom(residue, "N9");
    const bool has_n1 = has_atom(residue, "N1");
    const bool has_c6 = has_atom(residue, "C6");

    const bool is_purine = has_n9 || (has_n1 && has_c6);
    const bool is_pyrimidine = has_n1 && !has_c6;

    if (is_purine) {
        return determine_purine_type(residue);
    }
    if (is_pyrimidine) {
        return determine_pyrimidine_type(residue);
    }
    return '?';
}

} // namespace

HBondDetector::HBondDetector(const HBondDetectionParams& params) : params_(params) {}

std::vector<HBond> HBondDetector::detect_base_hbonds(const Residue& residue1, const Residue& residue2) const {
    auto result = detect_base_hbonds_detailed(residue1, residue2);
    return result.final_bonds;
}

HBondPipelineResult HBondDetector::detect_base_hbonds_detailed(const Residue& residue1,
                                                               const Residue& residue2) const {
    return detect_internal(residue1, residue2, true, MoleculeType::NUCLEIC_ACID, MoleculeType::NUCLEIC_ACID);
}

std::vector<HBond> HBondDetector::detect_all_hbonds_between(const Residue& residue1, const Residue& residue2,
                                                            MoleculeType mol1_type, MoleculeType mol2_type) const {
    auto result = detect_all_hbonds_detailed(residue1, residue2, mol1_type, mol2_type);
    return result.final_bonds;
}

HBondPipelineResult HBondDetector::detect_all_hbonds_detailed(const Residue& residue1, const Residue& residue2,
                                                              MoleculeType mol1_type, MoleculeType mol2_type) const {
    return detect_internal(residue1, residue2, false, mol1_type, mol2_type);
}

HBondPipelineResult HBondDetector::detect_internal(const Residue& residue1, const Residue& residue2,
                                                    bool base_atoms_only, MoleculeType mol1_type,
                                                    MoleculeType mol2_type) const {
    HBondPipelineResult result;

    // Step 1: Find candidate bonds - work in place using all_classified_bonds as working vector
    auto& bonds = result.all_classified_bonds;
    bonds = find_candidate_bonds(residue1, residue2, base_atoms_only, mol1_type, mol2_type);

    if (bonds.empty()) {
        return result;
    }

    // Step 2: Resolve atom-sharing conflicts (in place)
    resolve_atom_sharing_conflicts(bonds);

    // Step 3: Get base types for classification (for nucleic acids)
    char base1 = '?';
    char base2 = '?';
    if (mol1_type == MoleculeType::NUCLEIC_ACID) {
        base1 = get_base_type_for_hbond(residue1);
    }
    if (mol2_type == MoleculeType::NUCLEIC_ACID) {
        base2 = get_base_type_for_hbond(residue2);
    }

    // Step 4: Classify bonds (in place)
    classify_bonds(bonds, base1, base2);

    // Step 5: Calculate angles for all bonds (in place)
    calculate_angles(bonds, residue1, residue2);

    // Step 6: Apply post-validation filtering (marks bonds as INVALID but doesn't remove)
    apply_post_validation_filtering(bonds);

    // Step 7: Build final_bonds by copying only valid bonds (single allocation)
    result.final_bonds.reserve(bonds.size());
    for (const auto& bond : bonds) {
        if (bond.classification != HBondClassification::INVALID) {
            result.final_bonds.push_back(bond);
            if (bond.classification == HBondClassification::STANDARD) {
                result.standard_bond_count++;
                if (HBondRoleClassifier::is_good_hbond_distance(bond.distance)) {
                    result.good_bond_count++;
                }
            }
        }
    }

    return result;
}

void HBondDetector::count_potential_hbonds(const Residue& res1, const Residue& res2, int& num_base_hb,
                                           int& num_o2_hb) const {
    num_base_hb = 0;
    num_o2_hb = 0;

    const double hb_lower = params_.distances.min_distance;
    const double hb_dist1 = params_.distances.base_base_max;

    for (const auto& a1 : res1.atoms()) {
        for (const auto& a2 : res2.atoms()) {
            const double dist = (a1.position() - a2.position()).length();

            if (dist < hb_lower || dist > hb_dist1) {
                continue;
            }

            // Use legacy is_base_atom logic for pair validation counting
            const bool atom1_is_base = is_base_atom(a1.name());
            const bool atom2_is_base = is_base_atom(a2.name());
            const bool both_base = atom1_is_base && atom2_is_base;
            const bool not_o2prime = (a1.name() != "O2'" && a2.name() != "O2'");

            if (both_base && not_o2prime) {
                if (good_hb_atoms(a1.name(), a2.name(), params_.allowed_elements)) {
                    num_base_hb++;
                }
            }

            if (a1.name() == "O2'" || a2.name() == "O2'") {
                num_o2_hb++;
            }
        }
    }
}

std::vector<HBond> HBondDetector::find_candidate_bonds(const Residue& residue1, const Residue& residue2,
                                                        bool base_atoms_only, MoleculeType mol1_type,
                                                        MoleculeType mol2_type) const {
    std::vector<HBond> candidates;

    for (const auto& atom1 : residue1.atoms()) {
        for (const auto& atom2 : residue2.atoms()) {
            // Check if atoms can form H-bond based on elements
            if (!good_hb_atoms(atom1.name(), atom2.name(), params_.allowed_elements)) {
                continue;
            }

            // Skip non-base atoms if requested (for nucleic acid base-base detection)
            if (base_atoms_only) {
                const bool atom1_is_base = AtomClassifier::is_nucleobase_atom(atom1.name());
                const bool atom2_is_base = AtomClassifier::is_nucleobase_atom(atom2.name());
                if (!atom1_is_base || !atom2_is_base) {
                    continue;
                }
            }

            // Calculate distance
            const double dist = (atom1.position() - atom2.position()).length();

            // Determine context for distance threshold
            const HBondContext context = HBondGeometry::determine_context(atom1.name(), atom2.name(), mol1_type,
                                                                           mol2_type);
            const double max_dist = params_.distances.max_for_context(context);

            // Check distance is in valid range
            if (dist < params_.distances.min_distance || dist > max_dist) {
                continue;
            }

            // Create H-bond
            HBond hbond;
            hbond.donor_atom_name = atom1.name();
            hbond.acceptor_atom_name = atom2.name();
            hbond.distance = dist;
            hbond.context = context;
            hbond.classification = HBondClassification::UNKNOWN;
            hbond.conflict_state = ConflictState::NO_CONFLICT;

            candidates.push_back(hbond);
        }
    }

    return candidates;
}

void HBondDetector::resolve_atom_sharing_conflicts(std::vector<HBond>& bonds) const {
    if (bonds.empty()) {
        return;
    }

    // 3-phase algorithm from legacy
    resolve_conflicts_phase1(bonds);
    resolve_conflicts_phase2(bonds);
    resolve_conflicts_phase3(bonds);
}

void HBondDetector::resolve_conflicts_phase1(std::vector<HBond>& bonds) const {
    const size_t num_bonds = bonds.size();
    std::vector<bool> matched_idx(num_bonds, false);
    size_t num_iter = 0;
    size_t total_matched = 0;

    while (num_iter < num_bonds) {
        if (matched_idx[num_iter]) {
            ++num_iter;
            continue;
        }

        // Find shortest H-bonds sharing donor/acceptor atoms
        const size_t shortest_donor = find_shortest_hbond_sharing_atom(bonds, num_iter, matched_idx, true);
        const size_t shortest_acceptor = find_shortest_hbond_sharing_atom(bonds, num_iter, matched_idx, false);

        // If both point to same H-bond, it's a valid conflict
        const bool is_conflict = (shortest_donor == shortest_acceptor) && (shortest_donor < num_bonds);

        if (!is_conflict) {
            ++num_iter;
            continue;
        }

        // Mark this as the conflict winner
        bonds[shortest_donor].conflict_state = ConflictState::IS_CONFLICT_WINNER;

        // Mark all H-bonds sharing atoms with this conflict
        const size_t newly_matched = mark_sharing_hbonds_as_matched(bonds, shortest_donor, matched_idx);
        total_matched += newly_matched;

        // Restart iteration when new matches are found
        num_iter = 0;

        if (total_matched >= num_bonds) {
            break;
        }
    }
}

void HBondDetector::resolve_conflicts_phase2(std::vector<HBond>& bonds) const {
    // Calculate linkage markers for bonds that share atoms with conflict winners
    for (size_t k = 0; k < bonds.size(); ++k) {
        if (bonds[k].conflict_state != ConflictState::IS_CONFLICT_WINNER) {
            continue;
        }

        // Find non-conflicted bonds that share atoms with this winner
        for (size_t m = 0; m < bonds.size(); ++m) {
            if (m == k || bonds[m].conflict_state == ConflictState::IS_CONFLICT_WINNER) {
                continue;
            }

            const bool shares_donor = bonds[m].donor_atom_name == bonds[k].donor_atom_name;
            const bool shares_acceptor = bonds[m].acceptor_atom_name == bonds[k].acceptor_atom_name;

            if (shares_donor && shares_acceptor) {
                bonds[m].conflict_state = ConflictState::SHARES_BOTH_WITH_WINNER;
            } else if (shares_donor) {
                bonds[m].conflict_state = ConflictState::SHARES_DONOR_WITH_WINNER;
            } else if (shares_acceptor) {
                bonds[m].conflict_state = ConflictState::SHARES_ACCEPTOR_WITH_WINNER;
            }
        }
    }
}

void HBondDetector::resolve_conflicts_phase3(std::vector<HBond>& bonds) const {
    // Phase 3 in legacy sets linkage_type = idx2[k][0] + idx2[k][1] and negates distance
    // of bonds with linkage != 18 that are in [hb_lower, hb_dist2].
    // In our implementation, we DON'T modify conflict_state here because we need to
    // preserve the original conflict status for post-validation filtering.
    // Instead, classify_bonds will handle whether to classify based on distance range.
    (void)bonds; // Suppress unused parameter warning - logic is in classify_bonds
}

void HBondDetector::classify_bonds(std::vector<HBond>& bonds, char base1_type, char base2_type) const {
    // Legacy behavior: classify bonds that either:
    // 1. Are conflict winners (IS_CONFLICT_WINNER), OR
    // 2. Are in the extended distance range [hb_lower, hb_dist2] (Phase 3 promoted bonds)
    const double hb_lower = params_.distances.min_distance;
    const double hb_dist2 = params_.distances.conflict_filter_distance;

    for (auto& bond : bonds) {
        // Skip already invalidated bonds
        if (bond.classification == HBondClassification::INVALID) {
            continue;
        }

        // Determine if bond should be classified (matches legacy's negative distance check)
        bool should_classify = false;

        if (bond.conflict_state == ConflictState::IS_CONFLICT_WINNER) {
            // Conflict winners from Phase 1
            should_classify = true;
        } else if (bond.distance >= hb_lower && bond.distance <= hb_dist2) {
            // Phase 3: promote non-winners in extended distance range
            should_classify = true;
        }

        if (!should_classify) {
            bond.classification = HBondClassification::INVALID;
            continue;
        }

        // Use legacy classification for nucleotide base-base (matches baseline exactly)
        char legacy_type = BasePairValidator::donor_acceptor(
            base1_type, base2_type, bond.donor_atom_name, bond.acceptor_atom_name);

        if (legacy_type == '-') {
            bond.classification = HBondClassification::STANDARD;
        } else if (legacy_type == '*') {
            bond.classification = HBondClassification::NON_STANDARD;
        } else {
            bond.classification = HBondClassification::INVALID;
        }
    }
}

void HBondDetector::calculate_angles(std::vector<HBond>& bonds, const Residue& residue1,
                                     const Residue& residue2) const {
    for (auto& bond : bonds) {
        // Find donor and acceptor positions
        auto donor_opt = residue1.find_atom(bond.donor_atom_name);
        auto acceptor_opt = residue2.find_atom(bond.acceptor_atom_name);

        if (!donor_opt || !acceptor_opt) {
            continue;
        }

        const Vector3D donor_pos = donor_opt->position();
        const Vector3D acceptor_pos = acceptor_opt->position();

        // Find donor neighbor
        auto donor_neighbor_pos = HBondGeometry::find_neighbor_position(bond.donor_atom_name, residue1);
        if (donor_neighbor_pos) {
            bond.donor_angle = HBondGeometry::calculate_angle(*donor_neighbor_pos, donor_pos, acceptor_pos);
            bond.donor_neighbor_atom = HBondGeometry::get_neighbor_atom_name(bond.donor_atom_name);
        }

        // Find acceptor neighbor
        auto acceptor_neighbor_pos = HBondGeometry::find_neighbor_position(bond.acceptor_atom_name, residue2);
        if (acceptor_neighbor_pos) {
            bond.acceptor_angle = HBondGeometry::calculate_angle(donor_pos, acceptor_pos, *acceptor_neighbor_pos);
            bond.acceptor_neighbor_atom = HBondGeometry::get_neighbor_atom_name(bond.acceptor_atom_name);
        }

        // Calculate dihedral if both neighbors found
        if (donor_neighbor_pos && acceptor_neighbor_pos) {
            bond.dihedral_angle = HBondGeometry::calculate_dihedral(*donor_neighbor_pos, donor_pos, acceptor_pos,
                                                                    *acceptor_neighbor_pos);
            bond.dihedral_valid = true;
        }
    }
}

void HBondDetector::apply_post_validation_filtering(std::vector<HBond>& bonds) const {
    // Apply post-validation filtering (marks bonds as INVALID but doesn't remove them)
    const int num_good = HBondRoleClassifier::count_good_hbonds(bonds);

    if (num_good > 0) {
        for (auto& bond : bonds) {
            // Skip already invalid bonds
            if (bond.classification == HBondClassification::INVALID) {
                continue;
            }

            // Filter out H-bonds exceeding max distance
            if (bond.distance > params_.post_validation_max_distance) {
                bond.classification = HBondClassification::INVALID;
                continue;
            }

            // Filter out non-standard H-bonds outside valid range
            const bool is_nonstandard_invalid = bond.classification == HBondClassification::NON_STANDARD &&
                                                bond.conflict_state != ConflictState::IS_CONFLICT_WINNER &&
                                                (bond.distance < params_.nonstandard_min_distance ||
                                                 bond.distance > params_.nonstandard_max_distance);

            if (is_nonstandard_invalid) {
                bond.classification = HBondClassification::INVALID;
            }
        }
    }
}

char HBondDetector::get_base_type_for_hbond(const Residue& residue) {
    // Check one_letter_code from classification first
    const char code = residue.classification().one_letter_code;
    if (code != '?') {
        return code;
    }

    // Try base_type() for known types
    switch (residue.base_type()) {
        case BaseType::ADENINE:
            return 'A';
        case BaseType::CYTOSINE:
            return 'C';
        case BaseType::GUANINE:
            return 'G';
        case BaseType::THYMINE:
            return 'T';
        case BaseType::URACIL:
            return 'U';
        default:
            break;
    }

    // Fall back to atom-based detection
    return determine_base_type_from_atoms(residue);
}

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna

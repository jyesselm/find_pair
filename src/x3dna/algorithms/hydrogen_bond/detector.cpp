/**
 * @file detector.cpp
 * @brief Implementation of general-purpose H-bond detector
 */

#include <x3dna/algorithms/hydrogen_bond/detector.hpp>
#include <x3dna/algorithms/hydrogen_bond/geometry.hpp>
#include <x3dna/algorithms/hydrogen_bond/role_classifier.hpp>
#include <x3dna/algorithms/hydrogen_bond/hydrogen_bond_utils.hpp>
#include <x3dna/algorithms/hydrogen_bond/quality_scorer.hpp>
#include <x3dna/algorithms/hydrogen_bond/edge_classifier.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/core/typing/atom_classification.hpp>
#include <x3dna/core/typing/nucleotide_type.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

using namespace x3dna::core;
using namespace x3dna::core::typing;
using namespace x3dna::geometry;

namespace {

/**
 * @brief Check if two backbone atoms are covalently connected through the phosphodiester linkage
 *
 * For adjacent nucleotides in the backbone:
 * - Residue i: O3' connects to P of residue i+1
 * - The phosphate group has O1P, O2P, O5' connected to P
 *
 * So for sequential residues: O3'(i) is close to P(i+1), O1P(i+1), O2P(i+1), O5'(i+1)
 * And P(i), O1P(i), O2P(i), O5'(i) are close to O3'(i-1)
 *
 * @param atom1 First atom name
 * @param atom2 Second atom name
 * @return True if atoms are part of the same phosphodiester linkage
 */
[[nodiscard]] bool is_phosphodiester_pair(const std::string& atom1, const std::string& atom2) {
    // Atoms that connect the 3' end to the phosphate group
    static const std::set<std::string> three_prime = {"O3'"};
    // Atoms in the phosphate group that connect to the 5' side
    static const std::set<std::string> phosphate = {"P", "O1P", "O2P", "OP1", "OP2", "O5'"};

    // Check if one atom is O3' and the other is in the phosphate group
    const bool a1_3p = three_prime.count(atom1) > 0;
    const bool a2_3p = three_prime.count(atom2) > 0;
    const bool a1_phos = phosphate.count(atom1) > 0;
    const bool a2_phos = phosphate.count(atom2) > 0;

    // O3' to phosphate group atoms are covalently connected
    if ((a1_3p && a2_phos) || (a2_3p && a1_phos)) {
        return true;
    }

    // Within phosphate group - O1P/O2P to O5' are within 2 bonds
    // These can appear as short "H-bonds" but are just geometry
    static const std::set<std::string> phosphate_o = {"O1P", "O2P", "OP1", "OP2"};
    static const std::set<std::string> five_prime = {"O5'"};
    const bool a1_po = phosphate_o.count(atom1) > 0;
    const bool a2_po = phosphate_o.count(atom2) > 0;
    const bool a1_5p = five_prime.count(atom1) > 0;
    const bool a2_5p = five_prime.count(atom2) > 0;

    // O1P/O2P to O5' are connected through P
    if ((a1_po && a2_5p) || (a2_po && a1_5p)) {
        return true;
    }

    return false;
}

/**
 * @brief Convert HBondContext to HBondInteractionType for filtering
 */
[[nodiscard]] HBondInteractionType context_to_interaction_type(HBondContext ctx) {
    switch (ctx) {
        case HBondContext::BASE_BASE:
            return HBondInteractionType::BASE_BASE;
        case HBondContext::BASE_BACKBONE:
        case HBondContext::BACKBONE_BACKBONE:
            return HBondInteractionType::BASE_BACKBONE;
        case HBondContext::BASE_SUGAR:
        case HBondContext::SUGAR_SUGAR:
            return HBondInteractionType::BASE_SUGAR;
        case HBondContext::BASE_PROTEIN:
        case HBondContext::SUGAR_PROTEIN:
        case HBondContext::BACKBONE_PROTEIN:
            return HBondInteractionType::BASE_PROTEIN;
        case HBondContext::BASE_LIGAND:
            return HBondInteractionType::BASE_LIGAND;
        case HBondContext::PROTEIN_MAINCHAIN:
        case HBondContext::PROTEIN_SIDECHAIN:
            return HBondInteractionType::PROTEIN_PROTEIN;
        case HBondContext::PROTEIN_LIGAND:
        case HBondContext::LIGAND_LIGAND:
            return HBondInteractionType::PROTEIN_LIGAND;
        default:
            return HBondInteractionType::ANY;
    }
}

/**
 * @brief Check if a context passes the interaction type filter
 */
[[nodiscard]] bool passes_interaction_filter(HBondContext ctx, HBondInteractionType filter) {
    if (filter == HBondInteractionType::ANY) {
        return true;
    }
    // RNA_INTERNAL includes base-base, base-sugar, sugar-sugar, backbone interactions
    if (filter == HBondInteractionType::RNA_INTERNAL) {
        return ctx == HBondContext::BASE_BASE || ctx == HBondContext::BASE_SUGAR ||
               ctx == HBondContext::SUGAR_SUGAR || ctx == HBondContext::BASE_BACKBONE ||
               ctx == HBondContext::BACKBONE_BACKBONE;
    }
    return filter & context_to_interaction_type(ctx);
}

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
 * @brief Determine purine type (A or G) from atoms using AtomType
 */
char determine_purine_type(const Residue& residue) {
    const bool has_o6 = residue.has_atom_type(AtomType::O6);
    const bool has_n6 = residue.has_atom_type(AtomType::N6);
    return (has_o6 || !has_n6) ? 'G' : 'A';
}

/**
 * @brief Determine pyrimidine type (C, T, or U) from atoms using AtomType
 */
char determine_pyrimidine_type(const Residue& residue) {
    if (residue.has_atom_type(AtomType::N4)) {
        return 'C';
    }
    if (residue.has_atom_type(AtomType::C5M) || residue.has_atom_type(AtomType::C7)) {
        return 'T';
    }
    return 'U';
}

/**
 * @brief Determine base type from atoms for unknown residues using AtomType
 */
char determine_base_type_from_atoms(const Residue& residue) {
    const bool has_n9 = residue.has_atom_type(AtomType::N9);
    const bool has_n1 = residue.has_atom_type(AtomType::N1);
    const bool has_c6 = residue.has_atom_type(AtomType::C6);

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

    // Step 4b: Classify Leontis-Westhof edges for each H-bond
    for (auto& bond : bonds) {
        bond.donor_edge = EdgeClassifier::classify(bond.donor_atom_name, base1);
        bond.acceptor_edge = EdgeClassifier::classify(bond.acceptor_atom_name, base2);
    }

    // Step 5: Calculate angles for all bonds (in place)
    calculate_angles(bonds, residue1, residue2);

    // Step 6: Apply post-validation filtering (marks bonds as INVALID but doesn't remove)
    apply_post_validation_filtering(bonds);

    // Step 7: Apply optional angle-based filtering
    apply_angle_filtering(bonds);

    // Step 8: Apply optional quality scoring
    apply_quality_scoring(bonds);

    // Step 9: Build final_bonds by copying only valid bonds (single allocation)
    result.final_bonds.reserve(bonds.size());
    for (const auto& bond : bonds) {
        // Skip INVALID bonds
        if (bond.classification == HBondClassification::INVALID) {
            continue;
        }
        // Skip UNLIKELY_CHEMISTRY bonds unless exploration mode is enabled
        if (bond.classification == HBondClassification::UNLIKELY_CHEMISTRY &&
            !params_.include_unlikely_chemistry) {
            continue;
        }

        result.final_bonds.push_back(bond);
        if (bond.classification == HBondClassification::STANDARD) {
            result.standard_bond_count++;
            if (HBondRoleClassifier::is_good_hbond_distance(bond.distance)) {
                result.good_bond_count++;
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
            // Use AtomType for O2' check (O(1) comparison)
            const bool not_o2prime = !a1.is_o2_prime() && !a2.is_o2_prime();

            if (both_base && not_o2prime) {
                if (good_hb_atoms(a1.name(), a2.name(), params_.allowed_elements)) {
                    num_base_hb++;
                }
            }

            // Use AtomType for O2' check
            if (a1.is_o2_prime() || a2.is_o2_prime()) {
                num_o2_hb++;
            }
        }
    }
}

std::vector<HBond> HBondDetector::detect_intra_residue_hbonds(
    const Residue& residue, MoleculeType mol_type) const {

    std::vector<HBond> bonds;
    const auto& atoms = residue.atoms();
    const size_t n = atoms.size();

    // Need at least 2 atoms for an intra-residue H-bond
    if (n < 2) {
        return bonds;
    }

    // Check all atom pairs within the residue
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const auto& atom1 = atoms[i];
            const auto& atom2 = atoms[j];

            // Check if atoms can form H-bond based on elements
            if (!good_hb_atoms(atom1.name(), atom2.name(), params_.allowed_elements,
                               params_.include_backbone_backbone)) {
                continue;
            }

            // Calculate distance
            const double dist = (atom1.position() - atom2.position()).length();

            // Use a generic context for intra-residue (treat as base-sugar for nucleotides)
            const HBondContext context = HBondContext::BASE_SUGAR;

            // Check distance is in valid range
            const double max_dist = params_.distances.max_for_context(context);
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
            hbond.donor_res_id = residue.res_id();
            hbond.acceptor_res_id = residue.res_id();  // Same residue

            // Classify the bond
            char base_type = '?';
            if (mol_type == MoleculeType::NUCLEIC_ACID) {
                base_type = get_base_type_for_hbond(residue);
            }

            // Get atom roles and classify
            HBondAtomRole role1 = HBondRoleClassifier::get_nucleotide_atom_role(base_type, atom1.name());
            HBondAtomRole role2 = HBondRoleClassifier::get_nucleotide_atom_role(base_type, atom2.name());

            // Check for AA/DD (unlikely chemistry)
            const bool is_aa = (role1 == HBondAtomRole::ACCEPTOR && role2 == HBondAtomRole::ACCEPTOR);
            const bool is_dd = (role1 == HBondAtomRole::DONOR && role2 == HBondAtomRole::DONOR);

            if (is_aa || is_dd) {
                if (!params_.include_unlikely_chemistry) {
                    continue;  // Skip unlikely chemistry
                }
                hbond.classification = HBondClassification::UNLIKELY_CHEMISTRY;
            } else {
                hbond.classification = HBondClassification::STANDARD;
            }

            // Set Leontis-Westhof edges
            hbond.donor_edge = EdgeClassifier::classify(atom1.name(), base_type);
            hbond.acceptor_edge = EdgeClassifier::classify(atom2.name(), base_type);

            bonds.push_back(hbond);
        }
    }

    return bonds;
}

std::vector<HBond> HBondDetector::find_candidate_bonds(const Residue& residue1, const Residue& residue2,
                                                        bool base_atoms_only, MoleculeType mol1_type,
                                                        MoleculeType mol2_type) const {
    std::vector<HBond> candidates;

    for (const auto& atom1 : residue1.atoms()) {
        for (const auto& atom2 : residue2.atoms()) {
            // Check if atoms can form H-bond based on elements
            if (!good_hb_atoms(atom1.name(), atom2.name(), params_.allowed_elements,
                               params_.include_backbone_backbone)) {
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

            // Check interaction type filter
            if (!passes_interaction_filter(context, params_.interaction_filter)) {
                continue;
            }

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
            // Check if this is an AA or DD pair (chemically unlikely)
            // Get atom roles to distinguish AA/DD from ambiguous X-containing combinations
            HBondAtomRole role1 = HBondRoleClassifier::get_nucleotide_atom_role(base1_type, bond.donor_atom_name);
            HBondAtomRole role2 = HBondRoleClassifier::get_nucleotide_atom_role(base2_type, bond.acceptor_atom_name);

            // AA or DD pairs are chemically unlikely without protonation/tautomerization
            const bool is_aa = (role1 == HBondAtomRole::ACCEPTOR && role2 == HBondAtomRole::ACCEPTOR);
            const bool is_dd = (role1 == HBondAtomRole::DONOR && role2 == HBondAtomRole::DONOR);

            if (is_aa || is_dd) {
                // For backbone-backbone context, AA bonds (like O1P-O2P) are common and
                // can form via water-mediated H-bonds. Classify as NON_STANDARD instead
                // of UNLIKELY_CHEMISTRY to include them in DSSR-like detection.
                if (bond.context == HBondContext::BACKBONE_BACKBONE) {
                    bond.classification = HBondClassification::NON_STANDARD;
                } else {
                    bond.classification = HBondClassification::UNLIKELY_CHEMISTRY;
                }
            } else {
                bond.classification = HBondClassification::NON_STANDARD;
            }
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

void HBondDetector::apply_angle_filtering(std::vector<HBond>& bonds) const {
    // Only apply if enabled (off by default for legacy compatibility)
    if (!params_.enable_angle_filtering) {
        return;
    }

    for (auto& bond : bonds) {
        // Skip already invalid bonds
        if (bond.classification == HBondClassification::INVALID) {
            continue;
        }

        // Check donor angle if calculated (> 0 means angle was computed)
        if (bond.donor_angle > 0.0 && bond.donor_angle < params_.min_donor_angle) {
            bond.classification = HBondClassification::INVALID;
            continue;
        }

        // Check acceptor angle if calculated
        if (bond.acceptor_angle > 0.0 && bond.acceptor_angle < params_.min_acceptor_angle) {
            bond.classification = HBondClassification::INVALID;
        }
    }
}

void HBondDetector::apply_quality_scoring(std::vector<HBond>& bonds) const {
    // Only apply if enabled
    if (!params_.enable_quality_scoring) {
        return;
    }

    // Create scorer with default parameters
    HBondQualityScorer scorer;

    for (auto& bond : bonds) {
        // Skip already invalid bonds
        if (bond.classification == HBondClassification::INVALID) {
            continue;
        }

        // Score this bond
        bond.quality_score = scorer.score(bond);

        // Optionally filter out INVALID tier bonds
        if (params_.filter_invalid_scores &&
            bond.quality_score->tier == core::HBondQualityTier::INVALID) {
            bond.classification = HBondClassification::INVALID;
        }
    }
}

StructureHBondResult HBondDetector::detect_all_structure_hbonds(
    const Structure& structure,
    double max_residue_distance) const {

    StructureHBondResult result;

    // Get all residues
    auto residues = structure.all_residues();
    const size_t n_residues = residues.size();

    if (n_residues < 2) {
        return result;
    }

    // Pre-compute residue centers for distance filtering
    std::vector<Vector3D> centers(n_residues);
    for (size_t i = 0; i < n_residues; ++i) {
        // Use center of mass of heavy atoms
        Vector3D sum(0, 0, 0);
        size_t count = 0;
        for (const auto& atom : residues[i]->atoms()) {
            // Skip hydrogens
            if (atom.element() == "H" || atom.name()[0] == 'H') {
                continue;
            }
            sum = sum + atom.position();
            ++count;
        }
        if (count > 0) {
            centers[i] = sum * (1.0 / count);
        }
    }

    const double max_dist_sq = max_residue_distance * max_residue_distance;

    // Detect intra-residue H-bonds if enabled
    if (params_.include_intra_residue) {
        for (size_t i = 0; i < n_residues; ++i) {
            const MoleculeType mol_type = residues[i]->is_nucleotide() ? MoleculeType::NUCLEIC_ACID
                                        : residues[i]->is_protein() ? MoleculeType::PROTEIN
                                        : MoleculeType::LIGAND;

            auto intra_hbonds = detect_intra_residue_hbonds(*residues[i], mol_type);

            if (!intra_hbonds.empty()) {
                // Add to grouped result (same residue for both i and j)
                ResidueHBonds intra_result;
                intra_result.res_id_i = residues[i]->res_id();
                intra_result.res_id_j = residues[i]->res_id();  // Same residue
                intra_result.residue_idx_i = i;
                intra_result.residue_idx_j = i;
                intra_result.hbonds = std::move(intra_hbonds);

                // Add to flat list
                for (const auto& hb : intra_result.hbonds) {
                    result.all_hbonds.push_back(hb);
                }

                result.residue_pair_hbonds.push_back(std::move(intra_result));
                ++result.pairs_with_hbonds;
            }
        }
    }

    // Check all residue pairs
    for (size_t i = 0; i < n_residues; ++i) {
        for (size_t j = i + 1; j < n_residues; ++j) {
            // Early rejection based on residue center distance
            const double center_dist_sq = (centers[i] - centers[j]).length_squared();
            if (center_dist_sq > max_dist_sq) {
                continue;
            }

            ++result.total_residue_pairs_checked;

            // Determine molecule types
            const MoleculeType mol1_type = residues[i]->is_nucleotide() ? MoleculeType::NUCLEIC_ACID
                                         : residues[i]->is_protein() ? MoleculeType::PROTEIN
                                         : MoleculeType::LIGAND;
            const MoleculeType mol2_type = residues[j]->is_nucleotide() ? MoleculeType::NUCLEIC_ACID
                                         : residues[j]->is_protein() ? MoleculeType::PROTEIN
                                         : MoleculeType::LIGAND;

            // Detect H-bonds between this pair (all atoms, not just base atoms)
            auto hbonds = detect_all_hbonds_between(*residues[i], *residues[j], mol1_type, mol2_type);

            if (!hbonds.empty()) {
                // Check if residues are sequence-adjacent nucleotides
                // If so, filter out backbone-backbone bonds that are part of the phosphodiester linkage
                const bool both_nucleotides = mol1_type == MoleculeType::NUCLEIC_ACID &&
                                              mol2_type == MoleculeType::NUCLEIC_ACID;
                bool is_sequence_adjacent = false;

                if (both_nucleotides) {
                    // Check if same chain and consecutive sequence numbers
                    const auto& res1 = *residues[i];
                    const auto& res2 = *residues[j];
                    if (res1.chain_id() == res2.chain_id()) {
                        const int seq_diff = std::abs(res1.seq_num() - res2.seq_num());
                        is_sequence_adjacent = (seq_diff == 1);
                    }
                }

                // Filter out phosphodiester-linked backbone atoms for adjacent residues
                if (is_sequence_adjacent) {
                    hbonds.erase(
                        std::remove_if(hbonds.begin(), hbonds.end(),
                            [](const HBond& hb) {
                                return hb.context == HBondContext::BACKBONE_BACKBONE &&
                                       is_phosphodiester_pair(hb.donor_atom_name, hb.acceptor_atom_name);
                            }),
                        hbonds.end());
                }

                if (hbonds.empty()) {
                    continue;  // All bonds were filtered out
                }

                ++result.pairs_with_hbonds;

                // Set residue info on each H-bond
                for (auto& hb : hbonds) {
                    hb.donor_res_id = residues[i]->res_id();
                    hb.acceptor_res_id = residues[j]->res_id();
                    hb.donor_residue_idx = i;
                    hb.acceptor_residue_idx = j;
                }

                // Add to grouped result
                ResidueHBonds pair_result;
                pair_result.res_id_i = residues[i]->res_id();
                pair_result.res_id_j = residues[j]->res_id();
                pair_result.residue_idx_i = i;
                pair_result.residue_idx_j = j;
                pair_result.hbonds = std::move(hbonds);

                // Also add to flat list
                for (const auto& hb : pair_result.hbonds) {
                    result.all_hbonds.push_back(hb);
                }

                result.residue_pair_hbonds.push_back(std::move(pair_result));
            }
        }
    }

    return result;
}

void HBondDetector::apply_global_occupancy_filter(
    StructureHBondResult& result,
    int max_bonds_per_donor,
    int max_bonds_per_acceptor) const {

    if (result.all_hbonds.empty()) {
        return;
    }

    // Build unique atom identifiers: "res_id:atom_name"
    auto make_atom_id = [](const std::string& res_id, const std::string& atom_name) {
        return res_id + ":" + atom_name;
    };

    // Sort all H-bonds by distance only (shortest first)
    // Note: Our "donor" and "acceptor" labels are based on residue order, not chemistry
    // So we track each atom uniformly using max_bonds_per_acceptor as the limit
    std::vector<size_t> sorted_indices(result.all_hbonds.size());
    std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
    std::sort(sorted_indices.begin(), sorted_indices.end(),
        [&result](size_t a, size_t b) {
            return result.all_hbonds[a].distance < result.all_hbonds[b].distance;
        });

    // Track occupancy for each atom (uniform tracking, not role-based)
    // Since our donor/acceptor labels don't reflect actual chemistry,
    // we use the same limit for all atoms
    std::unordered_map<std::string, int> atom_occupancy;

    // Mark which bonds to keep
    std::vector<bool> keep_bond(result.all_hbonds.size(), false);

    // Use the more permissive limit for all atoms
    const int max_bonds_per_atom = std::max(max_bonds_per_donor, max_bonds_per_acceptor);

    // Greedy selection: process bonds in distance order
    for (size_t idx : sorted_indices) {
        const auto& hb = result.all_hbonds[idx];

        // Create unique atom identifiers for both atoms in the bond
        const std::string atom1_id = make_atom_id(hb.donor_res_id, hb.donor_atom_name);
        const std::string atom2_id = make_atom_id(hb.acceptor_res_id, hb.acceptor_atom_name);

        // Check if both atoms have capacity
        const int atom1_count = atom_occupancy[atom1_id];
        const int atom2_count = atom_occupancy[atom2_id];

        if (atom1_count < max_bonds_per_atom && atom2_count < max_bonds_per_atom) {
            // Keep this bond
            keep_bond[idx] = true;
            atom_occupancy[atom1_id]++;
            atom_occupancy[atom2_id]++;
        }
    }

    // Filter all_hbonds
    std::vector<core::HBond> filtered_hbonds;
    filtered_hbonds.reserve(result.all_hbonds.size());
    for (size_t i = 0; i < result.all_hbonds.size(); ++i) {
        if (keep_bond[i]) {
            filtered_hbonds.push_back(std::move(result.all_hbonds[i]));
        }
    }
    result.all_hbonds = std::move(filtered_hbonds);

    // Rebuild residue_pair_hbonds from filtered bonds
    // Group by (res_id_i, res_id_j) pair
    std::map<std::pair<std::string, std::string>, std::vector<core::HBond>> grouped;
    for (auto& hb : result.all_hbonds) {
        // Determine which residue is "i" and which is "j" based on original order
        std::string res_i = hb.donor_res_id;
        std::string res_j = hb.acceptor_res_id;
        // For intra-residue, both are the same
        if (res_i == res_j) {
            grouped[{res_i, res_j}].push_back(hb);
        } else {
            // Ensure consistent ordering (smaller res_id first)
            if (res_i > res_j) {
                std::swap(res_i, res_j);
            }
            grouped[{res_i, res_j}].push_back(hb);
        }
    }

    // Rebuild residue_pair_hbonds
    result.residue_pair_hbonds.clear();
    result.pairs_with_hbonds = 0;
    for (auto& [key, hbonds] : grouped) {
        ResidueHBonds pair_result;
        pair_result.res_id_i = key.first;
        pair_result.res_id_j = key.second;
        pair_result.hbonds = std::move(hbonds);
        result.residue_pair_hbonds.push_back(std::move(pair_result));
        ++result.pairs_with_hbonds;
    }
}

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna

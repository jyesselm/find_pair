/**
 * @file dssr_filter.cpp
 * @brief Implementation of DSSR-style H-bond filtering
 */

#include <x3dna/algorithms/hydrogen_bond/dssr_filter.hpp>
#include <x3dna/algorithms/hydrogen_bond/detector.hpp>
#include <x3dna/algorithms/hydrogen_bond/quality_scorer.hpp>
#include <algorithm>
#include <cctype>
#include <map>
#include <unordered_map>
#include <vector>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

namespace {

/**
 * @brief Check if atom is an amino group (NH2) that acts as H-bond donor only
 */
bool is_amino_group(const std::string& atom_name) {
    std::string name = atom_name;
    // Trim whitespace
    while (!name.empty() && std::isspace(static_cast<unsigned char>(name.front()))) {
        name.erase(0, 1);
    }
    while (!name.empty() && std::isspace(static_cast<unsigned char>(name.back()))) {
        name.pop_back();
    }

    // Amino groups on bases: N6 (adenine), N4 (cytosine), N2 (guanine)
    return name == "N6" || name == "N4" || name == "N2";
}

/**
 * @brief Get atom-specific H-bond capacity based on chemical properties
 *
 * Different atoms can participate in different numbers of H-bonds:
 * - Amino groups (N6, N4, N2): 2 bonds (have 2 hydrogens to donate)
 * - Carbonyl oxygens (O6, O4, O2): 2 bonds (have 2 lone pairs)
 * - Ring nitrogens (N1, N3, N7, N9): typically 1-2 bonds
 * - O2' (ribose hydroxyl): 3 bonds (1 H + 2 lone pairs)
 * - Phosphate oxygens: 2 bonds each
 */
int get_atom_capacity(const std::string& atom_name) {
    std::string name = atom_name;
    // Trim whitespace
    while (!name.empty() && std::isspace(static_cast<unsigned char>(name.front()))) {
        name.erase(0, 1);
    }
    while (!name.empty() && std::isspace(static_cast<unsigned char>(name.back()))) {
        name.pop_back();
    }

    // Amino groups: 2 bonds (2 H's)
    if (name == "N6" || name == "N4" || name == "N2") {
        return 2;
    }

    // Carbonyl oxygens: 2 bonds (2 lone pairs)
    if (name == "O6" || name == "O4" || name == "O2") {
        return 2;
    }

    // O2' (ribose 2' hydroxyl): 3 bonds (1 H + 2 lone pairs)
    if (name == "O2'" || name == "O2*") {
        return 3;
    }

    // Phosphate oxygens: 2 bonds each
    if (name == "OP1" || name == "OP2" || name == "O1P" || name == "O2P") {
        return 2;
    }

    // O5', O3': 2 bonds (2 lone pairs)
    if (name == "O5'" || name == "O3'" || name == "O5*" || name == "O3*") {
        return 2;
    }

    // Ring nitrogens: typically 1-2 bonds
    if (name == "N1" || name == "N3" || name == "N7" || name == "N9") {
        return 2;  // Can participate in multiple bonds (donor + acceptor roles)
    }

    // Default: 2 bonds
    return 2;
}

} // anonymous namespace

char DSSRStyleFilter::get_element(const std::string& atom_name) {
    // Trim whitespace
    std::string name = atom_name;
    while (!name.empty() && std::isspace(static_cast<unsigned char>(name.front()))) {
        name.erase(0, 1);
    }
    while (!name.empty() && std::isspace(static_cast<unsigned char>(name.back()))) {
        name.pop_back();
    }

    if (name.empty()) {
        return '?';
    }

    // Get first non-digit character (handles cases like "1H5'" -> H)
    for (char c : name) {
        if (std::isalpha(static_cast<unsigned char>(c))) {
            return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }
    }

    return '?';
}

bool DSSRStyleFilter::is_o2_prime(const std::string& atom_name) {
    // Check for O2' or O2* (alternate naming)
    std::string name = atom_name;
    // Trim whitespace
    while (!name.empty() && std::isspace(static_cast<unsigned char>(name.front()))) {
        name.erase(0, 1);
    }
    while (!name.empty() && std::isspace(static_cast<unsigned char>(name.back()))) {
        name.pop_back();
    }

    return name == "O2'" || name == "O2*";
}

double DSSRStyleFilter::get_distance_threshold(
    const std::string& atom1_name,
    const std::string& atom2_name,
    const DSSRFilterParams& params) {

    const char elem1 = get_element(atom1_name);
    const char elem2 = get_element(atom2_name);

    // If either atom contains N, use the N-containing threshold
    if (elem1 == 'N' || elem2 == 'N') {
        return params.n_containing_max_distance;
    }

    // Both are O (or other non-N elements)
    if (elem1 == 'O' && elem2 == 'O') {
        // Check if either is O2' (ribose hydroxyl)
        if (is_o2_prime(atom1_name) || is_o2_prime(atom2_name)) {
            return params.o2prime_oo_max_distance;
        }
        // Other O-O bonds get stricter threshold
        return params.other_oo_max_distance;
    }

    // Default to N-containing threshold for other cases (C-O, etc.)
    return params.n_containing_max_distance;
}

bool DSSRStyleFilter::is_chemically_unlikely_pair(
    const std::string& atom1_name,
    const std::string& atom2_name) {

    // Amino-amino pairs are chemically unlikely (both are donors)
    // N6-N4, N6-N2, N4-N2, N6-N6, N4-N4, N2-N2
    // These atoms have hydrogens to donate but no lone pairs to accept
    if (is_amino_group(atom1_name) && is_amino_group(atom2_name)) {
        return true;
    }

    // Note: We do NOT filter carbonyl-carbonyl pairs because DSSR reports them
    // O2-O2 (C-U), O6-O4, etc. are structurally important close contacts

    return false;
}

bool DSSRStyleFilter::should_keep(const core::HBond& hb, const DSSRFilterParams& params) {
    // Check minimum distance
    if (hb.distance < params.min_distance) {
        return false;
    }

    // Filter out chemically unlikely pairs by default
    if (params.filter_unlikely_pairs &&
        is_chemically_unlikely_pair(hb.donor_atom_name, hb.acceptor_atom_name)) {
        return false;
    }

    // Get appropriate threshold for this atom pair
    const double max_dist = get_distance_threshold(
        hb.donor_atom_name, hb.acceptor_atom_name, params);

    // Check if distance is within threshold
    return hb.distance <= max_dist;
}

std::vector<core::HBond> DSSRStyleFilter::filter(
    const std::vector<core::HBond>& hbonds,
    const DSSRFilterParams& params) {

    std::vector<core::HBond> result;
    result.reserve(hbonds.size());

    for (const auto& hb : hbonds) {
        if (should_keep(hb, params)) {
            result.push_back(hb);
        }
    }

    return result;
}

void DSSRStyleFilter::filter_in_place(StructureHBondResult& result, const DSSRFilterParams& params) {
    // Filter all_hbonds
    auto& all = result.all_hbonds;
    all.erase(
        std::remove_if(all.begin(), all.end(),
            [&params](const core::HBond& hb) {
                return !should_keep(hb, params);
            }),
        all.end());

    // Rebuild residue_pair_hbonds from filtered bonds
    // Group by (donor_res_id, acceptor_res_id) pair - preserve original order!
    // IMPORTANT: Do NOT swap res_ids, as donor_atom is defined to come from res_id_i
    std::map<std::pair<std::string, std::string>, std::vector<core::HBond>> grouped;
    for (const auto& hb : all) {
        // Keep original ordering: donor_res_id = res_id_i, acceptor_res_id = res_id_j
        grouped[{hb.donor_res_id, hb.acceptor_res_id}].push_back(hb);
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

void DSSRStyleFilter::apply_scored_occupancy_filter(StructureHBondResult& result, int max_bonds_per_atom) {
    if (result.all_hbonds.empty()) {
        return;
    }

    // Score all bonds first
    HBondQualityScorer scorer;
    for (auto& hb : result.all_hbonds) {
        if (!hb.quality_score.has_value()) {
            hb.quality_score = scorer.score(hb);
        }
    }

    // Helper to create atom ID
    auto make_atom_id = [](const std::string& res_id, const std::string& atom_name) {
        return res_id + ":" + atom_name;
    };

    // Helper to get score (use 0 if not available)
    auto get_score = [](const core::HBond& hb) -> double {
        if (hb.quality_score.has_value()) {
            return hb.quality_score->total_score;
        }
        // Fallback: use inverse distance as score (shorter = better)
        return 100.0 - (hb.distance * 20.0);
    };

    // Helper to get atom-specific capacity
    auto get_capacity = [max_bonds_per_atom](const std::string& atom_name) -> int {
        // If user specified a limit, use that as upper bound
        int chemical_capacity = get_atom_capacity(atom_name);
        if (max_bonds_per_atom > 0) {
            return std::min(chemical_capacity, max_bonds_per_atom);
        }
        return chemical_capacity;
    };

    // Create index + score pairs and sort by score descending
    std::vector<std::pair<size_t, double>> indexed_scores;
    indexed_scores.reserve(result.all_hbonds.size());
    for (size_t i = 0; i < result.all_hbonds.size(); ++i) {
        indexed_scores.emplace_back(i, get_score(result.all_hbonds[i]));
    }
    std::sort(indexed_scores.begin(), indexed_scores.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    // Track how many bonds each atom participates in
    std::unordered_map<std::string, int> atom_bond_count;
    // Track atom capacities (cache for efficiency)
    std::unordered_map<std::string, int> atom_capacities;

    // Mark which bonds to keep
    std::vector<bool> keep_bond(result.all_hbonds.size(), false);

    // Greedy selection: process bonds in score order (best first)
    for (const auto& [idx, score] : indexed_scores) {
        const auto& hb = result.all_hbonds[idx];

        // Get atom IDs for both atoms
        const std::string atom1_id = make_atom_id(hb.donor_res_id, hb.donor_atom_name);
        const std::string atom2_id = make_atom_id(hb.acceptor_res_id, hb.acceptor_atom_name);

        // Get or compute capacities
        if (atom_capacities.find(atom1_id) == atom_capacities.end()) {
            atom_capacities[atom1_id] = get_capacity(hb.donor_atom_name);
        }
        if (atom_capacities.find(atom2_id) == atom_capacities.end()) {
            atom_capacities[atom2_id] = get_capacity(hb.acceptor_atom_name);
        }

        // Check if both atoms have capacity (using atom-specific limits)
        if (atom_bond_count[atom1_id] < atom_capacities[atom1_id] &&
            atom_bond_count[atom2_id] < atom_capacities[atom2_id]) {
            keep_bond[idx] = true;
            atom_bond_count[atom1_id]++;
            atom_bond_count[atom2_id]++;
        }
    }

    // Filter all_hbonds
    std::vector<core::HBond> filtered;
    filtered.reserve(result.all_hbonds.size());
    for (size_t i = 0; i < result.all_hbonds.size(); ++i) {
        if (keep_bond[i]) {
            filtered.push_back(std::move(result.all_hbonds[i]));
        }
    }
    result.all_hbonds = std::move(filtered);

    // Rebuild residue_pair_hbonds
    std::map<std::pair<std::string, std::string>, std::vector<core::HBond>> grouped;
    for (const auto& hb : result.all_hbonds) {
        grouped[{hb.donor_res_id, hb.acceptor_res_id}].push_back(hb);
    }

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

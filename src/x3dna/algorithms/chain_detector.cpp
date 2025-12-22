/**
 * @file chain_detector.cpp
 * @brief Implementation of chain detection based on backbone connectivity
 */

#include <x3dna/algorithms/chain_detector.hpp>
#include <x3dna/core/typing/atom_type.hpp>
#include <algorithm>
#include <limits>
#include <functional>

namespace x3dna {
namespace algorithms {

using core::AtomType;

ChainDetector::ChainDetector(const Config& config) : config_(config) {}

// ============================================================================
// Main detection methods
// ============================================================================

std::vector<ConnectedChain> ChainDetector::detect_rna_chains(const core::Structure& structure) const {
    auto rna_residues = filter_rna_residues(structure);
    if (rna_residues.empty()) {
        return {};
    }

    sort_by_chain_and_num(rna_residues);

    // Build chains using RNA connectivity
    auto connectivity_func = [this](const core::Residue& r1, const core::Residue& r2) {
        return are_rna_residues_connected(r1, r2);
    };

    auto chains = build_chains(rna_residues, connectivity_func, true);

    // Merge adjacent chains if enabled
    if (config_.merge_adjacent_chains) {
        chains = merge_adjacent_chains(chains);
    }

    return chains;
}

std::vector<ConnectedChain> ChainDetector::detect_protein_chains(const core::Structure& structure) const {
    auto protein_residues = filter_protein_residues(structure);
    if (protein_residues.empty()) {
        return {};
    }

    sort_by_chain_and_num(protein_residues);

    // Build chains using protein connectivity
    auto connectivity_func = [this](const core::Residue& r1, const core::Residue& r2) {
        return are_protein_residues_connected(r1, r2);
    };

    auto chains = build_chains(protein_residues, connectivity_func, false);

    return chains;
}

std::vector<ConnectedChain> ChainDetector::detect_all_chains(const core::Structure& structure) const {
    auto rna_chains = detect_rna_chains(structure);
    auto protein_chains = detect_protein_chains(structure);

    // Combine results
    std::vector<ConnectedChain> all_chains;
    all_chains.reserve(rna_chains.size() + protein_chains.size());
    all_chains.insert(all_chains.end(), rna_chains.begin(), rna_chains.end());
    all_chains.insert(all_chains.end(), protein_chains.begin(), protein_chains.end());

    return all_chains;
}

// ============================================================================
// Connectivity checking
// ============================================================================

int ChainDetector::are_rna_residues_connected(const core::Residue& res1, const core::Residue& res2) const {
    auto bb1 = extract_backbone(res1);
    auto bb2 = extract_backbone(res2);

    // Check 5' to 3' connection: res1.O3' -> res2.P
    if (bb1.O3_prime && (bb2.P || bb2.PA)) {
        const auto& p_pos = bb2.P ? *bb2.P : *bb2.PA;
        double dist = bb1.O3_prime->distance_to(p_pos);
        if (dist < config_.rna_connectivity_cutoff) {
            return 1; // Forward connection
        }
    }

    // Check 3' to 5' connection: res2.O3' -> res1.P
    if (bb2.O3_prime && (bb1.P || bb1.PA)) {
        const auto& p_pos = bb1.P ? *bb1.P : *bb1.PA;
        double dist = bb2.O3_prime->distance_to(p_pos);
        if (dist < config_.rna_connectivity_cutoff) {
            return -1; // Reverse connection
        }
    }

    return 0; // Not connected
}

int ChainDetector::are_protein_residues_connected(const core::Residue& res1, const core::Residue& res2) const {
    auto bb1 = extract_backbone(res1);
    auto bb2 = extract_backbone(res2);

    // Check N-term to C-term connection: res1.C -> res2.N
    if (bb1.C && bb2.N) {
        double dist = bb1.C->distance_to(*bb2.N);
        if (dist < config_.protein_connectivity_cutoff) {
            return 1; // Forward connection
        }
    }

    // Check C-term to N-term connection: res2.C -> res1.N
    if (bb2.C && bb1.N) {
        double dist = bb2.C->distance_to(*bb1.N);
        if (dist < config_.protein_connectivity_cutoff) {
            return -1; // Reverse connection
        }
    }

    return 0; // Not connected
}

// ============================================================================
// Helper methods
// ============================================================================

BackboneConnectivity ChainDetector::extract_backbone(const core::Residue& residue) const {
    BackboneConnectivity bb;

    // RNA backbone atoms - use AtomType for O(1) lookup
    if (const auto* atom = residue.find_atom_by_type(AtomType::O3_PRIME)) {
        bb.O3_prime = atom->position();
    }
    if (const auto* atom = residue.find_atom_by_type(AtomType::P)) {
        bb.P = atom->position();
    }
    // PA not in AtomType enum - keep string lookup
    if (auto atom = residue.find_atom("PA")) {
        bb.PA = atom->position();
    }

    // Protein backbone atoms - use AtomType for O(1) lookup
    if (const auto* atom = residue.find_atom_by_type(AtomType::C)) {
        bb.C = atom->position();
    }
    if (const auto* atom = residue.find_atom_by_type(AtomType::N)) {
        bb.N = atom->position();
    }

    return bb;
}

std::vector<const core::Residue*> ChainDetector::filter_rna_residues(const core::Structure& structure) const {
    std::vector<const core::Residue*> rna_residues;

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            // Check if nucleotide
            if (residue.is_nucleotide()) {
                rna_residues.push_back(&residue);
                continue;
            }

            // Check for P or C1' atoms as fallback (for modified bases) - use AtomType
            if (residue.has_atom_type(AtomType::P) || residue.has_atom_type(AtomType::C1_PRIME)) {
                rna_residues.push_back(&residue);
            }
        }
    }

    return rna_residues;
}

std::vector<const core::Residue*> ChainDetector::filter_protein_residues(const core::Structure& structure) const {
    std::vector<const core::Residue*> protein_residues;

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            // Check if protein
            if (residue.is_protein()) {
                protein_residues.push_back(&residue);
                continue;
            }

            // Check for C and N atoms as fallback (for modified amino acids) - use AtomType
            if (residue.has_atom_type(AtomType::C) && residue.has_atom_type(AtomType::N)) {
                protein_residues.push_back(&residue);
            }
        }
    }

    return protein_residues;
}

void ChainDetector::sort_by_chain_and_num(std::vector<const core::Residue*>& residues) const {
    std::sort(residues.begin(), residues.end(), [](const core::Residue* a, const core::Residue* b) {
        // First by chain_id
        if (a->chain_id() != b->chain_id()) {
            return a->chain_id() < b->chain_id();
        }
        // Then by seq_num
        return a->seq_num() < b->seq_num();
    });
}

std::vector<ConnectedChain> ChainDetector::build_chains(std::vector<const core::Residue*>& residues,
                                                        const ConnectivityFunc& connectivity_func, bool is_rna) const {

    std::vector<ConnectedChain> chains;

    while (!residues.empty()) {
        // Start a new chain with first remaining residue
        ConnectedChain chain;
        chain.residues.push_back(residues.front());
        chain.chain_id = residues.front()->chain_id();
        chain.is_rna = is_rna;
        chain.is_protein = !is_rna;
        residues.erase(residues.begin());

        // Keep extending chain until no more connections found
        bool found_connection = true;
        while (found_connection) {
            found_connection = false;

            // Try to extend chain at both ends
            for (auto it = residues.begin(); it != residues.end(); ++it) {
                const auto* res = *it;

                // Only connect within same PDB chain
                if (res->chain_id() != chain.chain_id) {
                    continue;
                }

                // Check connection to end of chain (5' or N-term end)
                int conn_to_end = connectivity_func(*chain.residues.back(), *res);
                if (conn_to_end == 1) { // Forward connection
                    chain.residues.push_back(res);
                    residues.erase(it);
                    found_connection = true;
                    break;
                }

                // Check connection to start of chain (3' or C-term end)
                int conn_to_start = connectivity_func(*chain.residues.front(), *res);
                if (conn_to_start == -1) { // Reverse connection
                    chain.residues.insert(chain.residues.begin(), res);
                    residues.erase(it);
                    found_connection = true;
                    break;
                }
            }
        }

        chains.push_back(chain);
    }

    return chains;
}

std::vector<ConnectedChain> ChainDetector::merge_adjacent_chains(std::vector<ConnectedChain>& chains) const {
    if (chains.size() <= 1) {
        return chains;
    }

    bool merged_any = true;

    while (merged_any) {
        merged_any = false;
        std::vector<ConnectedChain> merged_chains;

        for (size_t i = 0; i < chains.size(); ++i) {
            // Handle last chain
            if (i == chains.size() - 1) {
                merged_chains.push_back(chains[i]);
                break;
            }

            const auto& chain1 = chains[i];
            const auto& chain2 = chains[i + 1];

            // Only merge chains from same PDB chain
            if (chain1.chain_id != chain2.chain_id) {
                merged_chains.push_back(chain1);
                continue;
            }

            // Check if chains are adjacent by sequence number
            const auto* res1_last = chain1.residues.back();
            const auto* res2_first = chain2.residues.front();

            int seq_diff = res2_first->seq_num() - res1_last->seq_num();
            bool should_merge = false;

            if (seq_diff == 1) {
                // Check physical distance
                double dist = get_residue_distance(*res1_last, *res2_first);

                // For RNA, try backbone distance first
                if (chain1.is_rna) {
                    auto bb1 = extract_backbone(*res1_last);
                    auto bb2 = extract_backbone(*res2_first);

                    // Try O3'-P distance
                    if (bb1.O3_prime && (bb2.P || bb2.PA)) {
                        const auto& p_pos = bb2.P ? *bb2.P : *bb2.PA;
                        double backbone_dist = bb1.O3_prime->distance_to(p_pos);
                        if (backbone_dist < 5.0) { // Relaxed cutoff for merging
                            should_merge = true;
                        }
                    }

                    // Fallback to sugar center distance
                    if (!should_merge && dist < config_.chain_merge_distance) {
                        should_merge = true;
                    }
                }
            }

            if (should_merge) {
                // Merge chain2 into chain1
                ConnectedChain merged = chain1;
                merged.residues.insert(merged.residues.end(), chain2.residues.begin(), chain2.residues.end());
                merged_chains.push_back(merged);
                merged_any = true;
                ++i; // Skip next chain since we merged it
            } else {
                merged_chains.push_back(chain1);
            }
        }

        chains = merged_chains;
    }

    return chains;
}

std::optional<geometry::Vector3D> ChainDetector::calculate_sugar_center(const core::Residue& residue) const {
    // Sugar atoms - use AtomType for O(1) lookup
    static constexpr AtomType sugar_atom_types[] = {
        AtomType::C1_PRIME, AtomType::C2_PRIME, AtomType::C3_PRIME,
        AtomType::C4_PRIME, AtomType::C5_PRIME, AtomType::O2_PRIME,
        AtomType::O3_PRIME, AtomType::O4_PRIME, AtomType::O5_PRIME
    };

    geometry::Vector3D sum(0, 0, 0);
    int count = 0;

    for (auto type : sugar_atom_types) {
        if (const auto* atom = residue.find_atom_by_type(type)) {
            sum += atom->position();
            ++count;
        }
    }

    if (count == 0) {
        return std::nullopt;
    }

    return sum / static_cast<double>(count);
}

double ChainDetector::get_residue_distance(const core::Residue& res1, const core::Residue& res2) const {
    // For RNA, use sugar center
    if (res1.is_nucleotide() && res2.is_nucleotide()) {
        auto center1 = calculate_sugar_center(res1);
        auto center2 = calculate_sugar_center(res2);
        if (center1 && center2) {
            return center1->distance_to(*center2);
        }
    }

    // For protein, use CA atoms - use AtomType for O(1) lookup
    if (res1.is_protein() && res2.is_protein()) {
        const auto* ca1 = res1.find_atom_by_type(AtomType::CA);
        const auto* ca2 = res2.find_atom_by_type(AtomType::CA);
        if (ca1 && ca2) {
            return ca1->position().distance_to(ca2->position());
        }
    }

    return std::numeric_limits<double>::infinity();
}

} // namespace algorithms
} // namespace x3dna

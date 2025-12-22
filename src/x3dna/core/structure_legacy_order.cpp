/**
 * @file structure_legacy_order.cpp
 * @brief Implementation of legacy order utilities
 *
 * Uses stored legacy_residue_idx values instead of recomputing from atoms.
 */

#include <x3dna/core/structure_legacy_order.hpp>
#include <x3dna/core/structure.hpp>
#include <algorithm>

namespace x3dna {
namespace core {

std::vector<const Residue*> get_residues_in_legacy_order(const Structure& structure) {
    // Collect all residues with their legacy indices
    std::vector<std::pair<int, const Residue*>> indexed_residues;

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            int legacy_idx = residue.legacy_residue_idx();
            if (legacy_idx > 0) {
                indexed_residues.push_back({legacy_idx, &residue});
            }
        }
    }

    // Sort by legacy_residue_idx
    std::sort(indexed_residues.begin(), indexed_residues.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    // Extract just the residue pointers
    std::vector<const Residue*> result;
    result.reserve(indexed_residues.size());
    for (const auto& [idx, res] : indexed_residues) {
        result.push_back(res);
    }

    return result;
}

const Residue* get_residue_by_legacy_idx(const Structure& structure, int legacy_idx) {
    return structure.get_residue_by_legacy_idx(legacy_idx);
}

int get_legacy_idx_for_residue(const Structure& structure, const Residue* residue) {
    return structure.get_legacy_idx_for_residue(residue);
}

} // namespace core
} // namespace x3dna

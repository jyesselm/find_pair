/**
 * @file residue_index_map.cpp
 * @brief Implementation of ResidueIndexMap
 */

#include <x3dna/algorithms/pair_identification/residue_index_map.hpp>
#include <x3dna/core/chain.hpp>
#include <limits>

namespace x3dna {
namespace algorithms {

void ResidueIndexMap::build(const core::Structure& structure) {
    clear();

    size_t modern_idx = 0;
    min_legacy_idx_ = std::numeric_limits<int>::max();
    max_legacy_idx_ = 0;

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            // Get legacy residue index from residue (set during PDB parsing)
            int legacy_idx = residue.legacy_residue_idx();
            if (legacy_idx > 0) {
                // Store mappings
                by_legacy_[legacy_idx] = &residue;
                by_modern_[modern_idx] = &residue;
                legacy_to_modern_[legacy_idx] = modern_idx;
                modern_to_legacy_[modern_idx] = legacy_idx;

                // Update range
                if (legacy_idx > max_legacy_idx_) {
                    max_legacy_idx_ = legacy_idx;
                }
                if (legacy_idx < min_legacy_idx_) {
                    min_legacy_idx_ = legacy_idx;
                }
            }
            modern_idx++;
        }
    }

    // Handle case where no valid indices were found
    if (by_legacy_.empty()) {
        min_legacy_idx_ = 0;
    }
}

void ResidueIndexMap::clear() {
    by_legacy_.clear();
    by_modern_.clear();
    legacy_to_modern_.clear();
    modern_to_legacy_.clear();
    max_legacy_idx_ = 0;
    min_legacy_idx_ = 0;
}

const core::Residue* ResidueIndexMap::get_by_legacy_idx(int legacy_idx) const {
    auto it = by_legacy_.find(legacy_idx);
    return (it != by_legacy_.end()) ? it->second : nullptr;
}

const core::Residue* ResidueIndexMap::get_by_modern_idx(size_t modern_idx) const {
    auto it = by_modern_.find(modern_idx);
    return (it != by_modern_.end()) ? it->second : nullptr;
}

bool ResidueIndexMap::has_legacy_idx(int legacy_idx) const {
    return by_legacy_.find(legacy_idx) != by_legacy_.end();
}

bool ResidueIndexMap::has_modern_idx(size_t modern_idx) const {
    return by_modern_.find(modern_idx) != by_modern_.end();
}

std::optional<size_t> ResidueIndexMap::to_modern(int legacy_idx) const {
    auto it = legacy_to_modern_.find(legacy_idx);
    if (it != legacy_to_modern_.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<int> ResidueIndexMap::to_legacy(size_t modern_idx) const {
    auto it = modern_to_legacy_.find(modern_idx);
    if (it != modern_to_legacy_.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::vector<int> ResidueIndexMap::legacy_indices() const {
    std::vector<int> result;
    result.reserve(by_legacy_.size());
    for (const auto& [legacy_idx, _] : by_legacy_) {
        result.push_back(legacy_idx);
    }
    return result;
}

} // namespace algorithms
} // namespace x3dna

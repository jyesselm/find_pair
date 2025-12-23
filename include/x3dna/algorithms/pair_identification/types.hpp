/**
 * @file types.hpp
 * @brief Internal types for base pair identification
 */

#pragma once

#include <x3dna/algorithms/pair_identification/base_pair_validator.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/io/json_writer.hpp>
#include <algorithm>
#include <map>
#include <vector>

namespace x3dna {
namespace algorithms {
namespace pair_identification {

/** @brief Results from Phase 1 validation of all pairs */
struct Phase1Results {
    std::map<std::pair<int, int>, ValidationResult> validation_results;
    std::map<std::pair<int, int>, int> bp_type_ids;

    [[nodiscard]] const ValidationResult* get_result(int idx1, int idx2) const {
        auto key = (idx1 < idx2) ? std::make_pair(idx1, idx2) : std::make_pair(idx2, idx1);
        auto it = validation_results.find(key);
        return (it != validation_results.end()) ? &it->second : nullptr;
    }

    [[nodiscard]] int get_bp_type_id(int idx1, int idx2) const {
        auto key = (idx1 < idx2) ? std::make_pair(idx1, idx2) : std::make_pair(idx2, idx1);
        auto it = bp_type_ids.find(key);
        return (it != bp_type_ids.end()) ? it->second : 0;
    }
};

/** @brief Mapping between legacy indices and residue pointers */
struct ResidueIndexMapping {
    std::map<int, const core::Residue*> by_legacy_idx;
    int max_legacy_idx = 0;

    [[nodiscard]] const core::Residue* get(int legacy_idx) const {
        auto it = by_legacy_idx.find(legacy_idx);
        return (it != by_legacy_idx.end()) ? it->second : nullptr;
    }

    [[nodiscard]] bool empty() const { return by_legacy_idx.empty(); }
};

/** @brief Context for partner search - groups related data to reduce parameters */
struct PartnerSearchContext {
    const std::vector<bool>& matched_indices;
    const ResidueIndexMapping& mapping;
    const Phase1Results& phase1;
    io::JsonWriter* writer;
};

/** @brief Mutable state during pair selection */
struct PairSelectionState {
    std::vector<bool> matched_indices;
    std::vector<core::BasePair> base_pairs;
    std::vector<std::pair<size_t, size_t>> selected_pairs_legacy_idx;
    std::vector<std::pair<int, int>> pairs_found_this_iteration;

    explicit PairSelectionState(int max_idx)
        : matched_indices(max_idx + 1, false) {}

    void mark_matched(int idx1, int idx2) {
        matched_indices[idx1] = true;
        matched_indices[idx2] = true;
    }

    [[nodiscard]] size_t count_matched() const {
        return std::count(matched_indices.begin(), matched_indices.end(), true);
    }
};

} // namespace pair_identification
} // namespace algorithms
} // namespace x3dna

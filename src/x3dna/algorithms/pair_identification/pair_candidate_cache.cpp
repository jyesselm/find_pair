/**
 * @file pair_candidate_cache.cpp
 * @brief Implementation of PairCandidateCache
 */

#include <x3dna/algorithms/pair_identification/pair_candidate_cache.hpp>

namespace x3dna {
namespace algorithms {

void PairCandidateCache::build(const core::Structure& structure, const BasePairValidator& validator,
                               const QualityScoreCalculator& quality_calc, NucleotideChecker is_nucleotide) {
    clear();

    // Build index map from structure
    index_map_.build(structure);

    if (index_map_.empty()) {
        return;
    }

    int max_idx = index_map_.max_legacy_idx();

    // PHASE 1: Validate ALL pairs (matches legacy check_pair loop)
    // Legacy: for (i = 1; i < num_residue; i++) { for (j = i + 1; j <= num_residue; j++) { ... } }
    // Note: Legacy uses i < num_residue, so i goes from 1 to num_residue-1 (inclusive)
    for (int legacy_idx1 = 1; legacy_idx1 <= max_idx - 1; ++legacy_idx1) {
        const core::Residue* res1 = index_map_.get_by_legacy_idx(legacy_idx1);
        if (!res1) {
            continue;
        }

        // Check if res1 is a nucleotide with a valid frame
        if (!is_nucleotide(*res1) || !res1->reference_frame().has_value()) {
            continue;
        }

        for (int legacy_idx2 = legacy_idx1 + 1; legacy_idx2 <= max_idx; ++legacy_idx2) {
            const core::Residue* res2 = index_map_.get_by_legacy_idx(legacy_idx2);
            if (!res2) {
                continue;
            }

            // Check if res2 is a nucleotide with a valid frame
            if (!is_nucleotide(*res2) || !res2->reference_frame().has_value()) {
                continue;
            }

            // Validate pair
            ValidationResult result = validator.validate(*res1, *res2);

            // Calculate adjusted quality score and bp_type_id
            double adjusted_score = quality_calc.calculate_selection_score(result, *res1, *res2);
            int bp_type_id = quality_calc.calculate_bp_type_id(*res1, *res2, result);

            // Store in cache (already normalized since legacy_idx1 < legacy_idx2)
            auto key = std::make_pair(legacy_idx1, legacy_idx2);
            CandidateInfo info{result, bp_type_id, adjusted_score};
            cache_[key] = info;

            // Track partners
            all_partners_[legacy_idx1].push_back(legacy_idx2);
            all_partners_[legacy_idx2].push_back(legacy_idx1);

            if (result.is_valid) {
                valid_partners_[legacy_idx1].push_back(legacy_idx2);
                valid_partners_[legacy_idx2].push_back(legacy_idx1);
            }
        }
    }
}

void PairCandidateCache::clear() {
    cache_.clear();
    valid_partners_.clear();
    all_partners_.clear();
    index_map_.clear();
}

size_t PairCandidateCache::valid_count() const {
    size_t count = 0;
    for (const auto& [key, info] : cache_) {
        if (info.is_valid()) {
            count++;
        }
    }
    return count;
}

std::optional<CandidateInfo> PairCandidateCache::get(int legacy_idx1, int legacy_idx2) const {
    auto key = normalize(legacy_idx1, legacy_idx2);
    auto it = cache_.find(key);
    if (it != cache_.end()) {
        return it->second;
    }
    return std::nullopt;
}

bool PairCandidateCache::contains(int legacy_idx1, int legacy_idx2) const {
    auto key = normalize(legacy_idx1, legacy_idx2);
    return cache_.find(key) != cache_.end();
}

std::vector<int> PairCandidateCache::valid_partners_for(int legacy_idx) const {
    auto it = valid_partners_.find(legacy_idx);
    if (it != valid_partners_.end()) {
        return it->second;
    }
    return {};
}

std::vector<std::pair<int, CandidateInfo>> PairCandidateCache::all_candidates_for(int legacy_idx) const {
    std::vector<std::pair<int, CandidateInfo>> result;

    auto it = all_partners_.find(legacy_idx);
    if (it == all_partners_.end()) {
        return result;
    }

    for (int partner_idx : it->second) {
        auto info = get(legacy_idx, partner_idx);
        if (info) {
            result.push_back({partner_idx, *info});
        }
    }
    return result;
}

void PairCandidateCache::for_each_valid(std::function<void(int, int, const CandidateInfo&)> callback) const {
    for (const auto& [key, info] : cache_) {
        if (info.is_valid()) {
            callback(key.first, key.second, info);
        }
    }
}

} // namespace algorithms
} // namespace x3dna

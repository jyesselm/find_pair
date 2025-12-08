/**
 * @file pair_candidate_cache.hpp
 * @brief Caches validation results for all candidate base pairs
 */

#pragma once

#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/quality_score_calculator.hpp>
#include <x3dna/algorithms/residue_index_map.hpp>
#include <x3dna/core/structure.hpp>
#include <map>
#include <vector>
#include <optional>
#include <functional>

namespace x3dna {
namespace algorithms {

/**
 * @struct CandidateInfo
 * @brief Information about a validated pair candidate
 */
struct CandidateInfo {
    ValidationResult validation;
    int bp_type_id;
    double adjusted_quality_score;

    bool is_valid() const { return validation.is_valid; }
};

/**
 * @class PairCandidateCache
 * @brief Caches validation results for all candidate base pairs
 *
 * Pre-computes and caches validation for all candidate pairs during Phase 1,
 * ensuring consistency between validation and selection phases.
 *
 * Usage:
 * @code
 * PairCandidateCache cache;
 * cache.build(structure, validator, quality_calc, is_nucleotide_func);
 *
 * // Get cached result for a specific pair
 * auto info = cache.get(legacy_idx1, legacy_idx2);
 * if (info && info->is_valid()) {
 *     // Use the cached validation result
 * }
 *
 * // Get all valid candidates for a residue
 * for (int partner_idx : cache.valid_partners_for(legacy_idx)) {
 *     // Process each valid partner
 * }
 * @endcode
 */
class PairCandidateCache {
public:
    using NucleotideChecker = std::function<bool(const core::Residue&)>;

    /**
     * @brief Build cache for all valid pairs in structure
     * @param structure Structure with frames calculated
     * @param validator Validator to use for pair checking
     * @param quality_calc Quality score calculator
     * @param is_nucleotide Function to check if residue is a nucleotide
     */
    void build(
        const core::Structure& structure,
        const BasePairValidator& validator,
        const QualityScoreCalculator& quality_calc,
        NucleotideChecker is_nucleotide
    );

    /**
     * @brief Clear all cached data
     */
    void clear();

    /**
     * @brief Check if cache is empty
     */
    bool empty() const { return cache_.empty(); }

    /**
     * @brief Get number of cached pairs
     */
    size_t size() const { return cache_.size(); }

    /**
     * @brief Get number of valid pairs
     */
    size_t valid_count() const;

    // ==================== Lookups ====================

    /**
     * @brief Get cached result for a pair (order-independent)
     * @param legacy_idx1 First residue legacy index
     * @param legacy_idx2 Second residue legacy index
     * @return CandidateInfo or nullopt if pair wasn't cached
     */
    std::optional<CandidateInfo> get(int legacy_idx1, int legacy_idx2) const;

    /**
     * @brief Check if pair exists in cache
     */
    bool contains(int legacy_idx1, int legacy_idx2) const;

    /**
     * @brief Get all valid partner indices for a residue
     * @param legacy_idx Residue legacy index
     * @return Vector of legacy indices that form valid pairs with this residue
     */
    std::vector<int> valid_partners_for(int legacy_idx) const;

    /**
     * @brief Get all candidates (valid or not) for a residue
     */
    std::vector<std::pair<int, CandidateInfo>> all_candidates_for(int legacy_idx) const;

    // ==================== Iteration ====================

    /**
     * @brief Get all cached pairs (for iteration)
     */
    const std::map<std::pair<int, int>, CandidateInfo>& all() const { return cache_; }

    /**
     * @brief Iterate over all valid pairs
     * @param callback Function called with (legacy_idx1, legacy_idx2, info) for each valid pair
     */
    void for_each_valid(std::function<void(int, int, const CandidateInfo&)> callback) const;

    // ==================== Index Map Access ====================

    /**
     * @brief Get the residue index map used during build
     */
    const ResidueIndexMap& index_map() const { return index_map_; }

    /**
     * @brief Get maximum legacy index
     */
    int max_legacy_idx() const { return index_map_.max_legacy_idx(); }

private:
    /**
     * @brief Normalize pair key (smaller index first)
     */
    static std::pair<int, int> normalize(int i, int j) {
        return (i < j) ? std::make_pair(i, j) : std::make_pair(j, i);
    }

    std::map<std::pair<int, int>, CandidateInfo> cache_;
    std::map<int, std::vector<int>> valid_partners_;  // legacy_idx -> valid partner indices
    std::map<int, std::vector<int>> all_partners_;    // legacy_idx -> all partner indices
    ResidueIndexMap index_map_;
};

}  // namespace algorithms
}  // namespace x3dna


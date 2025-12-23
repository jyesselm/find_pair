/**
 * @file base_pair_finder.hpp
 * @brief Base pair finding algorithm (matches legacy find_bestpair)
 */

#pragma once

#include <x3dna/common_types.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/algorithms/pair_identification/base_pair_validator.hpp>
#include <x3dna/algorithms/pair_identification/quality_score_calculator.hpp>
#include <x3dna/algorithms/pair_identification/pair_candidate_cache.hpp>
#include <x3dna/algorithms/pair_identification/pair_finding_observer.hpp>
#include <x3dna/algorithms/pair_identification/pair_selection_strategy.hpp>
#include <x3dna/io/json_writer.hpp>
#include <algorithm>
#include <memory>
#include <optional>
#include <vector>

namespace x3dna {
namespace algorithms {

// PairFindingStrategy is defined in common_types.hpp

/**
 * @class BasePairFinder
 * @brief Finds base pairs in a structure using various strategies
 *
 * This class serves as a facade over several specialized components:
 * - BasePairValidator: Validates individual base pairs
 * - QualityScoreCalculator: Calculates adjusted quality scores
 * - PairCandidateCache: Caches Phase 1 validation results
 * - IPairSelectionStrategy: Implements pair selection algorithm
 * - IPairFindingObserver: Records events during pair finding
 *
 * Implements the legacy find_bestpair algorithm which uses a greedy mutual
 * best match strategy: for each unpaired residue, find its best partner,
 * then check if that partner's best partner is the original residue.
 */
class BasePairFinder {
public:
    /**
     * @brief Constructor
     * @param params Validation parameters
     */
    explicit BasePairFinder(const ValidationParameters& params = ValidationParameters::defaults())
        : validator_(params), strategy_(PairFindingStrategy::BEST_PAIR) {}

    /**
     * @brief Find base pairs in a structure
     * @param structure Structure to search (residues must have frames calculated)
     * @return Vector of found base pairs
     */
    [[nodiscard]] std::vector<core::BasePair> find_pairs(core::Structure& structure);

    /**
     * @brief Find base pairs (const version)
     */
    [[nodiscard]] std::vector<core::BasePair> find_pairs(const core::Structure& structure) const;

    /**
     * @brief Find base pairs and record validation results to JSON
     * @param structure Structure to search
     * @param writer JsonWriter to record validation results (can be nullptr to skip recording)
     * @return Vector of found base pairs
     */
    [[nodiscard]] std::vector<core::BasePair> find_pairs_with_recording(core::Structure& structure,
                                                                        io::JsonWriter* writer = nullptr) const;

    /**
     * @brief Set finding strategy
     */
    void set_strategy(PairFindingStrategy strategy) {
        strategy_ = strategy;
    }

    /**
     * @brief Get finding strategy
     */
    [[nodiscard]] PairFindingStrategy strategy() const {
        return strategy_;
    }

    /**
     * @brief Set validation parameters
     */
    void set_parameters(const ValidationParameters& params) {
        validator_.set_parameters(params);
    }

    /**
     * @brief Get validation parameters
     */
    [[nodiscard]] const ValidationParameters& parameters() const {
        return validator_.parameters();
    }

    /**
     * @brief Check if residue is a nucleotide
     * @param residue Residue to check
     * @return True if residue is a nucleotide (standard or modified)
     */
    [[nodiscard]] static bool is_nucleotide(const core::Residue& residue);

private:
    // Core components
    BasePairValidator validator_;
    QualityScoreCalculator quality_calculator_;
    PairFindingStrategy strategy_;
    mutable PairCandidateCache cache_;

    // ============================================================================
    // Internal types - must be defined before methods that use them
    // ============================================================================

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

    // ============================================================================
    // Private methods
    // ============================================================================

    [[nodiscard]] std::vector<core::BasePair> find_best_pairs(core::Structure& structure,
                                                              io::JsonWriter* writer = nullptr) const;
    [[nodiscard]] std::vector<core::BasePair> find_all_pairs(const core::Structure& structure) const;
    [[nodiscard]] std::optional<std::pair<int, ValidationResult>> find_best_partner(
        int legacy_idx, const PartnerSearchContext& ctx) const;

    [[nodiscard]] double adjust_pair_quality(const std::vector<core::hydrogen_bond>& hbonds) const;
    [[nodiscard]] int calculate_bp_type_id(const core::Residue* res1, const core::Residue* res2,
                                           const ValidationResult& result, double quality_score) const;
    [[nodiscard]] double calculate_adjusted_score(const ValidationResult& result, int bp_type_id) const;

    void record_validation_results(int legacy_idx1, int legacy_idx2, const core::Residue* res1,
                                   const core::Residue* res2, const ValidationResult& result,
                                   io::JsonWriter* writer) const;

    [[nodiscard]] static size_t get_residue_index(const core::Structure& structure, const core::Residue& residue);
    [[nodiscard]] static bool can_participate_in_pairing(const core::Residue* res);
    [[nodiscard]] static bool is_matched(int legacy_idx, const std::vector<bool>& matched);

    [[nodiscard]] ResidueIndexMapping build_residue_index_mapping(const core::Structure& structure) const;
    [[nodiscard]] Phase1Results run_phase1_validation(const ResidueIndexMapping& mapping) const;
    [[nodiscard]] core::BasePair create_base_pair(int legacy_idx1, int legacy_idx2, const core::Residue* res1,
                                                  const core::Residue* res2, const ValidationResult& result) const;
    [[nodiscard]] bool try_select_mutual_pair(int legacy_idx1, int legacy_idx2,
                                               const core::Residue* res1, const core::Residue* res2,
                                               const ValidationResult& result,
                                               const PartnerSearchContext& ctx,
                                               PairSelectionState& state) const;
};

} // namespace algorithms
} // namespace x3dna

/**
 * @file pair_finding_observer.hpp
 * @brief Observer interface for base pair finding events
 */

#pragma once

#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/residue.hpp>
#include <vector>
#include <utility>

namespace x3dna {
namespace algorithms {

/**
 * @struct BestPartnerCandidate
 * @brief Information about a candidate during best partner selection
 */
struct BestPartnerCandidate {
    int partner_legacy_idx;
    double quality_score;
    int bp_type_id;
    bool is_valid;
};

/**
 * @class IPairFindingObserver
 * @brief Interface for observing base pair finding events
 *
 * Implementations can record, log, or process events during pair finding.
 * This decouples the recording logic from the core algorithm.
 */
class IPairFindingObserver {
public:
    virtual ~IPairFindingObserver() = default;

    // ==================== Phase 1: Validation Events ====================

    /**
     * @brief Called when a pair is validated during Phase 1
     * @param legacy_idx1 First residue legacy index
     * @param legacy_idx2 Second residue legacy index
     * @param res1 First residue
     * @param res2 Second residue
     * @param result Validation result
     * @param bp_type_id Base pair type ID (-1, 0, 1, or 2)
     */
    virtual void on_pair_validated(int legacy_idx1, int legacy_idx2, const core::Residue& res1,
                                   const core::Residue& res2, const ValidationResult& result, int bp_type_id) = 0;

    // ==================== Phase 2: Selection Events ====================

    /**
     * @brief Called when best partner candidates are evaluated for a residue
     * @param legacy_idx Residue being evaluated
     * @param candidates All candidates considered
     * @param best_partner_idx Index of selected best partner (0 if none)
     * @param best_score Best score found
     */
    virtual void on_best_partner_candidates(int legacy_idx, const std::vector<BestPartnerCandidate>& candidates,
                                            int best_partner_idx, double best_score) = 0;

    /**
     * @brief Called when mutual best partner check is performed
     * @param legacy_idx1 First residue legacy index
     * @param legacy_idx2 Second residue legacy index
     * @param best_j_for_i Best partner of idx1
     * @param best_i_for_j Best partner of idx2
     * @param is_mutual Whether they are mutual best partners
     * @param was_selected Whether the pair was selected
     */
    virtual void on_mutual_best_check(int legacy_idx1, int legacy_idx2, int best_j_for_i, int best_i_for_j,
                                      bool is_mutual, bool was_selected) = 0;

    /**
     * @brief Called after each iteration of best-pair selection
     * @param iteration_num Iteration number (1-based)
     * @param pairs_this_iteration Pairs selected in this iteration
     * @param matched_indices Array of matched status per legacy index
     * @param total_matched Total number of matched residues
     */
    virtual void on_iteration_complete(int iteration_num, const std::vector<std::pair<int, int>>& pairs_this_iteration,
                                       const std::vector<bool>& matched_indices, size_t total_matched) = 0;

    // ==================== Final Results ====================

    /**
     * @brief Called with final selection results
     * @param selected_pairs All selected pairs (legacy indices)
     */
    virtual void on_selection_complete(const std::vector<std::pair<int, int>>& selected_pairs) = 0;

    /**
     * @brief Called with final base pairs
     * @param base_pairs The final BasePair objects
     */
    virtual void on_pairs_finalized(const std::vector<core::BasePair>& base_pairs) = 0;
};

/**
 * @class NullPairFindingObserver
 * @brief No-op observer for when recording is not needed
 */
class NullPairFindingObserver : public IPairFindingObserver {
public:
    void on_pair_validated(int, int, const core::Residue&, const core::Residue&, const ValidationResult&,
                           int) override {}
    void on_best_partner_candidates(int, const std::vector<BestPartnerCandidate>&, int, double) override {}
    void on_mutual_best_check(int, int, int, int, bool, bool) override {}
    void on_iteration_complete(int, const std::vector<std::pair<int, int>>&, const std::vector<bool>&,
                               size_t) override {}
    void on_selection_complete(const std::vector<std::pair<int, int>>&) override {}
    void on_pairs_finalized(const std::vector<core::BasePair>&) override {}
};

} // namespace algorithms
} // namespace x3dna

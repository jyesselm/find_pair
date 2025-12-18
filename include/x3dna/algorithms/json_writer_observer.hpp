/**
 * @file json_writer_observer.hpp
 * @brief Observer implementation that records to JsonWriter
 */

#pragma once

#include <x3dna/algorithms/pair_finding_observer.hpp>
#include <x3dna/io/json_writer.hpp>

namespace x3dna {
namespace algorithms {

/**
 * @class JsonWriterObserver
 * @brief Pair finding observer that records events to a JsonWriter
 *
 * Bridges the Observer pattern to the existing JsonWriter implementation,
 * maintaining backward compatibility with the current recording format.
 */
class JsonWriterObserver : public IPairFindingObserver {
public:
    /**
     * @brief Construct observer with a JsonWriter
     * @param writer The JsonWriter to record to (must outlive this observer)
     * @param params Validation parameters (for recording)
     */
    explicit JsonWriterObserver(io::JsonWriter& writer, const ValidationParameters& params)
        : writer_(writer), params_(params) {}

    void on_pair_validated(int legacy_idx1, int legacy_idx2, const core::Residue& res1, const core::Residue& res2,
                           const ValidationResult& result, int bp_type_id) override;

    void on_best_partner_candidates(int legacy_idx, const std::vector<BestPartnerCandidate>& candidates,
                                    int best_partner_idx, double best_score) override;

    void on_mutual_best_check(int legacy_idx1, int legacy_idx2, int best_j_for_i, int best_i_for_j, bool is_mutual,
                              bool was_selected) override;

    void on_iteration_complete(int iteration_num, const std::vector<std::pair<int, int>>& pairs_this_iteration,
                               const std::vector<bool>& matched_indices, size_t total_matched) override;

    void on_selection_complete(const std::vector<std::pair<int, int>>& selected_pairs) override;

    void on_pairs_finalized(const std::vector<core::BasePair>& base_pairs) override;

private:
    io::JsonWriter& writer_;
    ValidationParameters params_;
};

} // namespace algorithms
} // namespace x3dna

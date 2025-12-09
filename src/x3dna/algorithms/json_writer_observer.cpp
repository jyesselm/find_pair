/**
 * @file json_writer_observer.cpp
 * @brief Implementation of JsonWriterObserver
 */

#include <x3dna/algorithms/json_writer_observer.hpp>
#include <x3dna/algorithms/quality_score_calculator.hpp>

namespace x3dna {
namespace algorithms {

void JsonWriterObserver::on_pair_validated(
    int legacy_idx1, int legacy_idx2,
    const core::Residue& /* res1 */, const core::Residue& /* res2 */,
    const ValidationResult& result, int bp_type_id
) {
    // Convert to 0-based indices for JsonWriter
    size_t base_i = static_cast<size_t>(legacy_idx1 - 1);
    size_t base_j = static_cast<size_t>(legacy_idx2 - 1);

    // Check if passes distance/angle checks (cdns) - matches legacy behavior
    bool passes_cdns = result.distance_check && result.d_v_check && 
                       result.plane_angle_check && result.dNN_check;

    // Only record when i < j to avoid recording both (i,j) and (j,i)
    // This cuts the output file size in half
    if (passes_cdns && legacy_idx1 < legacy_idx2) {
        // Prepare rtn_val array: [dorg, d_v, plane_angle, dNN, quality_score]
        // Note: The quality_score in rtn_val has already been adjusted and bp_type_id
        // bonus applied by the time on_pair_validated is called
        std::array<double, 5> rtn_val = {
            result.dorg, result.d_v, result.plane_angle, result.dNN, result.quality_score
        };

        // Record validation results (only for valid pairs to match legacy behavior)
        if (result.is_valid) {
            writer_.record_pair_validation(base_i, base_j, result.is_valid, bp_type_id,
                                           result.dir_x, result.dir_y, result.dir_z,
                                           rtn_val, params_);
        }
    }

    // Record distance checks only if also passes hydrogen bond check and i < j
    if (result.hbond_check && legacy_idx1 < legacy_idx2) {
        writer_.record_distance_checks(base_i, base_j, result.dorg, result.dNN,
                                       result.plane_angle, result.d_v, result.overlap_area);
    }

    // Record H-bond list if present and i < j
    if (!result.hbonds.empty() && legacy_idx1 < legacy_idx2) {
        writer_.record_hbond_list(base_i, base_j, result.hbonds);
    }
}

void JsonWriterObserver::on_best_partner_candidates(
    int legacy_idx,
    const std::vector<BestPartnerCandidate>& candidates,
    int best_partner_idx,
    double best_score
) {
    // Convert to tuple format expected by JsonWriter: (res_j, is_eligible, score, bp_type_id)
    std::vector<std::tuple<int, bool, double, int>> json_candidates;
    json_candidates.reserve(candidates.size());

    for (const auto& c : candidates) {
        json_candidates.emplace_back(c.partner_legacy_idx, c.is_valid, c.quality_score, c.bp_type_id);
    }

    writer_.record_best_partner_candidates(legacy_idx, json_candidates, best_partner_idx, best_score);
}

void JsonWriterObserver::on_mutual_best_check(
    int legacy_idx1, int legacy_idx2,
    int best_j_for_i, int best_i_for_j,
    bool is_mutual, bool was_selected
) {
    writer_.record_mutual_best_decision(legacy_idx1, legacy_idx2,
                                        best_j_for_i, best_i_for_j,
                                        is_mutual, was_selected);
}

void JsonWriterObserver::on_iteration_complete(
    int iteration_num,
    const std::vector<std::pair<int, int>>& pairs_this_iteration,
    const std::vector<bool>& matched_indices,
    size_t total_matched
) {
    int max_legacy_idx = static_cast<int>(matched_indices.size() - 1);
    writer_.record_iteration_state(iteration_num, total_matched, max_legacy_idx,
                                   matched_indices, pairs_this_iteration);
}

void JsonWriterObserver::on_selection_complete(
    const std::vector<std::pair<int, int>>& selected_pairs
) {
    // Convert int pairs to size_t pairs for JsonWriter
    std::vector<std::pair<size_t, size_t>> converted;
    converted.reserve(selected_pairs.size());
    for (const auto& p : selected_pairs) {
        converted.emplace_back(static_cast<size_t>(p.first), static_cast<size_t>(p.second));
    }
    writer_.record_find_bestpair_selection(converted);
}

void JsonWriterObserver::on_pairs_finalized(
    const std::vector<core::BasePair>& base_pairs
) {
    for (const auto& pair : base_pairs) {
        writer_.record_base_pair(pair);
    }
}

}  // namespace algorithms
}  // namespace x3dna


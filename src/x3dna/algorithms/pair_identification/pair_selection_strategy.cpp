/**
 * @file pair_selection_strategy.cpp
 * @brief Implementation of pair selection strategies
 */

#include <x3dna/algorithms/pair_identification/pair_selection_strategy.hpp>
#include <limits>

namespace x3dna {
namespace algorithms {

std::vector<std::pair<int, int>> MutualBestStrategy::select(SelectionContext& context, IPairFindingObserver* observer) {
    std::vector<std::pair<int, int>> selected_pairs;
    std::vector<std::pair<int, int>> pairs_found_this_iteration;

    size_t num_matched_prev = 0;
    size_t num_matched_curr = 0;
    int iteration_num = 0;

    do {
        iteration_num++;
        num_matched_prev = num_matched_curr;
        num_matched_curr = 0;
        pairs_found_this_iteration.clear();

        // Count current matches
        for (bool matched : context.matched_indices) {
            if (matched) {
                num_matched_curr++;
            }
        }

        // Try to find pairs for each unpaired residue
        // CRITICAL: Iterate sequentially from 1 to max_legacy_idx to match legacy iteration order
        for (int legacy_idx1 = 1; legacy_idx1 <= context.max_legacy_idx; ++legacy_idx1) {
            // Skip if already matched
            if (legacy_idx1 >= static_cast<int>(context.matched_indices.size()) ||
                context.matched_indices[legacy_idx1]) {
                continue;
            }

            // Skip if residue doesn't have valid candidates
            auto partners = context.cache.valid_partners_for(legacy_idx1);
            if (partners.empty()) {
                continue;
            }

            // Find best partner for this residue
            auto best_partner = find_best_partner(legacy_idx1, context, observer);
            if (!best_partner.has_value()) {
                continue;
            }

            int legacy_idx2 = best_partner->first;

            // Check if res_idx2's best partner is res_idx1 (mutual best match)
            auto partner_of_partner = find_best_partner(legacy_idx2, context, observer);
            bool is_mutual = (partner_of_partner.has_value() && partner_of_partner->first == legacy_idx1);

            // Notify observer of mutual best check
            if (observer) {
                int best_j_for_i = legacy_idx2;
                int best_i_for_j = partner_of_partner.has_value() ? partner_of_partner->first : 0;
                observer->on_mutual_best_check(legacy_idx1, legacy_idx2, best_j_for_i, best_i_for_j, is_mutual,
                                               is_mutual);
            }

            if (is_mutual) {
                // Verify the pair is valid in the cache
                auto info = context.cache.get(legacy_idx1, legacy_idx2);
                if (!info || !info->is_valid()) {
                    continue;
                }

                // Mutual best match found
                context.matched_indices[legacy_idx1] = true;
                context.matched_indices[legacy_idx2] = true;

                // Store with smaller index first for consistency
                if (legacy_idx1 < legacy_idx2) {
                    selected_pairs.emplace_back(legacy_idx1, legacy_idx2);
                    pairs_found_this_iteration.emplace_back(legacy_idx1, legacy_idx2);
                } else {
                    selected_pairs.emplace_back(legacy_idx2, legacy_idx1);
                    pairs_found_this_iteration.emplace_back(legacy_idx2, legacy_idx1);
                }
            }
        }

        // Notify observer of iteration complete
        if (observer) {
            size_t total_matched = 0;
            for (bool m : context.matched_indices) {
                if (m)
                    total_matched++;
            }
            observer->on_iteration_complete(iteration_num, pairs_found_this_iteration, context.matched_indices,
                                            total_matched);
        }

        // Recount matches after this iteration
        num_matched_curr = 0;
        for (bool matched : context.matched_indices) {
            if (matched) {
                num_matched_curr++;
            }
        }

    } while (num_matched_curr > num_matched_prev);

    return selected_pairs;
}

std::optional<std::pair<int, double>> MutualBestStrategy::find_best_partner(int legacy_idx,
                                                                            const SelectionContext& context,
                                                                            IPairFindingObserver* observer) const {
    std::optional<std::pair<int, double>> best_result;
    double best_score = std::numeric_limits<double>::max();

    // Collect candidates for observer
    std::vector<BestPartnerCandidate> candidates;

    // Get all valid partners from cache
    auto partners = context.cache.valid_partners_for(legacy_idx);

    for (int partner_idx : partners) {
        // Skip if partner is already matched
        if (partner_idx >= static_cast<int>(context.matched_indices.size()) || context.matched_indices[partner_idx]) {
            continue;
        }

        // Get cached info
        auto info = context.cache.get(legacy_idx, partner_idx);
        if (!info || !info->is_valid()) {
            continue;
        }

        // Record candidate for observer
        if (observer) {
            BestPartnerCandidate c;
            c.partner_legacy_idx = partner_idx;
            c.quality_score = info->adjusted_quality_score;
            c.bp_type_id = info->bp_type_id;
            c.is_valid = info->is_valid();
            candidates.push_back(c);
        }

        // Update best if this score is better (lower is better)
        if (info->adjusted_quality_score < best_score) {
            best_score = info->adjusted_quality_score;
            best_result = std::make_pair(partner_idx, best_score);
        }
    }

    // Notify observer
    if (observer && !candidates.empty()) {
        int best_partner = best_result.has_value() ? best_result->first : 0;
        observer->on_best_partner_candidates(legacy_idx, candidates, best_partner, best_score);
    }

    return best_result;
}

} // namespace algorithms
} // namespace x3dna

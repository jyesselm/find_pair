/**
 * @file pair_selection_strategy.hpp
 * @brief Strategy interface for base pair selection algorithms
 */

#pragma once

#include <x3dna/algorithms/pair_candidate_cache.hpp>
#include <x3dna/algorithms/pair_finding_observer.hpp>
#include <vector>
#include <utility>

namespace x3dna {
namespace algorithms {

/**
 * @struct SelectionContext
 * @brief Context provided to selection strategies
 */
struct SelectionContext {
    const PairCandidateCache& cache;
    std::vector<bool>& matched_indices;
    int max_legacy_idx;
};

/**
 * @class IPairSelectionStrategy
 * @brief Interface for pair selection algorithms
 *
 * Different strategies can implement different selection policies:
 * - MutualBestStrategy: Legacy behavior - only select mutual best partners
 * - BestAvailableStrategy: Select best partner without mutual check
 * - ScoreThresholdStrategy: Select all pairs above a quality threshold
 */
class IPairSelectionStrategy {
public:
    virtual ~IPairSelectionStrategy() = default;

    /**
     * @brief Select base pairs from validated candidates
     * @param context Selection context with cache and state
     * @param observer Observer for recording selection events (may be null)
     * @return Selected pairs as (legacy_idx1, legacy_idx2) with additional info
     */
    virtual std::vector<std::pair<int, int>> select(
        SelectionContext& context,
        IPairFindingObserver* observer
    ) = 0;

    /**
     * @brief Get the name of this strategy (for logging/debugging)
     */
    virtual std::string name() const = 0;
};

/**
 * @class MutualBestStrategy
 * @brief Legacy selection strategy - select only mutual best partners
 *
 * This is the default strategy that matches legacy X3DNA behavior:
 * - For each unmatched residue, find its best partner
 * - Check if that partner's best partner is the original residue
 * - Only select pairs that are mutual best partners
 * - Iterate until no new pairs can be found
 */
class MutualBestStrategy : public IPairSelectionStrategy {
public:
    std::vector<std::pair<int, int>> select(
        SelectionContext& context,
        IPairFindingObserver* observer
    ) override;

    std::string name() const override { return "MutualBest"; }

private:
    /**
     * @brief Find best partner for a residue
     * @param legacy_idx Legacy index of residue to find partner for
     * @param context Selection context
     * @param observer Observer for recording candidates
     * @return (partner_idx, adjusted_score) or nullopt if no valid partner
     */
    std::optional<std::pair<int, double>> find_best_partner(
        int legacy_idx,
        const SelectionContext& context,
        IPairFindingObserver* observer
    ) const;
};

}  // namespace algorithms
}  // namespace x3dna


/**
 * @file five_to_three_orderer.hpp
 * @brief Main five2three ordering algorithm orchestration
 *
 * Extracts the ensure_five_to_three algorithm from helix_organizer,
 * composing all the helix helper classes for proper 5'→3' ordering.
 */

#pragma once

#include <x3dna/core/base_pair.hpp>
#include <x3dna/algorithms/helix_organizer.hpp> // For types
#include <x3dna/algorithms/helix/strand_direction_checker.hpp>
#include <vector>

namespace x3dna::algorithms::helix {

// Re-use types from helix_organizer
using algorithms::BackboneData;
using algorithms::DirectionCounts;
using algorithms::HelixSegment;

/**
 * @struct FiveToThreeConfig
 * @brief Configuration for five-to-three ordering
 */
struct FiveToThreeConfig {
    double end_stack_xang = 125.0; ///< Max x-angle for stacked WC pairs
    double o3p_upper = 2.5;        ///< Max O3'-P distance for backbone linkage
};

/**
 * @class FiveToThreeOrderer
 * @brief Orchestrates the five2three algorithm for proper strand direction
 *
 * This class implements the main ensure_five_to_three algorithm that
 * ensures base pairs are ordered with proper 5'→3' strand direction.
 * It composes StrandDirectionChecker for the individual checks.
 */
class FiveToThreeOrderer {
public:
    explicit FiveToThreeOrderer(const FiveToThreeConfig& config = {})
        : config_(config), direction_checker_({config.end_stack_xang, config.o3p_upper}) {}

    /**
     * @brief Ensure 5'→3' direction for all helices
     *
     * Applies the full five2three algorithm to ensure proper strand
     * direction in each helix segment.
     *
     * @param pairs Vector of base pairs
     * @param backbone Backbone connectivity data
     * @param pair_order Ordered pair indices (may be modified)
     * @param helices Helix segments (may be modified)
     * @param swapped Output: strand swap flags for each pair
     */
    void ensure_five_to_three(const std::vector<core::BasePair>& pairs, const BackboneData& backbone,
                              std::vector<size_t>& pair_order, std::vector<HelixSegment>& helices,
                              std::vector<bool>& swapped) const;

    [[nodiscard]] const FiveToThreeConfig& config() const {
        return config_;
    }
    [[nodiscard]] const StrandDirectionChecker& direction_checker() const {
        return direction_checker_;
    }

private:
    FiveToThreeConfig config_;
    StrandDirectionChecker direction_checker_;
};

} // namespace x3dna::algorithms::helix

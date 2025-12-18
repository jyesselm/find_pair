/**
 * @file helix_context_calculator.hpp
 * @brief Neighbor context calculation for helix organization
 *
 * Extracts bp_context/locate_helix logic from helix_organizer for
 * determining pair neighbors and helix endpoints.
 */

#pragma once

#include <x3dna/core/base_pair.hpp>
#include <x3dna/algorithms/helix_organizer.hpp>  // For BackboneData, HelixSegment, PairContextInfo
#include <x3dna/algorithms/helix/pair_geometry_helper.hpp>
#include <x3dna/algorithms/helix/backbone_linkage_checker.hpp>
#include <vector>
#include <optional>
#include <utility>

namespace x3dna::algorithms::helix {

// Re-use types from helix_organizer
using algorithms::BackboneData;
using algorithms::HelixSegment;
using algorithms::PairContextInfo;

/**
 * @struct HelixContextConfig
 * @brief Configuration for context calculation
 */
struct HelixContextConfig {
    double helix_break = 7.8;      ///< Max distance (Å) between adjacent pairs in helix
    double neighbor_cutoff = 8.5;  ///< Cutoff for neighbor detection
};

/**
 * @struct PairContext
 * @brief Neighbor information for a base pair (internal use)
 */
struct PairContext {
    bool is_endpoint = true;         ///< True if at helix end (< 2 neighbors)
    std::optional<size_t> neighbor1; ///< Nearest neighbor
    std::optional<size_t> neighbor2; ///< 2nd neighbor (opposite z-side)
    double dist1 = 0.0;              ///< Distance to neighbor1
    double dist2 = 0.0;              ///< Distance to neighbor2
    bool has_backbone_link1 = false; ///< Backbone connected to neighbor1
    bool has_backbone_link2 = false; ///< Backbone connected to neighbor2
};

/**
 * @class HelixContextCalculator
 * @brief Calculates neighbor context and locates helices
 *
 * Implements the legacy bp_context and locate_helix algorithms for
 * determining which pairs are neighbors and organizing them into helices.
 */
class HelixContextCalculator {
public:
    explicit HelixContextCalculator(const HelixContextConfig& config = {})
        : config_(config), linkage_checker_({config.helix_break > 0 ? 2.5 : 2.5}) {}

    /**
     * @brief Calculate neighbor context for all pairs
     *
     * For each pair, finds the nearest neighbors and determines if
     * it's a helix endpoint.
     *
     * @param pairs Vector of base pairs
     * @param backbone Backbone connectivity data
     * @return Vector of PairContext for each pair
     */
    [[nodiscard]] std::vector<PairContext> calculate_context(
        const std::vector<core::BasePair>& pairs,
        const BackboneData& backbone) const;

    /**
     * @brief Find helix endpoints from context
     *
     * Endpoints are pairs with fewer than 2 neighbors within helix_break distance.
     *
     * @param context Calculated pair contexts
     * @return Vector of endpoint pair indices
     */
    [[nodiscard]] std::vector<size_t> find_endpoints(
        const std::vector<PairContext>& context) const;

    /**
     * @brief Locate and chain pairs into helices
     *
     * Starting from endpoints, traverses neighbors to build continuous
     * helix segments.
     *
     * @param context Calculated pair contexts
     * @param endpoints Endpoint pair indices
     * @param backbone Backbone connectivity data (unused, for API compatibility)
     * @param num_pairs Total number of pairs
     * @return Pair of (ordered pair indices, helix segments)
     */
    [[nodiscard]] std::pair<std::vector<size_t>, std::vector<HelixSegment>> locate_helices(
        const std::vector<PairContext>& context,
        const std::vector<size_t>& endpoints,
        const BackboneData& backbone,
        size_t num_pairs) const;

    /**
     * @brief Convert internal PairContext to public PairContextInfo
     */
    [[nodiscard]] static std::vector<PairContextInfo> to_public_context(
        const std::vector<PairContext>& context);

    [[nodiscard]] const HelixContextConfig& config() const { return config_; }

private:
    HelixContextConfig config_;
    BackboneLinkageChecker linkage_checker_;

    /**
     * @brief Check if two z-distances are on opposite sides
     */
    [[nodiscard]] static bool are_on_opposite_z_sides(double d1, double d2) {
        return d1 * d2 < 0;
    }
};

} // namespace x3dna::algorithms::helix

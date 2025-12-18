/**
 * @file backbone_linkage_checker.hpp
 * @brief Backbone linkage detection for helix organization
 *
 * Extracts O3'-P linkage detection logic from helix_organizer to enable
 * determination of 5'→3' strand direction.
 */

#pragma once

#include <x3dna/core/base_pair.hpp>
#include <x3dna/algorithms/helix_organizer.hpp>  // For LinkDirection, BackboneData

namespace x3dna::algorithms::helix {

// Re-use types from helix_organizer to avoid duplication
using algorithms::LinkDirection;
using algorithms::BackboneData;

/**
 * @struct BackboneLinkageConfig
 * @brief Configuration for backbone linkage detection
 */
struct BackboneLinkageConfig {
    double o3p_upper = 2.5;  ///< Maximum O3'-P distance for linkage (Angstroms)
};

/**
 * @class BackboneLinkageChecker
 * @brief Detects backbone connectivity between residues and base pairs
 *
 * This class handles the O3'-P distance-based linkage detection used by
 * the five2three algorithm to determine strand direction.
 */
class BackboneLinkageChecker {
public:
    explicit BackboneLinkageChecker(const BackboneLinkageConfig& config = {})
        : config_(config) {}

    /**
     * @brief Check linkage direction between two residues
     * @param res_i First residue index (1-based)
     * @param res_j Second residue index (1-based)
     * @param backbone Backbone atom data
     * @return LinkDirection indicating the type of linkage
     */
    [[nodiscard]] LinkDirection check_linkage(
        size_t res_i, size_t res_j,
        const BackboneData& backbone) const;

    /**
     * @brief Calculate O3'-O3' distance between residues
     * @param res_i First residue index (1-based)
     * @param res_j Second residue index (1-based)
     * @param backbone Backbone atom data
     * @return Distance in Angstroms, or -1.0 if atoms not found
     */
    [[nodiscard]] double o3_distance(
        size_t res_i, size_t res_j,
        const BackboneData& backbone) const;

    /**
     * @brief Check if two base pairs are connected via backbone
     * @param pair1 First base pair
     * @param pair2 Second base pair
     * @param backbone Backbone atom data
     * @return true if any backbone linkage exists between the pairs
     */
    [[nodiscard]] bool are_pairs_connected(
        const core::BasePair& pair1,
        const core::BasePair& pair2,
        const BackboneData& backbone) const;

    [[nodiscard]] const BackboneLinkageConfig& config() const { return config_; }

private:
    BackboneLinkageConfig config_;
};

/**
 * @brief Update direction count based on linkage type
 * @param link The linkage direction
 * @param forward Forward count to increment
 * @param reverse Reverse count to increment
 * @param none None count to increment
 */
inline void update_direction_count(LinkDirection link, int& forward, int& reverse, int& none) {
    if (link == LinkDirection::Forward) ++forward;
    else if (link == LinkDirection::Reverse) ++reverse;
    else ++none;
}

} // namespace x3dna::algorithms::helix

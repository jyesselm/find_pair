/**
 * @file strand_direction_checker.hpp
 * @brief Strand direction checking for five2three algorithm
 *
 * Extracts the five2three sub-functions from helix_organizer for
 * checking and correcting strand direction in base pair steps.
 */

#pragma once

#include <x3dna/core/base_pair.hpp>
#include <x3dna/algorithms/helix_organizer.hpp>  // For types
#include <x3dna/algorithms/helix/backbone_linkage_checker.hpp>
#include <x3dna/algorithms/helix/pair_geometry_helper.hpp>
#include <vector>

namespace x3dna::algorithms::helix {

// Re-use types from helix_organizer
using algorithms::BackboneData;
using algorithms::DirectionCounts;
using algorithms::HelixSegment;
using algorithms::StrandResidues;

/**
 * @struct StrandDirectionConfig
 * @brief Configuration for strand direction checking
 */
struct StrandDirectionConfig {
    double end_stack_xang = 125.0;  ///< Max x-angle for stacked WC pairs (degrees)
    double o3p_upper = 2.5;         ///< Max O3'-P distance for backbone linkage
};

/**
 * @class StrandDirectionChecker
 * @brief Checks and corrects strand direction in base pair steps
 *
 * Implements the five2three sub-functions for ensuring proper 5'→3'
 * strand direction in helix organization.
 */
class StrandDirectionChecker {
public:
    explicit StrandDirectionChecker(const StrandDirectionConfig& config = {})
        : config_(config), linkage_checker_({config.o3p_upper}) {}

    /**
     * @brief Set initial strand assignment for first pair in helix
     */
    void first_step(const std::vector<core::BasePair>& pairs,
                   const BackboneData& backbone,
                   std::vector<size_t>& pair_order,
                   const HelixSegment& helix,
                   std::vector<bool>& swapped) const;

    /**
     * @brief Check Watson-Crick base pair z-direction alignment
     * @return true if pair_n should be swapped
     */
    [[nodiscard]] bool wc_bporien(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                   bool swap_m, bool swap_n,
                                   const BackboneData& backbone) const;

    /**
     * @brief Check O3' distance patterns for swap indication
     */
    [[nodiscard]] bool check_o3dist(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                     bool swap_m, bool swap_n,
                                     const BackboneData& backbone) const;

    /**
     * @brief Check strand chain connectivity for swap indication
     */
    [[nodiscard]] bool check_schain(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                     bool swap_m, bool swap_n,
                                     const BackboneData& backbone) const;

    /**
     * @brief Check frame orientation alignment for swap indication
     */
    [[nodiscard]] bool check_others(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                     bool swap_m, bool swap_n,
                                     const BackboneData& backbone) const;

    /**
     * @brief Check if strand 1 direction is reversed
     */
    [[nodiscard]] bool chain1dir(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                  bool swap_m, bool swap_n,
                                  const BackboneData& backbone) const;

    /**
     * @brief Count backbone linkage directions and apply fixes
     * @note May modify swapped and pair_order
     */
    DirectionCounts check_direction(
        const std::vector<core::BasePair>& pairs,
        const BackboneData& backbone,
        std::vector<size_t>& pair_order,
        HelixSegment& helix,
        std::vector<bool>& swapped) const;

    /**
     * @brief Additional strand corrections based on direction counts
     */
    void check_strand2(const std::vector<core::BasePair>& pairs,
                      const BackboneData& backbone,
                      const std::vector<size_t>& pair_order,
                      HelixSegment& helix,
                      std::vector<bool>& swapped,
                      const DirectionCounts& direction) const;

    [[nodiscard]] const StrandDirectionConfig& config() const { return config_; }

private:
    StrandDirectionConfig config_;
    BackboneLinkageChecker linkage_checker_;

    /**
     * @brief Calculate angle between combined x-axes of two pairs
     */
    [[nodiscard]] double wcbp_xang(const core::BasePair& pair_m, const core::BasePair& pair_n) const;

    /**
     * @brief Calculate z-direction dot product for WC pairs
     */
    [[nodiscard]] double wcbp_zdir(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                    bool swap_m, bool swap_n) const;

    /**
     * @brief Check if pair has positive bpid (WC-like geometry)
     */
    [[nodiscard]] bool has_positive_bpid(const core::BasePair& pair) const;
};

} // namespace x3dna::algorithms::helix

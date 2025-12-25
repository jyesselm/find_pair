/**
 * @file slot_optimizer_params.hpp
 * @brief Parameters for slot-based H-bond optimization
 */

#pragma once

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {
namespace slot {

/**
 * @struct SlotOptimizerParams
 * @brief Configuration parameters for SlotOptimizer
 */
struct SlotOptimizerParams {
    // Distance thresholds
    double max_distance = 4.0;             ///< Maximum donor-acceptor distance (Å)
    double short_distance_threshold = 3.5; ///< Below this, skip alignment check

    // Alignment thresholds
    double min_alignment = 0.3;            ///< Minimum alignment score (0-2) for regular bonds
    double min_bifurcation_alignment = 0.5;///< Stricter alignment for bifurcated bonds

    // Bifurcation
    double min_bifurcation_angle = 43.0;   ///< Minimum angle between bonds sharing slot (degrees)

    // Baseline mode (legacy-compatible)
    bool baseline_mode = false;            ///< Use simpler distance-based selection
    double baseline_min_distance = 2.5;    ///< Minimum distance in baseline mode
    double baseline_max_distance = 3.5;    ///< Maximum distance in baseline mode

    /**
     * @brief Create default optimized mode parameters
     */
    [[nodiscard]] static SlotOptimizerParams optimized() {
        return SlotOptimizerParams{};
    }

    /**
     * @brief Create baseline (legacy-compatible) mode parameters
     */
    [[nodiscard]] static SlotOptimizerParams baseline() {
        SlotOptimizerParams params;
        params.baseline_mode = true;
        return params;
    }

    /**
     * @brief Create strict mode (higher alignment requirements)
     */
    [[nodiscard]] static SlotOptimizerParams strict() {
        SlotOptimizerParams params;
        params.min_alignment = 0.5;
        params.min_bifurcation_alignment = 0.7;
        return params;
    }
};

} // namespace slot
} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna

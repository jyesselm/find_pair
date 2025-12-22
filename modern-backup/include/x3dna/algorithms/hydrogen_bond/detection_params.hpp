/**
 * @file detection_params.hpp
 * @brief Parameters for hydrogen bond detection algorithm
 */

#pragma once

#include <string>
#include <x3dna/core/hbond_types.hpp>

namespace x3dna {
namespace algorithms {

/**
 * @brief Context-specific distance thresholds for H-bond detection
 */
struct HBondDistanceThresholds {
    double base_base_max = 4.0;         // Base-base H-bonds
    double base_backbone_max = 3.5;     // Base to backbone
    double backbone_backbone_max = 3.5; // Backbone to backbone
    double base_sugar_max = 3.5;        // Base to sugar
    double sugar_sugar_max = 3.5;       // Sugar to sugar

    double min_distance = 1.8;             // Minimum for all contexts
    double conflict_filter_distance = 4.5; // For conflict resolution phase 3

    [[nodiscard]] double max_for_context(core::HBondContext ctx) const;
};

/**
 * @brief Parameters for H-bond detection algorithm
 */
struct HBondDetectionParams {
    // Distance thresholds (context-aware)
    HBondDistanceThresholds distances;

    // Element filter (which elements can form H-bonds)
    std::string allowed_elements = ".O.N";

    // Quality thresholds for "good" H-bond determination
    double good_bond_min_distance = 2.5;
    double good_bond_max_distance = 3.5;

    // Post-validation filtering
    double post_validation_max_distance = 3.6;
    double nonstandard_min_distance = 2.6;
    double nonstandard_max_distance = 3.2;

    // Minimum H-bonds required for valid base pair
    int min_base_hbonds_required = 1;

    [[nodiscard]] static HBondDetectionParams legacy_compatible();
    [[nodiscard]] static HBondDetectionParams modern();
};

} // namespace algorithms
} // namespace x3dna

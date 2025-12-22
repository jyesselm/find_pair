/**
 * @file detection_params.cpp
 * @brief Implementation of H-bond detection parameters
 */

#include <x3dna/algorithms/hydrogen_bond/detection_params.hpp>

namespace x3dna {
namespace algorithms {

double HBondDistanceThresholds::max_for_context(core::HBondContext ctx) const {
    switch (ctx) {
        case core::HBondContext::BASE_BASE:
            return base_base_max;
        case core::HBondContext::BASE_BACKBONE:
            return base_backbone_max;
        case core::HBondContext::BACKBONE_BACKBONE:
            return backbone_backbone_max;
        case core::HBondContext::BASE_SUGAR:
            return base_sugar_max;
        case core::HBondContext::SUGAR_SUGAR:
            return sugar_sugar_max;
        case core::HBondContext::UNKNOWN:
        default:
            return base_base_max; // Default to most permissive
    }
}

HBondDetectionParams HBondDetectionParams::legacy_compatible() {
    HBondDetectionParams params;

    // Legacy used 4.0 for all contexts
    params.distances.base_base_max = 4.0;
    params.distances.base_backbone_max = 4.0;
    params.distances.backbone_backbone_max = 4.0;
    params.distances.base_sugar_max = 4.0;
    params.distances.sugar_sugar_max = 4.0;

    params.distances.min_distance = 1.8;
    params.distances.conflict_filter_distance = 4.5;

    params.allowed_elements = ".O.N";

    params.good_bond_min_distance = 2.5;
    params.good_bond_max_distance = 3.5;

    params.post_validation_max_distance = 3.6;
    params.nonstandard_min_distance = 2.6;
    params.nonstandard_max_distance = 3.2;

    params.min_base_hbonds_required = 1;

    return params;
}

HBondDetectionParams HBondDetectionParams::modern() {
    HBondDetectionParams params;

    // Modern uses context-specific thresholds
    params.distances.base_base_max = 4.0;
    params.distances.base_backbone_max = 3.5;
    params.distances.backbone_backbone_max = 3.5;
    params.distances.base_sugar_max = 3.5;
    params.distances.sugar_sugar_max = 3.5;

    params.distances.min_distance = 1.8;
    params.distances.conflict_filter_distance = 4.5;

    params.allowed_elements = ".O.N";

    params.good_bond_min_distance = 2.5;
    params.good_bond_max_distance = 3.5;

    params.post_validation_max_distance = 3.6;
    params.nonstandard_min_distance = 2.6;
    params.nonstandard_max_distance = 3.2;

    params.min_base_hbonds_required = 1;

    return params;
}

} // namespace algorithms
} // namespace x3dna

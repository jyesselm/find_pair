/**
 * @file hbond_parameters.hpp
 * @brief Unified H-bond parameter definitions with nested structures
 *
 * This file defines all H-bond related parameters in a single, organized structure.
 * Parameters are loaded from resources/config/hbond_parameters.json
 */

#pragma once

#include <string>

namespace x3dna {
namespace config {

// ============================================================================
// Range structure for min/max pairs
// ============================================================================

struct Range {
    double min = 0.0;
    double max = 0.0;
};

// ============================================================================
// Detection Parameters
// ============================================================================

/**
 * @brief Context-specific distance thresholds for H-bond detection
 */
struct HBondDistanceConfig {
    double min = 2.0;  // Must match hbond_parameters.json
    double base_base_max = 4.0;
    double base_backbone_max = 3.5;
    double backbone_backbone_max = 3.5;
    double base_sugar_max = 3.5;
    double sugar_sugar_max = 3.5;
    double protein_mainchain_max = 3.5;
    double protein_sidechain_max = 3.5;
    double base_protein_max = 3.5;
    double protein_ligand_max = 3.5;
    double base_ligand_max = 3.5;
    double conflict_filter = 4.5;
};

/**
 * @brief Element filter for H-bond donors/acceptors
 */
struct HBondElementConfig {
    std::string allowed = ".O.N.";
};

/**
 * @brief Distance thresholds for bond quality classification
 */
struct HBondThresholdsConfig {
    Range good_bond{2.5, 3.5};
    double post_validation_max = 3.6;
    Range nonstandard{2.6, 3.2};
};

/**
 * @brief Validation requirements
 */
struct HBondValidationConfig {
    int min_base_hbonds = 1;
};

/**
 * @brief Optional detection features
 */
struct HBondOptionsConfig {
    bool enable_angle_filtering = false;
    bool enable_quality_scoring = false;
    bool filter_invalid_scores = false;
    bool include_unlikely_chemistry = false;
    bool include_backbone_backbone = false;
    bool include_intra_residue = false;  // Detect H-bonds within same residue
};

/**
 * @brief All detection-related parameters
 */
struct HBondDetectionConfig {
    HBondDistanceConfig distance;
    HBondElementConfig elements;
    HBondThresholdsConfig thresholds;
    HBondValidationConfig validation;
    HBondOptionsConfig options;
};

// ============================================================================
// Geometry Parameters
// ============================================================================

/**
 * @brief Donor angle thresholds
 */
struct DonorAngleConfig {
    double min = 90.0;
    double ideal = 165.0;
};

/**
 * @brief Acceptor angle thresholds (different ideals for sp2/sp3)
 */
struct AcceptorAngleConfig {
    double min = 70.0;
    double ideal_sp2 = 130.0;
    double ideal_sp3 = 110.0;
};

/**
 * @brief All geometry-related parameters
 */
struct HBondGeometryConfig {
    DonorAngleConfig donor_angle;
    AcceptorAngleConfig acceptor_angle;
};

// ============================================================================
// Scoring Parameters
// ============================================================================

/**
 * @brief Distance scoring parameters (Gaussian)
 */
struct ScoringDistanceConfig {
    double ideal = 2.9;
    double sigma = 0.3;
    double min = 2.0;
    double max = 4.0;
};

/**
 * @brief Component weights for quality scoring
 */
struct ScoringWeightsConfig {
    double distance = 0.45;
    double donor_angle = 0.30;
    double acceptor_angle = 0.25;
};

/**
 * @brief Resolution-based score adjustment
 */
struct ScoringResolutionConfig {
    bool apply_penalty = true;
    double high_res_threshold = 2.0;
    double low_res_threshold = 3.5;
};

/**
 * @brief All scoring-related parameters
 */
struct HBondScoringConfig {
    ScoringDistanceConfig distance;
    ScoringWeightsConfig weights;
    ScoringResolutionConfig resolution;
};

// ============================================================================
// Quality Tier Parameters
// ============================================================================

/**
 * @brief Score thresholds for quality tier classification
 */
struct QualityTiersConfig {
    double excellent_min = 90.0;
    double standard_min = 70.0;
    double acceptable_min = 50.0;
    double questionable_min = 30.0;
};

// ============================================================================
// Top-Level Container
// ============================================================================

/**
 * @brief Complete H-bond parameter configuration
 *
 * Contains all parameters for H-bond detection, geometry validation,
 * quality scoring, and tier classification.
 *
 * Example usage:
 * @code
 *   auto params = HBondParametersLoader::load();
 *   double max_dist = params.detection.distance.base_base_max;
 *   bool is_excellent = score >= params.quality_tiers.excellent_min;
 * @endcode
 */
struct HBondParameters {
    HBondDetectionConfig detection;
    HBondGeometryConfig geometry;
    HBondScoringConfig scoring;
    QualityTiersConfig quality_tiers;

    /**
     * @brief Get default parameters (matches JSON defaults)
     */
    static HBondParameters defaults();
};

} // namespace config
} // namespace x3dna

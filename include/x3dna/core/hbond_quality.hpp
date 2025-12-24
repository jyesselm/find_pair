/**
 * @file hbond_quality.hpp
 * @brief H-bond quality scoring types and enums
 *
 * Provides quality assessment for hydrogen bonds based on geometric criteria:
 * - Distance from ideal (Gaussian scoring)
 * - Donor angle (D-H...A linearity)
 * - Acceptor angle (D...A-Y hybridization)
 *
 * Quality tiers map to DSSR's donAcc_type classification:
 * - EXCELLENT/STANDARD -> DSSR "standard"
 * - ACCEPTABLE -> DSSR "acceptable"
 * - QUESTIONABLE -> DSSR "questionable"
 */

#pragma once

#include <string>

namespace x3dna {
namespace core {

/**
 * @enum HBondQualityTier
 * @brief Quality tier for hydrogen bonds (extends DSSR's 3-tier system)
 *
 * Tiers are based on combined geometric score (0-100):
 * - EXCELLENT: 90-100 - Ideal geometry, unambiguous
 * - STANDARD: 70-89 - Good geometry, reliable
 * - ACCEPTABLE: 50-69 - Minor deviations, usable
 * - QUESTIONABLE: 30-49 - Marginal, use with caution
 * - INVALID: 0-29 - Clearly wrong geometry, should be filtered
 */
enum class HBondQualityTier {
    EXCELLENT,     // 90-100: Best geometry
    STANDARD,      // 70-89: Good geometry, reliable
    ACCEPTABLE,    // 50-69: Minor deviations
    QUESTIONABLE,  // 30-49: Marginal quality
    INVALID        // 0-29: Clearly wrong, filter out
};

/**
 * @brief Convert quality tier to string
 */
[[nodiscard]] inline const char* to_string(HBondQualityTier tier) {
    switch (tier) {
        case HBondQualityTier::EXCELLENT:
            return "EXCELLENT";
        case HBondQualityTier::STANDARD:
            return "STANDARD";
        case HBondQualityTier::ACCEPTABLE:
            return "ACCEPTABLE";
        case HBondQualityTier::QUESTIONABLE:
            return "QUESTIONABLE";
        case HBondQualityTier::INVALID:
            return "INVALID";
    }
    return "UNKNOWN";
}

/**
 * @brief Convert quality tier to DSSR-compatible string
 *
 * Maps our 5-tier system to DSSR's 3-tier donAcc_type:
 * - EXCELLENT/STANDARD -> "standard"
 * - ACCEPTABLE -> "acceptable"
 * - QUESTIONABLE/INVALID -> "questionable"
 */
[[nodiscard]] inline const char* to_dssr_string(HBondQualityTier tier) {
    switch (tier) {
        case HBondQualityTier::EXCELLENT:
        case HBondQualityTier::STANDARD:
            return "standard";
        case HBondQualityTier::ACCEPTABLE:
            return "acceptable";
        case HBondQualityTier::QUESTIONABLE:
        case HBondQualityTier::INVALID:
            return "questionable";
    }
    return "unknown";
}

/**
 * @struct HBondQualityScore
 * @brief Complete quality score for a hydrogen bond
 *
 * Contains:
 * - Total weighted score (0-100)
 * - Component scores for distance, donor angle, acceptor angle
 * - Quality tier classification
 * - Failure reason if invalid
 *
 * Scoring uses: 45% distance + 30% donor angle + 25% acceptor angle
 * (Dihedral not included - matching HBPLUS/DSSP/Chimera approach)
 */
struct HBondQualityScore {
    double total_score = 0.0;           // Combined weighted score (0-100)
    double distance_score = 0.0;        // Distance component (0-100)
    double donor_angle_score = 0.0;     // Donor angle component (0-100)
    double acceptor_angle_score = 0.0;  // Acceptor angle component (0-100)
    HBondQualityTier tier = HBondQualityTier::INVALID;
    std::string failure_reason;         // Empty if valid, explains rejection otherwise

    /**
     * @brief Check if this H-bond passes quality threshold
     * @return true if tier is ACCEPTABLE or better
     */
    [[nodiscard]] bool is_acceptable() const {
        return tier == HBondQualityTier::EXCELLENT ||
               tier == HBondQualityTier::STANDARD ||
               tier == HBondQualityTier::ACCEPTABLE;
    }

    /**
     * @brief Check if this H-bond is high quality
     * @return true if tier is STANDARD or EXCELLENT
     */
    [[nodiscard]] bool is_high_quality() const {
        return tier == HBondQualityTier::EXCELLENT ||
               tier == HBondQualityTier::STANDARD;
    }

    /**
     * @brief Check if this H-bond should be filtered out
     * @return true if tier is INVALID
     */
    [[nodiscard]] bool should_filter() const {
        return tier == HBondQualityTier::INVALID;
    }
};

/**
 * @struct QualityTierThresholds
 * @brief Configurable thresholds for quality tier classification
 */
struct QualityTierThresholds {
    double excellent_min = 90.0;
    double standard_min = 70.0;
    double acceptable_min = 50.0;
    double questionable_min = 30.0;

    /**
     * @brief Get default tier thresholds
     */
    static QualityTierThresholds defaults() {
        return QualityTierThresholds{};
    }
};

/**
 * @brief Convert numeric score to quality tier using default thresholds
 * @param score Score in range [0, 100]
 * @return Corresponding quality tier
 */
[[nodiscard]] inline HBondQualityTier score_to_tier(double score) {
    if (score >= 90.0) return HBondQualityTier::EXCELLENT;
    if (score >= 70.0) return HBondQualityTier::STANDARD;
    if (score >= 50.0) return HBondQualityTier::ACCEPTABLE;
    if (score >= 30.0) return HBondQualityTier::QUESTIONABLE;
    return HBondQualityTier::INVALID;
}

/**
 * @brief Convert numeric score to quality tier using custom thresholds
 * @param score Score in range [0, 100]
 * @param thresholds Custom tier thresholds
 * @return Corresponding quality tier
 */
[[nodiscard]] inline HBondQualityTier score_to_tier(double score, const QualityTierThresholds& thresholds) {
    if (score >= thresholds.excellent_min) return HBondQualityTier::EXCELLENT;
    if (score >= thresholds.standard_min) return HBondQualityTier::STANDARD;
    if (score >= thresholds.acceptable_min) return HBondQualityTier::ACCEPTABLE;
    if (score >= thresholds.questionable_min) return HBondQualityTier::QUESTIONABLE;
    return HBondQualityTier::INVALID;
}

} // namespace core
} // namespace x3dna

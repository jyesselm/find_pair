/**
 * @file quality_score_calculator.hpp
 * @brief Calculates adjusted quality scores for base pair selection
 */

#pragma once

#include <x3dna/core/residue.hpp>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/parameter_calculator.hpp>
#include <vector>
#include <string>

namespace x3dna {
namespace algorithms {

/**
 * @class QualityScoreCalculator
 * @brief Calculates adjusted quality scores for base pair selection
 *
 * Encapsulates the quality score adjustment logic from legacy code:
 * - adjust_pairQuality: Adjusts score based on H-bond quality
 * - calculate_bp_type_id: Determines Watson-Crick/Wobble pair type
 *
 * The adjusted quality score is used during pair selection to determine
 * the best partner for each residue.
 */
/**
 * @brief Named constants for quality score calculations
 *
 * These constants match the legacy X3DNA code behavior.
 * Extracting them to named constants improves code readability
 * and makes the algorithm parameters explicit.
 */
namespace quality_constants {
    // Quality score adjustments
    constexpr double WC_PAIR_BONUS = -2.0;          // Bonus for Watson-Crick pairs
    constexpr double GOOD_HBOND_ADJUSTMENT = -3.0;  // Adjustment for >= 2 good H-bonds

    // H-bond distance range for "good" H-bonds (in Angstroms)
    constexpr double GOOD_HBOND_MIN_DIST = 2.5;
    constexpr double GOOD_HBOND_MAX_DIST = 3.5;
    constexpr int MIN_GOOD_HBONDS_FOR_BONUS = 2;

    // Shear thresholds for pair type classification
    constexpr double WOBBLE_SHEAR_MIN = 1.8;        // Minimum shear for wobble pair
    constexpr double WOBBLE_SHEAR_MAX = 2.8;        // Maximum shear for wobble pair
    constexpr double WC_SHEAR_MAX = 1.8;            // Maximum shear for Watson-Crick pair

    // Parameter thresholds for bp_type_id calculation
    constexpr double STRETCH_THRESHOLD = 2.0;       // Maximum stretch for valid bp_type
    constexpr double OPENING_THRESHOLD = 60.0;      // Maximum opening angle (degrees)
}  // namespace quality_constants

class QualityScoreCalculator {
public:
    /**
     * @brief Calculate adjusted quality score for pair selection
     *
     * This applies:
     * 1. H-bond quality adjustment (adjust_pairQuality)
     * 2. bp_type_id == 2 bonus (-2.0 for Watson-Crick pairs)
     *
     * @param result Validation result with raw quality score
     * @param res1 First residue (for bp_type calculation)
     * @param res2 Second residue (for bp_type calculation)
     * @return Adjusted score (lower is better)
     */
    [[nodiscard]] double calculate_selection_score(const ValidationResult& result, const core::Residue& res1,
                                     const core::Residue& res2) const;

    /**
     * @brief Calculate H-bond quality adjustment (matches legacy adjust_pairQuality)
     *
     * Counts "good" hydrogen bonds (type '-' with distance in [2.5, 3.5] Å)
     * and returns adjustment:
     * - >= 2 good H-bonds: -3.0
     * - 1 good H-bond: -1.0
     * - 0 good H-bonds: 0.0
     *
     * @param hbonds Vector of hydrogen bonds
     * @return Adjustment value (negative, to be added to quality_score)
     */
    [[nodiscard]] double adjust_pair_quality(const std::vector<core::hydrogen_bond>& hbonds) const;

    /**
     * @brief Calculate bp_type_id (matches legacy check_wc_wobble_pair)
     *
     * Determines base pair type based on direction vectors and step parameters:
     * - -1: Unknown (direction check failed or thresholds exceeded)
     * - 0: Invalid pair
     * - 1: Wobble pair (shear in [1.8, 2.8])
     * - 2: Watson-Crick pair (shear <= 1.8 and pair in WC_LIST)
     *
     * @param res1 First residue
     * @param res2 Second residue
     * @param result Validation result (for direction vectors and validity)
     * @return bp_type_id (-1, 0, 1, or 2)
     */
    [[nodiscard]] int calculate_bp_type_id(const core::Residue& res1, const core::Residue& res2,
                             const ValidationResult& result) const;

private:
    ParameterCalculator param_calculator_;

    // Watson-Crick pair list
    static const std::vector<std::string> WC_LIST;
};

} // namespace algorithms
} // namespace x3dna

/**
 * @file quality_scorer.hpp
 * @brief H-bond quality scoring based on geometric criteria
 *
 * Scores hydrogen bonds using:
 * - Distance (Gaussian centered on ideal ~2.9Å)
 * - Donor angle (linear penalty from ideal 165°)
 * - Acceptor angle (depends on sp2/sp3 hybridization)
 *
 * Weights: 45% distance + 30% donor angle + 25% acceptor angle
 * (No dihedral - matching HBPLUS/DSSP/Chimera industry standard)
 */

#pragma once

#include <vector>
#include <x3dna/algorithms/hydrogen_bond/hbond.hpp>
#include <x3dna/algorithms/hydrogen_bond/hbond_quality.hpp>
#include <x3dna/algorithms/hydrogen_bond/hbond_types.hpp>

// Forward declaration
namespace x3dna {
namespace config {
struct HBondParameters;
}
}

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

/**
 * @struct HBondScoringParams
 * @brief Parameters for H-bond quality scoring
 */
struct HBondScoringParams {
    // Distance scoring (Gaussian)
    double ideal_distance = 2.9;       // Ideal D-A distance in Angstroms
    double distance_sigma = 0.3;       // Gaussian width

    // Distance hard cutoffs
    double min_distance = 2.0;         // Below = steric clash
    double max_distance = 4.0;         // Above = too far

    // Donor angle scoring (X-D...A)
    double ideal_donor_angle = 165.0;  // Ideal angle in degrees
    double min_donor_angle = 90.0;     // Below = impossible geometry

    // Acceptor angle scoring (D...A-Y)
    double ideal_acceptor_sp2 = 130.0; // For carbonyl O, aromatic N
    double ideal_acceptor_sp3 = 110.0; // For hydroxyl O, amino N
    double min_acceptor_angle = 70.0;  // Below = impossible geometry

    // Scoring weights (must sum to 1.0)
    double weight_distance = 0.45;
    double weight_donor_angle = 0.30;
    double weight_acceptor_angle = 0.25;

    // Resolution adjustment (optional)
    bool apply_resolution_penalty = true;
    double high_res_threshold = 2.0;   // No penalty below this
    double low_res_threshold = 3.5;    // Maximum penalty above this

    // Presets (load from config file)
    [[nodiscard]] static HBondScoringParams defaults();
    [[nodiscard]] static HBondScoringParams strict();
    [[nodiscard]] static HBondScoringParams lenient();

    /**
     * @brief Create scoring params from unified config
     * @param config Loaded HBondParameters configuration
     * @return Scoring params populated from config
     */
    [[nodiscard]] static HBondScoringParams from_config(const config::HBondParameters& config);
};

/**
 * @class HBondQualityScorer
 * @brief Scores hydrogen bonds based on geometric criteria
 *
 * This scorer is INFORMATIONAL ONLY - it does not affect pair detection.
 * Scores are assigned after H-bond detection to provide quality metrics.
 *
 * Usage:
 *   HBondQualityScorer scorer;
 *   scorer.score_all(hbonds);  // Scores in-place
 *   // or
 *   auto score = scorer.score(hbond);
 */
class HBondQualityScorer {
public:
    explicit HBondQualityScorer(const HBondScoringParams& params = HBondScoringParams::defaults());

    /**
     * @brief Score a single H-bond
     * @param hbond H-bond with geometry already calculated
     * @return Quality score with tier and component scores
     */
    [[nodiscard]] core::HBondQualityScore score(const core::HBond& hbond) const;

    /**
     * @brief Score all H-bonds in a vector (in-place)
     * @param hbonds Vector of H-bonds to score
     */
    void score_all(std::vector<core::HBond>& hbonds) const;

    /**
     * @brief Score with optional resolution adjustment
     * @param hbond H-bond to score
     * @param resolution Structure resolution in Angstroms (0 = unknown)
     * @return Quality score with resolution-aware adjustment
     */
    [[nodiscard]] core::HBondQualityScore score_with_resolution(
        const core::HBond& hbond, double resolution) const;

    /**
     * @brief Get the scoring parameters
     */
    [[nodiscard]] const HBondScoringParams& params() const { return params_; }

private:
    HBondScoringParams params_;

    /**
     * @brief Score distance component (Gaussian centered on ideal)
     * @param distance D-A distance in Angstroms
     * @return Score 0-100
     */
    [[nodiscard]] double score_distance(double distance) const;

    /**
     * @brief Score donor angle component (X-D...A)
     * @param angle Donor angle in degrees
     * @param failure_reason Output: set if hard failure
     * @return Score 0-100
     */
    [[nodiscard]] double score_donor_angle(double angle, std::string& failure_reason) const;

    /**
     * @brief Score acceptor angle component (D...A-Y)
     * @param angle Acceptor angle in degrees
     * @param acceptor_atom Acceptor atom name (for hybridization detection)
     * @param failure_reason Output: set if hard failure
     * @return Score 0-100
     */
    [[nodiscard]] double score_acceptor_angle(double angle, const std::string& acceptor_atom,
                                               std::string& failure_reason) const;

    /**
     * @brief Check for hard failures (impossible geometry)
     * @param hbond H-bond to check
     * @param failure_reason Output: reason for failure
     * @return true if hard failure detected
     */
    [[nodiscard]] bool check_hard_failures(const core::HBond& hbond, std::string& failure_reason) const;

    /**
     * @brief Apply soft failure caps (marginal geometry)
     * @param hbond H-bond to check
     * @param score Current score to cap
     * @return Capped score (max 40 if soft failure)
     */
    [[nodiscard]] double apply_soft_failure_caps(const core::HBond& hbond, double score) const;

    /**
     * @brief Determine if acceptor atom is sp2 hybridized
     * @param atom_name Acceptor atom name (trimmed)
     * @return true if sp2 (carbonyl O, aromatic N)
     */
    [[nodiscard]] static bool is_sp2_acceptor(const std::string& atom_name);

    /**
     * @brief Adjust score based on resolution
     * @param score Raw score
     * @param resolution Structure resolution in Angstroms
     * @return Adjusted score (penalized for low resolution)
     */
    [[nodiscard]] double adjust_for_resolution(double score, double resolution) const;
};

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna

/**
 * @file quality_scorer.cpp
 * @brief Implementation of H-bond quality scoring
 */

#include <x3dna/algorithms/hydrogen_bond/quality_scorer.hpp>
#include <cmath>
#include <algorithm>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

// === HBondScoringParams presets ===

HBondScoringParams HBondScoringParams::defaults() {
    return HBondScoringParams{};
}

HBondScoringParams HBondScoringParams::strict() {
    HBondScoringParams params;
    params.ideal_distance = 2.85;
    params.distance_sigma = 0.25;
    params.max_distance = 3.5;
    params.ideal_donor_angle = 170.0;
    params.min_donor_angle = 110.0;
    params.min_acceptor_angle = 90.0;
    return params;
}

HBondScoringParams HBondScoringParams::lenient() {
    HBondScoringParams params;
    params.ideal_distance = 2.9;
    params.distance_sigma = 0.4;
    params.max_distance = 4.0;
    params.ideal_donor_angle = 160.0;
    params.min_donor_angle = 80.0;
    params.min_acceptor_angle = 60.0;
    return params;
}

// === HBondQualityScorer implementation ===

HBondQualityScorer::HBondQualityScorer(const HBondScoringParams& params)
    : params_(params) {}

core::HBondQualityScore HBondQualityScorer::score(const core::HBond& hbond) const {
    core::HBondQualityScore result;

    // Check for hard failures first
    if (check_hard_failures(hbond, result.failure_reason)) {
        result.total_score = 0.0;
        result.tier = core::HBondQualityTier::INVALID;
        return result;
    }

    // Calculate component scores
    std::string donor_failure, acceptor_failure;
    result.distance_score = score_distance(hbond.distance);
    result.donor_angle_score = score_donor_angle(hbond.donor_angle, donor_failure);
    result.acceptor_angle_score = score_acceptor_angle(
        hbond.acceptor_angle, hbond.acceptor_atom_name, acceptor_failure);

    // Check for angle-based hard failures from scoring
    if (!donor_failure.empty()) {
        result.failure_reason = donor_failure;
        result.total_score = 0.0;
        result.tier = core::HBondQualityTier::INVALID;
        return result;
    }
    if (!acceptor_failure.empty()) {
        result.failure_reason = acceptor_failure;
        result.total_score = 0.0;
        result.tier = core::HBondQualityTier::INVALID;
        return result;
    }

    // Calculate weighted total
    result.total_score = params_.weight_distance * result.distance_score +
                         params_.weight_donor_angle * result.donor_angle_score +
                         params_.weight_acceptor_angle * result.acceptor_angle_score;

    // Apply soft failure caps
    result.total_score = apply_soft_failure_caps(hbond, result.total_score);

    // Determine tier
    result.tier = core::score_to_tier(result.total_score);

    return result;
}

void HBondQualityScorer::score_all(std::vector<core::HBond>& hbonds) const {
    for (auto& hbond : hbonds) {
        hbond.quality_score = score(hbond);
    }
}

core::HBondQualityScore HBondQualityScorer::score_with_resolution(
    const core::HBond& hbond, double resolution) const {

    auto result = score(hbond);

    // Apply resolution adjustment if enabled and not already invalid
    if (params_.apply_resolution_penalty && result.tier != core::HBondQualityTier::INVALID) {
        result.total_score = adjust_for_resolution(result.total_score, resolution);
        result.tier = core::score_to_tier(result.total_score);
    }

    return result;
}

double HBondQualityScorer::score_distance(double distance) const {
    // Hard cutoffs
    if (distance < params_.min_distance || distance > params_.max_distance) {
        return 0.0;
    }

    // Gaussian scoring centered on ideal
    double diff = distance - params_.ideal_distance;
    double score = 100.0 * std::exp(-0.5 * (diff / params_.distance_sigma) * (diff / params_.distance_sigma));

    return std::clamp(score, 0.0, 100.0);
}

double HBondQualityScorer::score_donor_angle(double angle, std::string& failure_reason) const {
    // Check if angle was not calculated (0.0 means no neighbor found)
    if (angle <= 0.0) {
        // No neighbor found - give neutral score
        return 60.0;
    }

    // Hard cutoff - impossible geometry
    if (angle < params_.min_donor_angle) {
        failure_reason = "Donor angle < " + std::to_string(static_cast<int>(params_.min_donor_angle)) + "° (impossible geometry)";
        return 0.0;
    }

    // Linear penalty from ideal
    // For theta >= 120°: score = max(0, 100 - 2.5 * |theta - 165|)
    if (angle >= 120.0) {
        double deviation = std::abs(angle - params_.ideal_donor_angle);
        double score = 100.0 - 2.5 * deviation;
        return std::clamp(score, 0.0, 100.0);
    }

    // Between min_donor_angle and 120°: reduced score
    // Highly strained geometry
    double fraction = (angle - params_.min_donor_angle) / (120.0 - params_.min_donor_angle);
    return 40.0 * fraction;  // Max 40 for strained angles
}

double HBondQualityScorer::score_acceptor_angle(double angle, const std::string& acceptor_atom,
                                                  std::string& failure_reason) const {
    // Check if angle was not calculated
    if (angle <= 0.0) {
        return 60.0;  // Neutral score
    }

    // Hard cutoff
    if (angle < params_.min_acceptor_angle) {
        failure_reason = "Acceptor angle < " + std::to_string(static_cast<int>(params_.min_acceptor_angle)) + "° (impossible geometry)";
        return 0.0;
    }

    // Determine ideal based on hybridization
    double ideal = is_sp2_acceptor(acceptor_atom) ? params_.ideal_acceptor_sp2 : params_.ideal_acceptor_sp3;

    // Linear penalty from ideal
    double deviation = std::abs(angle - ideal);
    double score = 100.0 - 2.0 * deviation;

    return std::clamp(score, 0.0, 100.0);
}

bool HBondQualityScorer::check_hard_failures(const core::HBond& hbond, std::string& failure_reason) const {
    // 1. Distance < min (steric clash)
    if (hbond.distance < params_.min_distance) {
        failure_reason = "Distance < " + std::to_string(params_.min_distance) + "Å (steric clash)";
        return true;
    }

    // 2. Distance > max (too far)
    if (hbond.distance > params_.max_distance) {
        failure_reason = "Distance > " + std::to_string(params_.max_distance) + "Å (too far)";
        return true;
    }

    // 3. Donor angle < min (impossible)
    if (hbond.donor_angle > 0.0 && hbond.donor_angle < params_.min_donor_angle) {
        failure_reason = "Donor angle " + std::to_string(static_cast<int>(hbond.donor_angle)) +
                         "° < " + std::to_string(static_cast<int>(params_.min_donor_angle)) + "° (impossible geometry)";
        return true;
    }

    // 4. Acceptor angle < min (impossible)
    if (hbond.acceptor_angle > 0.0 && hbond.acceptor_angle < params_.min_acceptor_angle) {
        failure_reason = "Acceptor angle " + std::to_string(static_cast<int>(hbond.acceptor_angle)) +
                         "° < " + std::to_string(static_cast<int>(params_.min_acceptor_angle)) + "° (impossible geometry)";
        return true;
    }

    return false;
}

double HBondQualityScorer::apply_soft_failure_caps(const core::HBond& hbond, double score) const {
    // Soft failure 1: Donor angle < 120° AND distance > 3.2Å
    if (hbond.donor_angle > 0.0 && hbond.donor_angle < 120.0 && hbond.distance > 3.2) {
        score = std::min(score, 40.0);
    }

    // Soft failure 2: Acceptor angle < 100° AND donor angle < 140°
    if (hbond.acceptor_angle > 0.0 && hbond.acceptor_angle < 100.0 &&
        hbond.donor_angle > 0.0 && hbond.donor_angle < 140.0) {
        score = std::min(score, 40.0);
    }

    return score;
}

bool HBondQualityScorer::is_sp2_acceptor(const std::string& atom_name) {
    // Trim whitespace
    std::string name = atom_name;
    name.erase(0, name.find_first_not_of(" \t"));
    name.erase(name.find_last_not_of(" \t") + 1);

    // Carbonyl oxygens (sp2)
    if (name == "O6" || name == "O4" || name == "O2") {
        return true;  // Guanine O6, Uracil/Thymine O4, Cytosine/Uracil O2
    }

    // Aromatic ring nitrogens (sp2)
    if (name == "N1" || name == "N3" || name == "N7") {
        return true;
    }

    // Phosphate oxygens (sp2-like)
    if (name == "O1P" || name == "O2P" || name == "OP1" || name == "OP2") {
        return true;
    }

    // Sugar/hydroxyl oxygens (sp3)
    // O2', O3', O4', O5' - these are sp3
    if (name.find('\'') != std::string::npos) {
        return false;  // Sugar atoms are sp3
    }

    // Default to sp2 (more common in nucleic acids)
    return true;
}

double HBondQualityScorer::adjust_for_resolution(double score, double resolution) const {
    if (!params_.apply_resolution_penalty) {
        return score;
    }

    // Unknown resolution (0 or negative) - apply conservative penalty
    if (resolution <= 0.0) {
        return score * 0.90;  // 10% penalty
    }

    // High resolution (< 2.0 Å) - no penalty
    if (resolution <= params_.high_res_threshold) {
        return score;
    }

    // Low resolution (> 3.5 Å) - maximum penalty
    if (resolution > params_.low_res_threshold) {
        return score * 0.95;  // 5% penalty
    }

    // Medium resolution (2.0-3.5 Å) - linear interpolation
    // 2.0 Å -> no penalty, 3.5 Å -> 5% penalty
    double fraction = (resolution - params_.high_res_threshold) /
                      (params_.low_res_threshold - params_.high_res_threshold);
    double penalty = 1.0 - (0.05 * fraction);  // 0% to 5% penalty

    return score * penalty;
}

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna

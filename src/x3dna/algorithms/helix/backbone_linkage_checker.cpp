/**
 * @file backbone_linkage_checker.cpp
 * @brief Implementation of backbone linkage detection
 */

#include <x3dna/algorithms/helix/backbone_linkage_checker.hpp>

namespace x3dna::algorithms::helix {

LinkDirection BackboneLinkageChecker::check_linkage(
    size_t res_i, size_t res_j,
    const BackboneData& backbone) const {

    auto it_i = backbone.find(res_i);
    auto it_j = backbone.find(res_j);

    if (it_i == backbone.end() || it_j == backbone.end()) {
        return LinkDirection::None;
    }

    const auto& atoms_i = it_i->second;
    const auto& atoms_j = it_j->second;

    // Check O3'[i] → P[j] direction (5'→3')
    if (atoms_i.O3_prime.has_value() && atoms_j.P.has_value()) {
        double dist = (atoms_i.O3_prime.value() - atoms_j.P.value()).length();
        if (dist <= config_.o3p_upper) {
            return LinkDirection::Forward;
        }
    }

    // Check O3'[j] → P[i] direction (reverse)
    if (atoms_j.O3_prime.has_value() && atoms_i.P.has_value()) {
        double dist = (atoms_j.O3_prime.value() - atoms_i.P.value()).length();
        if (dist <= config_.o3p_upper) {
            return LinkDirection::Reverse;
        }
    }

    return LinkDirection::None;
}

double BackboneLinkageChecker::o3_distance(
    size_t res_i, size_t res_j,
    const BackboneData& backbone) const {

    auto it_i = backbone.find(res_i);
    auto it_j = backbone.find(res_j);

    if (it_i == backbone.end() || it_j == backbone.end()) {
        return -1.0;
    }

    const auto& atoms_i = it_i->second;
    const auto& atoms_j = it_j->second;

    if (atoms_i.O3_prime.has_value() && atoms_j.O3_prime.has_value()) {
        return (atoms_i.O3_prime.value() - atoms_j.O3_prime.value()).length();
    }

    return -1.0;
}

bool BackboneLinkageChecker::are_pairs_connected(
    const core::BasePair& pair1,
    const core::BasePair& pair2,
    const BackboneData& backbone) const {

    if (backbone.empty()) {
        return true;  // Assume connected when no backbone data
    }

    // Get 1-based residue indices for both pairs
    size_t i1 = pair1.residue_idx1() + 1;
    size_t j1 = pair1.residue_idx2() + 1;
    size_t i2 = pair2.residue_idx1() + 1;
    size_t j2 = pair2.residue_idx2() + 1;

    // Check all possible strand linkages:
    // Strand 1: i1 -> i2 or i2 -> i1
    // Strand 2: j1 -> j2 or j2 -> j1
    // Cross: i1 -> j2 or j1 -> i2 (for strand swaps)
    return check_linkage(i1, i2, backbone) != LinkDirection::None ||
           check_linkage(j1, j2, backbone) != LinkDirection::None ||
           check_linkage(i1, j2, backbone) != LinkDirection::None ||
           check_linkage(j1, i2, backbone) != LinkDirection::None;
}

} // namespace x3dna::algorithms::helix

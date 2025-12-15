/**
 * @file helix_organizer.cpp
 * @brief Implementation of helix organization algorithm
 * 
 * This implements the legacy X3DNA five2three algorithm for ensuring
 * proper 5'→3' strand direction in base pair step calculations.
 */

#include <x3dna/algorithms/helix_organizer.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <iostream>  // For debug tracing

namespace x3dna::algorithms {

namespace {
    constexpr double PI = 3.14159265358979323846;
    
    // Convert dot product to angle in degrees
    double dot2ang(double d) {
        if (d > 1.0) d = 1.0;
        if (d < -1.0) d = -1.0;
        return std::acos(d) * 180.0 / PI;
    }
}

HelixOrganizer::HelixOrganizer(const Config& config) : config_(config) {}

// =============================================================================
// Geometry helpers
// =============================================================================

geometry::Vector3D HelixOrganizer::get_pair_origin(const core::BasePair& pair) const {
    // Legacy uses the average of both base origins (morg) for pair distance calculations
    // See refs_right_left() in cmn_fncs.c: morg[j] = (org[base1][j] + org[base2][j]) / 2
    if (!pair.frame1().has_value() && !pair.frame2().has_value()) {
        return geometry::Vector3D(0, 0, 0);
    }
    if (!pair.frame1().has_value()) {
        return pair.frame2().value().origin();
    }
    if (!pair.frame2().has_value()) {
        return pair.frame1().value().origin();
    }
    // Return average of both origins
    auto o1 = pair.frame1().value().origin();
    auto o2 = pair.frame2().value().origin();
    return geometry::Vector3D(
        (o1.x() + o2.x()) / 2.0,
        (o1.y() + o2.y()) / 2.0,
        (o1.z() + o2.z()) / 2.0
    );
}

geometry::Vector3D HelixOrganizer::get_pair_z_axis(const core::BasePair& pair) const {
    // Legacy uses average z-axis from both frames (or difference if they point opposite)
    // This matches bp_context: (d <= 0.0) ? ddxyz(z2, z1, zave) : sumxyz(z2, z1, zave)
    if (!pair.frame1().has_value() || !pair.frame2().has_value()) {
        if (pair.frame1().has_value()) {
            return pair.frame1().value().z_axis();
        }
        return geometry::Vector3D(0, 0, 1);
    }

    auto z1 = pair.frame1().value().z_axis();
    auto z2 = pair.frame2().value().z_axis();
    double d = z1.dot(z2);

    geometry::Vector3D zave;
    if (d <= 0.0) {
        // Opposite directions: take difference
        zave = z2 - z1;
    } else {
        // Same direction: take sum
        zave = z2 + z1;
    }
    zave.normalize();
    return zave;
}

geometry::Vector3D HelixOrganizer::get_frame_z(const core::BasePair& pair, bool swapped) const {
    if (swapped) {
        if (pair.frame2().has_value()) {
            return pair.frame2().value().z_axis();
        }
    } else {
        if (pair.frame1().has_value()) {
            return pair.frame1().value().z_axis();
        }
    }
    return geometry::Vector3D(0, 0, 1);
}

StrandResidues HelixOrganizer::get_strand_residues(
    const core::BasePair& pair, bool swapped) const {
    // BasePair stores 0-based indices, but backbone data uses 1-based (legacy) indices
    // Convert to 1-based for backbone lookup
    if (swapped) {
        return {pair.residue_idx2() + 1, pair.residue_idx1() + 1};
    }
    return {pair.residue_idx1() + 1, pair.residue_idx2() + 1};
}

// =============================================================================
// Backbone connectivity
// =============================================================================

LinkDirection HelixOrganizer::check_linkage(
    size_t res_i, size_t res_j, const BackboneData& backbone) const {

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

double HelixOrganizer::o3_distance(size_t i, size_t j, const BackboneData& backbone) const {
    auto it_i = backbone.find(i);
    auto it_j = backbone.find(j);

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

bool HelixOrganizer::are_pairs_backbone_connected(
    const core::BasePair& pair1, const core::BasePair& pair2,
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

// =============================================================================
// WC pair geometry checks
// =============================================================================

double HelixOrganizer::wcbp_xang(const core::BasePair& pair_m, const core::BasePair& pair_n) const {
    // Calculate angle between combined x-axes of the two pairs
    if (!pair_m.frame1().has_value() || !pair_m.frame2().has_value() ||
        !pair_n.frame1().has_value() || !pair_n.frame2().has_value()) {
        return 180.0;  // Return large angle if frames missing
    }
    
    // Sum of x-axes for each pair
    auto xm = pair_m.frame1().value().x_axis() + pair_m.frame2().value().x_axis();
    auto xn = pair_n.frame1().value().x_axis() + pair_n.frame2().value().x_axis();
    
    xm = xm.normalized();
    xn = xn.normalized();
    
    double dot = xm.dot(xn);
    return dot2ang(dot);
}

double HelixOrganizer::wcbp_zdir(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                  bool swap_m, bool swap_n) const {
    // Get z-direction vectors based on swap status
    geometry::Vector3D zm, zn;
    
    if (swap_m) {
        if (pair_m.frame2().has_value() && pair_m.frame1().has_value()) {
            zm = pair_m.frame1().value().z_axis() - pair_m.frame2().value().z_axis();
        } else {
            return 0.0;
        }
    } else {
        if (pair_m.frame1().has_value() && pair_m.frame2().has_value()) {
            zm = pair_m.frame2().value().z_axis() - pair_m.frame1().value().z_axis();
        } else {
            return 0.0;
        }
    }
    
    if (swap_n) {
        if (pair_n.frame2().has_value() && pair_n.frame1().has_value()) {
            zn = pair_n.frame1().value().z_axis() - pair_n.frame2().value().z_axis();
        } else {
            return 0.0;
        }
    } else {
        if (pair_n.frame1().has_value() && pair_n.frame2().has_value()) {
            zn = pair_n.frame2().value().z_axis() - pair_n.frame1().value().z_axis();
        } else {
            return 0.0;
        }
    }
    
    zm = zm.normalized();
    zn = zn.normalized();
    
    return zm.dot(zn);
}

// =============================================================================
// Five2three sub-functions
// =============================================================================

void HelixOrganizer::first_step(const std::vector<core::BasePair>& pairs,
                                const BackboneData& backbone,
                                std::vector<size_t>& pair_order,
                                const HelixSegment& helix,
                                std::vector<bool>& swapped) const {
    // For single-pair helices, nothing to do
    if (helix.end_idx <= helix.start_idx) {
        return;
    }

    // Look at ONLY the first step (matching legacy behavior)
    size_t pos = helix.start_idx;
    size_t first_pair = pair_order[pos];
    size_t second_pair = pair_order[pos + 1];

    auto res_m = get_strand_residues(pairs[first_pair], swapped[first_pair]);
    auto res_n = get_strand_residues(pairs[second_pair], swapped[second_pair]);

    // Check if strand 1 residues are linked via backbone
    auto link = check_linkage(res_m.strand1, res_n.strand1, backbone);

    if (link == LinkDirection::Reverse) {
        // Reverse linkage - swap the first pair
        swapped[first_pair] = !swapped[first_pair];
    } else if (link == LinkDirection::None) {
        // No linkage - reverse the entire helix segment and try again
        std::reverse(pair_order.begin() + helix.start_idx,
                     pair_order.begin() + helix.end_idx + 1);

        // Update pair indices after reversal
        first_pair = pair_order[pos];
        second_pair = pair_order[pos + 1];

        res_m = get_strand_residues(pairs[first_pair], swapped[first_pair]);
        res_n = get_strand_residues(pairs[second_pair], swapped[second_pair]);

        link = check_linkage(res_m.strand1, res_n.strand1, backbone);

        if (link == LinkDirection::Reverse) {
            swapped[first_pair] = !swapped[first_pair];
        } else if (link == LinkDirection::None) {
            // Still no linkage - undo the reversal
            std::reverse(pair_order.begin() + helix.start_idx,
                         pair_order.begin() + helix.end_idx + 1);
        }
    }
    // Forward linkage: already correct, no action needed
}

bool HelixOrganizer::wc_bporien(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                 bool swap_m, bool swap_n,
                                 const BackboneData& backbone) const {
    // Only check WC pairs (bp_type > 0 in legacy)
    // For now, check if both pairs have proper frames
    if (!pair_m.frame1().has_value() || !pair_m.frame2().has_value() ||
        !pair_n.frame1().has_value() || !pair_n.frame2().has_value()) {
        return false;
    }

    auto res_m = get_strand_residues(pair_m, swap_m);
    auto res_n = get_strand_residues(pair_n, swap_n);

    // If x-angle is too large or backbone is linked, don't swap
    if (wcbp_xang(pair_m, pair_n) > config_.end_stack_xang ||
        check_linkage(res_m.strand1, res_n.strand1, backbone) != LinkDirection::None ||
        check_linkage(res_m.strand2, res_n.strand2, backbone) != LinkDirection::None) {
        return false;
    }

    // Check z-direction alignment
    double zdir_normal = wcbp_zdir(pair_m, pair_n, swap_m, swap_n);
    double zdir_swapped = wcbp_zdir(pair_m, pair_n, swap_m, !swap_n);

    return zdir_normal < 0.0 && zdir_swapped > 0.0;  // Need to swap pair_n
}

bool HelixOrganizer::check_o3dist(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                   bool swap_m, bool swap_n,
                                   const BackboneData& backbone) const {
    auto res_m = get_strand_residues(pair_m, swap_m);
    auto res_n = get_strand_residues(pair_n, swap_n);

    double di1_i2 = o3_distance(res_m.strand1, res_n.strand1, backbone);
    double di1_j2 = o3_distance(res_m.strand1, res_n.strand2, backbone);
    double dj1_i2 = o3_distance(res_m.strand2, res_n.strand1, backbone);
    double dj1_j2 = o3_distance(res_m.strand2, res_n.strand2, backbone);

    // Debug trace: print all check_o3dist calls that return true
    bool result = (di1_i2 > 0.0 && di1_j2 > 0.0 && di1_i2 > di1_j2) &&
                  (dj1_i2 > 0.0 && dj1_j2 > 0.0 && dj1_j2 > dj1_i2);
    if (result) {
        std::cerr << "[O3DIST TRUE] res_m=(" << res_m.strand1 << "," << res_m.strand2
                  << ") res_n=(" << res_n.strand1 << "," << res_n.strand2 << ")"
                  << " di1_i2=" << di1_i2 << " di1_j2=" << di1_j2
                  << " dj1_i2=" << dj1_i2 << " dj1_j2=" << dj1_j2 << "\n";
    }

    return result;
}

bool HelixOrganizer::check_schain(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                   bool swap_m, bool swap_n,
                                   const BackboneData& backbone) const {
    auto res_m = get_strand_residues(pair_m, swap_m);
    auto res_n = get_strand_residues(pair_n, swap_n);

    // If no same-strand linkage but cross-strand linkage exists, swap is indicated
    bool no_same_strand = check_linkage(res_m.strand1, res_n.strand1, backbone) == LinkDirection::None &&
                          check_linkage(res_m.strand2, res_n.strand2, backbone) == LinkDirection::None;
    bool has_cross_strand = check_linkage(res_m.strand1, res_n.strand2, backbone) != LinkDirection::None ||
                            check_linkage(res_m.strand2, res_n.strand1, backbone) != LinkDirection::None;

    return no_same_strand && has_cross_strand;
}

bool HelixOrganizer::check_others(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                   bool swap_m, bool swap_n,
                                   const BackboneData& backbone) const {
    auto res_m = get_strand_residues(pair_m, swap_m);
    auto res_n = get_strand_residues(pair_n, swap_n);

    // If any backbone linkage exists, no swap needed
    if (check_linkage(res_m.strand1, res_n.strand1, backbone) != LinkDirection::None ||
        check_linkage(res_m.strand2, res_n.strand2, backbone) != LinkDirection::None ||
        check_linkage(res_m.strand1, res_n.strand2, backbone) != LinkDirection::None ||
        check_linkage(res_m.strand2, res_n.strand1, backbone) != LinkDirection::None) {
        return false;
    }
    
    // Check frame alignment patterns
    if (!pair_m.frame1().has_value() || !pair_m.frame2().has_value() ||
        !pair_n.frame1().has_value() || !pair_n.frame2().has_value()) {
        return false;
    }
    
    // Get frames based on swap status
    auto frame_m1 = swap_m ? pair_m.frame2().value() : pair_m.frame1().value();
    auto frame_m2 = swap_m ? pair_m.frame1().value() : pair_m.frame2().value();
    auto frame_n1 = swap_n ? pair_n.frame2().value() : pair_n.frame1().value();
    auto frame_n2 = swap_n ? pair_n.frame1().value() : pair_n.frame2().value();
    
    // Check axis alignment (similar to legacy dot product checks)
    double a1_x = frame_m1.x_axis().dot(frame_n1.x_axis());
    double a1_y = frame_m1.y_axis().dot(frame_n1.y_axis());
    double a1_z = frame_m1.z_axis().dot(frame_n1.z_axis());
    
    double a2_x = frame_m2.x_axis().dot(frame_n2.x_axis());
    double a2_y = frame_m2.y_axis().dot(frame_n2.y_axis());
    double a2_z = frame_m2.z_axis().dot(frame_n2.z_axis());
    
    bool aligned1 = (a1_x > 0.0 && a1_y > 0.0 && a1_z > 0.0);
    bool aligned2 = (a2_x > 0.0 && a2_y > 0.0 && a2_z > 0.0);
    
    if (aligned1 && aligned2) {
        return false;
    }
    
    // Check cross-alignment (m1 with n2, m2 with n1)
    double r1_x = frame_m1.x_axis().dot(frame_n2.x_axis());
    double r1_y = frame_m1.y_axis().dot(frame_n2.y_axis());
    double r1_z = frame_m1.z_axis().dot(frame_n2.z_axis());
    
    double r2_x = frame_m2.x_axis().dot(frame_n1.x_axis());
    double r2_y = frame_m2.y_axis().dot(frame_n1.y_axis());
    double r2_z = frame_m2.z_axis().dot(frame_n1.z_axis());
    
    bool cross1 = (r1_x > 0.0 && r1_y > 0.0 && r1_z > 0.0);
    bool cross2 = (r2_x > 0.0 && r2_y > 0.0 && r2_z > 0.0);
    
    if (!aligned1 && !aligned2) {
        if (cross1 || cross2) {
            return true;
        }
    }
    
    // Compare total angles for mixed cases
    if ((aligned1 || aligned2) && (cross1 || cross2)) {
        double sum_aligned = dot2ang(a1_x) + dot2ang(a1_y) + dot2ang(a1_z) +
                            dot2ang(a2_x) + dot2ang(a2_y) + dot2ang(a2_z);
        double sum_cross = dot2ang(r1_x) + dot2ang(r1_y) + dot2ang(r1_z) +
                          dot2ang(r2_x) + dot2ang(r2_y) + dot2ang(r2_z);
        
        if (sum_aligned > sum_cross) {
            return true;
        }
    }
    
    return false;
}

bool HelixOrganizer::chain1dir(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                bool swap_m, bool swap_n,
                                const BackboneData& backbone) const {
    auto res_m = get_strand_residues(pair_m, swap_m);
    auto res_n = get_strand_residues(pair_n, swap_n);

    // Reverse linkage on strand 1 indicates need to swap
    return check_linkage(res_m.strand1, res_n.strand1, backbone) == LinkDirection::Reverse;
}

DirectionCounts HelixOrganizer::check_direction(
    const std::vector<core::BasePair>& pairs,
    const BackboneData& backbone,
    std::vector<size_t>& pair_order,
    HelixSegment& helix,
    std::vector<bool>& swapped) const {

    DirectionCounts dir;

    // First pass: compute direction counts
    for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
        size_t idx_m = pair_order[pos];
        size_t idx_n = pair_order[pos + 1];

        auto res_m = get_strand_residues(pairs[idx_m], swapped[idx_m]);
        auto res_n = get_strand_residues(pairs[idx_n], swapped[idx_n]);

        // Strand 1 direction
        auto link1 = check_linkage(res_m.strand1, res_n.strand1, backbone);
        if (link1 == LinkDirection::Forward) dir.strand1_forward++;
        else if (link1 == LinkDirection::Reverse) dir.strand1_reverse++;
        else dir.strand1_none++;

        // Strand 2 direction
        auto link2 = check_linkage(res_m.strand2, res_n.strand2, backbone);
        if (link2 == LinkDirection::Forward) dir.strand2_forward++;
        else if (link2 == LinkDirection::Reverse) dir.strand2_reverse++;
        else dir.strand2_none++;
    }

    // Check for mixed direction (returns early in legacy)
    bool mixed = (dir.strand1_forward && dir.strand1_reverse) ||
                 (dir.strand2_forward && dir.strand2_reverse);
    if (mixed) {
        helix.has_mixed_direction = true;
        return dir;
    }

    // No linkages at all
    if (dir.strand1_forward + dir.strand1_reverse +
        dir.strand2_forward + dir.strand2_reverse == 0) {
        return dir;
    }

    // Get first and last pair info
    size_t first_pair_idx = pair_order[helix.start_idx];
    size_t last_pair_idx = pair_order[helix.end_idx];
    auto res_first = get_strand_residues(pairs[first_pair_idx], swapped[first_pair_idx]);
    auto res_last = get_strand_residues(pairs[last_pair_idx], swapped[last_pair_idx]);

    // Set break flag if there are gaps
    if (dir.strand1_none || dir.strand2_none) {
        helix.has_break = true;
    }

    // Anti-parallel case: strand1 forward, strand2 reverse
    if (dir.strand1_forward && !dir.strand1_reverse) {
        if (!dir.strand2_forward && dir.strand2_reverse) {
            // Normal anti-parallel
            if (res_first.strand1 > res_last.strand2) {
                // Flip all swapped values in helix AND reverse order
                for (size_t pos = helix.start_idx; pos <= helix.end_idx; ++pos) {
                    size_t idx = pair_order[pos];
                    swapped[idx] = !swapped[idx];
                }
                // Reverse helix order
                std::reverse(pair_order.begin() + helix.start_idx,
                             pair_order.begin() + helix.end_idx + 1);
            }
        } else if (dir.strand2_forward && !dir.strand2_reverse) {
            // Parallel case
            helix.is_parallel = true;
            if (res_first.strand1 > res_first.strand2) {
                // Flip all swapped values in helix
                for (size_t pos = helix.start_idx; pos <= helix.end_idx; ++pos) {
                    size_t idx = pair_order[pos];
                    swapped[idx] = !swapped[idx];
                }
            }
        }
    }

    return dir;
}

void HelixOrganizer::check_strand2(const std::vector<core::BasePair>& pairs,
                                   const BackboneData& backbone,
                                   const std::vector<size_t>& pair_order,
                                   HelixSegment& helix,
                                   std::vector<bool>& swapped,
                                   const DirectionCounts& direction) const {

    bool mixed_direction = (direction.strand1_forward && direction.strand1_reverse) ||
                          (direction.strand2_forward && direction.strand2_reverse);

    if (!mixed_direction) {
        // Normal case - check for cross-strand swaps
        if (direction.strand1_forward + direction.strand1_reverse +
            direction.strand2_forward + direction.strand2_reverse == 0) {
            return;  // No linkages to check
        }

        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];

            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];

            // Skip if WC orientation check would pass
            if (wc_bporien(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone)) {
                continue;
            }

            auto res_m = get_strand_residues(pair_m, swapped[idx_m]);
            auto res_n = get_strand_residues(pair_n, swapped[idx_n]);

            // Check for cross-strand linkage indicating swap needed
            bool no_same_strand = check_linkage(res_m.strand1, res_n.strand1, backbone) == LinkDirection::None &&
                                  check_linkage(res_m.strand2, res_n.strand2, backbone) == LinkDirection::None;
            auto cross_link = check_linkage(res_m.strand1, res_n.strand2, backbone);
            auto cross_link2 = check_linkage(res_m.strand2, res_n.strand1, backbone);

            if (no_same_strand &&
                (cross_link == LinkDirection::Forward ||
                 (cross_link != LinkDirection::None && cross_link2 != LinkDirection::None))) {
                swapped[idx_n] = !swapped[idx_n];
            }
        }
    } else {
        // Mixed direction case - more complex handling
        bool anti_p = (direction.strand1_forward > direction.strand1_reverse) &&
                     (direction.strand2_forward < direction.strand2_reverse);
        bool parallel = (direction.strand1_forward > direction.strand1_reverse) &&
                       (direction.strand2_forward > direction.strand2_reverse);

        helix.is_parallel = parallel;

        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];

            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];

            auto res_m = get_strand_residues(pair_m, swapped[idx_m]);
            auto res_n = get_strand_residues(pair_n, swapped[idx_n]);

            auto link_strand2 = check_linkage(res_m.strand2, res_n.strand2, backbone);
            if (check_linkage(res_m.strand1, res_n.strand1, backbone) == LinkDirection::None &&
                ((anti_p && link_strand2 == LinkDirection::Forward) ||
                 (parallel && link_strand2 == LinkDirection::Reverse))) {
                swapped[idx_n] = !swapped[idx_n];
            }

            // Re-get residues after potential swap
            res_n = get_strand_residues(pair_n, swapped[idx_n]);

            if (check_linkage(res_m.strand1, res_n.strand1, backbone) == LinkDirection::None &&
                check_linkage(res_m.strand2, res_n.strand2, backbone) == LinkDirection::None) {
                if ((anti_p && check_linkage(res_m.strand2, res_n.strand1, backbone) == LinkDirection::Forward) ||
                    (parallel && check_linkage(res_m.strand1, res_n.strand2, backbone) == LinkDirection::Reverse)) {
                    swapped[idx_m] = !swapped[idx_m];
                } else if ((anti_p && check_linkage(res_m.strand1, res_n.strand2, backbone) == LinkDirection::Forward) ||
                           (parallel && check_linkage(res_m.strand2, res_n.strand1, backbone) == LinkDirection::Reverse)) {
                    swapped[idx_n] = !swapped[idx_n];
                }
            }
        }
    }
}

// =============================================================================
// Context calculation (bp_context equivalent)
// =============================================================================

std::vector<HelixOrganizer::PairContext> HelixOrganizer::calculate_context(
    const std::vector<core::BasePair>& pairs,
    const BackboneData& backbone) const {

    size_t n = pairs.size();
    std::vector<PairContext> context(n);

    if (n < 2) return context;

    for (size_t i = 0; i < n; ++i) {
        auto org_i = get_pair_origin(pairs[i]);
        auto z_i = get_pair_z_axis(pairs[i]);

        std::vector<std::pair<double, size_t>> neighbors;

        for (size_t j = 0; j < n; ++j) {
            if (j == i) continue;

            auto org_j = get_pair_origin(pairs[j]);
            double dist = (org_j - org_i).length();

            if (dist <= config_.neighbor_cutoff) {
                neighbors.emplace_back(dist, j);
            }
        }

        std::sort(neighbors.begin(), neighbors.end());

        if (neighbors.empty()) {
            context[i].is_endpoint = true;
            continue;
        }

        context[i].neighbor1 = neighbors[0].second;
        context[i].dist1 = neighbors[0].first;

        // Check backbone connectivity to neighbor1
        context[i].has_backbone_link1 = are_pairs_backbone_connected(
            pairs[i], pairs[neighbors[0].second], backbone);

        if (context[i].dist1 > config_.helix_break) {
            context[i].is_endpoint = true;
            continue;
        }

        auto v1 = get_pair_origin(pairs[neighbors[0].second]) - org_i;
        double d1 = z_i.dot(v1);

        // Legacy lines 931-941: If 2nd and 3rd closest are both on opposite z-side,
        // swap them if 2nd has larger |z-distance| (prefer smaller |z-distance|)
        if (neighbors.size() >= 3 &&
            neighbors[1].first <= config_.helix_break &&
            neighbors[2].first <= config_.helix_break) {
            auto v2 = get_pair_origin(pairs[neighbors[1].second]) - org_i;
            auto v3 = get_pair_origin(pairs[neighbors[2].second]) - org_i;
            double d2 = z_i.dot(v2);
            double d3 = z_i.dot(v3);

            // Both on opposite z-side from n1, and 2nd has larger |z-dist|
            if (d1 * d2 < 0 && d1 * d3 < 0 && std::abs(d2) > std::abs(d3)) {
                // Swap 2nd and 3rd
                std::swap(neighbors[1], neighbors[2]);
            }
        }

        for (size_t k = 1; k < neighbors.size(); ++k) {
            if (neighbors[k].first > config_.helix_break) break;

            auto vk = get_pair_origin(pairs[neighbors[k].second]) - org_i;
            double dk = z_i.dot(vk);

            if (d1 * dk < 0) {
                context[i].neighbor2 = neighbors[k].second;
                context[i].dist2 = neighbors[k].first;
                // Check backbone connectivity to neighbor2
                context[i].has_backbone_link2 = are_pairs_backbone_connected(
                    pairs[i], pairs[neighbors[k].second], backbone);
                context[i].is_endpoint = false;
                break;
            }
        }

        if (!context[i].neighbor2.has_value()) {
            context[i].is_endpoint = true;

            // Legacy special case (find_pair.c lines 971-976):
            // Even for endpoints, try to find n2 through an indirect check:
            // If vector from n1 to 2nd closest is on opposite z-side AND within helix_break
            if (neighbors.size() >= 2) {
                size_t n2_idx = neighbors[1].second;
                auto org_n1 = get_pair_origin(pairs[neighbors[0].second]);
                auto org_n2 = get_pair_origin(pairs[n2_idx]);
                // Legacy ddxyz(n1, n2) = n1 - n2, so we compute org_n1 - org_n2
                auto v_n2_n1 = org_n1 - org_n2;
                double dist_n1_n2 = v_n2_n1.length();
                double d2 = z_i.dot(v_n2_n1);

                // If n2->n1 is on opposite z-side from i->n1 AND within helix_break
                if (d1 * d2 < 0 && dist_n1_n2 <= config_.helix_break) {
                    context[i].neighbor2 = n2_idx;
                    context[i].dist2 = neighbors[1].first;
                    context[i].has_backbone_link2 = are_pairs_backbone_connected(
                        pairs[i], pairs[n2_idx], backbone);
                    // Still an endpoint, but now has a neighbor2
                }
            }
        }
    }

    return context;
}

std::vector<size_t> HelixOrganizer::find_endpoints(
    const std::vector<PairContext>& context) const {

    std::vector<size_t> endpoints;

    for (size_t i = 0; i < context.size(); ++i) {
        if (context[i].is_endpoint) {
            endpoints.push_back(i);
        }
    }

    if (endpoints.empty() && !context.empty()) {
        endpoints.push_back(0);
    }

    return endpoints;
}

std::pair<std::vector<size_t>, std::vector<HelixSegment>> HelixOrganizer::locate_helices(
    const std::vector<PairContext>& context,
    const std::vector<size_t>& endpoints,
    const BackboneData& backbone,
    size_t num_pairs) const {

    std::vector<size_t> pair_order;
    std::vector<HelixSegment> helices;
    std::vector<bool> visited(num_pairs, false);

    (void)backbone;  // Not used in traverse - backbone connectivity checked in calculate_context
    pair_order.reserve(num_pairs);

    // Legacy locate_helix algorithm (find_pair.c lines 1029-1082):
    // For each endpoint, add endpoint + neighbors from end_list, then traverse
    for (size_t ep : endpoints) {
        // Skip if all items in this endpoint's "end_list" are already matched
        const auto& ep_ctx = context[ep];
        int matched_count = 0;
        int item_count = 1;  // The endpoint itself
        if (visited[ep]) matched_count++;
        if (ep_ctx.neighbor1.has_value()) {
            item_count++;
            if (visited[ep_ctx.neighbor1.value()]) matched_count++;
        }
        if (ep_ctx.neighbor2.has_value()) {
            item_count++;
            if (visited[ep_ctx.neighbor2.value()]) matched_count++;
        }
        if (matched_count == item_count) continue;

        HelixSegment helix;
        helix.start_idx = pair_order.size();

        // Legacy lines 1039-1045: Add endpoint and its neighbors from end_list
        // Add endpoint
        if (!visited[ep]) {
            pair_order.push_back(ep);
            visited[ep] = true;
        }
        // Add neighbor1
        if (ep_ctx.neighbor1.has_value() && !visited[ep_ctx.neighbor1.value()]) {
            pair_order.push_back(ep_ctx.neighbor1.value());
            visited[ep_ctx.neighbor1.value()] = true;
        }
        // Add neighbor2
        if (ep_ctx.neighbor2.has_value() && !visited[ep_ctx.neighbor2.value()]) {
            pair_order.push_back(ep_ctx.neighbor2.value());
            visited[ep_ctx.neighbor2.value()] = true;
        }

        // Legacy lines 1046-1068: Traverse from the last added pair
        // Continue adding neighbors until we can't anymore
        while (pair_order.size() >= helix.start_idx + 1) {
            size_t ip = pair_order.size() - 1;
            size_t current = pair_order[ip];
            const auto& ctx = context[current];

            // If this is an endpoint (bp_order[k][1] == 0 in legacy means endpoint)
            if (ctx.is_endpoint) {
                // Legacy lines 1050-1055: Only add neighbor1 if it exists, not matched, and no neighbor2
                if (ctx.neighbor1.has_value() && !visited[ctx.neighbor1.value()] &&
                    !ctx.neighbor2.has_value()) {
                    pair_order.push_back(ctx.neighbor1.value());
                    visited[ctx.neighbor1.value()] = true;
                }
                break;  // Stop traversal for this helix
            }

            // Not an endpoint - continue traversal
            bool n1_matched = !ctx.neighbor1.has_value() || visited[ctx.neighbor1.value()];
            bool n2_matched = !ctx.neighbor2.has_value() || visited[ctx.neighbor2.value()];

            // Legacy line 1057-1059: If both or neither neighbor matched, stop
            if ((n1_matched && n2_matched) || (!n1_matched && !n2_matched && ctx.neighbor1.has_value() && ctx.neighbor2.has_value())) {
                // Both matched or neither matched - stop
                // For "neither matched" case, only apply if both neighbors actually exist
                break;
            }

            // Get previous pair
            std::optional<size_t> prev;
            if (ip > helix.start_idx) {
                prev = pair_order[ip - 1];
            }

            // Legacy lines 1060-1067: Go to the other neighbor
            std::optional<size_t> next;
            if (ctx.neighbor1.has_value() && prev.has_value() && ctx.neighbor1.value() == prev.value()) {
                // Previous was neighbor1, add neighbor2
                if (ctx.neighbor2.has_value() && !visited[ctx.neighbor2.value()]) {
                    next = ctx.neighbor2;
                }
            } else if (ctx.neighbor2.has_value() && prev.has_value() && ctx.neighbor2.value() == prev.value()) {
                // Previous was neighbor2, add neighbor1
                if (ctx.neighbor1.has_value() && !visited[ctx.neighbor1.value()]) {
                    next = ctx.neighbor1;
                }
            }

            if (!next.has_value()) {
                break;
            }

            pair_order.push_back(next.value());
            visited[next.value()] = true;
        }

        helix.end_idx = pair_order.size() - 1;
        if (helix.end_idx >= helix.start_idx) {
            helices.push_back(helix);
        }
    }

    // Handle any leftover pairs not reached from endpoints
    for (size_t i = 0; i < num_pairs; ++i) {
        if (!visited[i]) {
            HelixSegment helix;
            helix.start_idx = pair_order.size();
            pair_order.push_back(i);
            helix.end_idx = pair_order.size() - 1;
            helices.push_back(helix);
        }
    }

    return {pair_order, helices};
}

// =============================================================================
// Main five2three algorithm
// =============================================================================

void HelixOrganizer::ensure_five_to_three(
    const std::vector<core::BasePair>& pairs,
    const BackboneData& backbone,
    std::vector<size_t>& pair_order,
    std::vector<HelixSegment>& helices,
    std::vector<bool>& swapped) const {
    
    swapped.resize(pairs.size(), false);
    
    if (backbone.empty()) {
        return;
    }
    
    // Process each helix
    for (auto& helix : helices) {
        if (helix.start_idx > helix.end_idx) continue;

        // STEP 1: first_step - set initial strand assignment
        first_step(pairs, backbone, pair_order, helix, swapped);

        // STEP 2: First pass through steps - check each consecutive pair
        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];

            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];

            bool rev_wc = wc_bporien(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            bool rev_o3d = check_o3dist(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            bool rev_csc = check_schain(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            bool rev_oth = check_others(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);

            // Apply swap based on checks
            if (rev_wc) {
                swapped[idx_n] = !swapped[idx_n];
            } else if (rev_o3d || rev_csc || rev_oth) {
                swapped[idx_n] = !swapped[idx_n];
            }

            // Check strand 1 direction
            bool rev_s1 = chain1dir(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            if (rev_s1) {
                swapped[idx_n] = !swapped[idx_n];
            }
        }

        // STEP 3: Second pass - re-check WC orientation
        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];

            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];

            bool rev_wc = wc_bporien(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            if (rev_wc) {
                swapped[idx_m] = !swapped[idx_m];
            }
        }

        // STEP 4: check_direction - count backbone linkage directions and apply fixes
        DirectionCounts direction = check_direction(pairs, backbone, pair_order, helix, swapped);

        // STEP 5: check_strand2 - additional corrections based on direction
        check_strand2(pairs, backbone, pair_order, helix, swapped, direction);

        // STEP 6: check_direction AGAIN (legacy line 1361 - at end of check_strand2)
        // This recomputes direction with updated swaps and may apply additional corrections
        check_direction(pairs, backbone, pair_order, helix, swapped);
    }
}

// =============================================================================
// Main organize function
// =============================================================================

HelixOrdering HelixOrganizer::organize(const std::vector<core::BasePair>& pairs,
                                        const BackboneData& backbone) const {
    HelixOrdering result;

    if (pairs.empty()) {
        return result;
    }

    if (pairs.size() == 1) {
        result.pair_order = {0};
        result.helices = {{0, 0, false, false, false}};
        result.strand_swapped = {false};
        result.helix_breaks = {};  // No breaks for single pair
        return result;
    }

    // Step 1: Calculate neighbor context (now includes backbone connectivity)
    auto context = calculate_context(pairs, backbone);

    // Step 2: Find helix endpoints
    auto endpoints = find_endpoints(context);

    // Step 3: Chain pairs into helices (respects backbone connectivity)
    auto [pair_order, helices] = locate_helices(context, endpoints, backbone, pairs.size());

    // Step 4: Ensure 5'→3' direction (full five2three algorithm)
    std::vector<bool> strand_swapped;
    ensure_five_to_three(pairs, backbone, pair_order, helices, strand_swapped);

    // Step 5: Identify helix breaks (positions without backbone connectivity)
    std::vector<bool> helix_breaks(pair_order.size() > 0 ? pair_order.size() - 1 : 0, false);
    if (!backbone.empty() && pair_order.size() > 1) {
        for (size_t i = 0; i + 1 < pair_order.size(); ++i) {
            size_t idx1 = pair_order[i];
            size_t idx2 = pair_order[i + 1];
            bool connected = are_pairs_backbone_connected(pairs[idx1], pairs[idx2], backbone);
            helix_breaks[i] = !connected;
        }
    }

    result.pair_order = std::move(pair_order);
    result.helices = std::move(helices);
    result.strand_swapped = std::move(strand_swapped);
    result.helix_breaks = std::move(helix_breaks);

    // Convert internal PairContext to public PairContextInfo for debugging
    result.context.reserve(context.size());
    for (const auto& ctx : context) {
        PairContextInfo info;
        info.is_endpoint = ctx.is_endpoint;
        info.neighbor1 = ctx.neighbor1;
        info.neighbor2 = ctx.neighbor2;
        result.context.push_back(info);
    }

    return result;
}

} // namespace x3dna::algorithms

/**
 * @file helix_organizer.cpp
 * @brief Implementation of helix organization algorithm
 * 
 * This implements the legacy X3DNA five2three algorithm for ensuring
 * proper 5'→3' strand direction in base pair step calculations.
 */

#include <x3dna/algorithms/helix_organizer.hpp>
#include <x3dna/config/config_manager.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace x3dna::algorithms {

namespace {
    constexpr double PI = 3.14159265358979323846;

    /**
     * @brief Result of frame axis alignment comparison
     */
    struct FrameAlignment {
        double dot_x;  ///< Dot product of x-axes
        double dot_y;  ///< Dot product of y-axes
        double dot_z;  ///< Dot product of z-axes

        [[nodiscard]] bool is_aligned() const {
            return dot_x > 0.0 && dot_y > 0.0 && dot_z > 0.0;
        }

        [[nodiscard]] double angle_sum() const;  // Forward declaration, defined after dot2ang
    };

    /**
     * @brief Check if five2three debugging is enabled
     *
     * Uses ConfigManager for debug settings, which reads from environment
     * variables on first access. This centralizes the debug configuration.
     */
    bool is_five2three_debug_enabled() {
        static bool initialized = false;
        static bool debug_enabled = false;
        if (!initialized) {
            // Initialize from ConfigManager (which reads env vars)
            auto& cfg = config::ConfigManager::instance();
            cfg.init_debug_from_environment();
            debug_enabled = cfg.debug_config().debug_five2three;
            initialized = true;
        }
        return debug_enabled;
    }

    // Convert dot product to angle in degrees
    double dot2ang(double d) {
        if (d > 1.0) d = 1.0;
        if (d < -1.0) d = -1.0;
        return std::acos(d) * 180.0 / PI;
    }

    // Implementation of FrameAlignment::angle_sum (after dot2ang is defined)
    double FrameAlignment::angle_sum() const {
        return dot2ang(dot_x) + dot2ang(dot_y) + dot2ang(dot_z);
    }

    /**
     * @brief Compute alignment between two reference frames
     * @param frame1 First reference frame
     * @param frame2 Second reference frame
     * @return FrameAlignment with dot products and alignment status
     */
    FrameAlignment compute_frame_alignment(const core::ReferenceFrame& frame1,
                                           const core::ReferenceFrame& frame2) {
        return FrameAlignment{
            frame1.x_axis().dot(frame2.x_axis()),
            frame1.y_axis().dot(frame2.y_axis()),
            frame1.z_axis().dot(frame2.z_axis())
        };
    }

    /**
     * @brief Update direction count based on linkage type
     * @param link The linkage direction
     * @param forward Forward count to increment
     * @param reverse Reverse count to increment
     * @param none None count to increment
     */
    void update_direction_count(LinkDirection link, int& forward, int& reverse, int& none) {
        if (link == LinkDirection::Forward) ++forward;
        else if (link == LinkDirection::Reverse) ++reverse;
        else ++none;
    }

    /**
     * @brief Check if two z-distances are on opposite sides of origin
     * @param d1 First z-distance
     * @param d2 Second z-distance
     * @return true if d1 and d2 have opposite signs
     */
    bool are_on_opposite_z_sides(double d1, double d2) {
        return d1 * d2 < 0;
    }

    /**
     * @brief Check if a pair has positive bpid (WC-like geometry)
     *
     * Legacy code in calculate_more_bppars sets bpid = -1 by default,
     * and only sets it positive if:
     *   1. dir_x > 0 && dir_y < 0 && dir_z < 0
     *   2. Then check_wc_wobble_pair validates shear/stretch/opening
     *
     * This function checks the first geometric condition.
     * The wc_bporien check requires base_pairs[m][3] > 0 && base_pairs[n][3] > 0.
     */
    bool has_positive_bpid(const x3dna::core::BasePair& pair) {
        bool debug = is_five2three_debug_enabled();

        // Legacy check: base_pairs[m][3] > 0 requires both:
        // 1. Pair is WC or wobble type (check_wc_wobble_pair only sets bpid for these)
        // 2. Geometric condition (dir_x > 0 && dir_y < 0 && dir_z < 0)
        // For non-WC/wobble pairs, bpid stays at -1 even if geometry matches.
        using x3dna::core::BasePairType;
        if (pair.type() != BasePairType::WATSON_CRICK &&
            pair.type() != BasePairType::WOBBLE) {
            if (debug) {
                std::cerr << "[has_positive_bpid] pair(" << pair.residue_idx1() << ","
                          << pair.residue_idx2() << ") type=" << pair.bp_type()
                          << " -> 0 (not WC/wobble)" << std::endl;
            }
            return false;
        }

        if (!pair.frame1().has_value() || !pair.frame2().has_value()) {
            return false;
        }

        // Copy to avoid dangling reference (optional returns by value)
        auto f1 = pair.frame1().value();
        auto f2 = pair.frame2().value();

        // Legacy: dir_x = dot(&orien[i][0], &orien[j][0])
        //         dir_y = dot(&orien[i][3], &orien[j][3])
        //         dir_z = dot(&orien[i][6], &orien[j][6])
        double dir_x = f1.x_axis().dot(f2.x_axis());
        double dir_y = f1.y_axis().dot(f2.y_axis());
        double dir_z = f1.z_axis().dot(f2.z_axis());

        bool result = (dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0);

        if (debug) {
            std::cerr << "[has_positive_bpid] pair(" << pair.residue_idx1() << "," << pair.residue_idx2() << ")"
                      << " type=" << pair.bp_type()
                      << " dir_x=" << dir_x << " dir_y=" << dir_y << " dir_z=" << dir_z
                      << " -> " << result << std::endl;
        }

        // Legacy condition: dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0
        return result;
    }
}

HelixOrganizer::HelixOrganizer(const Config& config) : config_(config) {}

// =============================================================================
// Geometry helpers
// =============================================================================

geometry::Vector3D HelixOrganizer::get_pair_origin(const core::BasePair& pair) const {
    // Legacy uses the average of both base origins (morg) for pair distance calculations
    // See refs_right_left() in cmn_fncs.c: morg[j] = (org[base1][j] + org[base2][j]) / 2
    // Note: organize() validates all pairs have both frames
    const auto& o1 = pair.frame1()->origin();
    const auto& o2 = pair.frame2()->origin();
    return geometry::Vector3D(
        (o1.x() + o2.x()) / 2.0,
        (o1.y() + o2.y()) / 2.0,
        (o1.z() + o2.z()) / 2.0
    );
}

geometry::Vector3D HelixOrganizer::get_pair_z_axis(const core::BasePair& pair) const {
    // Legacy uses average z-axis from both frames (or difference if they point opposite)
    // This matches bp_context: (d <= 0.0) ? ddxyz(z2, z1, zave) : sumxyz(z2, z1, zave)
    // Note: organize() validates all pairs have both frames
    auto z1 = pair.frame1()->z_axis();
    auto z2 = pair.frame2()->z_axis();
    double d = z1.dot(z2);

    geometry::Vector3D zave = (d <= 0.0) ? (z2 - z1) : (z2 + z1);
    zave.normalize();
    return zave;
}

geometry::Vector3D HelixOrganizer::get_frame_z(const core::BasePair& pair, bool swapped) const {
    // Note: organize() validates all pairs have both frames
    return swapped ? pair.frame2()->z_axis() : pair.frame1()->z_axis();
}

StrandResidues HelixOrganizer::get_strand_residues(
    const core::BasePair& pair, bool swapped) const {
    // BasePair stores 0-based indices normalized to (smaller, larger).
    // finding_order_swapped() indicates if original finding order was (larger, smaller).
    // Legacy code uses the ORIGINAL finding order for strand assignments.
    //
    // To match legacy:
    // - Start with original finding order (apply finding_order_swapped to restore it)
    // - Then apply the five2three swap flag
    // This is equivalent to XOR of finding_order_swapped and swapped.
    bool use_reversed_order = (pair.finding_order_swapped() != swapped);

    if (use_reversed_order) {
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
    // Note: organize() validates all pairs have both frames
    auto xm = (pair_m.frame1()->x_axis() + pair_m.frame2()->x_axis()).normalized();
    auto xn = (pair_n.frame1()->x_axis() + pair_n.frame2()->x_axis()).normalized();
    return dot2ang(xm.dot(xn));
}

double HelixOrganizer::wcbp_zdir(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                  bool swap_m, bool swap_n) const {
    // Get z-direction vectors based on swap status
    // Note: organize() validates all pairs have both frames
    // When swapped: z = z1 - z2, otherwise z = z2 - z1
    auto zm = swap_m ? (pair_m.frame1()->z_axis() - pair_m.frame2()->z_axis())
                     : (pair_m.frame2()->z_axis() - pair_m.frame1()->z_axis());
    auto zn = swap_n ? (pair_n.frame1()->z_axis() - pair_n.frame2()->z_axis())
                     : (pair_n.frame2()->z_axis() - pair_n.frame1()->z_axis());

    return zm.normalized().dot(zn.normalized());
}

// =============================================================================
// Five2three sub-functions
// =============================================================================

void HelixOrganizer::first_step(const std::vector<core::BasePair>& pairs,
                                const BackboneData& backbone,
                                std::vector<size_t>& pair_order,
                                const HelixSegment& helix,
                                std::vector<bool>& swapped) const {
    bool debug = is_five2three_debug_enabled();

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

    if (debug) {
        std::cerr << "[first_step] first_pair=" << first_pair << " ("
                  << pairs[first_pair].residue_idx1()+1 << "," << pairs[first_pair].residue_idx2()+1 << ")"
                  << " second_pair=" << second_pair << " ("
                  << pairs[second_pair].residue_idx1()+1 << "," << pairs[second_pair].residue_idx2()+1 << ")"
                  << " res_m.s1=" << res_m.strand1 << " res_n.s1=" << res_n.strand1
                  << " link=" << static_cast<int>(link) << std::endl;
    }

    if (link == LinkDirection::Reverse) {
        // Reverse linkage - swap the first pair
        swapped[first_pair] = !swapped[first_pair];
        if (debug) std::cerr << "[first_step] -> Reverse linkage, swapped first pair" << std::endl;
    } else if (link == LinkDirection::None) {
        if (debug) std::cerr << "[first_step] -> No linkage, reversing helix" << std::endl;
        // No linkage - reverse the entire helix segment and try again
        std::reverse(pair_order.begin() + helix.start_idx,
                     pair_order.begin() + helix.end_idx + 1);

        // Update pair indices after reversal
        first_pair = pair_order[pos];
        second_pair = pair_order[pos + 1];

        res_m = get_strand_residues(pairs[first_pair], swapped[first_pair]);
        res_n = get_strand_residues(pairs[second_pair], swapped[second_pair]);

        link = check_linkage(res_m.strand1, res_n.strand1, backbone);

        if (debug) {
            std::cerr << "[first_step] After reversal: first_pair=" << first_pair << " ("
                      << pairs[first_pair].residue_idx1()+1 << "," << pairs[first_pair].residue_idx2()+1 << ")"
                      << " res_m.s1=" << res_m.strand1 << " res_n.s1=" << res_n.strand1
                      << " link=" << static_cast<int>(link) << std::endl;
        }

        if (link == LinkDirection::Reverse) {
            swapped[first_pair] = !swapped[first_pair];
            if (debug) std::cerr << "[first_step] -> After reversal: Reverse linkage, swapped first pair" << std::endl;
        } else if (link == LinkDirection::None) {
            // Still no linkage - undo the reversal
            std::reverse(pair_order.begin() + helix.start_idx,
                         pair_order.begin() + helix.end_idx + 1);
            if (debug) std::cerr << "[first_step] -> Still no linkage, undoing reversal" << std::endl;
        } else {
            if (debug) std::cerr << "[first_step] -> After reversal: Forward linkage" << std::endl;
        }
    } else {
        if (debug) std::cerr << "[first_step] -> Forward linkage, no action" << std::endl;
    }
}

bool HelixOrganizer::wc_bporien(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                 bool swap_m, bool swap_n,
                                 const BackboneData& backbone) const {
    // Legacy check: only apply to Watson-Crick pairs (base_pairs[m][3] > 0)
    //
    // Legacy's bpid is set in calculate_more_bppars:
    //   1. Default bpid = -1
    //   2. If dir_x > 0 && dir_y < 0 && dir_z < 0, call check_wc_wobble_pair
    //   3. check_wc_wobble_pair sets bpid = 1 (wobble) or 2 (WC) based on geometry
    //
    // The wc_bporien check requires: base_pairs[m][3] > 0 && base_pairs[n][3] > 0
    // So we need to check the geometric condition, not just the bp_type string.
    bool debug = is_five2three_debug_enabled();

    if (!has_positive_bpid(pair_m) || !has_positive_bpid(pair_n)) {
        if (debug) {
            std::cerr << "[MODERN wc_bporien] SKIP: pair_m or pair_n has non-positive bpid" << std::endl;
        }
        return false;
    }

    // Check if pairs need swap based on z-direction alignment
    // Note: organize() validates all pairs have both frames
    auto res_m = get_strand_residues(pair_m, swap_m);
    auto res_n = get_strand_residues(pair_n, swap_n);

    double xang = wcbp_xang(pair_m, pair_n);
    auto link_s1 = check_linkage(res_m.strand1, res_n.strand1, backbone);
    auto link_s2 = check_linkage(res_m.strand2, res_n.strand2, backbone);

    if (debug) {
        std::cerr << "[MODERN wc_bporien] res_m=(" << pair_m.residue_idx1() << "," << pair_m.residue_idx2() << ")"
                  << " res_n=(" << pair_n.residue_idx1() << "," << pair_n.residue_idx2() << ")"
                  << " swap_m=" << swap_m << " swap_n=" << swap_n
                  << " xang=" << xang << " link_s1=" << static_cast<int>(link_s1)
                  << " link_s2=" << static_cast<int>(link_s2) << std::endl;
    }

    // If x-angle is too large or backbone is linked, don't swap
    if (xang > config_.end_stack_xang ||
        link_s1 != LinkDirection::None ||
        link_s2 != LinkDirection::None) {
        if (debug) {
            std::cerr << "[MODERN wc_bporien] -> false (early exit: xang>" << config_.end_stack_xang
                      << " or has linkage)" << std::endl;
        }
        return false;
    }

    // Check z-direction alignment
    double zdir_normal = wcbp_zdir(pair_m, pair_n, swap_m, swap_n);
    double zdir_swapped = wcbp_zdir(pair_m, pair_n, swap_m, !swap_n);

    bool result = zdir_normal < 0.0 && zdir_swapped > 0.0;
    if (debug) {
        std::cerr << "[MODERN wc_bporien] zdir_normal=" << zdir_normal
                  << " zdir_swapped=" << zdir_swapped << " -> " << result << std::endl;
    }

    return result;  // Need to swap pair_n
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

    return (di1_i2 > 0.0 && di1_j2 > 0.0 && di1_i2 > di1_j2) &&
           (dj1_i2 > 0.0 && dj1_j2 > 0.0 && dj1_j2 > dj1_i2);
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
    const bool has_any_linkage =
        check_linkage(res_m.strand1, res_n.strand1, backbone) != LinkDirection::None ||
        check_linkage(res_m.strand2, res_n.strand2, backbone) != LinkDirection::None ||
        check_linkage(res_m.strand1, res_n.strand2, backbone) != LinkDirection::None ||
        check_linkage(res_m.strand2, res_n.strand1, backbone) != LinkDirection::None;
    if (has_any_linkage) {
        return false;
    }

    // Get frames based on swap status
    // Note: organize() validates all pairs have both frames
    const auto frame_m1 = swap_m ? *pair_m.frame2() : *pair_m.frame1();
    const auto frame_m2 = swap_m ? *pair_m.frame1() : *pair_m.frame2();
    const auto frame_n1 = swap_n ? *pair_n.frame2() : *pair_n.frame1();
    const auto frame_n2 = swap_n ? *pair_n.frame1() : *pair_n.frame2();

    // Check same-strand alignment (m1 with n1, m2 with n2)
    const auto align1 = compute_frame_alignment(frame_m1, frame_n1);
    const auto align2 = compute_frame_alignment(frame_m2, frame_n2);

    if (align1.is_aligned() && align2.is_aligned()) {
        return false;
    }

    // Check cross-strand alignment (m1 with n2, m2 with n1)
    const auto cross1 = compute_frame_alignment(frame_m1, frame_n2);
    const auto cross2 = compute_frame_alignment(frame_m2, frame_n1);

    if (!align1.is_aligned() && !align2.is_aligned()) {
        // Neither same-strand pair is aligned - check if cross-aligned
        return cross1.is_aligned() || cross2.is_aligned();
    }

    // Legacy logic: Compare SPECIFIC frame pairs, not total sums
    // Compare aligned angle sum vs cross angle sum for the appropriate pairs
    if (align1.is_aligned() && cross1.is_aligned()) {
        return align1.angle_sum() > cross1.angle_sum();
    }
    if (align1.is_aligned() && cross2.is_aligned()) {
        return align1.angle_sum() > cross2.angle_sum();
    }
    if (align2.is_aligned() && cross1.is_aligned()) {
        return align2.angle_sum() > cross1.angle_sum();
    }
    if (align2.is_aligned() && cross2.is_aligned()) {
        return align2.angle_sum() > cross2.angle_sum();
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

    bool debug = is_five2three_debug_enabled();

    DirectionCounts dir;

    // First pass: compute direction counts
    for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
        size_t idx_m = pair_order[pos];
        size_t idx_n = pair_order[pos + 1];

        auto res_m = get_strand_residues(pairs[idx_m], swapped[idx_m]);
        auto res_n = get_strand_residues(pairs[idx_n], swapped[idx_n]);

        auto link1 = check_linkage(res_m.strand1, res_n.strand1, backbone);
        auto link2 = check_linkage(res_m.strand2, res_n.strand2, backbone);

        update_direction_count(link1, dir.strand1_forward, dir.strand1_reverse, dir.strand1_none);
        update_direction_count(link2, dir.strand2_forward, dir.strand2_reverse, dir.strand2_none);
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

    // Check strand1 direction pattern
    const bool is_strand1_forward_only = dir.strand1_forward && !dir.strand1_reverse;
    if (!is_strand1_forward_only) {
        return dir;
    }

    // Determine helix orientation type
    const bool is_anti_parallel = !dir.strand2_forward && dir.strand2_reverse;
    const bool is_parallel = dir.strand2_forward && !dir.strand2_reverse;

    // Helper to flip all swapped values in helix
    auto flip_helix_swaps = [&]() {
        for (size_t pos = helix.start_idx; pos <= helix.end_idx; ++pos) {
            size_t idx = pair_order[pos];
            swapped[idx] = !swapped[idx];
        }
    };

    if (is_anti_parallel) {
        const bool needs_flip_and_reverse = res_first.strand1 > res_last.strand2;
        if (debug) {
            std::cerr << "[check_direction] ANTI-PARALLEL: first.s1=" << res_first.strand1
                      << " last.s2=" << res_last.strand2
                      << " check=" << (needs_flip_and_reverse ? "YES" : "NO") << std::endl;
        }
        if (needs_flip_and_reverse) {
            if (debug) std::cerr << "[check_direction] -> Flipping all swaps and reversing" << std::endl;
            flip_helix_swaps();
            std::reverse(pair_order.begin() + helix.start_idx,
                         pair_order.begin() + helix.end_idx + 1);
        }
        return dir;
    }

    if (is_parallel) {
        helix.is_parallel = true;
        const bool needs_flip = res_first.strand1 > res_first.strand2;
        if (needs_flip) {
            flip_helix_swaps();
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
        const bool anti_p =
            (direction.strand1_forward > direction.strand1_reverse) &&
            (direction.strand2_forward < direction.strand2_reverse);
        const bool parallel =
            (direction.strand1_forward > direction.strand1_reverse) &&
            (direction.strand2_forward > direction.strand2_reverse);

        helix.is_parallel = parallel;

        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];

            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];

            auto res_m = get_strand_residues(pair_m, swapped[idx_m]);
            auto res_n = get_strand_residues(pair_n, swapped[idx_n]);

            // Check for strand2 linkage requiring swap
            const auto link_strand1 = check_linkage(res_m.strand1, res_n.strand1, backbone);
            const auto link_strand2 = check_linkage(res_m.strand2, res_n.strand2, backbone);

            const bool no_strand1_link = (link_strand1 == LinkDirection::None);
            const bool needs_strand2_swap =
                (anti_p && link_strand2 == LinkDirection::Forward) ||
                (parallel && link_strand2 == LinkDirection::Reverse);

            if (no_strand1_link && needs_strand2_swap) {
                swapped[idx_n] = !swapped[idx_n];
            }

            // Re-get residues after potential swap
            res_n = get_strand_residues(pair_n, swapped[idx_n]);

            // Check for cross-strand swaps
            const bool no_same_strand_links =
                check_linkage(res_m.strand1, res_n.strand1, backbone) == LinkDirection::None &&
                check_linkage(res_m.strand2, res_n.strand2, backbone) == LinkDirection::None;

            if (!no_same_strand_links) {
                continue;
            }

            // Determine which pair needs swapping based on cross-strand linkage
            const bool should_swap_m =
                (anti_p && check_linkage(res_m.strand2, res_n.strand1, backbone) == LinkDirection::Forward) ||
                (parallel && check_linkage(res_m.strand1, res_n.strand2, backbone) == LinkDirection::Reverse);

            const bool should_swap_n =
                (anti_p && check_linkage(res_m.strand1, res_n.strand2, backbone) == LinkDirection::Forward) ||
                (parallel && check_linkage(res_m.strand2, res_n.strand1, backbone) == LinkDirection::Reverse);

            if (should_swap_m) {
                swapped[idx_m] = !swapped[idx_m];
            } else if (should_swap_n) {
                swapped[idx_n] = !swapped[idx_n];
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

        // Legacy behavior: if no neighbors within helix_break, pair is an isolated endpoint
        // with NO stored neighbors (end_list only stores the endpoint itself)
        if (neighbors.empty() || neighbors[0].first > config_.helix_break) {
            context[i].is_endpoint = true;
            continue;  // Don't store neighbor1 - matches legacy line 963-965
        }

        context[i].neighbor1 = neighbors[0].second;
        context[i].dist1 = neighbors[0].first;

        // Check backbone connectivity to neighbor1
        context[i].has_backbone_link1 = are_pairs_backbone_connected(
            pairs[i], pairs[neighbors[0].second], backbone);

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
            const bool both_opposite_from_n1 = are_on_opposite_z_sides(d1, d2) && are_on_opposite_z_sides(d1, d3);
            const bool second_has_larger_z_dist = std::abs(d2) > std::abs(d3);
            if (both_opposite_from_n1 && second_has_larger_z_dist) {
                std::swap(neighbors[1], neighbors[2]);
            }
        }

        for (size_t k = 1; k < neighbors.size(); ++k) {
            if (neighbors[k].first > config_.helix_break) break;

            auto vk = get_pair_origin(pairs[neighbors[k].second]) - org_i;
            double dk = z_i.dot(vk);

            if (are_on_opposite_z_sides(d1, dk)) {
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
                if (are_on_opposite_z_sides(d1, d2) && dist_n1_n2 <= config_.helix_break) {
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

    bool debug = is_five2three_debug_enabled();

    std::vector<size_t> endpoints;

    for (size_t i = 0; i < context.size(); ++i) {
        if (context[i].is_endpoint) {
            endpoints.push_back(i);
        }
    }

    if (endpoints.empty() && !context.empty()) {
        endpoints.push_back(0);
    }

    if (debug) {
        std::cerr << "[find_endpoints] Endpoints found: ";
        for (size_t ep : endpoints) {
            std::cerr << ep << " ";
        }
        std::cerr << std::endl;
        for (size_t i = 0; i < context.size(); ++i) {
            std::cerr << "[context] pair " << i << ": ep=" << context[i].is_endpoint
                      << " n1=" << (context[i].neighbor1.has_value() ? std::to_string(context[i].neighbor1.value()) : "-")
                      << " n2=" << (context[i].neighbor2.has_value() ? std::to_string(context[i].neighbor2.value()) : "-")
                      << std::endl;
        }
    }

    return endpoints;
}

std::pair<std::vector<size_t>, std::vector<HelixSegment>> HelixOrganizer::locate_helices(
    const std::vector<PairContext>& context,
    const std::vector<size_t>& endpoints,
    const BackboneData& backbone,
    size_t num_pairs) const {

    bool debug = is_five2three_debug_enabled();

    std::vector<size_t> pair_order;
    std::vector<HelixSegment> helices;
    std::vector<bool> visited(num_pairs, false);

    (void)backbone;  // Not used in traverse - backbone connectivity checked in calculate_context
    pair_order.reserve(num_pairs);

    // Legacy locate_helix algorithm (find_pair.c lines 1029-1082):
    // For each endpoint, add endpoint + neighbors from end_list, then traverse
    for (size_t ep : endpoints) {
        if (debug) {
            std::cerr << "[locate_helices] Processing endpoint " << ep << std::endl;
        }
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
            if (debug) {
                std::cerr << "[locate_helices] Created helix " << helices.size()
                          << " (pos " << helix.start_idx << "-" << helix.end_idx << "): ";
                for (size_t p = helix.start_idx; p <= helix.end_idx; ++p) {
                    std::cerr << pair_order[p] << " ";
                }
                std::cerr << std::endl;
            }
        }
    }

    // Handle any leftover pairs not reached from endpoints
    // Legacy behavior (lines 1120-1128): Put ALL leftover pairs into ONE helix region
    // This handles "complicated structures" with isolated pairs
    std::vector<size_t> leftover;
    for (size_t i = 0; i < num_pairs; ++i) {
        if (!visited[i]) {
            leftover.push_back(i);
        }
    }

    if (!leftover.empty()) {
        // All leftover pairs go into a single helix (matching legacy behavior)
        HelixSegment helix;
        helix.start_idx = pair_order.size();
        for (size_t idx : leftover) {
            pair_order.push_back(idx);
        }
        helix.end_idx = pair_order.size() - 1;
        helices.push_back(helix);
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

    bool debug = is_five2three_debug_enabled();

    swapped.resize(pairs.size(), false);

    if (backbone.empty()) {
        return;
    }

    // Process each helix
    size_t helix_num = 0;
    for (auto& helix : helices) {
        ++helix_num;
        if (helix.start_idx > helix.end_idx) continue;

        if (debug) {
            std::cerr << "\n=== HELIX " << helix_num << " (pairs "
                      << helix.start_idx << "-" << helix.end_idx << ") ===" << std::endl;
        }

        // STEP 1: first_step - set initial strand assignment
        first_step(pairs, backbone, pair_order, helix, swapped);

        if (debug) {
            std::cerr << "[STEP1 first_step] After first_step:" << std::endl;
            for (size_t pos = helix.start_idx; pos <= helix.end_idx; ++pos) {
                size_t idx = pair_order[pos];
                std::cerr << "  pos " << pos << " pair_idx=" << idx
                          << " (" << pairs[idx].residue_idx1()+1 << "," << pairs[idx].residue_idx2()+1 << ")"
                          << " swap=" << swapped[idx] << std::endl;
            }
        }

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

            if (debug) {
                std::cerr << "[STEP2 pass1] pos " << pos << "->" << (pos+1)
                          << " m=(" << pair_m.residue_idx1()+1 << "," << pair_m.residue_idx2()+1 << ")"
                          << " n=(" << pair_n.residue_idx1()+1 << "," << pair_n.residue_idx2()+1 << ")"
                          << " swap_m=" << swapped[idx_m] << " swap_n=" << swapped[idx_n]
                          << " rev_wc=" << rev_wc << " rev_o3d=" << rev_o3d
                          << " rev_csc=" << rev_csc << " rev_oth=" << rev_oth << std::endl;
            }

            // Apply swap based on checks
            if (rev_wc) {
                swapped[idx_n] = !swapped[idx_n];
                if (debug) std::cerr << "  -> rev_wc: toggled swap_n to " << swapped[idx_n] << std::endl;
            } else if (rev_o3d || rev_csc || rev_oth) {
                swapped[idx_n] = !swapped[idx_n];
                if (debug) std::cerr << "  -> rev_o3d/csc/oth: toggled swap_n to " << swapped[idx_n] << std::endl;
            }

            // Check strand 1 direction
            bool rev_s1 = chain1dir(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            if (debug) std::cerr << "  chain1dir=" << rev_s1 << std::endl;

            if (rev_s1) {
                swapped[idx_n] = !swapped[idx_n];
                if (debug) std::cerr << "  -> rev_s1: toggled swap_n to " << swapped[idx_n] << std::endl;
            }
        }

        // STEP 3: Second pass - re-check WC orientation
        if (debug) std::cerr << "[STEP3 pass2] Second pass WC check:" << std::endl;
        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];

            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];

            bool rev_wc = wc_bporien(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            if (debug) {
                std::cerr << "  pos " << pos << " rev_wc=" << rev_wc
                          << " (swap_m=" << swapped[idx_m] << " swap_n=" << swapped[idx_n] << ")" << std::endl;
            }
            if (rev_wc) {
                swapped[idx_m] = !swapped[idx_m];
                if (debug) std::cerr << "  -> toggled swap_m to " << swapped[idx_m] << std::endl;
            }
        }

        if (debug) {
            std::cerr << "[After pass2] Swap state:" << std::endl;
            for (size_t pos = helix.start_idx; pos <= helix.end_idx; ++pos) {
                size_t idx = pair_order[pos];
                std::cerr << "  pos " << pos << " pair_idx=" << idx
                          << " (" << pairs[idx].residue_idx1()+1 << "," << pairs[idx].residue_idx2()+1 << ")"
                          << " swap=" << swapped[idx] << std::endl;
            }
        }

        // STEP 4: check_direction - count backbone linkage directions and apply fixes
        DirectionCounts direction = check_direction(pairs, backbone, pair_order, helix, swapped);

        if (debug) {
            std::cerr << "[STEP4 check_direction] s1_fwd=" << direction.strand1_forward
                      << " s1_rev=" << direction.strand1_reverse << " s1_none=" << direction.strand1_none
                      << " s2_fwd=" << direction.strand2_forward << " s2_rev=" << direction.strand2_reverse
                      << " s2_none=" << direction.strand2_none << std::endl;
        }

        // STEP 5: check_strand2 - additional corrections based on direction
        check_strand2(pairs, backbone, pair_order, helix, swapped, direction);

        if (debug) {
            std::cerr << "[After check_strand2] Swap state:" << std::endl;
            for (size_t pos = helix.start_idx; pos <= helix.end_idx; ++pos) {
                size_t idx = pair_order[pos];
                std::cerr << "  pos " << pos << " pair_idx=" << idx
                          << " (" << pairs[idx].residue_idx1()+1 << "," << pairs[idx].residue_idx2()+1 << ")"
                          << " swap=" << swapped[idx] << std::endl;
            }
        }

        // STEP 6: check_direction AGAIN (legacy line 1361 - at end of check_strand2)
        // This recomputes direction with updated swaps and may apply additional corrections
        check_direction(pairs, backbone, pair_order, helix, swapped);

        if (debug) {
            std::cerr << "[FINAL] Helix " << helix_num << " swap state:" << std::endl;
            for (size_t pos = helix.start_idx; pos <= helix.end_idx; ++pos) {
                size_t idx = pair_order[pos];
                std::cerr << "  pos " << pos << " pair_idx=" << idx
                          << " (" << pairs[idx].residue_idx1()+1 << "," << pairs[idx].residue_idx2()+1 << ")"
                          << " swap=" << swapped[idx] << std::endl;
            }
        }
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

    // Validate all pairs have both frames - required for helix organization
    for (size_t i = 0; i < pairs.size(); ++i) {
        if (!pairs[i].frame1().has_value() || !pairs[i].frame2().has_value()) {
            throw std::invalid_argument(
                "HelixOrganizer::organize: pair " + std::to_string(i) +
                " missing frame(s). All pairs must have valid frames.");
        }
    }

    if (pairs.size() == 1) {
        result.pair_order = {0};
        result.helices = {{0, 0, false, false, false, false, {}}};  // Added direction field
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

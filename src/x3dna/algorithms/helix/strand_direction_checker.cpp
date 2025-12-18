/**
 * @file strand_direction_checker.cpp
 * @brief Implementation of strand direction checking
 */

#include <x3dna/algorithms/helix/strand_direction_checker.hpp>
#include <x3dna/config/config_manager.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace x3dna::algorithms::helix {

namespace {
constexpr double PI = 3.14159265358979323846;

/**
 * @brief Check if five2three debugging is enabled
 */
bool is_debug_enabled() {
    static bool initialized = false;
    static bool debug_enabled = false;
    if (!initialized) {
        auto& cfg = config::ConfigManager::instance();
        cfg.init_debug_from_environment();
        debug_enabled = cfg.debug_config().debug_five2three;
        initialized = true;
    }
    return debug_enabled;
}

/**
 * @brief Convert dot product to angle in degrees
 */
double dot2ang(double d) {
    if (d > 1.0)
        d = 1.0;
    if (d < -1.0)
        d = -1.0;
    return std::acos(d) * 180.0 / PI;
}

/**
 * @brief Result of frame axis alignment comparison
 */
struct FrameAlignment {
    double dot_x;
    double dot_y;
    double dot_z;

    [[nodiscard]] bool is_aligned() const {
        return dot_x > 0.0 && dot_y > 0.0 && dot_z > 0.0;
    }

    [[nodiscard]] double angle_sum() const {
        return dot2ang(dot_x) + dot2ang(dot_y) + dot2ang(dot_z);
    }
};

/**
 * @brief Compute alignment between two reference frames
 */
FrameAlignment compute_frame_alignment(const core::ReferenceFrame& frame1, const core::ReferenceFrame& frame2) {
    return FrameAlignment{frame1.x_axis().dot(frame2.x_axis()), frame1.y_axis().dot(frame2.y_axis()),
                          frame1.z_axis().dot(frame2.z_axis())};
}
} // namespace

double StrandDirectionChecker::wcbp_xang(const core::BasePair& pair_m, const core::BasePair& pair_n) const {
    auto xm = (pair_m.frame1()->x_axis() + pair_m.frame2()->x_axis()).normalized();
    auto xn = (pair_n.frame1()->x_axis() + pair_n.frame2()->x_axis()).normalized();
    return dot2ang(xm.dot(xn));
}

double StrandDirectionChecker::wcbp_zdir(const core::BasePair& pair_m, const core::BasePair& pair_n, bool swap_m,
                                         bool swap_n) const {
    auto zm = swap_m ? (pair_m.frame1()->z_axis() - pair_m.frame2()->z_axis())
                     : (pair_m.frame2()->z_axis() - pair_m.frame1()->z_axis());
    auto zn = swap_n ? (pair_n.frame1()->z_axis() - pair_n.frame2()->z_axis())
                     : (pair_n.frame2()->z_axis() - pair_n.frame1()->z_axis());
    return zm.normalized().dot(zn.normalized());
}

bool StrandDirectionChecker::has_positive_bpid(const core::BasePair& pair) const {
    bool debug = is_debug_enabled();

    using core::BasePairType;
    if (pair.type() != BasePairType::WATSON_CRICK && pair.type() != BasePairType::WOBBLE) {
        if (debug) {
            std::cerr << "[has_positive_bpid] pair(" << pair.residue_idx1() << "," << pair.residue_idx2()
                      << ") type=" << pair.bp_type() << " -> 0 (not WC/wobble)" << std::endl;
        }
        return false;
    }

    if (!pair.frame1().has_value() || !pair.frame2().has_value()) {
        return false;
    }

    auto f1 = pair.frame1().value();
    auto f2 = pair.frame2().value();

    double dir_x = f1.x_axis().dot(f2.x_axis());
    double dir_y = f1.y_axis().dot(f2.y_axis());
    double dir_z = f1.z_axis().dot(f2.z_axis());

    bool result = (dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0);

    if (debug) {
        std::cerr << "[has_positive_bpid] pair(" << pair.residue_idx1() << "," << pair.residue_idx2() << ")"
                  << " type=" << pair.bp_type() << " dir_x=" << dir_x << " dir_y=" << dir_y << " dir_z=" << dir_z
                  << " -> " << result << std::endl;
    }

    return result;
}

void StrandDirectionChecker::first_step(const std::vector<core::BasePair>& pairs, const BackboneData& backbone,
                                        std::vector<size_t>& pair_order, const HelixSegment& helix,
                                        std::vector<bool>& swapped) const {
    bool debug = is_debug_enabled();

    if (helix.end_idx <= helix.start_idx) {
        return;
    }

    size_t pos = helix.start_idx;
    size_t first_pair = pair_order[pos];
    size_t second_pair = pair_order[pos + 1];

    auto res_m = PairGeometryHelper::get_strand_residues(pairs[first_pair], swapped[first_pair]);
    auto res_n = PairGeometryHelper::get_strand_residues(pairs[second_pair], swapped[second_pair]);

    auto link = linkage_checker_.check_linkage(res_m.strand1, res_n.strand1, backbone);

    if (debug) {
        std::cerr << "[first_step] first_pair=" << first_pair << " (" << pairs[first_pair].residue_idx1() + 1 << ","
                  << pairs[first_pair].residue_idx2() + 1 << ")"
                  << " second_pair=" << second_pair << " (" << pairs[second_pair].residue_idx1() + 1 << ","
                  << pairs[second_pair].residue_idx2() + 1 << ")"
                  << " res_m.s1=" << res_m.strand1 << " res_n.s1=" << res_n.strand1
                  << " link=" << static_cast<int>(link) << std::endl;
    }

    if (link == LinkDirection::Reverse) {
        swapped[first_pair] = !swapped[first_pair];
        if (debug)
            std::cerr << "[first_step] -> Reverse linkage, swapped first pair" << std::endl;
    } else if (link == LinkDirection::None) {
        if (debug)
            std::cerr << "[first_step] -> No linkage, reversing helix" << std::endl;
        std::reverse(pair_order.begin() + helix.start_idx, pair_order.begin() + helix.end_idx + 1);

        first_pair = pair_order[pos];
        second_pair = pair_order[pos + 1];

        res_m = PairGeometryHelper::get_strand_residues(pairs[first_pair], swapped[first_pair]);
        res_n = PairGeometryHelper::get_strand_residues(pairs[second_pair], swapped[second_pair]);

        link = linkage_checker_.check_linkage(res_m.strand1, res_n.strand1, backbone);

        if (debug) {
            std::cerr << "[first_step] After reversal: first_pair=" << first_pair << " ("
                      << pairs[first_pair].residue_idx1() + 1 << "," << pairs[first_pair].residue_idx2() + 1 << ")"
                      << " res_m.s1=" << res_m.strand1 << " res_n.s1=" << res_n.strand1
                      << " link=" << static_cast<int>(link) << std::endl;
        }

        if (link == LinkDirection::Reverse) {
            swapped[first_pair] = !swapped[first_pair];
            if (debug)
                std::cerr << "[first_step] -> After reversal: Reverse linkage, swapped first pair" << std::endl;
        } else if (link == LinkDirection::None) {
            std::reverse(pair_order.begin() + helix.start_idx, pair_order.begin() + helix.end_idx + 1);
            if (debug)
                std::cerr << "[first_step] -> Still no linkage, undoing reversal" << std::endl;
        } else {
            if (debug)
                std::cerr << "[first_step] -> After reversal: Forward linkage" << std::endl;
        }
    } else {
        if (debug)
            std::cerr << "[first_step] -> Forward linkage, no action" << std::endl;
    }
}

bool StrandDirectionChecker::wc_bporien(const core::BasePair& pair_m, const core::BasePair& pair_n, bool swap_m,
                                        bool swap_n, const BackboneData& backbone) const {
    bool debug = is_debug_enabled();

    if (!has_positive_bpid(pair_m) || !has_positive_bpid(pair_n)) {
        if (debug) {
            std::cerr << "[MODERN wc_bporien] SKIP: pair_m or pair_n has non-positive bpid" << std::endl;
        }
        return false;
    }

    auto res_m = PairGeometryHelper::get_strand_residues(pair_m, swap_m);
    auto res_n = PairGeometryHelper::get_strand_residues(pair_n, swap_n);

    double xang = wcbp_xang(pair_m, pair_n);
    auto link_s1 = linkage_checker_.check_linkage(res_m.strand1, res_n.strand1, backbone);
    auto link_s2 = linkage_checker_.check_linkage(res_m.strand2, res_n.strand2, backbone);

    if (debug) {
        std::cerr << "[MODERN wc_bporien] res_m=(" << pair_m.residue_idx1() << "," << pair_m.residue_idx2() << ")"
                  << " res_n=(" << pair_n.residue_idx1() << "," << pair_n.residue_idx2() << ")"
                  << " swap_m=" << swap_m << " swap_n=" << swap_n << " xang=" << xang
                  << " link_s1=" << static_cast<int>(link_s1) << " link_s2=" << static_cast<int>(link_s2) << std::endl;
    }

    if (xang > config_.end_stack_xang || link_s1 != LinkDirection::None || link_s2 != LinkDirection::None) {
        if (debug) {
            std::cerr << "[MODERN wc_bporien] -> false (early exit: xang>" << config_.end_stack_xang
                      << " or has linkage)" << std::endl;
        }
        return false;
    }

    double zdir_normal = wcbp_zdir(pair_m, pair_n, swap_m, swap_n);
    double zdir_swapped = wcbp_zdir(pair_m, pair_n, swap_m, !swap_n);

    bool result = zdir_normal < 0.0 && zdir_swapped > 0.0;
    if (debug) {
        std::cerr << "[MODERN wc_bporien] zdir_normal=" << zdir_normal << " zdir_swapped=" << zdir_swapped << " -> "
                  << result << std::endl;
    }

    return result;
}

bool StrandDirectionChecker::check_o3dist(const core::BasePair& pair_m, const core::BasePair& pair_n, bool swap_m,
                                          bool swap_n, const BackboneData& backbone) const {
    auto res_m = PairGeometryHelper::get_strand_residues(pair_m, swap_m);
    auto res_n = PairGeometryHelper::get_strand_residues(pair_n, swap_n);

    double di1_i2 = linkage_checker_.o3_distance(res_m.strand1, res_n.strand1, backbone);
    double di1_j2 = linkage_checker_.o3_distance(res_m.strand1, res_n.strand2, backbone);
    double dj1_i2 = linkage_checker_.o3_distance(res_m.strand2, res_n.strand1, backbone);
    double dj1_j2 = linkage_checker_.o3_distance(res_m.strand2, res_n.strand2, backbone);

    return (di1_i2 > 0.0 && di1_j2 > 0.0 && di1_i2 > di1_j2) && (dj1_i2 > 0.0 && dj1_j2 > 0.0 && dj1_j2 > dj1_i2);
}

bool StrandDirectionChecker::check_schain(const core::BasePair& pair_m, const core::BasePair& pair_n, bool swap_m,
                                          bool swap_n, const BackboneData& backbone) const {
    auto res_m = PairGeometryHelper::get_strand_residues(pair_m, swap_m);
    auto res_n = PairGeometryHelper::get_strand_residues(pair_n, swap_n);

    bool no_same_strand = linkage_checker_.check_linkage(res_m.strand1, res_n.strand1, backbone) ==
                              LinkDirection::None &&
                          linkage_checker_.check_linkage(res_m.strand2, res_n.strand2, backbone) == LinkDirection::None;
    bool has_cross_strand = linkage_checker_.check_linkage(res_m.strand1, res_n.strand2, backbone) !=
                                LinkDirection::None ||
                            linkage_checker_.check_linkage(res_m.strand2, res_n.strand1, backbone) !=
                                LinkDirection::None;

    return no_same_strand && has_cross_strand;
}

bool StrandDirectionChecker::check_others(const core::BasePair& pair_m, const core::BasePair& pair_n, bool swap_m,
                                          bool swap_n, const BackboneData& backbone) const {
    auto res_m = PairGeometryHelper::get_strand_residues(pair_m, swap_m);
    auto res_n = PairGeometryHelper::get_strand_residues(pair_n, swap_n);

    const bool has_any_linkage = linkage_checker_.check_linkage(res_m.strand1, res_n.strand1, backbone) !=
                                     LinkDirection::None ||
                                 linkage_checker_.check_linkage(res_m.strand2, res_n.strand2, backbone) !=
                                     LinkDirection::None ||
                                 linkage_checker_.check_linkage(res_m.strand1, res_n.strand2, backbone) !=
                                     LinkDirection::None ||
                                 linkage_checker_.check_linkage(res_m.strand2, res_n.strand1, backbone) !=
                                     LinkDirection::None;
    if (has_any_linkage) {
        return false;
    }

    const auto frame_m1 = swap_m ? *pair_m.frame2() : *pair_m.frame1();
    const auto frame_m2 = swap_m ? *pair_m.frame1() : *pair_m.frame2();
    const auto frame_n1 = swap_n ? *pair_n.frame2() : *pair_n.frame1();
    const auto frame_n2 = swap_n ? *pair_n.frame1() : *pair_n.frame2();

    const auto align1 = compute_frame_alignment(frame_m1, frame_n1);
    const auto align2 = compute_frame_alignment(frame_m2, frame_n2);

    if (align1.is_aligned() && align2.is_aligned()) {
        return false;
    }

    const auto cross1 = compute_frame_alignment(frame_m1, frame_n2);
    const auto cross2 = compute_frame_alignment(frame_m2, frame_n1);

    if (!align1.is_aligned() && !align2.is_aligned()) {
        return cross1.is_aligned() || cross2.is_aligned();
    }

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

bool StrandDirectionChecker::chain1dir(const core::BasePair& pair_m, const core::BasePair& pair_n, bool swap_m,
                                       bool swap_n, const BackboneData& backbone) const {
    auto res_m = PairGeometryHelper::get_strand_residues(pair_m, swap_m);
    auto res_n = PairGeometryHelper::get_strand_residues(pair_n, swap_n);
    return linkage_checker_.check_linkage(res_m.strand1, res_n.strand1, backbone) == LinkDirection::Reverse;
}

DirectionCounts StrandDirectionChecker::check_direction(const std::vector<core::BasePair>& pairs,
                                                        const BackboneData& backbone, std::vector<size_t>& pair_order,
                                                        HelixSegment& helix, std::vector<bool>& swapped) const {

    bool debug = is_debug_enabled();

    DirectionCounts dir;

    for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
        size_t idx_m = pair_order[pos];
        size_t idx_n = pair_order[pos + 1];

        auto res_m = PairGeometryHelper::get_strand_residues(pairs[idx_m], swapped[idx_m]);
        auto res_n = PairGeometryHelper::get_strand_residues(pairs[idx_n], swapped[idx_n]);

        auto link1 = linkage_checker_.check_linkage(res_m.strand1, res_n.strand1, backbone);
        auto link2 = linkage_checker_.check_linkage(res_m.strand2, res_n.strand2, backbone);

        update_direction_count(link1, dir.strand1_forward, dir.strand1_reverse, dir.strand1_none);
        update_direction_count(link2, dir.strand2_forward, dir.strand2_reverse, dir.strand2_none);
    }

    bool mixed = (dir.strand1_forward && dir.strand1_reverse) || (dir.strand2_forward && dir.strand2_reverse);
    if (mixed) {
        helix.has_mixed_direction = true;
        return dir;
    }

    if (dir.strand1_forward + dir.strand1_reverse + dir.strand2_forward + dir.strand2_reverse == 0) {
        return dir;
    }

    size_t first_pair_idx = pair_order[helix.start_idx];
    size_t last_pair_idx = pair_order[helix.end_idx];
    auto res_first = PairGeometryHelper::get_strand_residues(pairs[first_pair_idx], swapped[first_pair_idx]);
    auto res_last = PairGeometryHelper::get_strand_residues(pairs[last_pair_idx], swapped[last_pair_idx]);

    if (dir.strand1_none || dir.strand2_none) {
        helix.has_break = true;
    }

    const bool is_strand1_forward_only = dir.strand1_forward && !dir.strand1_reverse;
    if (!is_strand1_forward_only) {
        return dir;
    }

    const bool is_anti_parallel = !dir.strand2_forward && dir.strand2_reverse;
    const bool is_parallel = dir.strand2_forward && !dir.strand2_reverse;

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
                      << " last.s2=" << res_last.strand2 << " check=" << (needs_flip_and_reverse ? "YES" : "NO")
                      << std::endl;
        }
        if (needs_flip_and_reverse) {
            if (debug)
                std::cerr << "[check_direction] -> Flipping all swaps and reversing" << std::endl;
            flip_helix_swaps();
            std::reverse(pair_order.begin() + helix.start_idx, pair_order.begin() + helix.end_idx + 1);
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

void StrandDirectionChecker::check_strand2(const std::vector<core::BasePair>& pairs, const BackboneData& backbone,
                                           const std::vector<size_t>& pair_order, HelixSegment& helix,
                                           std::vector<bool>& swapped, const DirectionCounts& direction) const {

    bool mixed_direction = (direction.strand1_forward && direction.strand1_reverse) ||
                           (direction.strand2_forward && direction.strand2_reverse);

    if (!mixed_direction) {
        if (direction.strand1_forward + direction.strand1_reverse + direction.strand2_forward +
                direction.strand2_reverse ==
            0) {
            return;
        }

        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];

            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];

            if (wc_bporien(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone)) {
                continue;
            }

            auto res_m = PairGeometryHelper::get_strand_residues(pair_m, swapped[idx_m]);
            auto res_n = PairGeometryHelper::get_strand_residues(pair_n, swapped[idx_n]);

            bool no_same_strand = linkage_checker_.check_linkage(res_m.strand1, res_n.strand1, backbone) ==
                                      LinkDirection::None &&
                                  linkage_checker_.check_linkage(res_m.strand2, res_n.strand2, backbone) ==
                                      LinkDirection::None;
            auto cross_link = linkage_checker_.check_linkage(res_m.strand1, res_n.strand2, backbone);
            auto cross_link2 = linkage_checker_.check_linkage(res_m.strand2, res_n.strand1, backbone);

            if (no_same_strand && (cross_link == LinkDirection::Forward ||
                                   (cross_link != LinkDirection::None && cross_link2 != LinkDirection::None))) {
                swapped[idx_n] = !swapped[idx_n];
            }
        }
    } else {
        const bool anti_p = (direction.strand1_forward > direction.strand1_reverse) &&
                            (direction.strand2_forward < direction.strand2_reverse);
        const bool parallel = (direction.strand1_forward > direction.strand1_reverse) &&
                              (direction.strand2_forward > direction.strand2_reverse);

        helix.is_parallel = parallel;

        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];

            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];

            auto res_m = PairGeometryHelper::get_strand_residues(pair_m, swapped[idx_m]);
            auto res_n = PairGeometryHelper::get_strand_residues(pair_n, swapped[idx_n]);

            const auto link_strand1 = linkage_checker_.check_linkage(res_m.strand1, res_n.strand1, backbone);
            const auto link_strand2 = linkage_checker_.check_linkage(res_m.strand2, res_n.strand2, backbone);

            const bool no_strand1_link = (link_strand1 == LinkDirection::None);
            const bool needs_strand2_swap = (anti_p && link_strand2 == LinkDirection::Forward) ||
                                            (parallel && link_strand2 == LinkDirection::Reverse);

            if (no_strand1_link && needs_strand2_swap) {
                swapped[idx_n] = !swapped[idx_n];
            }

            res_n = PairGeometryHelper::get_strand_residues(pair_n, swapped[idx_n]);

            const bool no_same_strand_links = linkage_checker_.check_linkage(res_m.strand1, res_n.strand1, backbone) ==
                                                  LinkDirection::None &&
                                              linkage_checker_.check_linkage(res_m.strand2, res_n.strand2, backbone) ==
                                                  LinkDirection::None;

            if (!no_same_strand_links) {
                continue;
            }

            const bool should_swap_m = (anti_p && linkage_checker_.check_linkage(res_m.strand2, res_n.strand1,
                                                                                 backbone) == LinkDirection::Forward) ||
                                       (parallel && linkage_checker_.check_linkage(res_m.strand1, res_n.strand2,
                                                                                   backbone) == LinkDirection::Reverse);

            const bool should_swap_n = (anti_p && linkage_checker_.check_linkage(res_m.strand1, res_n.strand2,
                                                                                 backbone) == LinkDirection::Forward) ||
                                       (parallel && linkage_checker_.check_linkage(res_m.strand2, res_n.strand1,
                                                                                   backbone) == LinkDirection::Reverse);

            if (should_swap_m) {
                swapped[idx_m] = !swapped[idx_m];
            } else if (should_swap_n) {
                swapped[idx_n] = !swapped[idx_n];
            }
        }
    }
}

} // namespace x3dna::algorithms::helix

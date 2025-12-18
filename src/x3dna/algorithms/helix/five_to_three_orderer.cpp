/**
 * @file five_to_three_orderer.cpp
 * @brief Implementation of five-to-three ordering algorithm
 */

#include <x3dna/algorithms/helix/five_to_three_orderer.hpp>
#include <x3dna/config/config_manager.hpp>
#include <iostream>

namespace x3dna::algorithms::helix {

namespace {
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
} // namespace

void FiveToThreeOrderer::ensure_five_to_three(const std::vector<core::BasePair>& pairs, const BackboneData& backbone,
                                              std::vector<size_t>& pair_order, std::vector<HelixSegment>& helices,
                                              std::vector<bool>& swapped) const {

    bool debug = is_debug_enabled();

    swapped.resize(pairs.size(), false);

    if (backbone.empty()) {
        return;
    }

    // Process each helix
    size_t helix_num = 0;
    for (auto& helix : helices) {
        ++helix_num;
        if (helix.start_idx > helix.end_idx)
            continue;

        if (debug) {
            std::cerr << "\n=== HELIX " << helix_num << " (pairs " << helix.start_idx << "-" << helix.end_idx
                      << ") ===" << std::endl;
        }

        // STEP 1: first_step - set initial strand assignment
        direction_checker_.first_step(pairs, backbone, pair_order, helix, swapped);

        if (debug) {
            std::cerr << "[STEP1 first_step] After first_step:" << std::endl;
            for (size_t pos = helix.start_idx; pos <= helix.end_idx; ++pos) {
                size_t idx = pair_order[pos];
                std::cerr << "  pos " << pos << " pair_idx=" << idx << " (" << pairs[idx].residue_idx1() + 1 << ","
                          << pairs[idx].residue_idx2() + 1 << ")"
                          << " swap=" << swapped[idx] << std::endl;
            }
        }

        // STEP 2: First pass through steps - check each consecutive pair
        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];

            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];

            bool rev_wc = direction_checker_.wc_bporien(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            bool rev_o3d = direction_checker_.check_o3dist(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            bool rev_csc = direction_checker_.check_schain(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            bool rev_oth = direction_checker_.check_others(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);

            if (debug) {
                std::cerr << "[STEP2 pass1] pos " << pos << "->" << (pos + 1) << " m=(" << pair_m.residue_idx1() + 1
                          << "," << pair_m.residue_idx2() + 1 << ")"
                          << " n=(" << pair_n.residue_idx1() + 1 << "," << pair_n.residue_idx2() + 1 << ")"
                          << " swap_m=" << swapped[idx_m] << " swap_n=" << swapped[idx_n] << " rev_wc=" << rev_wc
                          << " rev_o3d=" << rev_o3d << " rev_csc=" << rev_csc << " rev_oth=" << rev_oth << std::endl;
            }

            // Apply swap based on checks
            if (rev_wc) {
                swapped[idx_n] = !swapped[idx_n];
                if (debug)
                    std::cerr << "  -> rev_wc: toggled swap_n to " << swapped[idx_n] << std::endl;
            } else if (rev_o3d || rev_csc || rev_oth) {
                swapped[idx_n] = !swapped[idx_n];
                if (debug)
                    std::cerr << "  -> rev_o3d/csc/oth: toggled swap_n to " << swapped[idx_n] << std::endl;
            }

            // Check strand 1 direction
            bool rev_s1 = direction_checker_.chain1dir(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            if (debug)
                std::cerr << "  chain1dir=" << rev_s1 << std::endl;

            if (rev_s1) {
                swapped[idx_n] = !swapped[idx_n];
                if (debug)
                    std::cerr << "  -> rev_s1: toggled swap_n to " << swapped[idx_n] << std::endl;
            }
        }

        // STEP 3: Second pass - re-check WC orientation
        if (debug)
            std::cerr << "[STEP3 pass2] Second pass WC check:" << std::endl;
        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];

            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];

            bool rev_wc = direction_checker_.wc_bporien(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            if (debug) {
                std::cerr << "  pos " << pos << " rev_wc=" << rev_wc << " (swap_m=" << swapped[idx_m]
                          << " swap_n=" << swapped[idx_n] << ")" << std::endl;
            }
            if (rev_wc) {
                swapped[idx_m] = !swapped[idx_m];
                if (debug)
                    std::cerr << "  -> toggled swap_m to " << swapped[idx_m] << std::endl;
            }
        }

        if (debug) {
            std::cerr << "[After pass2] Swap state:" << std::endl;
            for (size_t pos = helix.start_idx; pos <= helix.end_idx; ++pos) {
                size_t idx = pair_order[pos];
                std::cerr << "  pos " << pos << " pair_idx=" << idx << " (" << pairs[idx].residue_idx1() + 1 << ","
                          << pairs[idx].residue_idx2() + 1 << ")"
                          << " swap=" << swapped[idx] << std::endl;
            }
        }

        // STEP 4: check_direction - count backbone linkage directions and apply fixes
        DirectionCounts direction = direction_checker_.check_direction(pairs, backbone, pair_order, helix, swapped);

        if (debug) {
            std::cerr << "[STEP4 check_direction] s1_fwd=" << direction.strand1_forward
                      << " s1_rev=" << direction.strand1_reverse << " s1_none=" << direction.strand1_none
                      << " s2_fwd=" << direction.strand2_forward << " s2_rev=" << direction.strand2_reverse
                      << " s2_none=" << direction.strand2_none << std::endl;
        }

        // STEP 5: check_strand2 - additional corrections based on direction
        direction_checker_.check_strand2(pairs, backbone, pair_order, helix, swapped, direction);

        if (debug) {
            std::cerr << "[After check_strand2] Swap state:" << std::endl;
            for (size_t pos = helix.start_idx; pos <= helix.end_idx; ++pos) {
                size_t idx = pair_order[pos];
                std::cerr << "  pos " << pos << " pair_idx=" << idx << " (" << pairs[idx].residue_idx1() + 1 << ","
                          << pairs[idx].residue_idx2() + 1 << ")"
                          << " swap=" << swapped[idx] << std::endl;
            }
        }

        // STEP 6: check_direction AGAIN (legacy line 1361 - at end of check_strand2)
        // This recomputes direction with updated swaps and may apply additional corrections
        direction_checker_.check_direction(pairs, backbone, pair_order, helix, swapped);

        if (debug) {
            std::cerr << "[FINAL] Helix " << helix_num << " swap state:" << std::endl;
            for (size_t pos = helix.start_idx; pos <= helix.end_idx; ++pos) {
                size_t idx = pair_order[pos];
                std::cerr << "  pos " << pos << " pair_idx=" << idx << " (" << pairs[idx].residue_idx1() + 1 << ","
                          << pairs[idx].residue_idx2() + 1 << ")"
                          << " swap=" << swapped[idx] << std::endl;
            }
        }
    }
}

} // namespace x3dna::algorithms::helix

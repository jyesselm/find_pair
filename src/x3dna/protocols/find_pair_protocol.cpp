/**
 * @file find_pair_protocol.cpp
 * @brief FindPairProtocol implementation
 */

#include <x3dna/protocols/find_pair_protocol.hpp>
#include <x3dna/config/config_manager.hpp>
#include <iostream>

namespace x3dna {
namespace protocols {

FindPairProtocol::FindPairProtocol(const std::filesystem::path& template_path)
    : frame_calculator_(template_path)
    , pair_finder_(algorithms::ValidationParameters::defaults())
{
}

void FindPairProtocol::execute(core::Structure& structure) {
    // Check legacy mode from config if not explicitly set
    if (!legacy_mode_ && config_) {
        legacy_mode_ = config_->legacy_mode();
    }

    // Set legacy mode on frame calculator if needed
    // (Note: BaseFrameCalculator may have its own legacy mode support)
    // frame_calculator_.set_legacy_mode(legacy_mode_);

    // Step 1: Calculate frames for all residues
    calculate_frames(structure);

    // Step 2: Find base pairs
    find_pairs(structure);

    // Step 3: Detect helices (when available)
    // detect_helices(structure);

    // Step 4: Reorder pairs if needed
    // reorder_pairs(structure);
}

void FindPairProtocol::calculate_frames(core::Structure& structure) {
    // Detect RNA by checking for O2' atoms
    bool is_rna = false;
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == " O2'") {
                    is_rna = true;
                    break;
                }
            }
            if (is_rna) break;
        }
        if (is_rna) break;
    }

    frame_calculator_.set_is_rna(is_rna);
    frame_calculator_.calculate_all_frames(structure);

    // Note: Frame calculation recording to JSON is typically done manually
    // by iterating through residues and calling record_base_frame_calc().
    // This is handled at a higher level (e.g., in generate_modern_json.cpp).
    // The protocol focuses on orchestration, not low-level JSON recording.
}

void FindPairProtocol::find_pairs(core::Structure& structure) {
    // Set finding strategy based on options
    if (find_all_pairs_) {
        pair_finder_.set_strategy(algorithms::PairFindingStrategy::ALL_PAIRS);
    } else {
        pair_finder_.set_strategy(algorithms::PairFindingStrategy::BEST_PAIR);
    }

    // Update validation parameters from config if available
    if (config_) {
        const auto& thresholds = config_->thresholds();
        algorithms::ValidationParameters params;
        params.min_dorg = thresholds.min_dorg;
        params.max_dorg = thresholds.max_dorg;
        params.min_dv = thresholds.min_dv;
        params.max_dv = thresholds.max_dv;
        params.min_dNN = thresholds.min_dNN;
        params.max_dNN = thresholds.max_dNN;
        params.min_plane_angle = thresholds.min_plane_angle;
        params.max_plane_angle = thresholds.max_plane_angle;
        params.min_base_hb = thresholds.min_base_hb;
        params.hb_lower = thresholds.hb_lower;
        params.hb_dist1 = thresholds.hb_dist1;
        // Note: hb_dist2 is not part of ValidationParameters
        // It's used in hydrogen bond conflict resolution but not in validation
        params.hb_atoms = thresholds.hb_atoms;
        params.overlap_threshold = thresholds.overlap_threshold;
        pair_finder_.set_parameters(params);
    }

    // Find pairs (with JSON recording if writer provided)
    if (json_writer_) {
        base_pairs_ = pair_finder_.find_pairs_with_recording(structure, json_writer_);
    } else {
        base_pairs_ = pair_finder_.find_pairs(structure);
    }

    // In legacy mode, we may need to process pairs in legacy order
    // For now, the BasePairFinder handles this internally if needed
    // Future: May need to reorder pairs here for exact legacy match
}

void FindPairProtocol::detect_helices(core::Structure& structure) {
    // TODO: Implement helix detection when HelixDetector is available
    // For now, this is a placeholder
    (void)structure;  // Suppress unused parameter warning
}

void FindPairProtocol::reorder_pairs(core::Structure& structure) {
    // TODO: Implement pair reordering
    // This may be needed for exact legacy compatibility
    (void)structure;  // Suppress unused parameter warning
}

} // namespace protocols
} // namespace x3dna


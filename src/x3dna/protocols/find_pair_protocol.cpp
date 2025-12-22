/**
 * @file find_pair_protocol.cpp
 * @brief FindPairProtocol implementation
 */

#include <x3dna/protocols/find_pair_protocol.hpp>
#include <x3dna/config/config_manager.hpp>
#include <x3dna/io/json_writer.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/nucleotide_utils.hpp>
#include <iostream>
#include <filesystem>
#include <vector>

namespace x3dna {
namespace protocols {

FindPairProtocol::FindPairProtocol(const std::filesystem::path& template_path, const FindPairConfig& config)
    : frame_calculator_(template_path), pair_finder_(algorithms::ValidationParameters::defaults()), config_(config) {}

FindPairProtocol::FindPairProtocol(const std::filesystem::path& template_path, const std::filesystem::path& output_dir)
    : frame_calculator_(template_path), pair_finder_(algorithms::ValidationParameters::defaults()) {
    config_.output_dir = output_dir;
}

void FindPairProtocol::execute(core::Structure& structure) {
    // Step 1: Calculate frames for all residues
    calculate_frames(structure);

    // Step 2: Find base pairs (only if not just frames stage)
    if (config_.output_stage != "frames") {
        find_pairs(structure);
    }

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
                if (atom.name() == "O2'") {
                    is_rna = true;
                    break;
                }
            }
            if (is_rna)
                break;
        }
        if (is_rna)
            break;
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
    if (config_.find_all_pairs) {
        pair_finder_.set_strategy(algorithms::PairFindingStrategy::ALL_PAIRS);
    } else {
        pair_finder_.set_strategy(algorithms::PairFindingStrategy::BEST_PAIR);
    }

    // Find pairs (with JSON recording if writer provided)
    if (json_writer_) {
        base_pairs_ = pair_finder_.find_pairs_with_recording(structure, json_writer_);
    } else {
        base_pairs_ = pair_finder_.find_pairs(structure);
    }
}

size_t FindPairProtocol::write_frames_json(core::Structure& structure, const std::filesystem::path& pdb_file,
                                           const std::filesystem::path& output_dir) {
    io::JsonWriter writer(pdb_file);

    // Record residue indices (needed for frames)
    writer.record_residue_indices(structure);

    // Set legacy mode on frame calculator
    frame_calculator_.set_legacy_mode(config_.legacy_mode);

    // Record frame calculations for each residue in legacy order
    size_t frames_recorded = 0;

    // Get residues in legacy order (PDB file order)
    std::vector<core::Residue*> residues_in_order;
    for (const auto* residue_ptr : structure.residues_in_legacy_order()) {
        residues_in_order.push_back(const_cast<core::Residue*>(residue_ptr));
    }

    for (auto* residue : residues_in_order) {
        // Only process nucleotide residues
        bool is_nucleotide = residue->is_nucleotide();

        // Check for modified nucleotides that have ring atoms
        if (!is_nucleotide && residue->molecule_type() == core::typing::MoleculeType::UNKNOWN) {
            static const std::vector<std::string> common_ring_atoms = {"C4", "N3", "C2", "N1", "C6", "C5"};
            int ring_atom_count = 0;
            for (const auto& atom_name : common_ring_atoms) {
                for (const auto& atom : residue->atoms()) {
                    if (atom.name() == atom_name) {
                        ring_atom_count++;
                        break;
                    }
                }
            }
            if (ring_atom_count >= 3) {
                is_nucleotide = true;
            }
        }

        if (is_nucleotide) {
            // Get legacy_residue_idx from residue
            int legacy_residue_idx = residue->legacy_residue_idx();

            if (legacy_residue_idx <= 0) {
                continue;
            }

            // Calculate frame if not already calculated
            if (!residue->reference_frame().has_value()) {
                auto frame_result = frame_calculator_.calculate_frame(*residue);
                if (frame_result.is_valid) {
                    residue->set_reference_frame(frame_result.frame);
                } else {
                    continue;
                }
            }

            // Get frame calculation result (we need to recalculate to get matched atoms, etc.)
            auto frame_result = frame_calculator_.calculate_frame(*residue);
            if (!frame_result.is_valid) {
                continue;
            }

            char base_type = core::one_letter_code(*residue);
            size_t record_idx = static_cast<size_t>(legacy_residue_idx);

            // Record base_frame_calc
            writer.record_base_frame_calc(record_idx, base_type, frame_result.template_file, frame_result.rms_fit,
                                          frame_result.matched_atoms, residue->name(), residue->chain_id(),
                                          residue->seq_num(), residue->insertion());

            // Record ls_fitting
            writer.record_ls_fitting(record_idx, frame_result.num_matched, frame_result.rms_fit,
                                     frame_result.rotation_matrix, frame_result.translation, residue->name(),
                                     residue->chain_id(), residue->seq_num(), residue->insertion());

            // Record frame_calc
            std::vector<geometry::Vector3D> standard_coords, experimental_coords;
            writer.record_frame_calc(record_idx, base_type, frame_result.template_file, frame_result.rms_fit,
                                     standard_coords, experimental_coords, residue->name(), residue->chain_id(),
                                     residue->seq_num(), residue->insertion());

            frames_recorded++;
        }
    }

    // Write split JSON files
    writer.write_split_files(output_dir, true);

    return frames_recorded;
}

} // namespace protocols
} // namespace x3dna

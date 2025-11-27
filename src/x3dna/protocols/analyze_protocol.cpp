/**
 * @file analyze_protocol.cpp
 * @brief AnalyzeProtocol implementation
 */

#include <x3dna/protocols/analyze_protocol.hpp>
#include <x3dna/config/config_manager.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/structure_legacy_order.hpp>
#include <iostream>
#include <algorithm>

namespace x3dna {
namespace protocols {

AnalyzeProtocol::AnalyzeProtocol(const std::filesystem::path& template_path)
    : frame_calculator_(template_path)
{
}

void AnalyzeProtocol::execute(const std::filesystem::path& input_file) {
    // Step 1: Parse .inp file
    input_data_ = io::InputFileParser::parse(input_file);

    // Step 2: Load PDB structure
    core::Structure structure = load_structure(input_data_.pdb_file);

    // Step 3: Store base pairs from input file
    base_pairs_ = input_data_.base_pairs;

    // Step 4: Execute protocol on structure
    execute(structure);
}

void AnalyzeProtocol::execute(core::Structure& structure) {
    // Check legacy mode from config if not explicitly set
    if (!legacy_mode_ && config_) {
        legacy_mode_ = config_->legacy_mode();
    }

    // Step 1: Recalculate frames for all residues
    recalculate_frames(structure);

    // Step 2: Calculate step and helical parameters
    calculate_parameters(structure);
}

void AnalyzeProtocol::recalculate_frames(core::Structure& structure) {
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

    // Update base pairs with recalculated frames
    // For each base pair, get the frames from the corresponding residues
    // The base pair indices from .inp file are 0-based but represent legacy residue indices
    // So we convert back to 1-based to use get_residue_by_legacy_idx()
    for (auto& pair : base_pairs_) {
        // Convert from 0-based to 1-based legacy index
        int legacy_idx1 = static_cast<int>(pair.residue_idx1() + 1);
        int legacy_idx2 = static_cast<int>(pair.residue_idx2() + 1);

        const auto* res1 = core::get_residue_by_legacy_idx(structure, legacy_idx1);
        const auto* res2 = core::get_residue_by_legacy_idx(structure, legacy_idx2);

        if (res1 && res1->reference_frame().has_value()) {
            pair.set_frame1(res1->reference_frame().value());
        }

        if (res2 && res2->reference_frame().has_value()) {
            pair.set_frame2(res2->reference_frame().value());
        }
    }
}

void AnalyzeProtocol::calculate_parameters(core::Structure& /* structure */) {
    // Clear previous results
    step_parameters_.clear();
    helical_parameters_.clear();

    if (base_pairs_.size() < 2) {
        // Need at least 2 pairs to calculate step parameters
        return;
    }

    // Determine range of pairs to process (based on -S option)
    // step_start_ and step_size_ are 1-based (matching legacy)
    size_t start_idx = (step_start_ > 0) ? (step_start_ - 1) : 0;  // Convert to 0-based
    if (start_idx >= base_pairs_.size()) {
        return;  // Start index out of range
    }

    // Calculate parameters for consecutive pairs
    // Apply step_size_ by skipping pairs
    for (size_t i = start_idx; i + 1 < base_pairs_.size(); i += step_size_) {
        const auto& pair1 = base_pairs_[i];
        const auto& pair2 = base_pairs_[i + 1];

        // Verify both pairs have frames
        if (!pair1.frame1().has_value() || !pair1.frame2().has_value() ||
            !pair2.frame1().has_value() || !pair2.frame2().has_value()) {
            continue;  // Skip pairs without frames
        }

        // Calculate step parameters
        auto step_params = param_calculator_.calculate_step_parameters(pair1, pair2);
        step_parameters_.push_back(step_params);

        // Record to JSON if writer provided
        if (json_writer_) {
            // Note: JSON recording typically uses 1-based indices
            // Assuming json_writer has appropriate methods for step parameters
            // This is a placeholder - actual JSON recording may need specific methods
        }

        // Calculate helical parameters
        auto helical_params = param_calculator_.calculate_helical_parameters(pair1, pair2);
        helical_parameters_.push_back(helical_params);

        // Record to JSON if writer provided
        if (json_writer_) {
            // Placeholder for helical parameter JSON recording
        }
    }

    // Handle circular structure (last pair with first pair)
    if (circular_structure_ && base_pairs_.size() >= 2) {
        const auto& pair1 = base_pairs_.back();
        const auto& pair2 = base_pairs_.front();

        if (pair1.frame1().has_value() && pair1.frame2().has_value() &&
            pair2.frame1().has_value() && pair2.frame2().has_value()) {
            auto step_params = param_calculator_.calculate_step_parameters(pair1, pair2);
            step_parameters_.push_back(step_params);

            auto helical_params = param_calculator_.calculate_helical_parameters(pair1, pair2);
            helical_parameters_.push_back(helical_params);
        }
    }
}

core::Structure AnalyzeProtocol::load_structure(const std::filesystem::path& pdb_file) {
    // Configure parser based on input data
    // Note: input_data_.flags may contain information about HETATM inclusion
    // For now, we'll use default settings matching legacy behavior

    // Parse PDB file
    return pdb_parser_.parse_file(pdb_file);
}

} // namespace protocols
} // namespace x3dna


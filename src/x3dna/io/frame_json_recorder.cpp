/**
 * @file frame_json_recorder.cpp
 * @brief Implementation of frame JSON recorder
 */

#include <x3dna/io/frame_json_recorder.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/io/json_writer.hpp>
#include <x3dna/core/residue.hpp>

namespace x3dna {
namespace io {

FrameJsonRecorder::FrameJsonRecorder(algorithms::BaseFrameCalculator& calculator) : calculator_(calculator) {}

size_t FrameJsonRecorder::record_base_frame_calc(core::Structure& structure, JsonWriter& writer) {
    auto residues = structure.residues_in_legacy_order();
    size_t count = 0;

    for (const auto* residue_ptr : residues) {
        auto* residue = const_cast<core::Residue*>(residue_ptr);

        if (residue->residue_type() == core::ResidueType::AMINO_ACID) {
            continue;
        }

        algorithms::FrameCalculationResult frame_result = calculator_.calculate_frame(*residue);
        if (!frame_result.is_valid) {
            continue;
        }

        int legacy_residue_idx = 0;
        if (!residue->atoms().empty()) {
            legacy_residue_idx = residue->atoms()[0].legacy_residue_idx();
        }

        if (legacy_residue_idx <= 0) {
            continue;
        }

        size_t record_idx = static_cast<size_t>(legacy_residue_idx);
        char base_type = residue->one_letter_code();
        writer.record_base_frame_calc(record_idx, base_type, frame_result.template_file, frame_result.rms_fit,
                                      frame_result.matched_atoms, residue->name(), residue->chain_id(),
                                      residue->seq_num(), residue->insertion());
        count++;
    }

    return count;
}

size_t FrameJsonRecorder::record_ls_fitting(core::Structure& structure, JsonWriter& writer) {
    auto residues = structure.residues_in_legacy_order();
    size_t count = 0;

    for (const auto* residue_ptr : residues) {
        auto* residue = const_cast<core::Residue*>(residue_ptr);

        if (residue->residue_type() == core::ResidueType::AMINO_ACID) {
            continue;
        }

        algorithms::FrameCalculationResult frame_result = calculator_.calculate_frame(*residue);
        if (!frame_result.is_valid) {
            continue;
        }

        int legacy_residue_idx = 0;
        if (!residue->atoms().empty()) {
            legacy_residue_idx = residue->atoms()[0].legacy_residue_idx();
        }

        if (legacy_residue_idx <= 0) {
            continue;
        }

        size_t record_idx = static_cast<size_t>(legacy_residue_idx);
        writer.record_ls_fitting(record_idx, frame_result.num_matched, frame_result.rms_fit,
                                 frame_result.rotation_matrix, frame_result.translation, residue->name(),
                                 residue->chain_id(), residue->seq_num(), residue->insertion());
        count++;
    }

    return count;
}

size_t FrameJsonRecorder::record_frame_calc(core::Structure& structure, JsonWriter& writer) {
    auto residues = structure.residues_in_legacy_order();
    size_t count = 0;

    for (const auto* residue_ptr : residues) {
        auto* residue = const_cast<core::Residue*>(residue_ptr);

        if (residue->residue_type() == core::ResidueType::AMINO_ACID) {
            continue;
        }

        algorithms::FrameCalculationResult frame_result = calculator_.calculate_frame(*residue);
        if (!frame_result.is_valid) {
            continue;
        }

        int legacy_residue_idx = 0;
        if (!residue->atoms().empty()) {
            legacy_residue_idx = residue->atoms()[0].legacy_residue_idx();
        }

        if (legacy_residue_idx <= 0) {
            continue;
        }

        size_t record_idx = static_cast<size_t>(legacy_residue_idx);
        char base_type = residue->one_letter_code();
        writer.record_frame_calc(record_idx, base_type, frame_result.template_file, frame_result.rms_fit,
                                 frame_result.matched_standard_coords, frame_result.matched_experimental_coords,
                                 residue->name(), residue->chain_id(), residue->seq_num(), residue->insertion());
        count++;
    }

    return count;
}

size_t FrameJsonRecorder::record_all(core::Structure& structure, JsonWriter& writer) {
    size_t total = 0;
    total += record_base_frame_calc(structure, writer);
    total += record_ls_fitting(structure, writer);
    total += record_frame_calc(structure, writer);
    return total;
}

} // namespace io
} // namespace x3dna

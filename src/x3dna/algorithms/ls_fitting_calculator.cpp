/**
 * @file ls_fitting_calculator.cpp
 * @brief Implementation of least-squares fitting calculator
 */

#include <x3dna/algorithms/ls_fitting_calculator.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/io/json_writer.hpp>

namespace x3dna {
namespace algorithms {

LsFittingCalculator::LsFittingCalculator(const std::filesystem::path& template_path)
    : calculator_(new BaseFrameCalculator(template_path)) {}

LsFittingCalculator::~LsFittingCalculator() {
    delete calculator_;
}

size_t LsFittingCalculator::calculate_and_record(core::Structure& structure,
                                                 io::JsonWriter& writer) {
    // Auto-detect RNA vs DNA
    bool is_rna = detect_rna(structure);
    set_is_rna(is_rna);

    // Calculate and record ls_fitting in one pass
    // Get residues in legacy order (PDB file order)
    auto residues = structure.residues_in_legacy_order();
    size_t records_count = 0;

    for (const auto* residue_ptr : residues) {
        auto* residue = const_cast<core::Residue*>(residue_ptr);

        // Skip amino acids (calculate_frame handles other types including UNKNOWN)
        if (residue->residue_type() == core::ResidueType::AMINO_ACID) {
            continue;
        }

        // Calculate frame (stores frame on residue and returns full result)
        FrameCalculationResult frame_result = calculator_->calculate_frame(*residue);

        if (!frame_result.is_valid) {
            continue;
        }

        // Get legacy_residue_idx from atoms
        int legacy_residue_idx = 0;
        if (!residue->atoms().empty()) {
            legacy_residue_idx = residue->atoms()[0].legacy_residue_idx();
        }

        if (legacy_residue_idx <= 0) {
            continue;
        }

        size_t record_idx = static_cast<size_t>(legacy_residue_idx);

        // Record ONLY ls_fitting
        writer.record_ls_fitting(record_idx, frame_result.num_matched, frame_result.rms_fit,
                                frame_result.rotation_matrix, frame_result.translation,
                                residue->name(), residue->chain_id(), residue->seq_num(),
                                residue->insertion());

        records_count++;
    }

    return records_count;
}

void LsFittingCalculator::set_template_path(const std::filesystem::path& template_path) {
    calculator_->set_template_path(template_path);
}

std::filesystem::path LsFittingCalculator::template_path() const {
    return calculator_->template_path();
}

void LsFittingCalculator::set_is_rna(bool is_rna) {
    calculator_->set_is_rna(is_rna);
}

bool LsFittingCalculator::is_rna() const {
    return calculator_->is_rna();
}

void LsFittingCalculator::set_legacy_mode(bool legacy_mode) {
    calculator_->set_legacy_mode(legacy_mode);
}

bool LsFittingCalculator::legacy_mode() const {
    return calculator_->legacy_mode();
}

bool LsFittingCalculator::detect_rna(const core::Structure& structure) {
    return BaseFrameCalculator::detect_rna(structure);
}

} // namespace algorithms
} // namespace x3dna


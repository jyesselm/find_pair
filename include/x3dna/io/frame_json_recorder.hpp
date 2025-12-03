/**
 * @file frame_json_recorder.hpp
 * @brief Records frame calculation JSON (base_frame_calc, ls_fitting, frame_calc)
 */

#pragma once

#include <x3dna/core/structure.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/io/json_writer.hpp>
#include <vector>

namespace x3dna {
namespace io {

/**
 * @class FrameJsonRecorder
 * @brief Records frame calculation JSON using BaseFrameCalculator
 *
 * This class handles JSON recording for frame calculations.
 * It uses BaseFrameCalculator internally but only handles recording.
 *
 * Responsibilities:
 * - Iterate through residues in legacy order
 * - Use BaseFrameCalculator to calculate frames
 * - Record JSON via JsonWriter
 * - Handle different recording scenarios (base_frame_calc, ls_fitting, frame_calc)
 */
class FrameJsonRecorder {
public:
    /**
     * @brief Constructor
     * @param calculator BaseFrameCalculator to use for calculations
     */
    explicit FrameJsonRecorder(algorithms::BaseFrameCalculator& calculator);

    /**
     * @brief Record base_frame_calc JSON for all residues
     * @param structure Structure to process
     * @param writer JsonWriter to record results
     * @return Number of records written
     */
    size_t record_base_frame_calc(core::Structure& structure, JsonWriter& writer);

    /**
     * @brief Record ls_fitting JSON for all residues
     * @param structure Structure to process
     * @param writer JsonWriter to record results
     * @return Number of records written
     */
    size_t record_ls_fitting(core::Structure& structure, JsonWriter& writer);

    /**
     * @brief Record frame_calc JSON for all residues
     * @param structure Structure to process
     * @param writer JsonWriter to record results
     * @return Number of records written
     */
    size_t record_frame_calc(core::Structure& structure, JsonWriter& writer);

    /**
     * @brief Record all frame JSON types (base_frame_calc, ls_fitting, frame_calc)
     * @param structure Structure to process
     * @param writer JsonWriter to record results
     * @return Number of records written (total across all types)
     */
    size_t record_all(core::Structure& structure, JsonWriter& writer);

private:
    algorithms::BaseFrameCalculator& calculator_;

    /**
     * @brief Helper: iterate residues and call record_func for each valid frame
     * @param structure Structure to iterate
     * @param writer JsonWriter for recording
     * @param record_func Function to call for each valid frame
     * @return Number of records written
     */
    template <typename RecordFunc>
    size_t iterate_and_record(core::Structure& structure, JsonWriter& writer,
                              RecordFunc record_func) {
        // Get residues in legacy order (PDB file order)
        auto residues = structure.residues_in_legacy_order();
        size_t count = 0;

        for (const auto* residue_ptr : residues) {
            auto* residue = const_cast<core::Residue*>(residue_ptr);

            // Skip amino acids (calculate_frame handles other types including UNKNOWN)
            if (residue->residue_type() == core::ResidueType::AMINO_ACID) {
                continue;
            }

            // Calculate frame (stores frame on residue and returns full result)
            algorithms::FrameCalculationResult frame_result = calculator_.calculate_frame(*residue);

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
            record_func(record_idx, *residue, frame_result, writer);
            count++;
        }

        return count;
    }
};

} // namespace io
} // namespace x3dna

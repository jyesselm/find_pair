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
};

} // namespace io
} // namespace x3dna

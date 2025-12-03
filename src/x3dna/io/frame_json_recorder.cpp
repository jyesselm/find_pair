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

FrameJsonRecorder::FrameJsonRecorder(algorithms::BaseFrameCalculator& calculator)
    : calculator_(calculator) {}

size_t FrameJsonRecorder::record_base_frame_calc(core::Structure& structure, JsonWriter& writer) {
    return iterate_and_record(structure, writer,
                              [](size_t idx, const core::Residue& res,
                                 const algorithms::FrameCalculationResult& result, JsonWriter& w) {
                                  char base_type = res.one_letter_code();
                                  w.record_base_frame_calc(idx, base_type, result.template_file,
                                                           result.rms_fit, result.matched_atoms,
                                                           res.name(), res.chain_id(),
                                                           res.seq_num(), res.insertion());
                              });
}

size_t FrameJsonRecorder::record_ls_fitting(core::Structure& structure, JsonWriter& writer) {
    return iterate_and_record(structure, writer,
                              [](size_t idx, const core::Residue& res,
                                 const algorithms::FrameCalculationResult& result, JsonWriter& w) {
                                  w.record_ls_fitting(idx, result.num_matched, result.rms_fit,
                                                      result.rotation_matrix, result.translation,
                                                      res.name(), res.chain_id(), res.seq_num(),
                                                      res.insertion());
                              });
}

size_t FrameJsonRecorder::record_frame_calc(core::Structure& structure, JsonWriter& writer) {
    return iterate_and_record(
        structure, writer,
        [](size_t idx, const core::Residue& res, const algorithms::FrameCalculationResult& result,
           JsonWriter& w) {
            char base_type = res.one_letter_code();
            std::vector<geometry::Vector3D> standard_coords, experimental_coords;
            w.record_frame_calc(idx, base_type, result.template_file, result.rms_fit,
                                standard_coords, experimental_coords, res.name(), res.chain_id(),
                                res.seq_num(), res.insertion());
        });
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

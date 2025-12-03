/**
 * @file find_pair_protocol.hpp
 * @brief Protocol for finding base pairs (matches legacy find_pair workflow)
 */

#pragma once

#include <x3dna/protocols/protocol_base.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/structure_legacy_order.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/base_pair_finder.hpp>
#include <x3dna/algorithms/helix_detector.hpp>
#include <x3dna/io/json_writer.hpp>
#include <memory>
#include <vector>
#include <filesystem>

namespace x3dna {
namespace protocols {

/**
 * @class FindPairProtocol
 * @brief Orchestrates the find_pair workflow
 *
 * This protocol:
 * 1. Calculates reference frames for all residues
 * 2. Finds base pairs using specified strategy
 * 3. Detects helices (when HelixDetector is available)
 * 4. Records results to JSON (if writer provided)
 *
 * Supports legacy mode for exact compatibility with legacy code.
 */
class FindPairProtocol : public ProtocolBase {
public:
    /**
     * @brief Constructor
     * @param template_path Path to standard base template directory
     * @param output_dir Output directory for JSON files (optional, for Phase 2 refactor)
     */
    explicit FindPairProtocol(const std::filesystem::path& template_path = "data/templates",
                              const std::filesystem::path& output_dir = "");

    /**
     * @brief Execute the find_pair protocol
     * @param structure Structure to process
     */
    void execute(core::Structure& structure) override;

    // Options
    /**
     * @brief Set single strand mode
     */
    void set_single_strand_mode(bool value) {
        single_strand_mode_ = value;
    }

    /**
     * @brief Get single strand mode
     */
    bool single_strand_mode() const {
        return single_strand_mode_;
    }

    /**
     * @brief Set find all pairs mode
     */
    void set_find_all_pairs(bool value) {
        find_all_pairs_ = value;
    }

    /**
     * @brief Get find all pairs mode
     */
    bool find_all_pairs() const {
        return find_all_pairs_;
    }

    /**
     * @brief Set divide helices mode
     */
    void set_divide_helices(bool value) {
        divide_helices_ = value;
    }

    /**
     * @brief Get divide helices mode
     */
    bool divide_helices() const {
        return divide_helices_;
    }

    /**
     * @brief Set legacy mode (for exact compatibility)
     */
    void set_legacy_mode(bool value) {
        legacy_mode_ = value;
    }

    /**
     * @brief Get legacy mode
     */
    bool legacy_mode() const {
        return legacy_mode_;
    }

    /**
     * @brief Set fix indices from legacy JSON (for comparison with legacy)
     */
    void set_fix_indices_from_legacy_json(bool value, const std::string& json_file = "") {
        fix_indices_from_legacy_json_ = value;
        legacy_json_file_ = json_file;
    }

    /**
     * @brief Get fix indices from legacy JSON
     */
    bool fix_indices_from_legacy_json() const {
        return fix_indices_from_legacy_json_;
    }

    /**
     * @brief Set output directory for JSON files
     * @param output_dir Directory where JSON files will be written
     */
    void set_output_dir(const std::filesystem::path& output_dir) {
        output_dir_ = output_dir;
    }

    /**
     * @brief Get output directory
     */
    const std::filesystem::path& output_dir() const {
        return output_dir_;
    }

    /**
     * @brief Set output stage (controls which JSON to write)
     * @param stage Stage name: "frames", "distances", "hbonds", "validation", "selection", "all"
     */
    void set_output_stage(const std::string& stage) {
        output_stage_ = stage;
    }

    /**
     * @brief Get output stage
     */
    const std::string& output_stage() const {
        return output_stage_;
    }

    /**
     * @brief Set JSON writer for recording results
     * @deprecated Phase 2 refactor: Protocols should write their own JSON
     */
    void set_json_writer(io::JsonWriter* writer) {
        json_writer_ = writer;
    }

    /**
     * @brief Get found base pairs
     */
    const std::vector<core::BasePair>& base_pairs() const {
        return base_pairs_;
    }

    /**
     * @brief Get detected helices
     */
    const std::vector<algorithms::Helix>& helices() const {
        return helices_;
    }

    /**
     * @brief Write frames JSON (base_frame_calc, frame_calc, ls_fitting)
     * @param structure Structure with calculated frames
     * @param pdb_file PDB file path (for JsonWriter initialization)
     * @param output_dir Output directory for JSON files
     * @return Number of frames written
     */
    size_t write_frames_json(core::Structure& structure, const std::filesystem::path& pdb_file,
                             const std::filesystem::path& output_dir);

    /**
     * @brief Get frame calculator (for configuration)
     */
    algorithms::BaseFrameCalculator& frame_calculator() {
        return frame_calculator_;
    }

    /**
     * @brief Get pair finder (for configuration)
     */
    algorithms::BasePairFinder& pair_finder() {
        return pair_finder_;
    }

private:
    /**
     * @brief Calculate frames for all residues
     */
    void calculate_frames(core::Structure& structure);

    /**
     * @brief Find base pairs
     */
    void find_pairs(core::Structure& structure);

    /**
     * @brief Detect helices using HelixDetector
     */
    void detect_helices(core::Structure& structure);

    /**
     * @brief Reorder pairs using HelixDetector
     */
    void reorder_pairs(core::Structure& structure);

    // Algorithm components
    algorithms::BaseFrameCalculator frame_calculator_;
    algorithms::BasePairFinder pair_finder_;
    algorithms::HelixDetector helix_detector_;

    // Results - helices
    std::vector<algorithms::Helix> helices_;

    // Options
    bool single_strand_mode_ = false;
    bool find_all_pairs_ = false;
    bool divide_helices_ = false;
    bool legacy_mode_ = false;
    bool fix_indices_from_legacy_json_ = false;
    std::string legacy_json_file_ = ""; // Legacy compatibility mode

    // Results
    std::vector<core::BasePair> base_pairs_;

    // Output directory for JSON files (Phase 2 refactor)
    std::filesystem::path output_dir_;

    // Output stage (controls which JSON to write)
    std::string output_stage_ = "all";

    // JSON writer (optional, for recording)
    // TODO: Phase 4 - Remove JsonWriter once protocols write their own JSON
    io::JsonWriter* json_writer_ = nullptr;
};

} // namespace protocols
} // namespace x3dna

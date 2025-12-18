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
 * @struct FindPairConfig
 * @brief Configuration options for FindPairProtocol
 *
 * Consolidates all protocol options into a single struct for cleaner
 * initialization and reduced setter boilerplate.
 */
struct FindPairConfig {
    bool single_strand_mode = false;   ///< Process single strand only
    bool find_all_pairs = false;       ///< Find all potential pairs
    bool divide_helices = false;       ///< Divide structures into helices
    bool legacy_mode = false;          ///< Enable legacy compatibility mode
    std::filesystem::path output_dir;  ///< Output directory for JSON files
    std::string output_stage = "all";  ///< Output stage: "frames", "distances", "hbonds", etc.
};

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
     * @brief Constructor with config struct (preferred)
     * @param template_path Path to standard base template directory
     * @param config Protocol configuration options
     */
    explicit FindPairProtocol(const std::filesystem::path& template_path,
                              const FindPairConfig& config);

    /**
     * @brief Constructor (legacy - for backward compatibility)
     * @param template_path Path to standard base template directory
     * @param output_dir Output directory for JSON files (optional)
     */
    explicit FindPairProtocol(const std::filesystem::path& template_path = "data/templates",
                              const std::filesystem::path& output_dir = "");

    /**
     * @brief Execute the find_pair protocol
     * @param structure Structure to process
     */
    void execute(core::Structure& structure) override;

    // Configuration
    /**
     * @brief Get configuration (const)
     */
    [[nodiscard]] const FindPairConfig& config() const { return config_; }

    /**
     * @brief Get configuration (mutable, for modification)
     */
    FindPairConfig& config() { return config_; }

    // Individual setters (for backward compatibility)
    void set_single_strand_mode(bool value) { config_.single_strand_mode = value; }
    [[nodiscard]] bool single_strand_mode() const { return config_.single_strand_mode; }

    void set_find_all_pairs(bool value) { config_.find_all_pairs = value; }
    [[nodiscard]] bool find_all_pairs() const { return config_.find_all_pairs; }

    void set_divide_helices(bool value) { config_.divide_helices = value; }
    [[nodiscard]] bool divide_helices() const { return config_.divide_helices; }

    void set_legacy_mode(bool value) { config_.legacy_mode = value; }
    [[nodiscard]] bool legacy_mode() const { return config_.legacy_mode; }

    void set_output_dir(const std::filesystem::path& dir) { config_.output_dir = dir; }
    [[nodiscard]] const std::filesystem::path& output_dir() const { return config_.output_dir; }

    void set_output_stage(const std::string& stage) { config_.output_stage = stage; }
    [[nodiscard]] const std::string& output_stage() const { return config_.output_stage; }

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

    // Configuration (consolidates all options)
    FindPairConfig config_;

    // Results - helices
    std::vector<algorithms::Helix> helices_;

    // Results
    std::vector<core::BasePair> base_pairs_;

    // JSON writer (optional, for recording)
    // TODO: Phase 4 - Remove JsonWriter once protocols write their own JSON
    io::JsonWriter* json_writer_ = nullptr;
};

} // namespace protocols
} // namespace x3dna

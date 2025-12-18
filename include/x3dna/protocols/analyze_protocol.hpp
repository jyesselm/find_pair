/**
 * @file analyze_protocol.hpp
 * @brief Protocol for analyzing base pair step parameters (matches legacy analyze workflow)
 */

#pragma once

#include <x3dna/protocols/protocol_base.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/parameters.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/parameter_calculator.hpp>
#include <x3dna/io/input_file_parser.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/io/json_writer.hpp>
#include <memory>
#include <vector>
#include <filesystem>
#include <string>

namespace x3dna {
namespace protocols {

/**
 * @struct AnalyzeConfig
 * @brief Configuration options for AnalyzeProtocol
 *
 * Consolidates all protocol options into a single struct for cleaner
 * initialization and reduced setter boilerplate.
 */
struct AnalyzeConfig {
    bool calculate_torsions = false;   ///< Calculate torsion angles
    bool simple_parameters = false;    ///< Use simple parameter mode
    bool circular_structure = false;   ///< Treat as circular structure
    size_t step_start = 1;             ///< Step start index (1-based, for -S option)
    size_t step_size = 1;              ///< Step size (1-based, for -S option)
    bool legacy_mode = false;          ///< Enable legacy compatibility mode
};

/**
 * @class AnalyzeProtocol
 * @brief Orchestrates the analyze workflow
 *
 * This protocol:
 * 1. Reads .inp file (created by find_pair)
 * 2. Loads PDB file
 * 3. Recalculates reference frames for base pairs
 * 4. Calculates step parameters for consecutive pairs
 * 5. Calculates helical parameters
 * 6. Records results to JSON (if writer provided)
 *
 * Supports legacy mode for exact compatibility with legacy code.
 */
class AnalyzeProtocol : public ProtocolBase {
public:
    /**
     * @brief Constructor with config struct (preferred)
     * @param template_path Path to standard base template directory
     * @param config Protocol configuration options
     */
    explicit AnalyzeProtocol(const std::filesystem::path& template_path,
                             const AnalyzeConfig& config);

    /**
     * @brief Constructor (legacy - for backward compatibility)
     * @param template_path Path to standard base template directory
     */
    explicit AnalyzeProtocol(const std::filesystem::path& template_path = "data/templates");

    /**
     * @brief Execute the analyze protocol
     * @param input_file Path to .inp file to analyze
     */
    void execute(const std::filesystem::path& input_file);

    /**
     * @brief Execute the analyze protocol on a structure with base pairs
     * @param structure Structure to process (must have base pairs with frames)
     */
    void execute(core::Structure& structure) override;

    // Configuration
    /**
     * @brief Get configuration (const)
     */
    [[nodiscard]] const AnalyzeConfig& config() const { return config_; }

    /**
     * @brief Get configuration (mutable, for modification)
     */
    AnalyzeConfig& config() { return config_; }

    // Individual setters (for backward compatibility)
    void set_calculate_torsions(bool value) { config_.calculate_torsions = value; }
    [[nodiscard]] bool calculate_torsions() const { return config_.calculate_torsions; }

    void set_simple_parameters(bool value) { config_.simple_parameters = value; }
    [[nodiscard]] bool simple_parameters() const { return config_.simple_parameters; }

    void set_circular_structure(bool value) { config_.circular_structure = value; }
    [[nodiscard]] bool circular_structure() const { return config_.circular_structure; }

    void set_step_start(size_t start) { config_.step_start = start; }
    [[nodiscard]] size_t step_start() const { return config_.step_start; }

    void set_step_size(size_t size) { config_.step_size = size; }
    [[nodiscard]] size_t step_size() const { return config_.step_size; }

    void set_legacy_mode(bool value) { config_.legacy_mode = value; }
    [[nodiscard]] bool legacy_mode() const { return config_.legacy_mode; }

    /**
     * @brief Set JSON writer for recording results
     */
    void set_json_writer(io::JsonWriter* writer) {
        json_writer_ = writer;
    }

    /**
     * @brief Get calculated step parameters
     */
    const std::vector<core::BasePairStepParameters>& step_parameters() const {
        return step_parameters_;
    }

    /**
     * @brief Get calculated helical parameters
     */
    const std::vector<core::HelicalParameters>& helical_parameters() const {
        return helical_parameters_;
    }

    /**
     * @brief Get base pairs from input file
     */
    const std::vector<core::BasePair>& base_pairs() const {
        return base_pairs_;
    }

    /**
     * @brief Get frame calculator (for configuration)
     */
    algorithms::BaseFrameCalculator& frame_calculator() {
        return frame_calculator_;
    }

    /**
     * @brief Get parameter calculator (for configuration)
     */
    algorithms::ParameterCalculator& parameter_calculator() {
        return param_calculator_;
    }

private:
    /**
     * @brief Recalculate frames for all residues in base pairs
     */
    void recalculate_frames(core::Structure& structure);

    /**
     * @brief Calculate step and helical parameters
     */
    void calculate_parameters(core::Structure& structure);

    /**
     * @brief Load structure from PDB file
     */
    core::Structure load_structure(const std::filesystem::path& pdb_file);

    /**
     * @brief Convert atom indices to residue indices in base pairs
     *
     * Legacy input files may contain atom indices instead of residue indices.
     * This method detects atom indices (values larger than residue count) and
     * converts them to residue indices using the structure's atom-to-residue mapping.
     */
    void convert_atom_indices_to_residue_indices(const core::Structure& structure);

    // Algorithm components
    algorithms::BaseFrameCalculator frame_calculator_;
    algorithms::ParameterCalculator param_calculator_;
    io::PdbParser pdb_parser_;

    // Configuration (consolidates all options)
    AnalyzeConfig config_;

    // Results
    std::vector<core::BasePair> base_pairs_;
    std::vector<core::BasePairStepParameters> step_parameters_;
    std::vector<core::HelicalParameters> helical_parameters_;

    // JSON writer (optional, for recording)
    io::JsonWriter* json_writer_ = nullptr;

    // Input data (from .inp file)
    io::InputData input_data_;
};

} // namespace protocols
} // namespace x3dna

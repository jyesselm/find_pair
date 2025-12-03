/**
 * @file ls_fitting_calculator.hpp
 * @brief Least-squares fitting calculator for frame calculations
 */

#pragma once

#include <string>
#include <filesystem>
#include <x3dna/core/structure.hpp>

namespace x3dna {
namespace io {
class JsonWriter; // Forward declaration
} // namespace io

namespace algorithms {

// Forward declaration
class BaseFrameCalculator;

/**
 * @class LsFittingCalculator
 * @brief Calculates and records least-squares fitting data for nucleotide residues
 *
 * This class handles Stage 3: ls_fitting JSON generation.
 * It uses BaseFrameCalculator internally but only records ls_fitting data.
 */
class LsFittingCalculator {
public:
    /**
     * @brief Constructor
     * @param template_path Path to standard base template directory
     */
    explicit LsFittingCalculator(const std::filesystem::path& template_path = "data/templates");

    /**
     * @brief Calculate frames and record only ls_fitting JSON
     * 
     * Stage 3: Only records ls_fitting (least-squares fitting data)
     * 
     * @param structure Structure to process
     * @param writer JsonWriter to record results
     * @return Number of ls_fitting records calculated and recorded
     */
    size_t calculate_and_record(core::Structure& structure, io::JsonWriter& writer);

    /**
     * @brief Set template path
     * @param template_path Path to template directory
     */
    void set_template_path(const std::filesystem::path& template_path);

    /**
     * @brief Get template path
     * @return Current template path
     */
    std::filesystem::path template_path() const;

    /**
     * @brief Set whether to process RNA (includes C1' in matching)
     * @param is_rna True for RNA, false for DNA
     */
    void set_is_rna(bool is_rna);

    /**
     * @brief Get RNA flag
     * @return True if processing RNA
     */
    bool is_rna() const;

    /**
     * @brief Set legacy compatibility mode (excludes C4 atom from matching)
     * @param legacy_mode True to exclude C4 (matches legacy behavior), false to include C4
     */
    void set_legacy_mode(bool legacy_mode);

    /**
     * @brief Get legacy compatibility mode
     * @return True if legacy mode is enabled (C4 excluded)
     */
    bool legacy_mode() const;

    /**
     * @brief Detect if structure is RNA by checking for O2' atoms
     * @param structure Structure to check
     * @return True if RNA detected (O2' atoms found), false if DNA
     */
    static bool detect_rna(const core::Structure& structure);

    /**
     * @brief Destructor
     */
    ~LsFittingCalculator();

private:
    BaseFrameCalculator* calculator_; // Uses BaseFrameCalculator internally
};

} // namespace algorithms
} // namespace x3dna


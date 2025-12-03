/**
 * @file base_frame_calculator.hpp
 * @brief Base frame calculator for reference frame determination
 */

#pragma once

#include <string>
#include <filesystem>
#include <vector>
#include <optional>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/algorithms/standard_base_templates.hpp>
#include <x3dna/algorithms/ring_atom_matcher.hpp>
#include <x3dna/geometry/least_squares_fitter.hpp>

namespace x3dna {
namespace io {
class JsonWriter; // Forward declaration
} // namespace io
} // namespace x3dna

namespace x3dna {
namespace algorithms {

/**
 * @struct FrameCalculationResult
 * @brief Result of base frame calculation
 */
struct FrameCalculationResult {
    core::ReferenceFrame frame;             // Calculated reference frame
    double rms_fit = 0.0;                   // RMS fit quality
    std::vector<std::string> matched_atoms; // Names of matched atoms
    size_t num_matched = 0;                 // Number of matched atoms
    std::filesystem::path template_file;    // Path to template file used
    bool is_valid = false;                  // Whether calculation was successful

    // Additional metrics for debugging/comparison
    geometry::Matrix3D rotation_matrix; // Rotation matrix (3x3)
    geometry::Vector3D translation;     // Translation vector
};

/**
 * @class BaseFrameCalculator
 * @brief Calculates reference frames for nucleotide residues
 *
 * Uses least-squares fitting to align experimental ring atoms with standard
 * base template atoms. The resulting transformation defines the reference frame.
 */
class BaseFrameCalculator {
public:
    /**
     * @brief Constructor
     * @param template_path Path to standard base template directory
     */
    explicit BaseFrameCalculator(const std::filesystem::path& template_path = "data/templates");

    /**
     * @brief Calculate reference frame for a residue
     * @param residue Residue to calculate frame for (will be modified to store frame)
     * @return FrameCalculationResult with frame and metrics
     */
    FrameCalculationResult calculate_frame(core::Residue& residue);

    /**
     * @brief Calculate frame without modifying residue
     * @param residue Residue to calculate frame for
     * @return FrameCalculationResult with frame and metrics
     */
    FrameCalculationResult calculate_frame_const(const core::Residue& residue) const;

    /**
     * @brief Calculate frames for all residues in a structure
     * @param structure Structure to calculate frames for (residues will be modified)
     */
    void calculate_all_frames(core::Structure& structure);

    /**
     * @brief Set template path
     * @param template_path Path to template directory
     */
    void set_template_path(const std::filesystem::path& template_path);

    /**
     * @brief Get template path
     * @return Current template path
     */
    std::filesystem::path template_path() const {
        return templates_.template_path();
    }

    /**
     * @brief Set whether to process RNA (includes C1' in matching)
     * @param is_rna True for RNA, false for DNA
     */
    void set_is_rna(bool is_rna) {
        is_rna_ = is_rna;
    }

    /**
     * @brief Get RNA flag
     * @return True if processing RNA
     */
    bool is_rna() const {
        return is_rna_;
    }

    /**
     * @brief Set legacy compatibility mode (excludes C4 atom from matching)
     * @param legacy_mode True to exclude C4 (matches legacy behavior), false to include C4
     * (default, correct behavior)
     */
    void set_legacy_mode(bool legacy_mode) {
        legacy_mode_ = legacy_mode;
    }

    /**
     * @brief Get legacy compatibility mode
     * @return True if legacy mode is enabled (C4 excluded)
     */
    bool legacy_mode() const {
        return legacy_mode_;
    }

    /**
     * @brief Detect if structure is RNA by checking for O2' atoms
     * @param structure Structure to check
     * @return True if RNA detected (O2' atoms found), false if DNA
     */
    static bool detect_rna(const core::Structure& structure);

    /**
     * @brief Calculate frames and write JSON records
     * 
     * This is the standard workflow for frame calculation:
     * 1. Auto-detect RNA vs DNA
     * 2. Calculate all frames
     * 3. Record frame JSON (base_frame_calc, ls_fitting, frame_calc)
     * 
     * @param structure Structure to process
     * @param writer JsonWriter to record results
     * @return Number of frames calculated and recorded
     */
    size_t calculate_and_record_frames(core::Structure& structure, io::JsonWriter& writer);

    /**
     * @brief Calculate frames and record only ls_fitting JSON
     * 
     * Stage 3: Only records ls_fitting (least-squares fitting data)
     * 
     * @param structure Structure to process
     * @param writer JsonWriter to record results
     * @return Number of frames calculated and recorded
     */
    size_t calculate_and_record_ls_fitting(core::Structure& structure, io::JsonWriter& writer);

    /**
     * @brief Calculate frames and record only base_frame_calc and frame_calc JSON
     * 
     * Stage 4: Only records base_frame_calc and frame_calc (not ls_fitting)
     * 
     * @param structure Structure to process
     * @param writer JsonWriter to record results
     * @return Number of frames calculated and recorded
     */
    size_t calculate_and_record_frames_only(core::Structure& structure, io::JsonWriter& writer);

private:
    mutable StandardBaseTemplates
        templates_; // Mutable for caching (doesn't affect logical constness)
    bool is_rna_ = false;
    bool legacy_mode_ = false; // If true, exclude C4 atom to match legacy behavior

    /**
     * @brief Calculate frame for a single residue (implementation)
     * @param residue Residue to calculate frame for
     * @return FrameCalculationResult
     */
    FrameCalculationResult calculate_frame_impl(const core::Residue& residue) const;
};

} // namespace algorithms
} // namespace x3dna

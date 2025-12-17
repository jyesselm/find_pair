/**
 * @file residue_type_detector.hpp
 * @brief Detects nucleotide types via RMSD and atom analysis
 */

#pragma once

#include <vector>
#include <string>
#include <optional>
#include <x3dna/core/residue.hpp>
#include <x3dna/geometry/vector3d.hpp>

namespace x3dna {
namespace algorithms {

/**
 * @struct RmsdCheckResult
 * @brief Result of RMSD check for nucleotide type detection
 */
struct RmsdCheckResult {
    std::optional<double> rmsd;                                  // RMSD value if calculable
    bool found_purine_atoms;                                     // Whether any purine atoms (N7, C8, N9) were found
    std::vector<std::string> matched_atom_names;                 // Atom names that were matched in RMSD check
    std::vector<geometry::Vector3D> matched_experimental_coords; // Experimental coordinates from RMSD check
    std::vector<geometry::Vector3D> matched_standard_coords;     // Standard coordinates from RMSD check
};

/**
 * @struct TypeDetectionResult
 * @brief Result of residue type detection
 */
struct TypeDetectionResult {
    core::ResidueType detected_type; // Detected residue type
    bool used_fallback;              // Whether fallback logic was used
    std::optional<double> rmsd;      // RMSD value if calculated
    std::string detection_method;    // "registry", "rmsd", "atom_analysis", "standard"
};

/**
 * @class ResidueTypeDetector
 * @brief Determines nucleotide type via RMSD and atom analysis
 *
 * This class handles residue type detection using multiple methods:
 * 1. Registry check for modified nucleotides
 * 2. RMSD-based detection for non-standard nucleotides
 * 3. Atom-based detection (purine vs pyrimidine)
 * 4. Standard type for known nucleotides
 */
class ResidueTypeDetector {
public:
    /**
     * @brief Check nucleotide type by RMSD (matches legacy check_nt_type_by_rmsd)
     * @param residue Residue to check
     * @return RmsdCheckResult with RMSD value and matched atoms
     *
     * This performs least-squares fitting of ring atoms to standard geometry
     * to detect whether residue is a purine or pyrimidine.
     */
    [[nodiscard]] static RmsdCheckResult check_by_rmsd(const core::Residue& residue);

    /**
     * @brief Detect residue type using full detection logic
     * @param residue Residue to detect type for
     * @return TypeDetectionResult with detected type and method used
     */
    [[nodiscard]] static TypeDetectionResult detect_type(const core::Residue& residue);

private:
    /**
     * @brief Check if residue is in standard NT_LIST
     * @param res_name Residue name
     * @return true if in NT_LIST
     */
    [[nodiscard]] static bool is_in_nt_list(const std::string& res_name);

    /**
     * @brief Check if residue type is a purine
     * @param type Residue type
     * @return true if purine (A or G)
     */
    [[nodiscard]] static bool is_purine(core::ResidueType type);
};

} // namespace algorithms
} // namespace x3dna

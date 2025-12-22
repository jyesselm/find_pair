/**
 * @file base_pair_validator.hpp
 * @brief Base pair validation algorithm (matches legacy check_pair)
 */

#pragma once

#include <x3dna/core/residue.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>
#include <x3dna/algorithms/validation_constants.hpp>
#include <x3dna/algorithms/hydrogen_bond_finder.hpp>
#include <vector>
#include <optional>

namespace x3dna {
namespace algorithms {

/**
 * @struct ValidationParameters
 * @brief Parameters for base pair validation (matches legacy miscPars)
 */
struct ValidationParameters {
    // Distance constraints
    double min_dorg = 0.0;
    double max_dorg = 15.0;
    double min_dv = 0.0;
    double max_dv = 2.5;
    double min_dNN = 4.5;
    double max_dNN = 1e18; // XBIG (matches legacy)

    // Angle constraints
    double min_plane_angle = 0.0;
    double max_plane_angle = 65.0;

    // Hydrogen bond constraints
    int min_base_hb = 1;
    double hb_lower = 1.8;
    double hb_dist1 = 4.0;

    // H-bond atom list (default ".O.N" - matches legacy default)
    // Used to build hb_idx array for good_hbatoms() idx checks
    std::string hb_atoms = ".O.N";

    // Overlap threshold (matches legacy OVERLAP = 0.01)
    double overlap_threshold = 0.01; // OVERLAP constant from x3dna.h

    /**
     * @brief Set default values (matches set_default_misc_pars)
     */
    [[nodiscard]] static ValidationParameters defaults() {
        return ValidationParameters{};
    }
};

/**
 * @struct ValidationResult
 * @brief Result of base pair validation
 */
struct ValidationResult {
    bool is_valid = false;
    core::BasePairType bp_type = core::BasePairType::UNKNOWN;

    // Direction vectors (dot products of frame axes)
    double dir_x = 0.0; // x-axis dot product
    double dir_y = 0.0; // y-axis dot product
    double dir_z = 0.0; // z-axis dot product (most important)

    // Calculated values (matches rtn_val array)
    double dorg = 0.0;          // Distance between origins
    double d_v = 0.0;           // Vertical distance (abs(dot(dorg, zave)))
    double plane_angle = 0.0;   // Angle between z-axes (0-90 degrees)
    double dNN = 0.0;           // Distance between N1/N9 atoms
    double quality_score = 0.0; // Quality score (rtn_val[5])

    // Overlap area
    double overlap_area = 0.0;

    // Hydrogen bonds
    std::vector<core::hydrogen_bond> hbonds;
    int num_base_hb = 0;
    int num_o2_hb = 0;

    // Validation checks
    bool distance_check = false;
    bool d_v_check = false;
    bool plane_angle_check = false;
    bool dNN_check = false;
    bool overlap_check = false;
    bool hbond_check = false;
};

/**
 * @class BasePairValidator
 * @brief Validates base pairs using legacy check_pair algorithm
 *
 * This class implements the exact validation logic from the legacy check_pair function,
 * including distance checks, angle checks, overlap detection, and hydrogen bond validation.
 */
class BasePairValidator {
public:
    /**
     * @brief Constructor with default parameters
     */
    explicit BasePairValidator(const ValidationParameters& params = ValidationParameters::defaults())
        : params_(params) {}

    /**
     * @brief Validate a potential base pair
     * @param res1 First residue
     * @param res2 Second residue
     * @return ValidationResult with all validation details
     */
    [[nodiscard]] ValidationResult validate(const core::Residue& res1, const core::Residue& res2) const;

    /**
     * @brief Set validation parameters
     */
    void set_parameters(const ValidationParameters& params) {
        params_ = params;
    }

    /**
     * @brief Get validation parameters
     */
    [[nodiscard]] const ValidationParameters& parameters() const {
        return params_;
    }

    /**
     * @brief Calculate overlap area between two residues (matches legacy get_oarea)
     * @param res1 First residue
     * @param res2 Second residue
     * @param oave Average origin
     * @param zave Average z-axis
     * @return Overlap area in Angstrom²
     */
    [[nodiscard]] double calculate_overlap_area(const core::Residue& res1, const core::Residue& res2,
                                                const geometry::Vector3D& oave, const geometry::Vector3D& zave) const;

    /**
     * @brief Determine H-bond type based on donor-acceptor relationship (matches legacy
     * donor_acceptor)
     * @param base1 First base (A, C, G, I, T, U)
     * @param base2 Second base (A, C, G, I, T, U)
     * @param atom1 First atom name
     * @param atom2 Second atom name
     * @return '-' if valid standard H-bond, '*' otherwise
     */
    [[nodiscard]] static char donor_acceptor(char base1, char base2, const std::string& atom1,
                                             const std::string& atom2);

private:
    ValidationParameters params_;

    /**
     * @brief Pattern match function (matches legacy str_pmatch)
     * Checks if pattern matches string (where '.' in pattern matches any char)
     */
    [[nodiscard]] static bool pattern_match(const std::string& str, const std::string& pattern);

    /**
     * @brief Calculate direction vectors (dir_x, dir_y, dir_z)
     */
    static void calculate_direction_vectors(const core::ReferenceFrame& frame1, const core::ReferenceFrame& frame2,
                                            double& dir_x, double& dir_y, double& dir_z);

    /**
     * @brief Calculate average z-axis and origin
     */
    static void get_bp_zoave(const core::ReferenceFrame& frame1, const core::ReferenceFrame& frame2,
                             geometry::Vector3D& oave, geometry::Vector3D& zave);

    /**
     * @brief Calculate angle between two z-axes (0-90 degrees)
     */
    [[nodiscard]] static double z1_z2_angle_in_0_to_90(const geometry::Vector3D& z1, const geometry::Vector3D& z2);

    /**
     * @brief Check if value is in range
     */
    [[nodiscard]] static bool in_range(double value, double min_val, double max_val) {
        return value >= min_val && value <= max_val;
    }

    /**
     * @brief Find N1/N9 atoms for dNN calculation
     */
    [[nodiscard]] static std::optional<geometry::Vector3D> find_n1_n9_position(const core::Residue& residue);

    /**
     * @brief Find hydrogen bonds between two residues (with validation)
     * Used for adjust_pairQuality - matches hb_numlist behavior
     */
    [[nodiscard]] std::vector<core::hydrogen_bond> find_hydrogen_bonds(const core::Residue& res1,
                                                                       const core::Residue& res2) const;

    /**
     * @brief Count hydrogen bonds simply (before validation) - matches legacy check_pair behavior
     */
    void count_hydrogen_bonds_simple(const core::Residue& res1, const core::Residue& res2, int& num_base_hb,
                                     int& num_o2_hb) const;
};

} // namespace algorithms
} // namespace x3dna

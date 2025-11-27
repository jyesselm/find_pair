/**
 * @file parameter_calculator.hpp
 * @brief ParameterCalculator class for calculating base pair step and helical parameters
 */

#pragma once

#include <vector>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/parameters.hpp>
#include <x3dna/core/reference_frame.hpp>

namespace x3dna {
namespace algorithms {

/**
 * @class ParameterCalculator
 * @brief Calculates base pair step parameters and helical parameters
 *
 * This class implements the bpstep_par algorithm from the legacy code,
 * which calculates the 6 base pair step parameters (Shift, Slide, Rise, Tilt, Roll, Twist)
 * from two base reference frames.
 */
class ParameterCalculator {
public:
    /**
     * @brief Calculate step parameters from two consecutive base pairs
     * @param pair1 First base pair
     * @param pair2 Second base pair
     * @return BasePairStepParameters structure
     */
    core::BasePairStepParameters calculate_step_parameters(const core::BasePair& pair1,
                                                           const core::BasePair& pair2) const;

    /**
     * @brief Calculate step parameters for a single base pair (for bp_type_id calculation)
     * @param pair Base pair
     * @return BasePairStepParameters structure
     * 
     * This matches legacy: bpstep_par(r2, org[j], r1, org[i], ...)
     * where r2/org[j] is frame2 and r1/org[i] is frame1
     */
    core::BasePairStepParameters calculate_step_parameters_for_pair(const core::BasePair& pair);

    /**
     * @brief Calculate step parameters from two reference frames
     * @param frame1 First reference frame
     * @param frame2 Second reference frame
     * @return BasePairStepParameters structure
     */
    core::BasePairStepParameters calculate_step_parameters(const core::ReferenceFrame& frame1,
                                                           const core::ReferenceFrame& frame2) const;

    /**
     * @brief Calculate helical parameters from two base pairs
     * @param pair1 First base pair
     * @param pair2 Second base pair
     * @return HelicalParameters structure
     */
    core::HelicalParameters calculate_helical_parameters(const core::BasePair& pair1,
                                                         const core::BasePair& pair2);

    /**
     * @brief Calculate helical parameters from two reference frames
     * @param frame1 First reference frame
     * @param frame2 Second reference frame
     * @return HelicalParameters structure
     */
    core::HelicalParameters calculate_helical_parameters_impl(const core::ReferenceFrame& frame1,
                                                              const core::ReferenceFrame& frame2);

    /**
     * @brief Calculate step parameters for all consecutive pairs
     * @param pairs Vector of base pairs
     * @return Vector of step parameters (one per consecutive pair)
     */
    std::vector<core::BasePairStepParameters> calculate_all_step_parameters(
        const std::vector<core::BasePair>& pairs);

    /**
     * @brief Calculate midstep reference frame
     * @param frame1 First reference frame
     * @param frame2 Second reference frame
     * @return Midstep reference frame
     */
    core::ReferenceFrame calculate_midstep_frame(const core::ReferenceFrame& frame1,
                                                  const core::ReferenceFrame& frame2);

private:
    /**
     * @brief Core algorithm matching legacy bpstep_par
     * @param r1 Rotation matrix of first frame
     * @param o1 Origin of first frame
     * @param r2 Rotation matrix of second frame
     * @param o2 Origin of second frame
     * @param params Output: step parameters
     * @param midstep_frame Output: midstep reference frame
     */
    void bpstep_par_impl(const geometry::Matrix3D& r1, const geometry::Vector3D& o1,
                         const geometry::Matrix3D& r2, const geometry::Vector3D& o2,
                         core::BasePairStepParameters& params,
                         core::ReferenceFrame& midstep_frame) const;

    // Geometry utility functions (matching legacy implementations)
    static double magang(const geometry::Vector3D& va, const geometry::Vector3D& vb);
    static double vec_ang(const geometry::Vector3D& va, const geometry::Vector3D& vb,
                          const geometry::Vector3D& vref);
    static geometry::Matrix3D arb_rotation(const geometry::Vector3D& axis, double angle_deg);
    static geometry::Vector3D get_vector(const geometry::Vector3D& va,
                                         const geometry::Vector3D& vref, double deg_ang);
    static geometry::Matrix3D x_y_z_2_mtx(const geometry::Vector3D& x, const geometry::Vector3D& y,
                                          const geometry::Vector3D& z);
    static double deg2rad(double ang);
    static double rad2deg(double ang);
    static constexpr double XEPS = 1.0e-10; // Small epsilon for comparisons
    static constexpr double PI = 3.14159265358979323846;
};

} // namespace algorithms
} // namespace x3dna


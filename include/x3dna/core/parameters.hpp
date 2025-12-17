/**
 * @file parameters.hpp
 * @brief Parameter structures for base pair step and helical parameters
 */

#pragma once

#include <array>
#include <cmath>
#include <nlohmann/json.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>

namespace x3dna {
namespace core {

/**
 * @struct BasePairStepParameters
 * @brief Base pair step parameters (6 parameters: Shift, Slide, Rise, Tilt, Roll, Twist)
 *
 * These parameters describe the relative orientation and position of two consecutive
 * base pairs in a nucleic acid structure. The parameters are calculated using the
 * standard 3DNA convention.
 */
struct BasePairStepParameters {
    double shift = 0.0; // x-displacement (Angstroms)
    double slide = 0.0; // y-displacement (Angstroms)
    double rise = 0.0;  // z-displacement (Angstroms)
    double tilt = 0.0;  // rotation about x-axis (degrees)
    double roll = 0.0;  // rotation about y-axis (degrees)
    double twist = 0.0; // rotation about z-axis (degrees)

    // Midstep reference frame (calculated during parameter computation)
    std::optional<ReferenceFrame> midstep_frame;

    /**
     * @brief Default constructor
     */
    BasePairStepParameters() = default;

    /**
     * @brief Constructor with all parameters
     */
    BasePairStepParameters(double s, double sl, double r, double ti, double ro, double tw)
        : shift(s), slide(sl), rise(r), tilt(ti), roll(ro), twist(tw) {}

    /**
     * @brief Convert parameters to array format
     * @return Array of 6 parameters: [shift, slide, rise, tilt, roll, twist]
     */
    [[nodiscard]] std::array<double, 6> as_array() const {
        return {shift, slide, rise, tilt, roll, twist};
    }

    /**
     * @brief Create from array format
     * @param arr Array of 6 parameters: [shift, slide, rise, tilt, roll, twist]
     */
    [[nodiscard]] static BasePairStepParameters from_array(const std::array<double, 6>& arr) {
        BasePairStepParameters params;
        params.shift = arr[0];
        params.slide = arr[1];
        params.rise = arr[2];
        params.tilt = arr[3];
        params.roll = arr[4];
        params.twist = arr[5];
        return params;
    }

    /**
     * @brief Equality comparison
     */
    [[nodiscard]] bool operator==(const BasePairStepParameters& other) const {
        return shift == other.shift && slide == other.slide && rise == other.rise && tilt == other.tilt &&
               roll == other.roll && twist == other.twist;
    }

    /**
     * @brief Approximate equality comparison (within tolerance)
     * @param other Other parameters to compare
     * @param tolerance Tolerance for comparison (default: 1e-6)
     */
    [[nodiscard]] bool approximately_equal(const BasePairStepParameters& other, double tolerance = 1e-6) const {
        return std::abs(shift - other.shift) < tolerance && std::abs(slide - other.slide) < tolerance &&
               std::abs(rise - other.rise) < tolerance && std::abs(tilt - other.tilt) < tolerance &&
               std::abs(roll - other.roll) < tolerance && std::abs(twist - other.twist) < tolerance;
    }

    /**
     * @brief Convert to modern JSON format
     */
    [[nodiscard]] nlohmann::json to_json() const {
        nlohmann::json j;
        j["shift"] = shift;
        j["slide"] = slide;
        j["rise"] = rise;
        j["tilt"] = tilt;
        j["roll"] = roll;
        j["twist"] = twist;
        if (midstep_frame.has_value()) {
            j["midstep_frame"] = midstep_frame->to_json();
        }
        return j;
    }

    /**
     * @brief Create from modern JSON format
     */
    [[nodiscard]] static BasePairStepParameters from_json(const nlohmann::json& j) {
        BasePairStepParameters params;
        params.shift = j.value("shift", 0.0);
        params.slide = j.value("slide", 0.0);
        params.rise = j.value("rise", 0.0);
        params.tilt = j.value("tilt", 0.0);
        params.roll = j.value("roll", 0.0);
        params.twist = j.value("twist", 0.0);
        if (j.contains("midstep_frame")) {
            params.midstep_frame = ReferenceFrame::from_json(j["midstep_frame"]);
        }
        return params;
    }

    /**
     * @brief Convert to legacy JSON format (bpstep_params record)
     * Format: {"type": "bpstep_params", "bp_idx1": ..., "bp_idx2": ...,
     *          "params": {"Shift": ..., "Slide": ..., ...},
     *          "mst_org": [...], "mst_orien": [[...], [...], [...]]}
     */
    [[nodiscard]] nlohmann::json to_json_legacy(size_t bp_idx1 = 0, size_t bp_idx2 = 0) const {
        nlohmann::json j;
        j["type"] = "bpstep_params";
        j["bp_idx1"] = static_cast<long>(bp_idx1);
        j["bp_idx2"] = static_cast<long>(bp_idx2);
        j["params"]["Shift"] = shift;
        j["params"]["Slide"] = slide;
        j["params"]["Rise"] = rise;
        j["params"]["Tilt"] = tilt;
        j["params"]["Roll"] = roll;
        j["params"]["Twist"] = twist;

        if (midstep_frame.has_value()) {
            j["mst_org"] = midstep_frame->origin().to_json();
            j["mst_orien"] = midstep_frame->rotation().to_json_legacy();
        }

        return j;
    }

    /**
     * @brief Create from legacy JSON format (bpstep_params record)
     */
    [[nodiscard]] static BasePairStepParameters from_json_legacy(const nlohmann::json& j) {
        BasePairStepParameters params;

        if (j.contains("params") && j["params"].is_object()) {
            const auto& p = j["params"];
            params.shift = p.value("Shift", 0.0);
            params.slide = p.value("Slide", 0.0);
            params.rise = p.value("Rise", 0.0);
            params.tilt = p.value("Tilt", 0.0);
            params.roll = p.value("Roll", 0.0);
            params.twist = p.value("Twist", 0.0);
        }

        // Parse midstep frame if present
        if (j.contains("mst_org") && j.contains("mst_orien")) {
            nlohmann::json frame_json;
            frame_json["org"] = j["mst_org"];
            frame_json["orien"] = j["mst_orien"];
            params.midstep_frame = ReferenceFrame::from_json_legacy(frame_json);
        }

        return params;
    }
};

/**
 * @struct HelicalParameters
 * @brief Helical parameters (6 parameters: x_displacement, y_displacement, rise, inclination, tip,
 * twist)
 *
 * These parameters describe the helical geometry of a base pair step in a nucleic
 * acid structure. The parameters are calculated using the standard 3DNA convention.
 */
struct HelicalParameters {
    double x_displacement = 0.0; // x-displacement (Angstroms)
    double y_displacement = 0.0; // y-displacement (Angstroms)
    double rise = 0.0;           // z-displacement (Angstroms)
    double inclination = 0.0;    // inclination angle (degrees)
    double tip = 0.0;            // tip angle (degrees)
    double twist = 0.0;          // twist angle (degrees)

    // Helical midstep reference frame (calculated during parameter computation)
    std::optional<ReferenceFrame> midstep_frame;

    /**
     * @brief Default constructor
     */
    HelicalParameters() = default;

    /**
     * @brief Constructor with all parameters
     */
    HelicalParameters(double xd, double yd, double r, double inc, double t, double tw)
        : x_displacement(xd), y_displacement(yd), rise(r), inclination(inc), tip(t), twist(tw) {}

    /**
     * @brief Convert parameters to array format
     * @return Array of 6 parameters: [x_displacement, y_displacement, rise, inclination, tip,
     * twist]
     */
    [[nodiscard]] std::array<double, 6> as_array() const {
        return {x_displacement, y_displacement, rise, inclination, tip, twist};
    }

    /**
     * @brief Create from array format
     * @param arr Array of 6 parameters: [x_displacement, y_displacement, rise, inclination, tip,
     * twist]
     */
    [[nodiscard]] static HelicalParameters from_array(const std::array<double, 6>& arr) {
        HelicalParameters params;
        params.x_displacement = arr[0];
        params.y_displacement = arr[1];
        params.rise = arr[2];
        params.inclination = arr[3];
        params.tip = arr[4];
        params.twist = arr[5];
        return params;
    }

    /**
     * @brief Equality comparison
     */
    [[nodiscard]] bool operator==(const HelicalParameters& other) const {
        return x_displacement == other.x_displacement && y_displacement == other.y_displacement && rise == other.rise &&
               inclination == other.inclination && tip == other.tip && twist == other.twist;
    }

    /**
     * @brief Approximate equality comparison (within tolerance)
     * @param other Other parameters to compare
     * @param tolerance Tolerance for comparison (default: 1e-6)
     */
    [[nodiscard]] bool approximately_equal(const HelicalParameters& other, double tolerance = 1e-6) const {
        return std::abs(x_displacement - other.x_displacement) < tolerance &&
               std::abs(y_displacement - other.y_displacement) < tolerance && std::abs(rise - other.rise) < tolerance &&
               std::abs(inclination - other.inclination) < tolerance && std::abs(tip - other.tip) < tolerance &&
               std::abs(twist - other.twist) < tolerance;
    }

    /**
     * @brief Convert to modern JSON format
     */
    [[nodiscard]] nlohmann::json to_json() const {
        nlohmann::json j;
        j["x_displacement"] = x_displacement;
        j["y_displacement"] = y_displacement;
        j["rise"] = rise;
        j["inclination"] = inclination;
        j["tip"] = tip;
        j["twist"] = twist;
        if (midstep_frame.has_value()) {
            j["midstep_frame"] = midstep_frame->to_json();
        }
        return j;
    }

    /**
     * @brief Create from modern JSON format
     */
    [[nodiscard]] static HelicalParameters from_json(const nlohmann::json& j) {
        HelicalParameters params;
        params.x_displacement = j.value("x_displacement", 0.0);
        params.y_displacement = j.value("y_displacement", 0.0);
        params.rise = j.value("rise", 0.0);
        params.inclination = j.value("inclination", 0.0);
        params.tip = j.value("tip", 0.0);
        params.twist = j.value("twist", 0.0);
        if (j.contains("midstep_frame")) {
            params.midstep_frame = ReferenceFrame::from_json(j["midstep_frame"]);
        }
        return params;
    }

    /**
     * @brief Convert to legacy JSON format (helical_params record)
     * Format: {"type": "helical_params", "bp_idx1": ..., "bp_idx2": ...,
     *          "params": [x_displacement, y_displacement, rise, inclination, tip, twist],
     *          "mst_orgH": [...], "mst_orienH": [[...], [...], [...]]}
     */
    [[nodiscard]] nlohmann::json to_json_legacy(size_t bp_idx1 = 0, size_t bp_idx2 = 0) const {
        nlohmann::json j;
        j["type"] = "helical_params";
        j["bp_idx1"] = static_cast<long>(bp_idx1);
        j["bp_idx2"] = static_cast<long>(bp_idx2);
        j["params"] = as_array();

        if (midstep_frame.has_value()) {
            j["mst_orgH"] = midstep_frame->origin().to_json();
            j["mst_orienH"] = midstep_frame->rotation().to_json_legacy();
        }

        return j;
    }

    /**
     * @brief Create from legacy JSON format (helical_params record)
     */
    [[nodiscard]] static HelicalParameters from_json_legacy(const nlohmann::json& j) {
        HelicalParameters params;

        if (j.contains("params") && j["params"].is_array() && j["params"].size() >= 6) {
            const auto& p = j["params"];
            params.x_displacement = p[0].get<double>();
            params.y_displacement = p[1].get<double>();
            params.rise = p[2].get<double>();
            params.inclination = p[3].get<double>();
            params.tip = p[4].get<double>();
            params.twist = p[5].get<double>();
        }

        // Parse helical midstep frame if present
        if (j.contains("mst_orgH") && j.contains("mst_orienH")) {
            nlohmann::json frame_json;
            frame_json["org"] = j["mst_orgH"];
            frame_json["orien"] = j["mst_orienH"];
            params.midstep_frame = ReferenceFrame::from_json_legacy(frame_json);
        }

        return params;
    }
};

} // namespace core
} // namespace x3dna

/**
 * @file reference_frame.hpp
 * @brief ReferenceFrame class representing a local coordinate frame
 */

#pragma once

#include <array>
#include <nlohmann/json.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>

namespace x3dna {
namespace core {

/**
 * @class ReferenceFrame
 * @brief Represents a local coordinate frame with rotation matrix and origin
 *
 * The reference frame is used to represent the orientation and position of
 * nucleic acid bases. The rotation matrix (orien) defines the orientation,
 * and the origin (org) defines the position.
 */
class ReferenceFrame {
public:
    /**
     * @brief Default constructor (identity rotation, zero origin)
     */
    ReferenceFrame() : rotation_(geometry::Matrix3D::identity()), origin_(geometry::Vector3D(0, 0, 0)) {}

    /**
     * @brief Constructor with rotation matrix and origin
     * @param rotation 3x3 rotation matrix (orientation)
     * @param origin 3D origin vector (position)
     */
    ReferenceFrame(const geometry::Matrix3D& rotation, const geometry::Vector3D& origin)
        : rotation_(rotation), origin_(origin) {}

    /**
     * @brief Constructor from 9-element rotation array and 3-element origin array
     * @param rotation_array Array of 9 elements (row-major: [r00, r01, r02, r10, r11, r12, r20,
     * r21, r22])
     * @param origin_array Array of 3 elements [x, y, z]
     */
    ReferenceFrame(const std::array<double, 9>& rotation_array, const std::array<double, 3>& origin_array)
        : rotation_(geometry::Matrix3D(rotation_array[0], rotation_array[1], rotation_array[2], rotation_array[3],
                                       rotation_array[4], rotation_array[5], rotation_array[6], rotation_array[7],
                                       rotation_array[8])),
          origin_(geometry::Vector3D(origin_array[0], origin_array[1], origin_array[2])) {}

    // Getters
    [[nodiscard]] const geometry::Matrix3D& rotation() const {
        return rotation_;
    }
    [[nodiscard]] const geometry::Vector3D& origin() const {
        return origin_;
    }

    /**
     * @brief Get x-axis (first column of rotation matrix)
     */
    [[nodiscard]] geometry::Vector3D x_axis() const {
        return geometry::Vector3D(rotation_.at(0, 0), rotation_.at(1, 0), rotation_.at(2, 0));
    }

    /**
     * @brief Get y-axis (second column of rotation matrix)
     */
    [[nodiscard]] geometry::Vector3D y_axis() const {
        return geometry::Vector3D(rotation_.at(0, 1), rotation_.at(1, 1), rotation_.at(2, 1));
    }

    /**
     * @brief Get z-axis (third column of rotation matrix)
     *
     * The z-axis is the normal to the base plane in nucleic acids.
     */
    [[nodiscard]] geometry::Vector3D z_axis() const {
        return geometry::Vector3D(rotation_.at(0, 2), rotation_.at(1, 2), rotation_.at(2, 2));
    }

    /**
     * @brief Calculate direction dot product with another frame
     *
     * This is used to validate base pairs - the z-axes should point
     * in opposite directions (dot product should be negative).
     *
     * @param other Another reference frame
     * @return Dot product of z-axes
     */
    [[nodiscard]] double direction_dot_product(const ReferenceFrame& other) const {
        return z_axis().dot(other.z_axis());
    }

    /**
     * @brief Transform a point from local to global coordinates
     * @param local_point Point in local frame
     * @return Point in global coordinates
     */
    [[nodiscard]] geometry::Vector3D transform(const geometry::Vector3D& local_point) const {
        return rotation_ * local_point + origin_;
    }

    /**
     * @brief Transform a point from global to local coordinates
     * @param global_point Point in global coordinates
     * @return Point in local frame
     */
    [[nodiscard]] geometry::Vector3D inverse_transform(const geometry::Vector3D& global_point) const {
        geometry::Vector3D translated = global_point - origin_;
        return rotation_.transpose() * translated;
    }

    /**
     * @brief Get rotation matrix as 9-element array (row-major)
     */
    [[nodiscard]] std::array<double, 9> rotation_as_array() const {
        return {rotation_.at(0, 0), rotation_.at(0, 1), rotation_.at(0, 2), rotation_.at(1, 0), rotation_.at(1, 1),
                rotation_.at(1, 2), rotation_.at(2, 0), rotation_.at(2, 1), rotation_.at(2, 2)};
    }

    /**
     * @brief Get origin as 3-element array
     */
    [[nodiscard]] std::array<double, 3> origin_as_array() const {
        return {origin_.x(), origin_.y(), origin_.z()};
    }

    /**
     * @brief Convert to legacy JSON format
     *
     * Legacy format uses "orien" (3x3 nested array) and "org" (3-element array)
     */
    [[nodiscard]] nlohmann::json to_json_legacy() const {
        nlohmann::json j;

        // orien: 3x3 nested array (row-major)
        j["orien"] = nlohmann::json::array();
        for (int i = 0; i < 3; ++i) {
            nlohmann::json row = nlohmann::json::array();
            for (int j_col = 0; j_col < 3; ++j_col) {
                row.push_back(rotation_.at(i, j_col));
            }
            j["orien"].push_back(row);
        }

        // org: 3-element array
        j["org"] = {origin_.x(), origin_.y(), origin_.z()};

        return j;
    }

    /**
     * @brief Create ReferenceFrame from legacy JSON format
     * @param j JSON object with "orien" and "org" fields
     * @return ReferenceFrame object
     */
    [[nodiscard]] static ReferenceFrame from_json_legacy(const nlohmann::json& j) {
        // Parse orien (3x3 nested array)
        std::array<double, 9> rotation_array;
        if (j.contains("orien") && j["orien"].is_array() && j["orien"].size() == 3) {
            int idx = 0;
            for (int i = 0; i < 3; ++i) {
                if (j["orien"][i].is_array() && j["orien"][i].size() == 3) {
                    for (int j_col = 0; j_col < 3; ++j_col) {
                        rotation_array[idx++] = j["orien"][i][j_col].get<double>();
                    }
                } else {
                    throw std::invalid_argument("Invalid orien format in JSON");
                }
            }
        } else {
            throw std::invalid_argument("Missing or invalid orien in JSON");
        }

        // Parse org (3-element array)
        std::array<double, 3> origin_array;
        if (j.contains("org") && j["org"].is_array() && j["org"].size() == 3) {
            origin_array[0] = j["org"][0].get<double>();
            origin_array[1] = j["org"][1].get<double>();
            origin_array[2] = j["org"][2].get<double>();
        } else {
            throw std::invalid_argument("Missing or invalid org in JSON");
        }

        return ReferenceFrame(rotation_array, origin_array);
    }

    /**
     * @brief Convert to modern JSON format
     */
    [[nodiscard]] nlohmann::json to_json() const {
        nlohmann::json j;
        j["rotation"] = rotation_.to_json_legacy(); // Use nested array format
        j["origin"] = origin_.to_json();
        return j;
    }

    /**
     * @brief Create ReferenceFrame from modern JSON format
     */
    [[nodiscard]] static ReferenceFrame from_json(const nlohmann::json& j) {
        geometry::Matrix3D rotation = geometry::Matrix3D::from_json_legacy(j["rotation"]);
        geometry::Vector3D origin = geometry::Vector3D::from_json(j["origin"]);
        return ReferenceFrame(rotation, origin);
    }

private:
    geometry::Matrix3D rotation_; // 3x3 rotation matrix (orientation)
    geometry::Vector3D origin_;   // 3D origin vector (position)
};

} // namespace core
} // namespace x3dna

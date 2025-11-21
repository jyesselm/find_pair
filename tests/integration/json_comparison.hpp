/**
 * @file json_comparison.hpp
 * @brief Utilities for comparing JSON data with legacy format
 */

#pragma once

#include <nlohmann/json.hpp>
#include <cmath>
#include <gtest/gtest.h>

namespace x3dna::test {

constexpr double DEFAULT_TOLERANCE = 0.001;

/**
 * @brief Check if two double values are approximately equal
 */
inline bool approximately_equal(double a, double b, double tolerance = DEFAULT_TOLERANCE) {
    return std::abs(a - b) < tolerance;
}

/**
 * @brief Compare two JSON arrays of doubles
 */
inline bool compare_double_array(const nlohmann::json& arr1, const nlohmann::json& arr2,
                                 double tolerance = DEFAULT_TOLERANCE) {
    if (!arr1.is_array() || !arr2.is_array()) {
        return false;
    }

    if (arr1.size() != arr2.size()) {
        return false;
    }

    for (size_t i = 0; i < arr1.size(); ++i) {
        if (!approximately_equal(arr1[i].get<double>(), arr2[i].get<double>(), tolerance)) {
            return false;
        }
    }

    return true;
}

/**
 * @brief Compare rotation matrices (3x3 nested arrays)
 * @note Will be implemented when Matrix3D is available
 */
// bool compare_rotation_matrix(const nlohmann::json& legacy_matrix,
//                              const Matrix3D& our_matrix,
//                              double tolerance = DEFAULT_TOLERANCE);

/**
 * @brief Compare vectors (3-element arrays)
 * @note Will be implemented when Vector3D is available
 */
// bool compare_vector(const nlohmann::json& legacy_vector,
//                    const Vector3D& our_vector,
//                    double tolerance = DEFAULT_TOLERANCE);

} // namespace x3dna::test

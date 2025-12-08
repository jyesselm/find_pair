/**
 * @file matrix3d.hpp
 * @brief 3x3 matrix class for geometric transformations
 */

#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <nlohmann/json.hpp>
#include "vector3d.hpp"

namespace x3dna {
namespace geometry {

/**
 * @class Matrix3D
 * @brief Represents a 3x3 matrix for rotations and transformations
 *
 * Matrix is stored in row-major order: [r11, r12, r13, r21, r22, r23, r31, r32, r33]
 */
class Matrix3D {
public:
    /**
     * @brief Default constructor (creates identity matrix)
     */
    Matrix3D() {
        data_ = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    }

    /**
     * @brief Constructor from array (row-major order)
     */
    explicit Matrix3D(const std::array<double, 9>& values) : data_(values) {}

    /**
     * @brief Constructor from individual elements (row-major)
     */
    Matrix3D(double m00, double m01, double m02, double m10, double m11, double m12, double m20, double m21,
             double m22) {
        data_ = {m00, m01, m02, m10, m11, m12, m20, m21, m22};
    }

    /**
     * @brief Get element at row i, column j (0-indexed)
     */
    double at(size_t i, size_t j) const {
        if (i >= 3 || j >= 3) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data_[i * 3 + j];
    }

    /**
     * @brief Set element at row i, column j (0-indexed)
     */
    void set(size_t i, size_t j, double value) {
        if (i >= 3 || j >= 3) {
            throw std::out_of_range("Matrix index out of range");
        }
        data_[i * 3 + j] = value;
    }

    /**
     * @brief Get row i (0-indexed)
     */
    Vector3D row(size_t i) const {
        if (i >= 3) {
            throw std::out_of_range("Row index out of range");
        }
        return Vector3D(data_[i * 3], data_[i * 3 + 1], data_[i * 3 + 2]);
    }

    /**
     * @brief Get column j (0-indexed)
     */
    Vector3D column(size_t j) const {
        if (j >= 3) {
            throw std::out_of_range("Column index out of range");
        }
        return Vector3D(data_[j], data_[3 + j], data_[6 + j]);
    }

    /**
     * @brief Set row i (0-indexed)
     */
    void set_row(size_t i, const Vector3D& vec) {
        if (i >= 3) {
            throw std::out_of_range("Row index out of range");
        }
        data_[i * 3] = vec.x();
        data_[i * 3 + 1] = vec.y();
        data_[i * 3 + 2] = vec.z();
    }

    /**
     * @brief Set column j (0-indexed)
     */
    void set_column(size_t j, const Vector3D& vec) {
        if (j >= 3) {
            throw std::out_of_range("Column index out of range");
        }
        data_[j] = vec.x();
        data_[3 + j] = vec.y();
        data_[6 + j] = vec.z();
    }

    /**
     * @brief Matrix-vector multiplication
     */
    Vector3D operator*(const Vector3D& vec) const {
        return Vector3D(data_[0] * vec.x() + data_[1] * vec.y() + data_[2] * vec.z(),
                        data_[3] * vec.x() + data_[4] * vec.y() + data_[5] * vec.z(),
                        data_[6] * vec.x() + data_[7] * vec.y() + data_[8] * vec.z());
    }

    /**
     * @brief Matrix-matrix multiplication
     */
    Matrix3D operator*(const Matrix3D& other) const {
        Matrix3D result;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                double sum = 0.0;
                for (size_t k = 0; k < 3; ++k) {
                    sum += at(i, k) * other.at(k, j);
                }
                result.set(i, j, sum);
            }
        }
        return result;
    }

    /**
     * @brief Matrix addition
     */
    Matrix3D operator+(const Matrix3D& other) const {
        Matrix3D result;
        for (size_t i = 0; i < 9; ++i) {
            result.data_[i] = data_[i] + other.data_[i];
        }
        return result;
    }

    /**
     * @brief Matrix subtraction
     */
    Matrix3D operator-(const Matrix3D& other) const {
        Matrix3D result;
        for (size_t i = 0; i < 9; ++i) {
            result.data_[i] = data_[i] - other.data_[i];
        }
        return result;
    }

    /**
     * @brief Scalar multiplication
     */
    Matrix3D operator*(double scalar) const {
        Matrix3D result;
        for (size_t i = 0; i < 9; ++i) {
            result.data_[i] = data_[i] * scalar;
        }
        return result;
    }

    /**
     * @brief Scalar division
     */
    Matrix3D operator/(double scalar) const {
        return *this * (1.0 / scalar);
    }

    /**
     * @brief Transpose matrix
     */
    Matrix3D transpose() const {
        Matrix3D result;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                result.set(i, j, at(j, i));
            }
        }
        return result;
    }

    /**
     * @brief Calculate determinant
     */
    double determinant() const {
        return data_[0] * (data_[4] * data_[8] - data_[5] * data_[7]) -
               data_[1] * (data_[3] * data_[8] - data_[5] * data_[6]) +
               data_[2] * (data_[3] * data_[7] - data_[4] * data_[6]);
    }

    /**
     * @brief Calculate inverse matrix
     * @throws std::runtime_error if matrix is singular (determinant == 0)
     */
    Matrix3D inverse() const {
        double det = determinant();
        if (std::abs(det) < 1e-9) {
            throw std::runtime_error("Matrix is singular (determinant is zero)");
        }

        Matrix3D result;
        double inv_det = 1.0 / det;

        // Calculate adjugate matrix (transpose of cofactor matrix)
        result.set(0, 0, (data_[4] * data_[8] - data_[5] * data_[7]) * inv_det);
        result.set(0, 1, (data_[2] * data_[7] - data_[1] * data_[8]) * inv_det);
        result.set(0, 2, (data_[1] * data_[5] - data_[2] * data_[4]) * inv_det);
        result.set(1, 0, (data_[5] * data_[6] - data_[3] * data_[8]) * inv_det);
        result.set(1, 1, (data_[0] * data_[8] - data_[2] * data_[6]) * inv_det);
        result.set(1, 2, (data_[2] * data_[3] - data_[0] * data_[5]) * inv_det);
        result.set(2, 0, (data_[3] * data_[7] - data_[4] * data_[6]) * inv_det);
        result.set(2, 1, (data_[1] * data_[6] - data_[0] * data_[7]) * inv_det);
        result.set(2, 2, (data_[0] * data_[4] - data_[1] * data_[3]) * inv_det);

        return result;
    }

    /**
     * @brief Check if matrix is approximately equal to another
     */
    bool approximately_equals(const Matrix3D& other, double tolerance = 1e-6) const {
        for (size_t i = 0; i < 9; ++i) {
            if (std::abs(data_[i] - other.data_[i]) >= tolerance) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Convert to array (row-major order)
     */
    std::array<double, 9> as_array() const {
        return data_;
    }

    /**
     * @brief Convert to JSON array (row-major: [r11, r12, r13, r21, r22, r23, r31, r32, r33])
     */
    nlohmann::json to_json() const {
        nlohmann::json result = nlohmann::json::array();
        for (double val : data_) {
            result.push_back(val);
        }
        return result;
    }

    /**
     * @brief Convert to legacy JSON format (3x3 nested array: [[r11, r12, r13], [r21, r22, r23],
     * [r31, r32, r33]])
     */
    nlohmann::json to_json_legacy() const {
        nlohmann::json result = nlohmann::json::array();
        for (size_t i = 0; i < 3; ++i) {
            nlohmann::json row = nlohmann::json::array();
            for (size_t j = 0; j < 3; ++j) {
                row.push_back(at(i, j));
            }
            result.push_back(row);
        }
        return result;
    }

    /**
     * @brief Create from JSON array (row-major)
     */
    static Matrix3D from_json(const nlohmann::json& json) {
        if (!json.is_array() || json.size() != 9) {
            throw std::invalid_argument("JSON must be array of 9 numbers");
        }
        std::array<double, 9> arr;
        for (size_t i = 0; i < 9; ++i) {
            arr[i] = json[i].get<double>();
        }
        return Matrix3D(arr);
    }

    /**
     * @brief Create from legacy JSON format (3x3 nested array)
     */
    static Matrix3D from_json_legacy(const nlohmann::json& json) {
        if (!json.is_array() || json.size() != 3) {
            throw std::invalid_argument("JSON must be array of 3 arrays");
        }
        Matrix3D result;
        for (size_t i = 0; i < 3; ++i) {
            if (!json[i].is_array() || json[i].size() != 3) {
                throw std::invalid_argument("Each row must be array of 3 numbers");
            }
            for (size_t j = 0; j < 3; ++j) {
                result.set(i, j, json[i][j].get<double>());
            }
        }
        return result;
    }

    /**
     * @brief Create identity matrix
     */
    static Matrix3D identity() {
        return Matrix3D();
    }

    /**
     * @brief Create rotation matrix around X axis (in radians)
     */
    static Matrix3D rotation_x(double angle) {
        double c = std::cos(angle);
        double s = std::sin(angle);
        return Matrix3D(1.0, 0.0, 0.0, 0.0, c, -s, 0.0, s, c);
    }

    /**
     * @brief Create rotation matrix around Y axis (in radians)
     */
    static Matrix3D rotation_y(double angle) {
        double c = std::cos(angle);
        double s = std::sin(angle);
        return Matrix3D(c, 0.0, s, 0.0, 1.0, 0.0, -s, 0.0, c);
    }

    /**
     * @brief Create rotation matrix around Z axis (in radians)
     */
    static Matrix3D rotation_z(double angle) {
        double c = std::cos(angle);
        double s = std::sin(angle);
        return Matrix3D(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0);
    }

private:
    std::array<double, 9> data_; // Row-major: [r11, r12, r13, r21, r22, r23, r31, r32, r33]
};

// Scalar multiplication from left
inline Matrix3D operator*(double scalar, const Matrix3D& mat) {
    return mat * scalar;
}

} // namespace geometry
} // namespace x3dna

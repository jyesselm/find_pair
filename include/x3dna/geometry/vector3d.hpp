/**
 * @file vector3d.hpp
 * @brief 3D vector class for geometric calculations
 */

#pragma once

#include <array>
#include <cmath>
#include <nlohmann/json.hpp>

namespace x3dna {
namespace geometry {

/**
 * @class Vector3D
 * @brief Represents a 3D vector with x, y, z components
 */
class Vector3D {
public:
    /**
     * @brief Default constructor (creates zero vector)
     */
    Vector3D() : x_(0.0), y_(0.0), z_(0.0) {}

    /**
     * @brief Constructor with x, y, z components
     */
    Vector3D(double x, double y, double z) : x_(x), y_(y), z_(z) {}

    /**
     * @brief Constructor from array
     */
    explicit Vector3D(const std::array<double, 3>& arr) : x_(arr[0]), y_(arr[1]), z_(arr[2]) {}

    // Getters
    double x() const {
        return x_;
    }
    double y() const {
        return y_;
    }
    double z() const {
        return z_;
    }

    // Setters
    void set_x(double x) {
        x_ = x;
    }
    void set_y(double y) {
        y_ = y;
    }
    void set_z(double z) {
        z_ = z;
    }
    void set(double x, double y, double z) {
        x_ = x;
        y_ = y;
        z_ = z;
    }

    // Arithmetic operations
    Vector3D operator+(const Vector3D& other) const {
        return Vector3D(x_ + other.x_, y_ + other.y_, z_ + other.z_);
    }

    Vector3D operator-(const Vector3D& other) const {
        return Vector3D(x_ - other.x_, y_ - other.y_, z_ - other.z_);
    }

    Vector3D operator*(double scalar) const {
        return Vector3D(x_ * scalar, y_ * scalar, z_ * scalar);
    }

    Vector3D operator/(double scalar) const {
        return Vector3D(x_ / scalar, y_ / scalar, z_ / scalar);
    }

    Vector3D operator-() const {
        return Vector3D(-x_, -y_, -z_);
    }

    // Compound assignment operators
    Vector3D& operator+=(const Vector3D& other) {
        x_ += other.x_;
        y_ += other.y_;
        z_ += other.z_;
        return *this;
    }

    Vector3D& operator-=(const Vector3D& other) {
        x_ -= other.x_;
        y_ -= other.y_;
        z_ -= other.z_;
        return *this;
    }

    Vector3D& operator*=(double scalar) {
        x_ *= scalar;
        y_ *= scalar;
        z_ *= scalar;
        return *this;
    }

    Vector3D& operator/=(double scalar) {
        x_ /= scalar;
        y_ /= scalar;
        z_ /= scalar;
        return *this;
    }

    // Comparison operators
    bool operator==(const Vector3D& other) const {
        constexpr double epsilon = 1e-9;
        return std::abs(x_ - other.x_) < epsilon && std::abs(y_ - other.y_) < epsilon &&
               std::abs(z_ - other.z_) < epsilon;
    }

    bool operator!=(const Vector3D& other) const {
        return !(*this == other);
    }

    // Vector operations
    /**
     * @brief Calculate dot product with another vector
     */
    double dot(const Vector3D& other) const {
        return x_ * other.x_ + y_ * other.y_ + z_ * other.z_;
    }

    /**
     * @brief Calculate cross product with another vector
     */
    Vector3D cross(const Vector3D& other) const {
        return Vector3D(y_ * other.z_ - z_ * other.y_, z_ * other.x_ - x_ * other.z_,
                        x_ * other.y_ - y_ * other.x_);
    }

    /**
     * @brief Calculate length (magnitude) of vector
     */
    double length() const {
        return std::sqrt(x_ * x_ + y_ * y_ + z_ * z_);
    }

    /**
     * @brief Calculate squared length (faster, avoids sqrt)
     */
    double length_squared() const {
        return x_ * x_ + y_ * y_ + z_ * z_;
    }

    /**
     * @brief Return normalized vector (unit vector)
     * @note Returns zero vector if length is zero
     */
    Vector3D normalized() const {
        double len = length();
        if (len < 1e-9) {
            return Vector3D(0, 0, 0);
        }
        return *this / len;
    }

    /**
     * @brief Normalize this vector in place
     * @return True if normalization succeeded, false if vector is zero
     */
    bool normalize() {
        double len = length();
        if (len < 1e-9) {
            return false;
        }
        *this /= len;
        return true;
    }

    /**
     * @brief Calculate distance to another vector
     */
    double distance_to(const Vector3D& other) const {
        return (*this - other).length();
    }

    /**
     * @brief Calculate squared distance to another vector
     */
    double distance_squared_to(const Vector3D& other) const {
        return (*this - other).length_squared();
    }

    /**
     * @brief Convert to array [x, y, z]
     */
    std::array<double, 3> to_array() const {
        return {x_, y_, z_};
    }

    /**
     * @brief Convert to JSON array [x, y, z]
     */
    nlohmann::json to_json() const {
        nlohmann::json result = nlohmann::json::array();
        result.push_back(x_);
        result.push_back(y_);
        result.push_back(z_);
        return result;
    }

    /**
     * @brief Create from JSON array [x, y, z]
     */
    static Vector3D from_json(const nlohmann::json& json) {
        if (!json.is_array() || json.size() != 3) {
            throw std::invalid_argument("JSON must be array of 3 numbers");
        }
        return Vector3D(json[0].get<double>(), json[1].get<double>(), json[2].get<double>());
    }

private:
    double x_;
    double y_;
    double z_;
};

// Scalar multiplication from left
inline Vector3D operator*(double scalar, const Vector3D& vec) {
    return vec * scalar;
}

} // namespace geometry
} // namespace x3dna

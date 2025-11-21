/**
 * @file least_squares_fitter.hpp
 * @brief Least squares fitting for aligning point sets (quaternion-based algorithm)
 *
 * Implements the same algorithm as the original ls_fitting() function:
 * 1. Compute covariance matrix
 * 2. Build 4x4 quaternion matrix
 * 3. Find largest eigenvalue/eigenvector (quaternion)
 * 4. Extract rotation matrix from quaternion
 * 5. Compute translation
 * 6. Calculate RMS
 */

#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <nlohmann/json.hpp>
#include "vector3d.hpp"
#include "matrix3d.hpp"

namespace x3dna {
namespace geometry {

/**
 * @class LeastSquaresFitter
 * @brief Fits one set of points to another using least squares
 */
class LeastSquaresFitter {
public:
    /**
     * @struct FitResult
     * @brief Result of least squares fitting
     */
    struct FitResult {
        Matrix3D rotation;    // Rotation matrix
        Vector3D translation; // Translation vector
        double rms;           // Root mean square deviation

        /**
         * @brief Convert to legacy JSON format (matches ls_fitting record)
         */
        nlohmann::json to_json_legacy() const {
            nlohmann::json result;
            result["rotation_matrix"] = rotation.to_json_legacy();
            result["translation"] = translation.to_json();
            result["rms"] = rms;
            return result;
        }
    };

    /**
     * @brief Fit points1 to points2 using least squares
     * @param points1 Source points (standard/template)
     * @param points2 Target points (experimental)
     * @return FitResult with rotation, translation, and RMS
     * @throws std::invalid_argument if fewer than 3 points or sizes don't match
     */
    FitResult fit(const std::vector<Vector3D>& points1, const std::vector<Vector3D>& points2) {
        if (points1.size() != points2.size()) {
            throw std::invalid_argument("Point sets must have same size");
        }
        if (points1.size() < 3) {
            throw std::invalid_argument("Need at least 3 points for fitting");
        }

        size_t n = points1.size();

        // Compute covariance matrix
        Matrix3D cov = compute_covariance_matrix(points1, points2, n);

        // Build 4x4 quaternion matrix N
        Matrix4D N = build_quaternion_matrix(cov);

        // Find largest eigenvalue and corresponding eigenvector (quaternion)
        Vector4D quaternion = find_largest_eigenvector(N);

        // Extract rotation matrix from quaternion
        Matrix3D rotation = quaternion_to_rotation_matrix(quaternion);

        // Compute centroids
        Vector3D centroid1 = compute_centroid(points1);
        Vector3D centroid2 = compute_centroid(points2);

        // Compute translation: t = centroid2 - R * centroid1
        Vector3D translation = centroid2 - (rotation * centroid1);

        // Compute RMS
        double rms = compute_rms(points1, points2, rotation, translation, n);

        return {rotation, translation, rms};
    }

private:
    // Helper types for 4x4 matrix and 4D vector
    using Matrix4D = std::array<std::array<double, 4>, 4>;
    using Vector4D = std::array<double, 4>;

    /**
     * @brief Compute covariance matrix between two point sets
     */
    Matrix3D compute_covariance_matrix(const std::vector<Vector3D>& points1,
                                       const std::vector<Vector3D>& points2, size_t n) {
        // Compute centroids
        Vector3D centroid1 = compute_centroid(points1);
        Vector3D centroid2 = compute_centroid(points2);

        // Compute covariance: U = (1/(n-1)) * sum((p1 - c1) * (p2 - c2)^T)
        Matrix3D cov;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                double sum = 0.0;
                for (size_t k = 0; k < n; ++k) {
                    Vector3D d1 = points1[k] - centroid1;
                    Vector3D d2 = points2[k] - centroid2;
                    double d1_component;
                    if (i == 0) {
                        d1_component = d1.x();
                    } else if (i == 1) {
                        d1_component = d1.y();
                    } else {
                        d1_component = d1.z();
                    }
                    double d2_component;
                    if (j == 0) {
                        d2_component = d2.x();
                    } else if (j == 1) {
                        d2_component = d2.y();
                    } else {
                        d2_component = d2.z();
                    }
                    sum += d1_component * d2_component;
                }
                cov.set(i, j, sum / (n - 1));
            }
        }
        return cov;
    }

    /**
     * @brief Build 4x4 quaternion matrix from 3x3 covariance matrix
     */
    Matrix4D build_quaternion_matrix(const Matrix3D& U) {
        Matrix4D N;

        // Diagonal elements
        double u11 = U.at(0, 0), u22 = U.at(1, 1), u33 = U.at(2, 2);
        N[0][0] = u11 + u22 + u33;
        N[1][1] = u11 - u22 - u33;
        N[2][2] = -u11 + u22 - u33;
        N[3][3] = -u11 - u22 + u33;

        // Off-diagonal elements
        double u12 = U.at(0, 1), u21 = U.at(1, 0);
        double u13 = U.at(0, 2), u31 = U.at(2, 0);
        double u23 = U.at(1, 2), u32 = U.at(2, 1);

        N[0][1] = u23 - u32;
        N[1][0] = N[0][1];
        N[0][2] = u31 - u13;
        N[2][0] = N[0][2];
        N[0][3] = u12 - u21;
        N[3][0] = N[0][3];
        N[1][2] = u12 + u21;
        N[2][1] = N[1][2];
        N[1][3] = u31 + u13;
        N[3][1] = N[1][3];
        N[2][3] = u23 + u32;
        N[3][2] = N[2][3];

        return N;
    }

    /**
     * @brief Find largest eigenvalue and corresponding eigenvector using Jacobi method
     * @note Implements simplified Jacobi for 4x4 symmetric matrix (matching original algorithm)
     */
    Vector4D find_largest_eigenvector(const Matrix4D& N) {
        constexpr double XEPS = 1.0e-7;
        constexpr size_t max_iterations = 100;

        // Working copies
        Matrix4D a = N;
        Matrix4D v;

        // Initialize V as identity
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                if (i == j) {
                    v[i][j] = 1.0;
                } else {
                    v[i][j] = 0.0;
                }
            }
        }

        // Diagonal elements (eigenvalues)
        std::array<double, 4> d;
        std::array<double, 4> b;
        std::array<double, 4> z;

        for (size_t i = 0; i < 4; ++i) {
            b[i] = d[i] = a[i][i];
            z[i] = 0.0;
        }

        // Jacobi iterations
        for (size_t iter = 0; iter < max_iterations; ++iter) {
            // Sum of off-diagonal elements
            double sm = 0.0;
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = i + 1; j < 4; ++j) {
                    sm += std::abs(a[i][j]);
                }
            }

            if (sm < XEPS) {
                // Converged - sort eigenvalues and return largest eigenvector
                sort_eigenvalues(d, v);
                Vector4D result;
                for (size_t i = 0; i < 4; ++i) {
                    result[i] = v[i][3]; // Last column (largest eigenvalue)
                }
                return result;
            }

            double tresh;
            if (iter < 4) {
                tresh = 0.2 * sm / 16.0;
            } else {
                tresh = 0.0;
            }

            for (size_t ip = 0; ip < 3; ++ip) {
                for (size_t iq = ip + 1; iq < 4; ++iq) {
                    double g = 100.0 * std::abs(a[ip][iq]);

                    // Skip if element is negligible
                    if (iter > 4 && (std::abs(d[ip]) + g) == std::abs(d[ip]) &&
                        (std::abs(d[iq]) + g) == std::abs(d[iq])) {
                        a[ip][iq] = 0.0;
                        continue;
                    }

                    // Skip if below threshold
                    if (std::abs(a[ip][iq]) <= tresh) {
                        continue;
                    }

                    // Perform Jacobi rotation
                    perform_jacobi_rotation(a, v, d, z, ip, iq, g);
                }
            }

            for (size_t i = 0; i < 4; ++i) {
                b[i] += z[i];
                d[i] = b[i];
                z[i] = 0.0;
            }
        }

        // If we get here, sort and return
        sort_eigenvalues(d, v);
        Vector4D result;
        for (size_t i = 0; i < 4; ++i) {
            result[i] = v[i][3]; // Last column (largest eigenvalue)
        }
        return result;
    }

    /**
     * @brief Rotate matrix elements (helper for Jacobi)
     */
    void rotate_matrix(Matrix4D& a, size_t i, size_t j, size_t k, size_t l, double s, double tau) {
        double g = a[i][j];
        double h = a[k][l];
        a[i][j] = g - s * (h + g * tau);
        a[k][l] = h + s * (g - h * tau);
    }

    /**
     * @brief Compute rotation angle for Jacobi method
     * Reduces nesting by extracting angle calculation
     */
    double compute_rotation_angle(const Matrix4D& a, size_t ip, size_t iq, double h, double g) {
        if ((std::abs(h) + g) == std::abs(h)) {
            return a[ip][iq] / h;
        }

        double theta = 0.5 * h / a[ip][iq];
        double t = 1.0 / (std::abs(theta) + std::sqrt(1.0 + theta * theta));
        if (theta < 0.0) {
            t = -t;
        }
        return t;
    }

    /**
     * @brief Perform a single Jacobi rotation step
     * Reduces nesting by extracting complex rotation logic
     */
    void perform_jacobi_rotation(Matrix4D& a, Matrix4D& v, std::array<double, 4>& d,
                                 std::array<double, 4>& z, size_t ip, size_t iq, double g) {
        double h = d[iq] - d[ip];
        double t = compute_rotation_angle(a, ip, iq, h, g);

        double c = 1.0 / std::sqrt(1.0 + t * t);
        double s = t * c;
        double tau = s / (1.0 + c);
        h = t * a[ip][iq];

        z[ip] -= h;
        z[iq] += h;
        d[ip] -= h;
        d[iq] += h;
        a[ip][iq] = 0.0;

        // Rotate matrix a
        for (size_t j = 0; j < ip; ++j) {
            rotate_matrix(a, j, ip, j, iq, s, tau);
        }
        for (size_t j = ip + 1; j < iq; ++j) {
            rotate_matrix(a, ip, j, j, iq, s, tau);
        }
        for (size_t j = iq + 1; j < 4; ++j) {
            rotate_matrix(a, ip, j, iq, j, s, tau);
        }

        // Rotate eigenvector matrix v
        for (size_t j = 0; j < 4; ++j) {
            rotate_matrix(v, j, ip, j, iq, s, tau);
        }
    }

    /**
     * @brief Sort eigenvalues in ascending order and reorder eigenvectors
     */
    void sort_eigenvalues(std::array<double, 4>& d, Matrix4D& v) {
        // Simple selection sort
        for (size_t i = 0; i < 3; ++i) {
            size_t k = i;
            double p = d[i];
            for (size_t j = i + 1; j < 4; ++j) {
                if (d[j] < p) {
                    p = d[j];
                    k = j;
                }
            }
            if (k != i) {
                std::swap(d[i], d[k]);
                // Swap columns in v
                for (size_t j = 0; j < 4; ++j) {
                    std::swap(v[j][i], v[j][k]);
                }
            }
        }
    }

    /**
     * @brief Convert quaternion to rotation matrix
     * @param q Quaternion [q0, q1, q2, q3] where q0 is scalar part
     */
    Matrix3D quaternion_to_rotation_matrix(const Vector4D& q) {
        // Build quaternion outer product matrix N (as in original code)
        Matrix4D N;
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                N[i][j] = q[i] * q[j];
            }
        }

        // Extract rotation matrix (matching original code)
        Matrix3D R;
        R.set(0, 0, N[0][0] + N[1][1] - N[2][2] - N[3][3]);
        R.set(0, 1, 2.0 * (N[1][2] - N[0][3]));
        R.set(0, 2, 2.0 * (N[1][3] + N[0][2]));
        R.set(1, 0, 2.0 * (N[2][1] + N[0][3]));
        R.set(1, 1, N[0][0] - N[1][1] + N[2][2] - N[3][3]);
        R.set(1, 2, 2.0 * (N[2][3] - N[0][1]));
        R.set(2, 0, 2.0 * (N[3][1] - N[0][2]));
        R.set(2, 1, 2.0 * (N[3][2] + N[0][1]));
        R.set(2, 2, N[0][0] - N[1][1] - N[2][2] + N[3][3]);

        return R;
    }

    /**
     * @brief Compute centroid of point set
     */
    Vector3D compute_centroid(const std::vector<Vector3D>& points) {
        Vector3D sum(0, 0, 0);
        for (const auto& p : points) {
            sum += p;
        }
        return sum / static_cast<double>(points.size());
    }

    /**
     * @brief Compute RMS between transformed points1 and points2
     */
    double compute_rms(const std::vector<Vector3D>& points1, const std::vector<Vector3D>& points2,
                       const Matrix3D& rotation, const Vector3D& translation, size_t n) {
        double sum_sq_diff = 0.0;
        for (size_t i = 0; i < n; ++i) {
            Vector3D transformed = rotation * points1[i] + translation;
            Vector3D diff = transformed - points2[i];
            sum_sq_diff += diff.length_squared();
        }
        return std::sqrt(sum_sq_diff / n);
    }
};

} // namespace geometry
} // namespace x3dna

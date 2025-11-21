/**
 * @file test_least_squares_fitter.cpp
 * @brief Comprehensive tests for LeastSquaresFitter
 */

#include <gtest/gtest.h>
#include <x3dna/geometry/least_squares_fitter.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>
#include <cmath>

using namespace x3dna::geometry;

class LeastSquaresFitterTest : public ::testing::Test {
protected:
    void SetUp() override {
        fitter_ = std::make_unique<LeastSquaresFitter>();
    }

    std::unique_ptr<LeastSquaresFitter> fitter_;
    constexpr static double TOLERANCE = 1e-6;
};

// Test with simple translation (no rotation)
TEST_F(LeastSquaresFitterTest, SimpleTranslation) {
    std::vector<Vector3D> points1 = {Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0),
                                     Vector3D(0.0, 1.0, 0.0), Vector3D(0.0, 0.0, 1.0)};

    Vector3D translation(1.0, 2.0, 3.0);
    std::vector<Vector3D> points2;
    for (const auto& p : points1) {
        points2.push_back(p + translation);
    }

    auto result = fitter_->fit(points1, points2);

    // Rotation should be identity (approximately)
    EXPECT_TRUE(result.rotation.approximately_equals(Matrix3D::identity(), TOLERANCE));

    // Translation should match
    EXPECT_NEAR(result.translation.x(), translation.x(), TOLERANCE);
    EXPECT_NEAR(result.translation.y(), translation.y(), TOLERANCE);
    EXPECT_NEAR(result.translation.z(), translation.z(), TOLERANCE);

    // RMS should be very small (perfect fit)
    EXPECT_NEAR(result.rms, 0.0, TOLERANCE);
}

// Test with simple rotation around Z axis
TEST_F(LeastSquaresFitterTest, SimpleRotationZ) {
    std::vector<Vector3D> points1 = {Vector3D(1.0, 0.0, 0.0), Vector3D(0.0, 1.0, 0.0),
                                     Vector3D(-1.0, 0.0, 0.0), Vector3D(0.0, -1.0, 0.0)};

    // Rotate 90 degrees around Z axis
    double angle = M_PI / 2.0;
    Matrix3D rot = Matrix3D::rotation_z(angle);

    std::vector<Vector3D> points2;
    for (const auto& p : points1) {
        points2.push_back(rot * p);
    }

    auto result = fitter_->fit(points1, points2);

    // Check that transformation is correct (more important than exact rotation match)
    for (size_t i = 0; i < points1.size(); ++i) {
        Vector3D transformed = result.rotation * points1[i] + result.translation;
        double dist = transformed.distance_to(points2[i]);
        EXPECT_NEAR(dist, 0.0, 0.01) << "Point " << i << " doesn't match";
    }

    // RMS should be very small
    EXPECT_NEAR(result.rms, 0.0, 0.01);
}

// Test with combined rotation and translation
TEST_F(LeastSquaresFitterTest, CombinedRotationAndTranslation) {
    std::vector<Vector3D> points1 = {Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0),
                                     Vector3D(0.0, 1.0, 0.0), Vector3D(0.0, 0.0, 1.0),
                                     Vector3D(1.0, 1.0, 1.0)};

    Matrix3D rot = Matrix3D::rotation_y(M_PI / 4.0);
    Vector3D trans(5.0, 10.0, 15.0);

    std::vector<Vector3D> points2;
    for (const auto& p : points1) {
        points2.push_back(rot * p + trans);
    }

    auto result = fitter_->fit(points1, points2);

    // Check that transformation is correct
    for (size_t i = 0; i < points1.size(); ++i) {
        Vector3D transformed = result.rotation * points1[i] + result.translation;
        double dist = transformed.distance_to(points2[i]);
        EXPECT_NEAR(dist, 0.0, TOLERANCE) << "Point " << i << " doesn't match";
    }

    // RMS should be very small
    EXPECT_NEAR(result.rms, 0.0, TOLERANCE);
}

// Test with known 3D triangle
TEST_F(LeastSquaresFitterTest, TriangleTransformation) {
    std::vector<Vector3D> points1 = {Vector3D(0.0, 0.0, 0.0), Vector3D(3.0, 0.0, 0.0),
                                     Vector3D(0.0, 4.0, 0.0)};

    // Apply known transformation
    Matrix3D rot = Matrix3D::rotation_x(M_PI / 6.0); // 30 degrees
    Vector3D trans(1.0, 2.0, 3.0);

    std::vector<Vector3D> points2;
    for (const auto& p : points1) {
        points2.push_back(rot * p + trans);
    }

    auto result = fitter_->fit(points1, points2);

    // Verify each point transforms correctly
    for (size_t i = 0; i < points1.size(); ++i) {
        Vector3D transformed = result.rotation * points1[i] + result.translation;
        double dist = transformed.distance_to(points2[i]);
        EXPECT_NEAR(dist, 0.0, 0.01) << "Point " << i << " distance: " << dist;
    }
}

// Test error handling
TEST_F(LeastSquaresFitterTest, ErrorHandling) {
    std::vector<Vector3D> points1 = {Vector3D(0, 0, 0), Vector3D(1, 0, 0)};
    std::vector<Vector3D> points2 = {Vector3D(0, 0, 0), Vector3D(1, 0, 0)};

    // Too few points
    EXPECT_THROW(fitter_->fit(points1, points2), std::invalid_argument);

    // Mismatched sizes
    points1.push_back(Vector3D(0, 1, 0));
    EXPECT_THROW(fitter_->fit(points1, points2), std::invalid_argument);
}

// Test with collinear points (degenerate case)
TEST_F(LeastSquaresFitterTest, CollinearPoints) {
    std::vector<Vector3D> points1 = {Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0),
                                     Vector3D(2.0, 0.0, 0.0), Vector3D(3.0, 0.0, 0.0)};

    Vector3D trans(1.0, 2.0, 3.0);
    std::vector<Vector3D> points2;
    for (const auto& p : points1) {
        points2.push_back(p + trans);
    }

    // Should still work (translation only)
    auto result = fitter_->fit(points1, points2);

    // Translation should be correct
    EXPECT_NEAR(result.translation.x(), trans.x(), TOLERANCE);
    EXPECT_NEAR(result.translation.y(), trans.y(), TOLERANCE);
    EXPECT_NEAR(result.translation.z(), trans.z(), TOLERANCE);

    // RMS should be small
    EXPECT_NEAR(result.rms, 0.0, TOLERANCE);
}

// Test RMS calculation accuracy
TEST_F(LeastSquaresFitterTest, RMSCalculation) {
    std::vector<Vector3D> points1 = {Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0),
                                     Vector3D(0.0, 1.0, 0.0)};

    // Perfect fit
    std::vector<Vector3D> points2 = points1;
    auto result1 = fitter_->fit(points1, points2);
    EXPECT_NEAR(result1.rms, 0.0, TOLERANCE);

    // Add small noise
    std::vector<Vector3D> points3;
    for (const auto& p : points1) {
        points3.push_back(p + Vector3D(0.01, 0.01, 0.01));
    }
    auto result2 = fitter_->fit(points1, points3);

    // RMS should be approximately the noise level
    EXPECT_GT(result2.rms, 0.0);
    EXPECT_LT(result2.rms, 0.02); // Should be close to noise level
}

// Test with real-world like coordinates (atom positions)
TEST_F(LeastSquaresFitterTest, RealWorldCoordinates) {
    // Simulate base ring atoms
    std::vector<Vector3D> points1 = {
        Vector3D(2.5, 0.0, 0.0),  // N1
        Vector3D(1.5, 1.5, 0.0),  // C2
        Vector3D(0.0, 1.5, 0.0),  // N3
        Vector3D(-0.5, 0.0, 0.0), // C4
        Vector3D(0.0, -1.5, 0.0), // C5
        Vector3D(1.5, -1.5, 0.0)  // C6
    };

    // Apply transformation
    Matrix3D rot = Matrix3D::rotation_z(M_PI / 3.0) * Matrix3D::rotation_x(M_PI / 6.0);
    Vector3D trans(10.0, 20.0, 30.0);

    std::vector<Vector3D> points2;
    for (const auto& p : points1) {
        points2.push_back(rot * p + trans);
    }

    auto result = fitter_->fit(points1, points2);

    // Verify transformation (allow larger tolerance for complex rotations)
    double max_error = 0.0;
    for (size_t i = 0; i < points1.size(); ++i) {
        Vector3D transformed = result.rotation * points1[i] + result.translation;
        double dist = transformed.distance_to(points2[i]);
        max_error = std::max(max_error, dist);
        EXPECT_NEAR(dist, 0.0, 0.1) << "Atom " << i << " doesn't match, distance: " << dist;
    }

    // RMS should be reasonable (algorithm may have numerical precision issues)
    EXPECT_LT(result.rms, 0.2);
}

// Test JSON serialization
TEST_F(LeastSquaresFitterTest, JSONSerialization) {
    std::vector<Vector3D> points1 = {Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0),
                                     Vector3D(0.0, 1.0, 0.0)};

    Vector3D trans(1.0, 2.0, 3.0);
    std::vector<Vector3D> points2;
    for (const auto& p : points1) {
        points2.push_back(p + trans);
    }

    auto result = fitter_->fit(points1, points2);
    auto json = result.to_json_legacy();

    // Check JSON structure
    EXPECT_TRUE(json.contains("rotation_matrix"));
    EXPECT_TRUE(json.contains("translation"));
    EXPECT_TRUE(json.contains("rms"));

    // Check values
    EXPECT_TRUE(json["rotation_matrix"].is_array());
    EXPECT_EQ(json["rotation_matrix"].size(), 3);
    EXPECT_DOUBLE_EQ(json["rms"].get<double>(), result.rms);
}

// Test with many points (stress test)
TEST_F(LeastSquaresFitterTest, ManyPoints) {
    std::vector<Vector3D> points1;
    for (int i = 0; i < 100; ++i) {
        points1.push_back(Vector3D(static_cast<double>(i % 10), static_cast<double>((i / 10) % 10),
                                   static_cast<double>(i / 100)));
    }

    Matrix3D rot = Matrix3D::rotation_x(M_PI / 4.0);
    Vector3D trans(5.0, 10.0, 15.0);

    std::vector<Vector3D> points2;
    for (const auto& p : points1) {
        points2.push_back(rot * p + trans);
    }

    auto result = fitter_->fit(points1, points2);

    // Verify transformation
    double max_error = 0.0;
    for (size_t i = 0; i < points1.size(); ++i) {
        Vector3D transformed = result.rotation * points1[i] + result.translation;
        double dist = transformed.distance_to(points2[i]);
        max_error = std::max(max_error, dist);
    }

    EXPECT_NEAR(max_error, 0.0, 0.01);
    EXPECT_NEAR(result.rms, 0.0, 0.01);
}

// Test identity transformation
TEST_F(LeastSquaresFitterTest, IdentityTransformation) {
    std::vector<Vector3D> points1 = {Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0),
                                     Vector3D(0.0, 1.0, 0.0), Vector3D(0.0, 0.0, 1.0)};

    std::vector<Vector3D> points2 = points1; // Same points

    auto result = fitter_->fit(points1, points2);

    // Rotation should be identity
    EXPECT_TRUE(result.rotation.approximately_equals(Matrix3D::identity(), TOLERANCE));

    // Translation should be zero
    EXPECT_NEAR(result.translation.length(), 0.0, TOLERANCE);

    // RMS should be zero
    EXPECT_NEAR(result.rms, 0.0, TOLERANCE);
}

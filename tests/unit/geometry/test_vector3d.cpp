/**
 * @file test_vector3d.cpp
 * @brief Unit tests for Vector3D class
 */

#include <gtest/gtest.h>
#include <x3dna/geometry/vector3d.hpp>
#include <cmath>

using namespace x3dna::geometry;

class Vector3DTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

// Constructor tests
TEST_F(Vector3DTest, DefaultConstructor) {
    Vector3D v;
    EXPECT_DOUBLE_EQ(v.x(), 0.0);
    EXPECT_DOUBLE_EQ(v.y(), 0.0);
    EXPECT_DOUBLE_EQ(v.z(), 0.0);
}

TEST_F(Vector3DTest, ConstructorWithValues) {
    Vector3D v(1.0, 2.0, 3.0);
    EXPECT_DOUBLE_EQ(v.x(), 1.0);
    EXPECT_DOUBLE_EQ(v.y(), 2.0);
    EXPECT_DOUBLE_EQ(v.z(), 3.0);
}

TEST_F(Vector3DTest, ConstructorFromArray) {
    std::array<double, 3> arr = {4.0, 5.0, 6.0};
    Vector3D v(arr);
    EXPECT_DOUBLE_EQ(v.x(), 4.0);
    EXPECT_DOUBLE_EQ(v.y(), 5.0);
    EXPECT_DOUBLE_EQ(v.z(), 6.0);
}

// Getter/Setter tests
TEST_F(Vector3DTest, GettersSetters) {
    Vector3D v;
    v.set_x(1.0);
    v.set_y(2.0);
    v.set_z(3.0);
    EXPECT_DOUBLE_EQ(v.x(), 1.0);
    EXPECT_DOUBLE_EQ(v.y(), 2.0);
    EXPECT_DOUBLE_EQ(v.z(), 3.0);

    v.set(4.0, 5.0, 6.0);
    EXPECT_DOUBLE_EQ(v.x(), 4.0);
    EXPECT_DOUBLE_EQ(v.y(), 5.0);
    EXPECT_DOUBLE_EQ(v.z(), 6.0);
}

// Arithmetic operations
TEST_F(Vector3DTest, Addition) {
    Vector3D v1(1.0, 2.0, 3.0);
    Vector3D v2(4.0, 5.0, 6.0);
    Vector3D result = v1 + v2;
    EXPECT_DOUBLE_EQ(result.x(), 5.0);
    EXPECT_DOUBLE_EQ(result.y(), 7.0);
    EXPECT_DOUBLE_EQ(result.z(), 9.0);
}

TEST_F(Vector3DTest, Subtraction) {
    Vector3D v1(5.0, 7.0, 9.0);
    Vector3D v2(1.0, 2.0, 3.0);
    Vector3D result = v1 - v2;
    EXPECT_DOUBLE_EQ(result.x(), 4.0);
    EXPECT_DOUBLE_EQ(result.y(), 5.0);
    EXPECT_DOUBLE_EQ(result.z(), 6.0);
}

TEST_F(Vector3DTest, ScalarMultiplication) {
    Vector3D v(1.0, 2.0, 3.0);
    Vector3D result = v * 2.0;
    EXPECT_DOUBLE_EQ(result.x(), 2.0);
    EXPECT_DOUBLE_EQ(result.y(), 4.0);
    EXPECT_DOUBLE_EQ(result.z(), 6.0);

    // Test left multiplication
    Vector3D result2 = 2.0 * v;
    EXPECT_DOUBLE_EQ(result2.x(), 2.0);
    EXPECT_DOUBLE_EQ(result2.y(), 4.0);
    EXPECT_DOUBLE_EQ(result2.z(), 6.0);
}

TEST_F(Vector3DTest, ScalarDivision) {
    Vector3D v(4.0, 6.0, 8.0);
    Vector3D result = v / 2.0;
    EXPECT_DOUBLE_EQ(result.x(), 2.0);
    EXPECT_DOUBLE_EQ(result.y(), 3.0);
    EXPECT_DOUBLE_EQ(result.z(), 4.0);
}

TEST_F(Vector3DTest, Negation) {
    Vector3D v(1.0, 2.0, 3.0);
    Vector3D result = -v;
    EXPECT_DOUBLE_EQ(result.x(), -1.0);
    EXPECT_DOUBLE_EQ(result.y(), -2.0);
    EXPECT_DOUBLE_EQ(result.z(), -3.0);
}

// Dot product
TEST_F(Vector3DTest, DotProduct) {
    Vector3D v1(1.0, 2.0, 3.0);
    Vector3D v2(4.0, 5.0, 6.0);
    double result = v1.dot(v2);
    EXPECT_DOUBLE_EQ(result, 32.0); // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
}

// Cross product
TEST_F(Vector3DTest, CrossProduct) {
    Vector3D v1(1.0, 0.0, 0.0);
    Vector3D v2(0.0, 1.0, 0.0);
    Vector3D result = v1.cross(v2);
    EXPECT_NEAR(result.x(), 0.0, 1e-9);
    EXPECT_NEAR(result.y(), 0.0, 1e-9);
    EXPECT_NEAR(result.z(), 1.0, 1e-9);

    // Test orthogonality
    EXPECT_NEAR(result.dot(v1), 0.0, 1e-9);
    EXPECT_NEAR(result.dot(v2), 0.0, 1e-9);
}

// Length
TEST_F(Vector3DTest, Length) {
    Vector3D v(3.0, 4.0, 0.0);
    EXPECT_DOUBLE_EQ(v.length(), 5.0); // 3-4-5 triangle

    Vector3D v2(1.0, 1.0, 1.0);
    EXPECT_NEAR(v2.length(), std::sqrt(3.0), 1e-9);
}

TEST_F(Vector3DTest, LengthSquared) {
    Vector3D v(3.0, 4.0, 0.0);
    EXPECT_DOUBLE_EQ(v.length_squared(), 25.0);
}

// Normalization
TEST_F(Vector3DTest, Normalize) {
    Vector3D v(3.0, 4.0, 0.0);
    Vector3D normalized = v.normalized();
    EXPECT_NEAR(normalized.length(), 1.0, 1e-9);
    EXPECT_NEAR(normalized.x(), 0.6, 1e-9);
    EXPECT_NEAR(normalized.y(), 0.8, 1e-9);
    EXPECT_NEAR(normalized.z(), 0.0, 1e-9);

    // Original vector should be unchanged
    EXPECT_DOUBLE_EQ(v.length(), 5.0);
}

TEST_F(Vector3DTest, NormalizeInPlace) {
    Vector3D v(3.0, 4.0, 0.0);
    bool success = v.normalize();
    EXPECT_TRUE(success);
    EXPECT_NEAR(v.length(), 1.0, 1e-9);

    // Test zero vector
    Vector3D zero;
    bool fail = zero.normalize();
    EXPECT_FALSE(fail);
}

// Distance
TEST_F(Vector3DTest, Distance) {
    Vector3D v1(0.0, 0.0, 0.0);
    Vector3D v2(3.0, 4.0, 0.0);
    EXPECT_DOUBLE_EQ(v1.distance_to(v2), 5.0);
}

TEST_F(Vector3DTest, DistanceSquared) {
    Vector3D v1(0.0, 0.0, 0.0);
    Vector3D v2(3.0, 4.0, 0.0);
    EXPECT_DOUBLE_EQ(v1.distance_squared_to(v2), 25.0);
}

// Comparison
TEST_F(Vector3DTest, Equality) {
    Vector3D v1(1.0, 2.0, 3.0);
    Vector3D v2(1.0, 2.0, 3.0);
    Vector3D v3(1.0, 2.0, 3.0001);

    EXPECT_TRUE(v1 == v2);
    EXPECT_FALSE(v1 == v3);
    EXPECT_TRUE(v1 != v3);
}

// JSON serialization
TEST_F(Vector3DTest, JSONSerialization) {
    Vector3D v(1.0, 2.0, 3.0);
    auto json = v.to_json();

    EXPECT_TRUE(json.is_array());
    EXPECT_EQ(json.size(), 3);
    EXPECT_DOUBLE_EQ(json[0].get<double>(), 1.0);
    EXPECT_DOUBLE_EQ(json[1].get<double>(), 2.0);
    EXPECT_DOUBLE_EQ(json[2].get<double>(), 3.0);

    Vector3D restored = Vector3D::from_json(json);
    EXPECT_TRUE(v == restored);
}

// Edge cases
TEST_F(Vector3DTest, ZeroVector) {
    Vector3D zero;
    EXPECT_DOUBLE_EQ(zero.length(), 0.0);
    EXPECT_DOUBLE_EQ(zero.dot(zero), 0.0);

    Vector3D normalized = zero.normalized();
    EXPECT_DOUBLE_EQ(normalized.length(), 0.0);
}

TEST_F(Vector3DTest, CompoundAssignment) {
    Vector3D v1(1.0, 2.0, 3.0);
    Vector3D v2(4.0, 5.0, 6.0);

    v1 += v2;
    EXPECT_DOUBLE_EQ(v1.x(), 5.0);
    EXPECT_DOUBLE_EQ(v1.y(), 7.0);
    EXPECT_DOUBLE_EQ(v1.z(), 9.0);

    v1 -= v2;
    EXPECT_DOUBLE_EQ(v1.x(), 1.0);
    EXPECT_DOUBLE_EQ(v1.y(), 2.0);
    EXPECT_DOUBLE_EQ(v1.z(), 3.0);

    v1 *= 2.0;
    EXPECT_DOUBLE_EQ(v1.x(), 2.0);
    EXPECT_DOUBLE_EQ(v1.y(), 4.0);
    EXPECT_DOUBLE_EQ(v1.z(), 6.0);

    v1 /= 2.0;
    EXPECT_DOUBLE_EQ(v1.x(), 1.0);
    EXPECT_DOUBLE_EQ(v1.y(), 2.0);
    EXPECT_DOUBLE_EQ(v1.z(), 3.0);
}

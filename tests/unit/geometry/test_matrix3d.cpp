/**
 * @file test_matrix3d.cpp
 * @brief Unit tests for Matrix3D class
 */

#include <gtest/gtest.h>
#include <x3dna/geometry/matrix3d.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <cmath>

using namespace x3dna::geometry;

class Matrix3DTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

// Constructor tests
TEST_F(Matrix3DTest, DefaultConstructor) {
    Matrix3D m;
    EXPECT_DOUBLE_EQ(m.at(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(m.at(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(m.at(2, 2), 1.0);
    EXPECT_DOUBLE_EQ(m.at(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(m.at(1, 0), 0.0);
}

TEST_F(Matrix3DTest, ConstructorFromArray) {
    std::array<double, 9> arr = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    Matrix3D m(arr);
    EXPECT_DOUBLE_EQ(m.at(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(m.at(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(m.at(1, 0), 4.0);
    EXPECT_DOUBLE_EQ(m.at(2, 2), 9.0);
}

TEST_F(Matrix3DTest, ConstructorFromElements) {
    Matrix3D m(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    EXPECT_DOUBLE_EQ(m.at(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(m.at(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(m.at(1, 0), 4.0);
}

// Accessor tests
TEST_F(Matrix3DTest, AtAndSet) {
    Matrix3D m;
    m.set(1, 2, 5.0);
    EXPECT_DOUBLE_EQ(m.at(1, 2), 5.0);

    EXPECT_THROW(m.at(3, 0), std::out_of_range);
    EXPECT_THROW(m.set(0, 3, 1.0), std::out_of_range);
}

TEST_F(Matrix3DTest, RowAndColumn) {
    Matrix3D m(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

    Vector3D row0 = m.row(0);
    EXPECT_DOUBLE_EQ(row0.x(), 1.0);
    EXPECT_DOUBLE_EQ(row0.y(), 2.0);
    EXPECT_DOUBLE_EQ(row0.z(), 3.0);

    Vector3D col0 = m.column(0);
    EXPECT_DOUBLE_EQ(col0.x(), 1.0);
    EXPECT_DOUBLE_EQ(col0.y(), 4.0);
    EXPECT_DOUBLE_EQ(col0.z(), 7.0);
}

TEST_F(Matrix3DTest, SetRowAndColumn) {
    Matrix3D m;
    Vector3D vec(10.0, 20.0, 30.0);

    m.set_row(1, vec);
    Vector3D row1 = m.row(1);
    EXPECT_DOUBLE_EQ(row1.x(), 10.0);
    EXPECT_DOUBLE_EQ(row1.y(), 20.0);
    EXPECT_DOUBLE_EQ(row1.z(), 30.0);

    // Reset matrix and test column setting
    m = Matrix3D();
    m.set_column(2, vec);
    Vector3D col2 = m.column(2);
    EXPECT_DOUBLE_EQ(col2.x(), 10.0);
    EXPECT_DOUBLE_EQ(col2.y(), 20.0);
    EXPECT_DOUBLE_EQ(col2.z(), 30.0);
}

// Matrix-vector multiplication
TEST_F(Matrix3DTest, MatrixVectorMultiplication) {
    Matrix3D m(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Vector3D v(1.0, 2.0, 3.0);

    Vector3D result = m * v;
    // First row: 1*1 + 2*2 + 3*3 = 1 + 4 + 9 = 14
    EXPECT_DOUBLE_EQ(result.x(), 14.0);
    // Second row: 4*1 + 5*2 + 6*3 = 4 + 10 + 18 = 32
    EXPECT_DOUBLE_EQ(result.y(), 32.0);
    // Third row: 7*1 + 8*2 + 9*3 = 7 + 16 + 27 = 50
    EXPECT_DOUBLE_EQ(result.z(), 50.0);
}

TEST_F(Matrix3DTest, IdentityMatrixVectorMultiplication) {
    Matrix3D I = Matrix3D::identity();
    Vector3D v(1.0, 2.0, 3.0);
    Vector3D result = I * v;
    EXPECT_DOUBLE_EQ(result.x(), 1.0);
    EXPECT_DOUBLE_EQ(result.y(), 2.0);
    EXPECT_DOUBLE_EQ(result.z(), 3.0);
}

// Matrix-matrix multiplication
TEST_F(Matrix3DTest, MatrixMatrixMultiplication) {
    Matrix3D m1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Matrix3D m2(9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0);

    Matrix3D result = m1 * m2;
    // First element: (1*9 + 2*6 + 3*3) = 9 + 12 + 9 = 30
    EXPECT_DOUBLE_EQ(result.at(0, 0), 30.0);
    // Second element: (1*8 + 2*5 + 3*2) = 8 + 10 + 6 = 24
    EXPECT_DOUBLE_EQ(result.at(0, 1), 24.0);
}

TEST_F(Matrix3DTest, IdentityMatrixMultiplication) {
    Matrix3D I = Matrix3D::identity();
    Matrix3D m(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

    Matrix3D result1 = I * m;
    EXPECT_TRUE(result1.approximately_equals(m, 1e-9));

    Matrix3D result2 = m * I;
    EXPECT_TRUE(result2.approximately_equals(m, 1e-9));
}

// Transpose
TEST_F(Matrix3DTest, Transpose) {
    Matrix3D m(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Matrix3D transposed = m.transpose();

    EXPECT_DOUBLE_EQ(transposed.at(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(transposed.at(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(transposed.at(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(transposed.at(2, 0), 3.0);

    // Transpose twice should give original
    Matrix3D double_transposed = transposed.transpose();
    EXPECT_TRUE(double_transposed.approximately_equals(m, 1e-9));
}

// Determinant
TEST_F(Matrix3DTest, Determinant) {
    Matrix3D I = Matrix3D::identity();
    EXPECT_DOUBLE_EQ(I.determinant(), 1.0);

    Matrix3D m(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    // det = 1*(5*9-6*8) - 2*(4*9-6*7) + 3*(4*8-5*7)
    //     = 1*(45-48) - 2*(36-42) + 3*(32-35)
    //     = 1*(-3) - 2*(-6) + 3*(-3)
    //     = -3 + 12 - 9 = 0
    EXPECT_NEAR(m.determinant(), 0.0, 1e-9);
}

// Inverse
TEST_F(Matrix3DTest, Inverse) {
    Matrix3D I = Matrix3D::identity();
    Matrix3D I_inv = I.inverse();
    EXPECT_TRUE(I_inv.approximately_equals(I, 1e-9));

    // Test with a known invertible matrix
    Matrix3D m(1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0);
    Matrix3D m_inv = m.inverse();
    Matrix3D product = m * m_inv;
    EXPECT_TRUE(product.approximately_equals(Matrix3D::identity(), 1e-6));
}

TEST_F(Matrix3DTest, InverseSingularMatrix) {
    // Matrix with zero determinant
    Matrix3D singular(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    EXPECT_THROW(singular.inverse(), std::runtime_error);
}

// Rotation matrices
TEST_F(Matrix3DTest, RotationX) {
    double angle = M_PI / 2.0; // 90 degrees
    Matrix3D rot = Matrix3D::rotation_x(angle);

    Vector3D v(0.0, 1.0, 0.0);
    Vector3D result = rot * v;
    // Rotating (0,1,0) around X by 90° should give (0,0,1)
    EXPECT_NEAR(result.x(), 0.0, 1e-9);
    EXPECT_NEAR(result.y(), 0.0, 1e-9);
    EXPECT_NEAR(result.z(), 1.0, 1e-9);
}

TEST_F(Matrix3DTest, RotationY) {
    double angle = M_PI / 2.0; // 90 degrees
    Matrix3D rot = Matrix3D::rotation_y(angle);

    Vector3D v(1.0, 0.0, 0.0);
    Vector3D result = rot * v;
    // Rotating (1,0,0) around Y by 90° should give (0,0,-1)
    EXPECT_NEAR(result.x(), 0.0, 1e-9);
    EXPECT_NEAR(result.y(), 0.0, 1e-9);
    EXPECT_NEAR(result.z(), -1.0, 1e-9);
}

TEST_F(Matrix3DTest, RotationZ) {
    double angle = M_PI / 2.0; // 90 degrees
    Matrix3D rot = Matrix3D::rotation_z(angle);

    Vector3D v(1.0, 0.0, 0.0);
    Vector3D result = rot * v;
    // Rotating (1,0,0) around Z by 90° should give (0,1,0)
    EXPECT_NEAR(result.x(), 0.0, 1e-9);
    EXPECT_NEAR(result.y(), 1.0, 1e-9);
    EXPECT_NEAR(result.z(), 0.0, 1e-9);
}

// JSON serialization
TEST_F(Matrix3DTest, JSONSerialization) {
    Matrix3D m(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

    auto json = m.to_json();
    EXPECT_TRUE(json.is_array());
    EXPECT_EQ(json.size(), 9);
    EXPECT_DOUBLE_EQ(json[0].get<double>(), 1.0);

    Matrix3D restored = Matrix3D::from_json(json);
    EXPECT_TRUE(restored.approximately_equals(m, 1e-9));
}

TEST_F(Matrix3DTest, LegacyJSONSerialization) {
    Matrix3D m(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

    auto json = m.to_json_legacy();
    EXPECT_TRUE(json.is_array());
    EXPECT_EQ(json.size(), 3);
    EXPECT_TRUE(json[0].is_array());
    EXPECT_EQ(json[0].size(), 3);
    EXPECT_DOUBLE_EQ(json[0][0].get<double>(), 1.0);
    EXPECT_DOUBLE_EQ(json[1][1].get<double>(), 5.0);

    Matrix3D restored = Matrix3D::from_json_legacy(json);
    EXPECT_TRUE(restored.approximately_equals(m, 1e-9));
}

// Arithmetic operations
TEST_F(Matrix3DTest, Addition) {
    Matrix3D m1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Matrix3D m2(9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0);
    Matrix3D result = m1 + m2;
    EXPECT_DOUBLE_EQ(result.at(0, 0), 10.0);
    EXPECT_DOUBLE_EQ(result.at(1, 1), 10.0);
}

TEST_F(Matrix3DTest, Subtraction) {
    Matrix3D m1(5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0);
    Matrix3D m2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    Matrix3D result = m1 - m2;
    EXPECT_DOUBLE_EQ(result.at(0, 0), 4.0);
}

TEST_F(Matrix3DTest, ScalarMultiplication) {
    Matrix3D m(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Matrix3D result = m * 2.0;
    EXPECT_DOUBLE_EQ(result.at(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(result.at(1, 1), 10.0);

    Matrix3D result2 = 2.0 * m;
    EXPECT_TRUE(result2.approximately_equals(result, 1e-9));
}

// Approximately equals
TEST_F(Matrix3DTest, ApproximatelyEquals) {
    Matrix3D m1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Matrix3D m2(1.0000001, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

    EXPECT_TRUE(m1.approximately_equals(m2, 1e-6));
    EXPECT_FALSE(m1.approximately_equals(m2, 1e-9));
}

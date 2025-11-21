/**
 * @file test_reference_frame.cpp
 * @brief Unit tests for ReferenceFrame class
 */

#include <gtest/gtest.h>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>

using namespace x3dna::core;
using namespace x3dna::geometry;

class ReferenceFrameTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create identity frame at origin
        identity_frame_ = ReferenceFrame();

        // Create a test frame with rotation and translation
        Matrix3D rotation = Matrix3D::rotation_z(M_PI / 4.0); // 45 degree rotation around z
        Vector3D origin(1.0, 2.0, 3.0);
        test_frame_ = ReferenceFrame(rotation, origin);

        // Create frame from arrays
        std::array<double, 9> rot_arr = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
        std::array<double, 3> org_arr = {5.0, 6.0, 7.0};
        array_frame_ = ReferenceFrame(rot_arr, org_arr);
    }

    ReferenceFrame identity_frame_;
    ReferenceFrame test_frame_;
    ReferenceFrame array_frame_;
};

// Constructor tests
TEST_F(ReferenceFrameTest, DefaultConstructor) {
    Matrix3D expected_rot = Matrix3D::identity();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(identity_frame_.rotation().at(i, j), expected_rot.at(i, j), 1e-9);
        }
    }
    EXPECT_EQ(identity_frame_.origin(), Vector3D(0, 0, 0));
}

TEST_F(ReferenceFrameTest, MatrixVectorConstructor) {
    Matrix3D rot = Matrix3D::rotation_z(M_PI / 2.0);
    Vector3D org(10.0, 20.0, 30.0);
    ReferenceFrame frame(rot, org);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(frame.rotation().at(i, j), rot.at(i, j), 1e-9);
        }
    }
    EXPECT_EQ(frame.origin(), org);
}

TEST_F(ReferenceFrameTest, ArrayConstructor) {
    EXPECT_EQ(array_frame_.origin(), Vector3D(5.0, 6.0, 7.0));
    // Check rotation is identity
    Matrix3D expected_rot = Matrix3D::identity();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(array_frame_.rotation().at(i, j), expected_rot.at(i, j), 1e-9);
        }
    }
}

// Axis access tests
TEST_F(ReferenceFrameTest, XAxis) {
    Vector3D x_axis = identity_frame_.x_axis();
    EXPECT_NEAR(x_axis.x(), 1.0, 1e-9);
    EXPECT_NEAR(x_axis.y(), 0.0, 1e-9);
    EXPECT_NEAR(x_axis.z(), 0.0, 1e-9);
}

TEST_F(ReferenceFrameTest, YAxis) {
    Vector3D y_axis = identity_frame_.y_axis();
    EXPECT_NEAR(y_axis.x(), 0.0, 1e-9);
    EXPECT_NEAR(y_axis.y(), 1.0, 1e-9);
    EXPECT_NEAR(y_axis.z(), 0.0, 1e-9);
}

TEST_F(ReferenceFrameTest, ZAxis) {
    Vector3D z_axis = identity_frame_.z_axis();
    EXPECT_NEAR(z_axis.x(), 0.0, 1e-9);
    EXPECT_NEAR(z_axis.y(), 0.0, 1e-9);
    EXPECT_NEAR(z_axis.z(), 1.0, 1e-9);
}

// Direction dot product tests
TEST_F(ReferenceFrameTest, DirectionDotProductSame) {
    double dot = identity_frame_.direction_dot_product(identity_frame_);
    EXPECT_NEAR(dot, 1.0, 1e-9); // Same direction
}

TEST_F(ReferenceFrameTest, DirectionDotProductOpposite) {
    // Create frame with flipped z-axis
    Matrix3D flipped_rot = Matrix3D::identity();
    flipped_rot.set(2, 2, -1.0); // Flip z-axis
    ReferenceFrame flipped_frame(flipped_rot, Vector3D(0, 0, 0));

    double dot = identity_frame_.direction_dot_product(flipped_frame);
    EXPECT_NEAR(dot, -1.0, 1e-9); // Opposite direction
}

// Transform tests
TEST_F(ReferenceFrameTest, Transform) {
    Vector3D local(1.0, 0.0, 0.0);
    Vector3D global = identity_frame_.transform(local);
    EXPECT_NEAR(global.x(), 1.0, 1e-9);
    EXPECT_NEAR(global.y(), 0.0, 1e-9);
    EXPECT_NEAR(global.z(), 0.0, 1e-9);
}

TEST_F(ReferenceFrameTest, TransformWithTranslation) {
    Vector3D local(0.0, 0.0, 0.0);
    Vector3D global = test_frame_.transform(local);
    EXPECT_NEAR(global.x(), 1.0, 1e-9);
    EXPECT_NEAR(global.y(), 2.0, 1e-9);
    EXPECT_NEAR(global.z(), 3.0, 1e-9);
}

TEST_F(ReferenceFrameTest, InverseTransform) {
    Vector3D global(1.0, 2.0, 3.0);
    Vector3D local = identity_frame_.inverse_transform(global);
    EXPECT_NEAR(local.x(), 1.0, 1e-9);
    EXPECT_NEAR(local.y(), 2.0, 1e-9);
    EXPECT_NEAR(local.z(), 3.0, 1e-9);
}

TEST_F(ReferenceFrameTest, TransformRoundTrip) {
    Vector3D original(5.0, 10.0, 15.0);
    Vector3D transformed = test_frame_.transform(original);
    Vector3D back = test_frame_.inverse_transform(transformed);
    EXPECT_NEAR(back.x(), original.x(), 1e-9);
    EXPECT_NEAR(back.y(), original.y(), 1e-9);
    EXPECT_NEAR(back.z(), original.z(), 1e-9);
}

// Array conversion tests
TEST_F(ReferenceFrameTest, RotationAsArray) {
    auto arr = identity_frame_.rotation_as_array();
    EXPECT_EQ(arr.size(), 9);
    EXPECT_NEAR(arr[0], 1.0, 1e-9);
    EXPECT_NEAR(arr[4], 1.0, 1e-9);
    EXPECT_NEAR(arr[8], 1.0, 1e-9);
}

TEST_F(ReferenceFrameTest, OriginAsArray) {
    auto arr = test_frame_.origin_as_array();
    EXPECT_EQ(arr.size(), 3);
    EXPECT_NEAR(arr[0], 1.0, 1e-9);
    EXPECT_NEAR(arr[1], 2.0, 1e-9);
    EXPECT_NEAR(arr[2], 3.0, 1e-9);
}

// JSON serialization tests - Legacy format
TEST_F(ReferenceFrameTest, ToJsonLegacy) {
    auto json = test_frame_.to_json_legacy();

    EXPECT_TRUE(json.contains("orien"));
    EXPECT_TRUE(json.contains("org"));
    EXPECT_TRUE(json["orien"].is_array());
    EXPECT_EQ(json["orien"].size(), 3);
    EXPECT_TRUE(json["org"].is_array());
    EXPECT_EQ(json["org"].size(), 3);
}

TEST_F(ReferenceFrameTest, FromJsonLegacy) {
    nlohmann::json j = {{"orien", {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}},
                        {"org", {10.0, 20.0, 30.0}}};

    ReferenceFrame frame = ReferenceFrame::from_json_legacy(j);

    EXPECT_NEAR(frame.origin().x(), 10.0, 1e-9);
    EXPECT_NEAR(frame.origin().y(), 20.0, 1e-9);
    EXPECT_NEAR(frame.origin().z(), 30.0, 1e-9);

    // Check rotation is identity
    Matrix3D expected_rot = Matrix3D::identity();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(frame.rotation().at(i, j), expected_rot.at(i, j), 1e-9);
        }
    }
}

TEST_F(ReferenceFrameTest, JsonLegacyRoundTrip) {
    auto json = test_frame_.to_json_legacy();
    ReferenceFrame frame = ReferenceFrame::from_json_legacy(json);

    // Compare rotations
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(frame.rotation().at(i, j), test_frame_.rotation().at(i, j), 1e-9);
        }
    }

    // Compare origins
    EXPECT_NEAR(frame.origin().x(), test_frame_.origin().x(), 1e-9);
    EXPECT_NEAR(frame.origin().y(), test_frame_.origin().y(), 1e-9);
    EXPECT_NEAR(frame.origin().z(), test_frame_.origin().z(), 1e-9);
}

// JSON serialization tests - Modern format
TEST_F(ReferenceFrameTest, ToJsonModern) {
    auto json = test_frame_.to_json();

    EXPECT_TRUE(json.contains("rotation"));
    EXPECT_TRUE(json.contains("origin"));
    EXPECT_TRUE(json["rotation"].is_array());
    EXPECT_TRUE(json["origin"].is_array());
}

TEST_F(ReferenceFrameTest, FromJsonModern) {
    nlohmann::json j = {{"rotation", {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}},
                        {"origin", {10.0, 20.0, 30.0}}};

    ReferenceFrame frame = ReferenceFrame::from_json(j);

    EXPECT_NEAR(frame.origin().x(), 10.0, 1e-9);
    EXPECT_NEAR(frame.origin().y(), 20.0, 1e-9);
    EXPECT_NEAR(frame.origin().z(), 30.0, 1e-9);
}

TEST_F(ReferenceFrameTest, JsonModernRoundTrip) {
    auto json = test_frame_.to_json();
    ReferenceFrame frame = ReferenceFrame::from_json(json);

    // Compare rotations
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(frame.rotation().at(i, j), test_frame_.rotation().at(i, j), 1e-9);
        }
    }

    // Compare origins
    EXPECT_NEAR(frame.origin().x(), test_frame_.origin().x(), 1e-9);
    EXPECT_NEAR(frame.origin().y(), test_frame_.origin().y(), 1e-9);
    EXPECT_NEAR(frame.origin().z(), test_frame_.origin().z(), 1e-9);
}

/**
 * @file test_base_pair.cpp
 * @brief Unit tests for BasePair class
 */

#include <gtest/gtest.h>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>

using namespace x3dna::core;
using namespace x3dna::geometry;

class BasePairTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create reference frames
        // Frame 1: identity rotation at origin
        Matrix3D rot1 = Matrix3D::identity();
        Vector3D org1(0.0, 0.0, 0.0);
        frame1_ = ReferenceFrame(rot1, org1);

        // Frame 2: rotated 180 degrees around z-axis, translated
        // For 180 deg rotation around z, z-axis stays the same direction
        // So we'll rotate around y-axis to get opposite z-axes
        Matrix3D rot2 = Matrix3D::rotation_y(M_PI); // This flips z-axis
        Vector3D org2(10.0, 10.0, 10.0);
        frame2_ = ReferenceFrame(rot2, org2);

        // Create base pair
        bp_ = BasePair(1, 2, BasePairType::WATSON_CRICK);
        bp_.set_bp_type("CG");
        bp_.set_frame1(frame1_);
        bp_.set_frame2(frame2_);

        // Add hydrogen bond
        hydrogen_bond hbond;
        hbond.donor_atom = " N3 ";
        hbond.acceptor_atom = " O6 ";
        hbond.distance = 2.88;
        hbond.type = '-';
        bp_.add_hydrogen_bond(hbond);
    }

    ReferenceFrame frame1_;
    ReferenceFrame frame2_;
    BasePair bp_;
};

// Constructor tests
TEST_F(BasePairTest, DefaultConstructor) {
    BasePair bp;
    EXPECT_EQ(bp.residue_idx1(), 0);
    EXPECT_EQ(bp.residue_idx2(), 0);
    EXPECT_EQ(bp.type(), BasePairType::UNKNOWN);
    EXPECT_FALSE(bp.frame1().has_value());
    EXPECT_FALSE(bp.frame2().has_value());
}

TEST_F(BasePairTest, IndexTypeConstructor) {
    BasePair bp(5, 10, BasePairType::WOBBLE);
    EXPECT_EQ(bp.residue_idx1(), 5);
    EXPECT_EQ(bp.residue_idx2(), 10);
    EXPECT_EQ(bp.type(), BasePairType::WOBBLE);
}

// Reference frame tests
TEST_F(BasePairTest, SetGetFrames) {
    BasePair bp(1, 2, BasePairType::WATSON_CRICK);
    bp.set_frame1(frame1_);
    bp.set_frame2(frame2_);

    ASSERT_TRUE(bp.frame1().has_value());
    ASSERT_TRUE(bp.frame2().has_value());
    EXPECT_EQ(bp.frame1()->origin(), frame1_.origin());
    EXPECT_EQ(bp.frame2()->origin(), frame2_.origin());
}

// Distance/angle tests
TEST_F(BasePairTest, OriginDistance) {
    double dist = bp_.origin_distance();
    double expected = frame1_.origin().distance_to(frame2_.origin());
    EXPECT_NEAR(dist, expected, 1e-9);
    EXPECT_NEAR(dist, std::sqrt(300.0), 1e-9); // sqrt(10^2 + 10^2 + 10^2)
}

TEST_F(BasePairTest, PlaneAngle) {
    double angle = bp_.plane_angle();
    // Angle between z-axes
    // With rotation_y(M_PI), z-axis is flipped, so angle should be close to PI
    EXPECT_GT(angle, 0.0);
    EXPECT_LE(angle, M_PI);
    // Z-axes should be nearly opposite (angle close to PI)
    EXPECT_NEAR(angle, M_PI, 0.2); // Allow tolerance for floating point
}

TEST_F(BasePairTest, DirectionDotProduct) {
    double dot = bp_.direction_dot_product();
    // Z-axes should point in opposite directions for valid base pair
    // With rotation_y(M_PI), z-axes are opposite, so dot product should be negative
    EXPECT_LT(dot, 0.0); // Negative dot product
    // Should be close to -1 for opposite directions (allowing tolerance)
    EXPECT_NEAR(dot, -1.0, 0.2);
}

// Hydrogen bond tests
TEST_F(BasePairTest, AddHydrogenBond) {
    BasePair bp(1, 2, BasePairType::WATSON_CRICK);
    EXPECT_EQ(bp.hydrogen_bonds().size(), 0);

    hydrogen_bond hbond;
    hbond.donor_atom = " N3 ";
    hbond.acceptor_atom = " O6 ";
    hbond.distance = 2.88;
    hbond.type = '-';
    bp.add_hydrogen_bond(hbond);

    EXPECT_EQ(bp.hydrogen_bonds().size(), 1);
    EXPECT_EQ(bp.hydrogen_bonds()[0].donor_atom, " N3 ");
    EXPECT_EQ(bp.hydrogen_bonds()[0].acceptor_atom, " O6 ");
}

// JSON serialization tests - Legacy format
TEST_F(BasePairTest, ToJsonLegacy) {
    auto json = bp_.to_json_legacy();

    EXPECT_EQ(json["type"], "base_pair");
    EXPECT_EQ(json["base_i"], 1);
    EXPECT_EQ(json["base_j"], 2);
    EXPECT_EQ(json["bp_type"], "CG");
    EXPECT_TRUE(json.contains("orien_i"));
    EXPECT_TRUE(json.contains("orien_j"));
    EXPECT_TRUE(json.contains("org_i"));
    EXPECT_TRUE(json.contains("org_j"));
    EXPECT_TRUE(json.contains("dir_xyz"));
}

TEST_F(BasePairTest, FromJsonLegacy) {
    nlohmann::json j = {{"type", "base_pair"},
                        {"base_i", 1},
                        {"base_j", 24},
                        {"bp_type", "CG"},
                        {"orien_i", {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}},
                        {"orien_j", {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}},
                        {"org_i", {0.0, 0.0, 0.0}},
                        {"org_j", {10.0, 10.0, 10.0}},
                        {"hbonds",
                         {{{"donor_atom", " N3 "}, {"acceptor_atom", " O6 "}, {"distance", 2.88}, {"type", "-"}}}}};

    BasePair bp = BasePair::from_json_legacy(j);

    EXPECT_EQ(bp.residue_idx1(), 1);
    EXPECT_EQ(bp.residue_idx2(), 24);
    EXPECT_EQ(bp.bp_type(), "CG");
    EXPECT_EQ(bp.type(), BasePairType::WATSON_CRICK);
    EXPECT_EQ(bp.hydrogen_bonds().size(), 1);
}

TEST_F(BasePairTest, JsonLegacyRoundTrip) {
    auto json = bp_.to_json_legacy();
    BasePair bp = BasePair::from_json_legacy(json);

    EXPECT_EQ(bp.residue_idx1(), bp_.residue_idx1());
    EXPECT_EQ(bp.residue_idx2(), bp_.residue_idx2());
    EXPECT_EQ(bp.bp_type(), bp_.bp_type());
    EXPECT_EQ(bp.type(), bp_.type());
    EXPECT_EQ(bp.hydrogen_bonds().size(), bp_.hydrogen_bonds().size());
}

// JSON serialization tests - Modern format
TEST_F(BasePairTest, ToJsonModern) {
    auto json = bp_.to_json();

    EXPECT_EQ(json["residue_idx1"], 1);
    EXPECT_EQ(json["residue_idx2"], 2);
    EXPECT_EQ(json["bp_type"], "CG");
    EXPECT_TRUE(json.contains("frame1"));
    EXPECT_TRUE(json.contains("frame2"));
    EXPECT_TRUE(json.contains("hydrogen_bonds"));
}

TEST_F(BasePairTest, FromJsonModern) {
    nlohmann::json j = {{"residue_idx1", 5},
                        {"residue_idx2", 10},
                        {"bp_type", "AT"},
                        {"frame1",
                         {{"rotation", {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}},
                          {"origin", {0.0, 0.0, 0.0}}}},
                        {"hydrogen_bonds", {}}};

    BasePair bp = BasePair::from_json(j);

    EXPECT_EQ(bp.residue_idx1(), 5);
    EXPECT_EQ(bp.residue_idx2(), 10);
    EXPECT_EQ(bp.bp_type(), "AT");
    EXPECT_EQ(bp.type(), BasePairType::WATSON_CRICK);
}

TEST_F(BasePairTest, JsonModernRoundTrip) {
    auto json = bp_.to_json();
    BasePair bp = BasePair::from_json(json);

    EXPECT_EQ(bp.residue_idx1(), bp_.residue_idx1());
    EXPECT_EQ(bp.residue_idx2(), bp_.residue_idx2());
    EXPECT_EQ(bp.bp_type(), bp_.bp_type());
}

// Base pair type tests
TEST_F(BasePairTest, BasePairTypeDetection) {
    BasePair at(1, 2, BasePairType::UNKNOWN);
    at.set_bp_type("AT");
    EXPECT_EQ(at.type(), BasePairType::WATSON_CRICK);

    BasePair gc(1, 2, BasePairType::UNKNOWN);
    gc.set_bp_type("GC");
    EXPECT_EQ(gc.type(), BasePairType::WATSON_CRICK);

    BasePair gt(1, 2, BasePairType::UNKNOWN);
    gt.set_bp_type("GT");
    EXPECT_EQ(gt.type(), BasePairType::WOBBLE);
}

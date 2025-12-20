/**
 * @file test_base_pair_validator.cpp
 * @brief Unit tests for BasePairValidator
 */

#include <gtest/gtest.h>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>

using namespace x3dna::algorithms;
using namespace x3dna::core;
using namespace x3dna::geometry;

class BasePairValidatorTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create two residues with reference frames
        // Residue 1: Adenine at origin
        Matrix3D rot1 = Matrix3D::identity();
        Vector3D org1(0.0, 0.0, 0.0);
        frame1_ = ReferenceFrame(rot1, org1);

        // Residue 2: Thymine at distance (Watson-Crick pair)
        Matrix3D rot2 = Matrix3D::identity();
        // Flip z-axis for Watson-Crick pairing
        rot2.set_column(2, -rot2.column(2));
        Vector3D org2(10.0, 0.0, 0.0);
        frame2_ = ReferenceFrame(rot2, org2);

        // Create residues
        res1_ = Residue("  A", 1, "A");
        res1_.set_reference_frame(frame1_);

        res2_ = Residue("  T", 2, "A");
        res2_.set_reference_frame(frame2_);

        // Add N1/N9 atoms for dNN calculation
        res1_.add_atom(Atom(" N9 ", Vector3D(0.0, 0.0, 0.0)));
        res2_.add_atom(Atom(" N1 ", Vector3D(10.0, 0.0, 0.0)));
    }

    ReferenceFrame frame1_;
    ReferenceFrame frame2_;
    Residue res1_;
    Residue res2_;
    BasePairValidator validator_;
};

// Test basic validation
TEST_F(BasePairValidatorTest, BasicValidation) {
    ValidationResult result = validator_.validate(res1_, res2_);

    // Should calculate basic values
    EXPECT_GT(result.dorg, 0.0);
    EXPECT_GE(result.plane_angle, 0.0);
    EXPECT_LE(result.plane_angle, 90.0);
}

// Test direction vectors calculation
TEST_F(BasePairValidatorTest, DirectionVectors) {
    ValidationResult result = validator_.validate(res1_, res2_);

    // dir_z should be negative for Watson-Crick pairs (opposite z-axes)
    EXPECT_LT(result.dir_z, 0.0);
}

// Test distance checks
TEST_F(BasePairValidatorTest, DistanceChecks) {
    ValidationResult result = validator_.validate(res1_, res2_);

    // Check that distance values are calculated
    EXPECT_GE(result.dorg, 0.0);
    EXPECT_GE(result.dNN, 0.0);
    EXPECT_GE(result.d_v, 0.0);
}

// Test parameter modification
TEST_F(BasePairValidatorTest, ParameterModification) {
    ValidationParameters params = ValidationParameters::defaults();
    params.max_dorg = 5.0; // Very restrictive

    validator_.set_parameters(params);
    ValidationResult result = validator_.validate(res1_, res2_);

    // With restrictive parameters, pair should be invalid
    // (residues are 10.0 apart, max_dorg is 5.0)
    EXPECT_FALSE(result.is_valid);
}

// Test same residue (should return invalid)
TEST_F(BasePairValidatorTest, SameResidue) {
    ValidationResult result = validator_.validate(res1_, res1_);

    // Should return early and be invalid
    EXPECT_FALSE(result.is_valid);
}

// Test residues without frames
TEST_F(BasePairValidatorTest, NoFrames) {
    Residue res_no_frame("  A", 3, "A");
    ValidationResult result = validator_.validate(res1_, res_no_frame);

    // Should return invalid if no frame
    EXPECT_FALSE(result.is_valid);
}

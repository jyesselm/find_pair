/**
 * @file test_base_pair_validator_poly.cpp
 * @brief Unit tests for BasePairValidator with polymorphic types
 */

#include <gtest/gtest.h>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/residue/residue.hpp>
#include <filesystem>

using namespace x3dna::algorithms;
using namespace x3dna::io;
using namespace x3dna::core::poly;

class BasePairValidatorPolyTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Skip tests if templates directory doesn't exist
        if (!std::filesystem::exists("data/templates")) {
            GTEST_SKIP() << "Templates directory not found: data/templates";
        }
    }

    PdbParser parser;

    // Two adenines that can form a base pair (A-A)
    const std::string pair_pdb = R"(HEADER    BASE PAIR TEST
ATOM      1  P     A A   1       0.000   0.000   0.000  1.00 20.00           P
ATOM      2  O5'   A A   1       1.000   0.000   0.000  1.00 20.00           O
ATOM      3  C5'   A A   1       2.000   0.000   0.000  1.00 20.00           C
ATOM      4  C4'   A A   1       3.000   0.000   0.000  1.00 20.00           C
ATOM      5  O4'   A A   1       4.000   0.000   0.000  1.00 20.00           O
ATOM      6  C3'   A A   1       5.000   0.000   0.000  1.00 20.00           C
ATOM      7  O3'   A A   1       6.000   0.000   0.000  1.00 20.00           O
ATOM      8  C2'   A A   1       7.000   0.000   0.000  1.00 20.00           C
ATOM      9  O2'   A A   1       8.000   0.000   0.000  1.00 20.00           O
ATOM     10  C1'   A A   1       9.000   0.000   0.000  1.00 20.00           C
ATOM     11  N9    A A   1      -1.289   4.551   0.000  1.00 20.00           N
ATOM     12  C8    A A   1       0.023   4.962   0.000  1.00 20.00           C
ATOM     13  N7    A A   1       0.870   3.969   0.000  1.00 20.00           N
ATOM     14  C5    A A   1       0.071   2.833   0.000  1.00 20.00           C
ATOM     15  C6    A A   1       0.424   1.460   0.000  1.00 20.00           C
ATOM     16  N6    A A   1       1.689   1.024   0.000  1.00 20.00           N
ATOM     17  N1    A A   1      -0.700   0.641   0.000  1.00 20.00           N
ATOM     18  C2    A A   1      -1.999   1.087   0.000  1.00 20.00           C
ATOM     19  N3    A A   1      -2.342   2.364   0.001  1.00 20.00           N
ATOM     20  C4    A A   1      -1.265   3.177   0.000  1.00 20.00           C
ATOM     21  P     U A   2      10.000   0.000   0.000  1.00 20.00           P
ATOM     22  O5'   U A   2      11.000   0.000   0.000  1.00 20.00           O
ATOM     23  C5'   U A   2      12.000   0.000   0.000  1.00 20.00           C
ATOM     24  C4'   U A   2      13.000   0.000   0.000  1.00 20.00           C
ATOM     25  O4'   U A   2      14.000   0.000   0.000  1.00 20.00           O
ATOM     26  C3'   U A   2      15.000   0.000   0.000  1.00 20.00           C
ATOM     27  O3'   U A   2      16.000   0.000   0.000  1.00 20.00           O
ATOM     28  C2'   U A   2      17.000   0.000   0.000  1.00 20.00           C
ATOM     29  O2'   U A   2      18.000   0.000   0.000  1.00 20.00           O
ATOM     30  C1'   U A   2      19.000   0.000   0.000  1.00 20.00           C
ATOM     31  N1    U A   2      -0.700   7.100   0.000  1.00 20.00           N
ATOM     32  C2    U A   2      -1.999   7.546   0.000  1.00 20.00           C
ATOM     33  O2    U A   2      -2.955   6.746   0.000  1.00 20.00           O
ATOM     34  N3    U A   2      -2.245   8.903   0.000  1.00 20.00           N
ATOM     35  C4    U A   2      -1.189   9.815   0.000  1.00 20.00           C
ATOM     36  O4    U A   2      -1.381  11.019   0.000  1.00 20.00           O
ATOM     37  C5    U A   2       0.125   9.225   0.000  1.00 20.00           C
ATOM     38  C6    U A   2       0.424   7.919   0.000  1.00 20.00           C
)";
};

TEST_F(BasePairValidatorPolyTest, ValidateReturnsInvalidForNonNucleotides) {
    // Create a structure with protein residues (no nucleotides)
    const std::string protein_pdb = R"(HEADER    PROTEIN STRUCTURE
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       3.000   0.000   0.000  1.00 20.00           O
ATOM      5  N   GLY A   2       4.000   0.000   0.000  1.00 20.00           N
ATOM      6  CA  GLY A   2       5.000   0.000   0.000  1.00 20.00           C
ATOM      7  C   GLY A   2       6.000   0.000   0.000  1.00 20.00           C
ATOM      8  O   GLY A   2       7.000   0.000   0.000  1.00 20.00           O
)";

    Structure structure = parser.parse_string_poly(protein_pdb);
    ASSERT_GE(structure[0].size(), 2u);

    BasePairValidator validator;
    ValidationResult result = validator.validate(structure[0][0], structure[0][1]);

    // Should be invalid since proteins can't form base pairs
    EXPECT_FALSE(result.is_valid);
}

TEST_F(BasePairValidatorPolyTest, ValidateReturnsInvalidWithoutFrames) {
    Structure structure = parser.parse_string_poly(pair_pdb);
    ASSERT_GE(structure[0].size(), 2u);

    BasePairValidator validator;
    ValidationResult result = validator.validate(structure[0][0], structure[0][1]);

    // Should be invalid since no reference frames are set
    EXPECT_FALSE(result.is_valid);
}

TEST_F(BasePairValidatorPolyTest, ValidateWithFramesCalculated) {
    Structure structure = parser.parse_string_poly(pair_pdb);
    ASSERT_GE(structure[0].size(), 2u);

    // Calculate frames first
    BaseFrameCalculator frame_calc("data/templates");
    frame_calc.set_is_rna(true);
    frame_calc.calculate_all_frames(structure);

    // Verify frames were calculated
    auto* nuc1 = dynamic_cast<INucleotide*>(&structure[0][0]);
    auto* nuc2 = dynamic_cast<INucleotide*>(&structure[0][1]);
    ASSERT_NE(nuc1, nullptr);
    ASSERT_NE(nuc2, nullptr);
    EXPECT_TRUE(nuc1->reference_frame().has_value());
    EXPECT_TRUE(nuc2->reference_frame().has_value());

    // Now validate
    BasePairValidator validator;
    ValidationResult result = validator.validate(structure[0][0], structure[0][1]);

    // Result should have calculated values (may or may not be valid pair)
    // The key is that it doesn't crash and calculates something
    EXPECT_GE(result.dorg, 0.0);
    EXPECT_GE(result.dNN, 0.0);
    EXPECT_GE(result.plane_angle, 0.0);
    EXPECT_LE(result.plane_angle, 90.0);
}

TEST_F(BasePairValidatorPolyTest, FindN1N9PositionForPurine) {
    Structure structure = parser.parse_string_poly(pair_pdb);
    ASSERT_GE(structure[0].size(), 1u);

    // First residue is adenine (purine) - should find N9
    auto pos = BasePairValidator::find_n1_n9_position(structure[0][0]);
    EXPECT_TRUE(pos.has_value());

    // Check that the position is near where we placed N9 in the PDB
    if (pos.has_value()) {
        EXPECT_NEAR(pos->x(), -1.289, 0.1);
        EXPECT_NEAR(pos->y(), 4.551, 0.1);
        EXPECT_NEAR(pos->z(), 0.0, 0.1);
    }
}

TEST_F(BasePairValidatorPolyTest, FindN1N9PositionForPyrimidine) {
    Structure structure = parser.parse_string_poly(pair_pdb);
    ASSERT_GE(structure[0].size(), 2u);

    // Second residue is uracil (pyrimidine) - should find N1
    auto pos = BasePairValidator::find_n1_n9_position(structure[0][1]);
    EXPECT_TRUE(pos.has_value());

    // Check that the position is near where we placed N1 in the PDB
    if (pos.has_value()) {
        EXPECT_NEAR(pos->x(), -0.700, 0.1);
        EXPECT_NEAR(pos->y(), 7.100, 0.1);
        EXPECT_NEAR(pos->z(), 0.0, 0.1);
    }
}

TEST_F(BasePairValidatorPolyTest, CalculateOverlapAreaWithFrames) {
    Structure structure = parser.parse_string_poly(pair_pdb);
    ASSERT_GE(structure[0].size(), 2u);

    // Calculate frames first
    BaseFrameCalculator frame_calc("data/templates");
    frame_calc.set_is_rna(true);
    frame_calc.calculate_all_frames(structure);

    // Get frames for oave/zave calculation
    auto* nuc1 = dynamic_cast<INucleotide*>(&structure[0][0]);
    auto* nuc2 = dynamic_cast<INucleotide*>(&structure[0][1]);
    ASSERT_NE(nuc1, nullptr);
    ASSERT_NE(nuc2, nullptr);
    ASSERT_TRUE(nuc1->reference_frame().has_value());
    ASSERT_TRUE(nuc2->reference_frame().has_value());

    auto frame1 = nuc1->reference_frame().value();
    auto frame2 = nuc2->reference_frame().value();

    x3dna::geometry::Vector3D oave = (frame1.origin() + frame2.origin()) * 0.5;
    x3dna::geometry::Vector3D zave = frame1.z_axis() + frame2.z_axis();
    zave = zave / zave.length();

    // Calculate overlap
    BasePairValidator validator;
    double overlap = validator.calculate_overlap_area(structure[0][0], structure[0][1], oave, zave);

    // Overlap should be a non-negative value
    EXPECT_GE(overlap, 0.0);
}

TEST_F(BasePairValidatorPolyTest, ValidateSameResidueFails) {
    Structure structure = parser.parse_string_poly(pair_pdb);
    ASSERT_GE(structure[0].size(), 1u);

    BasePairValidator validator;
    ValidationResult result = validator.validate(structure[0][0], structure[0][0]);

    // Validating a residue against itself should fail
    EXPECT_FALSE(result.is_valid);
}

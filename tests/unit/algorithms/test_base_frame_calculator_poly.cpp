/**
 * @file test_base_frame_calculator_poly.cpp
 * @brief Unit tests for BaseFrameCalculator with polymorphic types
 */

#include <gtest/gtest.h>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/structure/residue.hpp>
#include <filesystem>

using namespace x3dna::algorithms;
using namespace x3dna::io;
using namespace x3dna::core::structure;

class BaseFrameCalculatorPolyTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Skip tests if templates directory doesn't exist
        if (!std::filesystem::exists("data/templates")) {
            GTEST_SKIP() << "Templates directory not found: data/templates";
        }
    }

    PdbParser parser;

    // RNA with proper atoms for frame calculation
    const std::string rna_pdb = R"(HEADER    RNA STRUCTURE
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
)";
};

TEST_F(BaseFrameCalculatorPolyTest, DetectRNAFromPolyStructure) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    bool is_rna = BaseFrameCalculator::detect_rna(structure);
    EXPECT_TRUE(is_rna);
}

TEST_F(BaseFrameCalculatorPolyTest, CalculateFrameForNucleotide) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    BaseFrameCalculator calculator("data/templates");
    calculator.set_is_rna(true);

    // Calculate frame for first residue
    auto result = calculator.calculate_frame(structure[0][0]);

    EXPECT_TRUE(result.is_valid);
    EXPECT_GT(result.num_matched, 0u);
}

TEST_F(BaseFrameCalculatorPolyTest, CalculateAllFramesOnPolyStructure) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    BaseFrameCalculator calculator("data/templates");
    calculator.set_is_rna(true);

    // Calculate frames for all nucleotides
    calculator.calculate_all_frames(structure);

    // Check that frame was set on the nucleotide
    auto* nuc = dynamic_cast<INucleotide*>(&structure[0][0]);
    ASSERT_NE(nuc, nullptr);
    EXPECT_TRUE(nuc->reference_frame().has_value());
}

TEST_F(BaseFrameCalculatorPolyTest, FrameNotSetOnProtein) {
    const std::string protein_pdb = R"(HEADER    PROTEIN STRUCTURE
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       3.000   0.000   0.000  1.00 20.00           O
)";

    Structure structure = parser.parse_string_poly(protein_pdb);

    BaseFrameCalculator calculator("data/templates");
    calculator.calculate_all_frames(structure);

    // Protein should not have a frame
    EXPECT_TRUE(structure[0][0].is_protein());
    EXPECT_FALSE(structure[0][0].is_nucleotide());
}

TEST_F(BaseFrameCalculatorPolyTest, CalculateFrameConstDoesNotModify) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    BaseFrameCalculator calculator("data/templates");
    calculator.set_is_rna(true);

    // Check no frame before
    auto* nuc = dynamic_cast<INucleotide*>(&structure[0][0]);
    ASSERT_NE(nuc, nullptr);
    EXPECT_FALSE(nuc->reference_frame().has_value());

    // Calculate frame (const version)
    auto result = calculator.calculate_frame_const(structure[0][0]);

    // Frame should be valid in result
    EXPECT_TRUE(result.is_valid);

    // But not set on the residue
    EXPECT_FALSE(nuc->reference_frame().has_value());
}


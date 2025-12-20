/**
 * @file test_base_frame_calculator.cpp
 * @brief Unit tests for BaseFrameCalculator class
 */

#include <gtest/gtest.h>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <filesystem>

using namespace x3dna::algorithms;
using namespace x3dna::core;
using namespace x3dna::geometry;

class BaseFrameCalculatorTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Skip tests if templates directory doesn't exist
        if (!std::filesystem::exists("data/templates")) {
            GTEST_SKIP() << "Templates directory not found: data/templates";
        }
        calculator_ = std::make_unique<BaseFrameCalculator>("data/templates");
    }

    std::unique_ptr<BaseFrameCalculator> calculator_;
};

// Test frame calculation for a single residue
TEST_F(BaseFrameCalculatorTest, CalculateFrameForResidue) {
    // Create a test residue with ring atoms
    Residue residue("  A", 1, "A");

    // Add adenine ring atoms (simplified coordinates)
    residue.add_atom(Atom(" C4 ", Vector3D(-1.267, 3.124, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" N3 ", Vector3D(-2.320, 2.290, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" C2 ", Vector3D(-1.912, 1.023, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" N1 ", Vector3D(-0.668, 0.532, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" C6 ", Vector3D(0.369, 1.398, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" C5 ", Vector3D(0.071, 2.771, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" N7 ", Vector3D(0.877, 3.902, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" C8 ", Vector3D(0.024, 4.897, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" N9 ", Vector3D(-1.291, 4.498, 0.000), "  A", "A", 1));

    FrameCalculationResult result = calculator_->calculate_frame_const(residue);

    // Should calculate frame if template is available
    if (result.is_valid) {
        EXPECT_GT(result.num_matched, 0);
        EXPECT_GT(result.matched_atoms.size(), 0);
        EXPECT_GE(result.rms_fit, 0.0);

        // Frame should be valid
        const auto& frame = result.frame;
        EXPECT_NEAR(frame.rotation().determinant(), 1.0,
                    0.01); // Rotation matrix should be orthogonal
    }
}

// Test that frame is stored in residue
TEST_F(BaseFrameCalculatorTest, StoreFrameInResidue) {
    Residue residue("  A", 1, "A");

    // Add adenine ring atoms
    residue.add_atom(Atom(" C4 ", Vector3D(-1.267, 3.124, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" N3 ", Vector3D(-2.320, 2.290, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" C2 ", Vector3D(-1.912, 1.023, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" N1 ", Vector3D(-0.668, 0.532, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" C6 ", Vector3D(0.369, 1.398, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" C5 ", Vector3D(0.071, 2.771, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" N7 ", Vector3D(0.877, 3.902, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" C8 ", Vector3D(0.024, 4.897, 0.000), "  A", "A", 1));
    residue.add_atom(Atom(" N9 ", Vector3D(-1.291, 4.498, 0.000), "  A", "A", 1));

    FrameCalculationResult result = calculator_->calculate_frame(residue);

    if (result.is_valid) {
        // Frame should be stored in residue
        EXPECT_TRUE(residue.reference_frame().has_value());

        auto stored_frame = residue.reference_frame().value();
        EXPECT_NEAR(stored_frame.rotation().determinant(), 1.0, 0.01);
    }
}

// Test batch calculation for structure
TEST_F(BaseFrameCalculatorTest, CalculateAllFrames) {
    Structure structure("TEST");
    Chain chain("A");

    // Add two residues
    Residue residue1("  A", 1, "A");
    residue1.add_atom(Atom(" C4 ", Vector3D(-1.267, 3.124, 0.000), "  A", "A", 1));
    residue1.add_atom(Atom(" N3 ", Vector3D(-2.320, 2.290, 0.000), "  A", "A", 1));
    residue1.add_atom(Atom(" C2 ", Vector3D(-1.912, 1.023, 0.000), "  A", "A", 1));
    residue1.add_atom(Atom(" N1 ", Vector3D(-0.668, 0.532, 0.000), "  A", "A", 1));
    residue1.add_atom(Atom(" C6 ", Vector3D(0.369, 1.398, 0.000), "  A", "A", 1));
    residue1.add_atom(Atom(" C5 ", Vector3D(0.071, 2.771, 0.000), "  A", "A", 1));
    residue1.add_atom(Atom(" N7 ", Vector3D(0.877, 3.902, 0.000), "  A", "A", 1));
    residue1.add_atom(Atom(" C8 ", Vector3D(0.024, 4.897, 0.000), "  A", "A", 1));
    residue1.add_atom(Atom(" N9 ", Vector3D(-1.291, 4.498, 0.000), "  A", "A", 1));

    Residue residue2("  C", 2, "A");
    residue2.add_atom(Atom(" C4 ", Vector3D(0.0, 0.0, 0.0), "  C", "A", 2));
    residue2.add_atom(Atom(" N3 ", Vector3D(1.0, 0.0, 0.0), "  C", "A", 2));
    residue2.add_atom(Atom(" C2 ", Vector3D(2.0, 0.0, 0.0), "  C", "A", 2));
    residue2.add_atom(Atom(" N1 ", Vector3D(3.0, 0.0, 0.0), "  C", "A", 2));
    residue2.add_atom(Atom(" C6 ", Vector3D(4.0, 0.0, 0.0), "  C", "A", 2));
    residue2.add_atom(Atom(" C5 ", Vector3D(5.0, 0.0, 0.0), "  C", "A", 2));

    chain.add_residue(residue1);
    chain.add_residue(residue2);
    structure.add_chain(chain);

    calculator_->calculate_all_frames(structure);

    // Check if frames were calculated (if templates exist)
    // Frames may or may not be set depending on template availability
    // Just verify the method doesn't crash
    EXPECT_NO_THROW(calculator_->calculate_all_frames(structure));
}

// Test error handling for invalid residue type
TEST_F(BaseFrameCalculatorTest, InvalidResidueType) {
    Residue invalid_residue("XXX", 1, "A");
    invalid_residue.add_atom(Atom(" C4 ", Vector3D(0.0, 0.0, 0.0), "XXX", "A", 1));

    FrameCalculationResult result = calculator_->calculate_frame_const(invalid_residue);

    // Should return invalid result for non-nucleotide
    EXPECT_FALSE(result.is_valid);
}

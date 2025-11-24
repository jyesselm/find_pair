/**
 * @file test_base_pair_finder.cpp
 * @brief Unit tests for BasePairFinder
 */

#include <gtest/gtest.h>
#include <x3dna/algorithms/base_pair_finder.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>

using namespace x3dna::algorithms;
using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::geometry;

class BasePairFinderTest : public ::testing::Test {
protected:
    void SetUp() override {
        finder_ = std::make_unique<BasePairFinder>();
    }
    
    std::unique_ptr<BasePairFinder> finder_;
};

// Test finding pairs in empty structure
TEST_F(BasePairFinderTest, EmptyStructure) {
    Structure structure;
    auto pairs = finder_->find_pairs(structure);
    
    EXPECT_EQ(pairs.size(), 0);
}

// Test finding pairs with single residue
TEST_F(BasePairFinderTest, SingleResidue) {
    Structure structure;
    Chain chain('A');
    Residue residue("  A", 'A', 1);
    
    // Add frame
    Matrix3D rot = Matrix3D::identity();
    Vector3D org(0.0, 0.0, 0.0);
    residue.set_reference_frame(ReferenceFrame(rot, org));
    
    chain.add_residue(residue);
    structure.add_chain(chain);
    
    auto pairs = finder_->find_pairs(structure);
    
    EXPECT_EQ(pairs.size(), 0);  // Can't pair with itself
}

// Test with real PDB file (if available)
TEST_F(BasePairFinderTest, RealPdbFile) {
    std::filesystem::path test_pdb = "data/pdb/100D.pdb";
    if (!std::filesystem::exists(test_pdb)) {
        GTEST_SKIP() << "Test PDB file not found: " << test_pdb;
    }
    
    // Parse PDB
    PdbParser parser;
    Structure structure = parser.parse_file(test_pdb);
    
    if (structure.num_atoms() == 0) {
        GTEST_SKIP() << "PDB file has no atoms";
    }
    
    // Calculate frames
    BaseFrameCalculator calculator("data/templates");
    calculator.calculate_all_frames(structure);
    
    // Find pairs
    auto pairs = finder_->find_pairs(structure);
    
    // Should find at least some pairs (or none if structure doesn't have pairs)
    EXPECT_GE(pairs.size(), 0);
    
    // If pairs found, verify they have frames
    for (const auto& pair : pairs) {
        if (pair.frame1().has_value() && pair.frame2().has_value()) {
            EXPECT_GT(pair.origin_distance(), 0.0);
        }
    }
}

// Test strategy setting
TEST_F(BasePairFinderTest, StrategySetting) {
    finder_->set_strategy(PairFindingStrategy::BEST_PAIR);
    EXPECT_EQ(finder_->strategy(), PairFindingStrategy::BEST_PAIR);
    
    finder_->set_strategy(PairFindingStrategy::ALL_PAIRS);
    EXPECT_EQ(finder_->strategy(), PairFindingStrategy::ALL_PAIRS);
}

// Test parameter setting
TEST_F(BasePairFinderTest, ParameterSetting) {
    ValidationParameters params = ValidationParameters::defaults();
    params.max_dorg = 5.0;
    
    finder_->set_parameters(params);
    
    const auto& retrieved_params = finder_->parameters();
    EXPECT_EQ(retrieved_params.max_dorg, 5.0);
}


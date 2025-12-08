/**
 * @file test_pair_candidate_cache.cpp
 * @brief Unit tests for PairCandidateCache
 */

#include <gtest/gtest.h>
#include <x3dna/algorithms/pair_candidate_cache.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/base_pair_finder.hpp>
#include <filesystem>

using namespace x3dna::algorithms;
using namespace x3dna::core;
using namespace x3dna::io;

class PairCandidateCacheTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Load a real PDB file for testing
        pdb_path_ = std::filesystem::path(X3DNA_SOURCE_DIR) / "data" / "pdb" / "100D.pdb";
        
        if (std::filesystem::exists(pdb_path_)) {
            PdbParser parser;
            structure_ = parser.parse_file(pdb_path_);
            
            // Calculate frames for all residues
            std::filesystem::path template_path = 
                std::filesystem::path(X3DNA_SOURCE_DIR) / "data" / "templates";
            BaseFrameCalculator frame_calc(template_path);
            frame_calc.calculate_all_frames(structure_);
            
            has_structure_ = true;
        }
    }

    std::filesystem::path pdb_path_;
    Structure structure_;
    bool has_structure_ = false;
};

TEST_F(PairCandidateCacheTest, BuildFromStructure) {
    if (!has_structure_) {
        GTEST_SKIP() << "PDB file not available";
    }

    PairCandidateCache cache;
    BasePairValidator validator;
    QualityScoreCalculator quality_calc;

    cache.build(structure_, validator, quality_calc, BasePairFinder::is_nucleotide);

    EXPECT_FALSE(cache.empty());
    EXPECT_GT(cache.size(), 0);
    EXPECT_GT(cache.max_legacy_idx(), 0);
}

TEST_F(PairCandidateCacheTest, GetPairInfo) {
    if (!has_structure_) {
        GTEST_SKIP() << "PDB file not available";
    }

    PairCandidateCache cache;
    BasePairValidator validator;
    QualityScoreCalculator quality_calc;

    cache.build(structure_, validator, quality_calc, BasePairFinder::is_nucleotide);

    // Get info for a pair that should exist
    // 100D typically has pairs like (1, 24), (2, 23), etc.
    auto info = cache.get(1, 24);
    
    // The pair might or might not be valid, but it should exist
    // if both residues are nucleotides
    if (info) {
        // Just verify the structure is populated
        EXPECT_TRUE(info->bp_type_id >= -1 && info->bp_type_id <= 2);
    }
}

TEST_F(PairCandidateCacheTest, GetOrderIndependent) {
    if (!has_structure_) {
        GTEST_SKIP() << "PDB file not available";
    }

    PairCandidateCache cache;
    BasePairValidator validator;
    QualityScoreCalculator quality_calc;

    cache.build(structure_, validator, quality_calc, BasePairFinder::is_nucleotide);

    // Get same pair in both orders - should return same result
    auto info1 = cache.get(1, 24);
    auto info2 = cache.get(24, 1);

    EXPECT_EQ(info1.has_value(), info2.has_value());
    
    if (info1 && info2) {
        EXPECT_EQ(info1->is_valid(), info2->is_valid());
        EXPECT_EQ(info1->bp_type_id, info2->bp_type_id);
        EXPECT_DOUBLE_EQ(info1->adjusted_quality_score, info2->adjusted_quality_score);
    }
}

TEST_F(PairCandidateCacheTest, ValidPartnersFor) {
    if (!has_structure_) {
        GTEST_SKIP() << "PDB file not available";
    }

    PairCandidateCache cache;
    BasePairValidator validator;
    QualityScoreCalculator quality_calc;

    cache.build(structure_, validator, quality_calc, BasePairFinder::is_nucleotide);

    // Get valid partners for residue 1
    auto partners = cache.valid_partners_for(1);
    
    // Each partner should have a valid pair with residue 1
    for (int partner_idx : partners) {
        auto info = cache.get(1, partner_idx);
        ASSERT_TRUE(info.has_value());
        EXPECT_TRUE(info->is_valid());
    }
}

TEST_F(PairCandidateCacheTest, ValidCount) {
    if (!has_structure_) {
        GTEST_SKIP() << "PDB file not available";
    }

    PairCandidateCache cache;
    BasePairValidator validator;
    QualityScoreCalculator quality_calc;

    cache.build(structure_, validator, quality_calc, BasePairFinder::is_nucleotide);

    size_t valid_count = cache.valid_count();
    
    // Count manually
    size_t manual_count = 0;
    for (const auto& [key, info] : cache.all()) {
        if (info.is_valid()) {
            manual_count++;
        }
    }

    EXPECT_EQ(valid_count, manual_count);
}

TEST_F(PairCandidateCacheTest, ForEachValid) {
    if (!has_structure_) {
        GTEST_SKIP() << "PDB file not available";
    }

    PairCandidateCache cache;
    BasePairValidator validator;
    QualityScoreCalculator quality_calc;

    cache.build(structure_, validator, quality_calc, BasePairFinder::is_nucleotide);

    size_t callback_count = 0;
    cache.for_each_valid([&](int idx1, int idx2, const CandidateInfo& info) {
        EXPECT_TRUE(info.is_valid());
        EXPECT_LT(idx1, idx2);  // Should be normalized
        callback_count++;
    });

    EXPECT_EQ(callback_count, cache.valid_count());
}

TEST_F(PairCandidateCacheTest, Clear) {
    if (!has_structure_) {
        GTEST_SKIP() << "PDB file not available";
    }

    PairCandidateCache cache;
    BasePairValidator validator;
    QualityScoreCalculator quality_calc;

    cache.build(structure_, validator, quality_calc, BasePairFinder::is_nucleotide);
    
    EXPECT_FALSE(cache.empty());
    
    cache.clear();
    
    EXPECT_TRUE(cache.empty());
    EXPECT_EQ(cache.size(), 0);
    EXPECT_EQ(cache.valid_count(), 0);
}

TEST_F(PairCandidateCacheTest, IndexMapAccess) {
    if (!has_structure_) {
        GTEST_SKIP() << "PDB file not available";
    }

    PairCandidateCache cache;
    BasePairValidator validator;
    QualityScoreCalculator quality_calc;

    cache.build(structure_, validator, quality_calc, BasePairFinder::is_nucleotide);

    // Index map should be populated
    const auto& index_map = cache.index_map();
    EXPECT_FALSE(index_map.empty());
    EXPECT_EQ(cache.max_legacy_idx(), index_map.max_legacy_idx());
}

TEST_F(PairCandidateCacheTest, EmptyStructure) {
    Structure empty_structure;
    
    PairCandidateCache cache;
    BasePairValidator validator;
    QualityScoreCalculator quality_calc;

    cache.build(empty_structure, validator, quality_calc, BasePairFinder::is_nucleotide);

    EXPECT_TRUE(cache.empty());
    EXPECT_EQ(cache.size(), 0);
    EXPECT_EQ(cache.valid_count(), 0);
}


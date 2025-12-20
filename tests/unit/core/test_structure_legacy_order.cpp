/**
 * @file test_structure_legacy_order.cpp
 * @brief Test that residue ordering matches legacy (PDB file order)
 *
 * This test verifies that:
 * 1. Residue counts match legacy (includes HETATMs, waters, etc.)
 * 2. Residue indices match legacy (same residue at same index)
 * 3. Ordering is preserved (PDB file order, not sorted)
 */

#include <gtest/gtest.h>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <filesystem>

using namespace x3dna::core;
using namespace x3dna::io;

class StructureLegacyOrderTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Use 3G8T as test case (known to have 1070 residues in legacy)
        pdb_file_ = "data/pdb/3G8T.pdb";

        // Check if test data exists
        if (!std::filesystem::exists(pdb_file_)) {
            GTEST_SKIP() << "Test data file not found: " << pdb_file_;
        }
    }

    std::filesystem::path pdb_file_;
};

/**
 * @test Verify total residue count matches legacy
 *
 * Legacy counts ALL residues including HETATMs and waters.
 * Modern should match this count when include_hetatm and include_waters are true.
 */
TEST_F(StructureLegacyOrderTest, ResidueCountMatchesLegacy) {
    // Parse with legacy-compatible settings
    PdbParser parser;
    parser.set_include_hetatm(true); // Include HETATM records (legacy includes all)
    parser.set_include_waters(true); // Include water molecules (legacy includes all)
    Structure structure = parser.parse_file(pdb_file_);

    // Get residues in legacy order (using Structure's built-in method)
    auto residues = structure.residues_in_legacy_order();

    // Legacy counts 1070 residues for 3G8T
    // This is the known count from legacy's residue_idx() function
    EXPECT_EQ(residues.size(), 1070) << "Residue count should match legacy (1070 for 3G8T)";
}

/**
 * @test Verify specific residue indices match legacy
 *
 * Test that known residue indices point to the correct residues.
 * Legacy index 946 = C (chain S, seq 113)
 * Legacy index 947 = U (chain S, seq 114)
 */
TEST_F(StructureLegacyOrderTest, SpecificResidueIndicesMatchLegacy) {
    // Parse with legacy-compatible settings
    PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(true);
    Structure structure = parser.parse_file(pdb_file_);

    // Test known residue indices from legacy
    const Residue* res946 = structure.get_residue_by_legacy_idx(946);
    const Residue* res947 = structure.get_residue_by_legacy_idx(947);

    ASSERT_NE(res946, nullptr) << "Residue 946 should exist";
    ASSERT_NE(res947, nullptr) << "Residue 947 should exist";

    // Legacy index 946 = C (chain S, seq 113)
    EXPECT_EQ(res946->name(), "  C") << "Residue 946 should be C";
    EXPECT_EQ(res946->chain_id(), "S") << "Residue 946 should be in chain S";
    EXPECT_EQ(res946->seq_num(), 113) << "Residue 946 should have seq 113";

    // Legacy index 947 = U (chain S, seq 114)
    EXPECT_EQ(res947->name(), "  U") << "Residue 947 should be U";
    EXPECT_EQ(res947->chain_id(), "S") << "Residue 947 should be in chain S";
    EXPECT_EQ(res947->seq_num(), 114) << "Residue 947 should have seq 114";
}

/**
 * @test Verify residue ordering is consistent
 *
 * Test that residues appear in the same order when retrieved multiple times.
 * This ensures the ordering function is deterministic.
 */
TEST_F(StructureLegacyOrderTest, ResidueOrderingIsConsistent) {
    // Parse with legacy-compatible settings
    PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(true);
    Structure structure = parser.parse_file(pdb_file_);

    // Get residues in legacy order twice (using Structure's built-in method)
    auto residues1 = structure.residues_in_legacy_order();
    auto residues2 = structure.residues_in_legacy_order();

    ASSERT_EQ(residues1.size(), residues2.size()) << "Both calls should return same number of residues";

    // Verify all residues match
    for (size_t i = 0; i < residues1.size(); i++) {
        EXPECT_EQ(residues1[i], residues2[i]) << "Residue at index " << i << " should be the same in both calls";

        if (residues1[i] != residues2[i]) {
            break; // Stop after first mismatch for cleaner output
        }
    }
}

/**
 * @test Verify get_legacy_idx_for_residue works correctly
 *
 * Test that we can get the legacy index for a known residue.
 */
TEST_F(StructureLegacyOrderTest, GetLegacyIdxForResidue) {
    // Parse with legacy-compatible settings
    PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(true);
    Structure structure = parser.parse_file(pdb_file_);

    // Get residue at known index
    const Residue* res946 = structure.get_residue_by_legacy_idx(946);
    ASSERT_NE(res946, nullptr) << "Residue 946 should exist";

    // Get legacy index for this residue
    int legacy_idx = structure.get_legacy_idx_for_residue(res946);

    EXPECT_EQ(legacy_idx, 946) << "get_legacy_idx_for_residue should return 946 for residue 946";
}

/**
 * @test Verify edge cases
 *
 * Test that invalid indices return nullptr and zero.
 */
TEST_F(StructureLegacyOrderTest, EdgeCases) {
    // Parse with legacy-compatible settings
    PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(true);
    Structure structure = parser.parse_file(pdb_file_);

    // Test invalid indices
    const Residue* res0 = structure.get_residue_by_legacy_idx(0);
    EXPECT_EQ(res0, nullptr) << "Index 0 should return nullptr";

    const Residue* res_negative = structure.get_residue_by_legacy_idx(-1);
    EXPECT_EQ(res_negative, nullptr) << "Negative index should return nullptr";

    const Residue* res_too_large = structure.get_residue_by_legacy_idx(10000);
    EXPECT_EQ(res_too_large, nullptr) << "Index too large should return nullptr";

    // Test get_legacy_idx_for_residue with nullptr
    int idx_null = structure.get_legacy_idx_for_residue(nullptr);
    EXPECT_EQ(idx_null, 0) << "nullptr residue should return 0";
}

/**
 * @test Verify parser settings affect residue count
 *
 * Test that excluding HETATMs or waters changes the residue count.
 * This ensures the test is actually checking the right thing.
 */
TEST_F(StructureLegacyOrderTest, ParserSettingsAffectCount) {
    // Parse WITHOUT HETATMs and waters
    PdbParser parser_excluded;
    parser_excluded.set_include_hetatm(false);
    parser_excluded.set_include_waters(false);
    Structure structure_excluded = parser_excluded.parse_file(pdb_file_);
    auto residues_excluded = structure_excluded.residues_in_legacy_order();

    // Parse WITH HETATMs and waters (legacy-compatible)
    PdbParser parser_included;
    parser_included.set_include_hetatm(true);
    parser_included.set_include_waters(true);
    Structure structure_included = parser_included.parse_file(pdb_file_);
    auto residues_included = structure_included.residues_in_legacy_order();

    // Included should have more residues
    EXPECT_GT(residues_included.size(), residues_excluded.size())
        << "Including HETATMs and waters should increase residue count";

    // Included should match legacy count (1070)
    EXPECT_EQ(residues_included.size(), 1070) << "With HETATMs and waters included, count should match legacy";
}

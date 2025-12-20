/**
 * @file test_residue_index_map.cpp
 * @brief Unit tests for ResidueIndexMap
 */

#include <gtest/gtest.h>
#include <x3dna/algorithms/residue_index_map.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>

using namespace x3dna::algorithms;
using namespace x3dna::core;
using namespace x3dna::geometry;

class ResidueIndexMapTest : public ::testing::Test {
protected:
    // Helper to create an atom with legacy residue index
    Atom create_atom(const std::string& name, const Vector3D& pos, int legacy_idx) {
        Atom atom(name, pos);
        atom.set_legacy_residue_idx(legacy_idx);
        return atom;
    }

    // Create a simple structure for testing
    Structure create_test_structure() {
        Structure structure;

        // Create chain A with 3 residues
        Chain chain_a("A");

        // Residue 1 (legacy idx 1)
        Residue res1("A", 1, "A");
        res1.add_atom(create_atom(" N9 ", Vector3D(0, 0, 0), 1));
        chain_a.add_residue(res1);

        // Residue 2 (legacy idx 2)
        Residue res2("G", 2, "A");
        res2.add_atom(create_atom(" N9 ", Vector3D(1, 0, 0), 2));
        chain_a.add_residue(res2);

        // Residue 3 (legacy idx 3)
        Residue res3("C", 3, "A");
        res3.add_atom(create_atom(" N1 ", Vector3D(2, 0, 0), 3));
        chain_a.add_residue(res3);

        structure.add_chain(chain_a);

        // Create chain B with 2 residues (legacy idx continues from 4)
        Chain chain_b("B");

        // Residue 4 (legacy idx 4)
        Residue res4("U", 1, "B");
        res4.add_atom(create_atom(" N1 ", Vector3D(3, 0, 0), 4));
        chain_b.add_residue(res4);

        // Residue 5 (legacy idx 5)
        Residue res5("T", 2, "B");
        res5.add_atom(create_atom(" N1 ", Vector3D(4, 0, 0), 5));
        chain_b.add_residue(res5);

        structure.add_chain(chain_b);

        return structure;
    }
};

TEST_F(ResidueIndexMapTest, BuildFromStructure) {
    Structure structure = create_test_structure();
    ResidueIndexMap index_map;

    index_map.build(structure);

    EXPECT_EQ(index_map.size(), 5);
    EXPECT_EQ(index_map.min_legacy_idx(), 1);
    EXPECT_EQ(index_map.max_legacy_idx(), 5);
    EXPECT_FALSE(index_map.empty());
}

TEST_F(ResidueIndexMapTest, GetByLegacyIdx) {
    Structure structure = create_test_structure();
    ResidueIndexMap index_map;
    index_map.build(structure);

    // Valid indices
    const Residue* res1 = index_map.get_by_legacy_idx(1);
    ASSERT_NE(res1, nullptr);
    EXPECT_EQ(res1->name(), "A");

    const Residue* res3 = index_map.get_by_legacy_idx(3);
    ASSERT_NE(res3, nullptr);
    EXPECT_EQ(res3->name(), "C");

    const Residue* res5 = index_map.get_by_legacy_idx(5);
    ASSERT_NE(res5, nullptr);
    EXPECT_EQ(res5->name(), "T");

    // Invalid indices
    EXPECT_EQ(index_map.get_by_legacy_idx(0), nullptr);
    EXPECT_EQ(index_map.get_by_legacy_idx(6), nullptr);
    EXPECT_EQ(index_map.get_by_legacy_idx(-1), nullptr);
}

TEST_F(ResidueIndexMapTest, GetByModernIdx) {
    Structure structure = create_test_structure();
    ResidueIndexMap index_map;
    index_map.build(structure);

    // Valid indices (0-based)
    const Residue* res0 = index_map.get_by_modern_idx(0);
    ASSERT_NE(res0, nullptr);
    EXPECT_EQ(res0->name(), "A");

    const Residue* res4 = index_map.get_by_modern_idx(4);
    ASSERT_NE(res4, nullptr);
    EXPECT_EQ(res4->name(), "T");

    // Invalid indices
    EXPECT_EQ(index_map.get_by_modern_idx(5), nullptr);
    EXPECT_EQ(index_map.get_by_modern_idx(100), nullptr);
}

TEST_F(ResidueIndexMapTest, HasIndices) {
    Structure structure = create_test_structure();
    ResidueIndexMap index_map;
    index_map.build(structure);

    // Legacy indices
    EXPECT_TRUE(index_map.has_legacy_idx(1));
    EXPECT_TRUE(index_map.has_legacy_idx(5));
    EXPECT_FALSE(index_map.has_legacy_idx(0));
    EXPECT_FALSE(index_map.has_legacy_idx(6));

    // Modern indices
    EXPECT_TRUE(index_map.has_modern_idx(0));
    EXPECT_TRUE(index_map.has_modern_idx(4));
    EXPECT_FALSE(index_map.has_modern_idx(5));
}

TEST_F(ResidueIndexMapTest, ToModern) {
    Structure structure = create_test_structure();
    ResidueIndexMap index_map;
    index_map.build(structure);

    // Valid conversions
    auto modern1 = index_map.to_modern(1);
    ASSERT_TRUE(modern1.has_value());
    EXPECT_EQ(*modern1, 0);

    auto modern5 = index_map.to_modern(5);
    ASSERT_TRUE(modern5.has_value());
    EXPECT_EQ(*modern5, 4);

    // Invalid conversions
    EXPECT_FALSE(index_map.to_modern(0).has_value());
    EXPECT_FALSE(index_map.to_modern(6).has_value());
}

TEST_F(ResidueIndexMapTest, ToLegacy) {
    Structure structure = create_test_structure();
    ResidueIndexMap index_map;
    index_map.build(structure);

    // Valid conversions
    auto legacy0 = index_map.to_legacy(0);
    ASSERT_TRUE(legacy0.has_value());
    EXPECT_EQ(*legacy0, 1);

    auto legacy4 = index_map.to_legacy(4);
    ASSERT_TRUE(legacy4.has_value());
    EXPECT_EQ(*legacy4, 5);

    // Invalid conversions
    EXPECT_FALSE(index_map.to_legacy(5).has_value());
    EXPECT_FALSE(index_map.to_legacy(100).has_value());
}

TEST_F(ResidueIndexMapTest, LegacyIndices) {
    Structure structure = create_test_structure();
    ResidueIndexMap index_map;
    index_map.build(structure);

    std::vector<int> indices = index_map.legacy_indices();

    EXPECT_EQ(indices.size(), 5);
    // Should be in ascending order (from std::map)
    EXPECT_EQ(indices[0], 1);
    EXPECT_EQ(indices[1], 2);
    EXPECT_EQ(indices[2], 3);
    EXPECT_EQ(indices[3], 4);
    EXPECT_EQ(indices[4], 5);
}

TEST_F(ResidueIndexMapTest, NucleotideLegacyIndices) {
    Structure structure = create_test_structure();
    ResidueIndexMap index_map;
    index_map.build(structure);

    // Filter for purines only (A and G)
    auto is_purine = [](const Residue& res) {
        return res.name() == "A" || res.name() == "G";
    };

    std::vector<int> purine_indices = index_map.nucleotide_legacy_indices(is_purine);

    EXPECT_EQ(purine_indices.size(), 2);
    // Should contain indices 1 (A) and 2 (G)
    EXPECT_TRUE(std::find(purine_indices.begin(), purine_indices.end(), 1) != purine_indices.end());
    EXPECT_TRUE(std::find(purine_indices.begin(), purine_indices.end(), 2) != purine_indices.end());
}

TEST_F(ResidueIndexMapTest, Clear) {
    Structure structure = create_test_structure();
    ResidueIndexMap index_map;
    index_map.build(structure);

    EXPECT_FALSE(index_map.empty());
    EXPECT_EQ(index_map.size(), 5);

    index_map.clear();

    EXPECT_TRUE(index_map.empty());
    EXPECT_EQ(index_map.size(), 0);
    EXPECT_EQ(index_map.max_legacy_idx(), 0);
    EXPECT_EQ(index_map.min_legacy_idx(), 0);
}

TEST_F(ResidueIndexMapTest, EmptyStructure) {
    Structure structure; // Empty
    ResidueIndexMap index_map;
    index_map.build(structure);

    EXPECT_TRUE(index_map.empty());
    EXPECT_EQ(index_map.size(), 0);
    EXPECT_EQ(index_map.get_by_legacy_idx(1), nullptr);
}

TEST_F(ResidueIndexMapTest, AllIteration) {
    Structure structure = create_test_structure();
    ResidueIndexMap index_map;
    index_map.build(structure);

    const auto& all = index_map.all();
    EXPECT_EQ(all.size(), 5);

    // Verify iteration order (ascending by legacy index)
    int expected_legacy = 1;
    for (const auto& [legacy_idx, residue] : all) {
        EXPECT_EQ(legacy_idx, expected_legacy);
        EXPECT_NE(residue, nullptr);
        expected_legacy++;
    }
}

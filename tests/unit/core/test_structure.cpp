/**
 * @file test_structure.cpp
 * @brief Unit tests for Structure class
 */

#include <gtest/gtest.h>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>

using namespace x3dna::core;
using namespace x3dna::geometry;

class StructureTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create structure with two chains
        structure_ = Structure("TEST");

        // Chain A
        Chain chain_a('A');
        Residue c1("  C", 1, 'A');
        c1.add_atom(Atom(" C1'", Vector3D(1, 2, 3), "  C", 'A', 1));
        chain_a.add_residue(c1);

        Residue g1("  G", 2, 'A');
        g1.add_atom(Atom(" C1'", Vector3D(4, 5, 6), "  G", 'A', 2));
        chain_a.add_residue(g1);
        structure_.add_chain(chain_a);

        // Chain B
        Chain chain_b('B');
        Residue a1("  A", 1, 'B');
        a1.add_atom(Atom(" C1'", Vector3D(7, 8, 9), "  A", 'B', 1));
        chain_b.add_residue(a1);
        structure_.add_chain(chain_b);
    }

    Structure structure_;
};

// Constructor tests
TEST_F(StructureTest, DefaultConstructor) {
    Structure structure;
    EXPECT_EQ(structure.pdb_id(), "");
    EXPECT_EQ(structure.num_chains(), 0);
    EXPECT_EQ(structure.num_residues(), 0);
    EXPECT_EQ(structure.num_atoms(), 0);
}

TEST_F(StructureTest, PdbIdConstructor) {
    Structure structure("157D");
    EXPECT_EQ(structure.pdb_id(), "157D");
    EXPECT_EQ(structure.num_chains(), 0);
}

// Chain management tests
TEST_F(StructureTest, AddChain) {
    Structure structure("TEST");
    EXPECT_EQ(structure.num_chains(), 0);

    Chain chain('A');
    structure.add_chain(chain);
    EXPECT_EQ(structure.num_chains(), 1);
}

TEST_F(StructureTest, NumChains) {
    EXPECT_EQ(structure_.num_chains(), 2);
}

TEST_F(StructureTest, NumResidues) {
    EXPECT_EQ(structure_.num_residues(), 3); // 2 in chain A, 1 in chain B
}

TEST_F(StructureTest, NumAtoms) {
    EXPECT_EQ(structure_.num_atoms(), 3); // 1 atom per residue
}

TEST_F(StructureTest, FindChain) {
    auto chain = structure_.find_chain('A');
    ASSERT_TRUE(chain.has_value());
    EXPECT_EQ(chain->chain_id(), 'A');
    EXPECT_EQ(chain->num_residues(), 2);
}

TEST_F(StructureTest, FindChainNotFound) {
    auto chain = structure_.find_chain('Z');
    EXPECT_FALSE(chain.has_value());
}

// Residue access tests
TEST_F(StructureTest, AllResidues) {
    auto residues = structure_.all_residues();
    EXPECT_EQ(residues.size(), 3);
}

TEST_F(StructureTest, Nucleotides) {
    auto nts = structure_.nucleotides();
    EXPECT_EQ(nts.size(), 3); // All are nucleotides
    for (const auto* nt : nts) {
        EXPECT_TRUE(nt->is_nucleotide());
    }
}

// JSON serialization tests - Legacy format
TEST_F(StructureTest, ToJsonLegacy) {
    auto json = structure_.to_json_legacy();

    EXPECT_EQ(json["pdb_id"], "TEST");
    EXPECT_EQ(json["num_atoms"], 3);
    EXPECT_EQ(json["num_residues"], 3);
    EXPECT_EQ(json["num_chains"], 2);
    EXPECT_TRUE(json.contains("atoms"));
    EXPECT_TRUE(json["atoms"].is_array());
    EXPECT_EQ(json["atoms"].size(), 3);
}

TEST_F(StructureTest, FromJsonLegacy) {
    nlohmann::json j = {{"pdb_id", "157D"},
                        {"num_atoms", 2},
                        {"num_residues", 2},
                        {"num_chains", 1},
                        {"atoms",
                         {{{"atom_name", " C1'"},
                           {"residue_name", "  C"},
                           {"chain_id", "A"},
                           {"residue_seq", 1},
                           {"xyz", {1.0, 2.0, 3.0}}},
                          {{"atom_name", " N1 "},
                           {"residue_name", "  C"},
                           {"chain_id", "A"},
                           {"residue_seq", 1},
                           {"xyz", {2.0, 3.0, 4.0}}}}}};

    Structure structure = Structure::from_json_legacy(j);

    EXPECT_EQ(structure.pdb_id(), "157D");
    EXPECT_EQ(structure.num_chains(), 1);
    EXPECT_EQ(structure.num_residues(), 1); // Both atoms in same residue
    EXPECT_EQ(structure.num_atoms(), 2);
}

TEST_F(StructureTest, JsonLegacyRoundTrip) {
    auto json = structure_.to_json_legacy();
    Structure structure = Structure::from_json_legacy(json);

    EXPECT_EQ(structure.pdb_id(), structure_.pdb_id());
    EXPECT_EQ(structure.num_chains(), structure_.num_chains());
    EXPECT_EQ(structure.num_residues(), structure_.num_residues());
    EXPECT_EQ(structure.num_atoms(), structure_.num_atoms());
}

// JSON serialization tests - Modern format
TEST_F(StructureTest, ToJsonModern) {
    auto json = structure_.to_json();

    EXPECT_EQ(json["pdb_id"], "TEST");
    EXPECT_TRUE(json.contains("chains"));
    EXPECT_TRUE(json["chains"].is_array());
    EXPECT_EQ(json["chains"].size(), 2);
}

TEST_F(StructureTest, FromJsonModern) {
    nlohmann::json j = {{"pdb_id", "100D"},
                        {"chains",
                         {{{"chain_id", "A"},
                           {"residues", {{{"name", "  C"}, {"seq_num", 1}, {"chain_id", "A"}, {"atoms", {}}}}}}}}};

    Structure structure = Structure::from_json(j);

    EXPECT_EQ(structure.pdb_id(), "100D");
    EXPECT_EQ(structure.num_chains(), 1);
    EXPECT_EQ(structure.num_residues(), 1);
}

TEST_F(StructureTest, JsonModernRoundTrip) {
    auto json = structure_.to_json();
    Structure structure = Structure::from_json(json);

    EXPECT_EQ(structure.pdb_id(), structure_.pdb_id());
    EXPECT_EQ(structure.num_chains(), structure_.num_chains());
    EXPECT_EQ(structure.num_residues(), structure_.num_residues());
}

// Edge cases
TEST_F(StructureTest, EmptyStructure) {
    Structure structure("EMPTY");
    EXPECT_EQ(structure.num_chains(), 0);
    EXPECT_EQ(structure.num_residues(), 0);
    EXPECT_EQ(structure.num_atoms(), 0);
    EXPECT_TRUE(structure.all_residues().empty());
    EXPECT_TRUE(structure.nucleotides().empty());
}

TEST_F(StructureTest, MultipleResiduesPerChain) {
    Structure structure("MULTI");
    Chain chain('A');

    for (int i = 1; i <= 5; ++i) {
        Residue residue("  C", i, 'A');
        residue.add_atom(Atom(" C1'", Vector3D(i, i, i)));
        chain.add_residue(residue);
    }

    structure.add_chain(chain);
    EXPECT_EQ(structure.num_residues(), 5);
    EXPECT_EQ(structure.num_atoms(), 5);
}

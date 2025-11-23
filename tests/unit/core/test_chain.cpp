/**
 * @file test_chain.cpp
 * @brief Unit tests for Chain class
 */

#include <gtest/gtest.h>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>

using namespace x3dna::core;
using namespace x3dna::geometry;

class ChainTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create chain A with nucleotides
        chain_a_ = Chain('A');

        Residue c1("  C", 1, 'A');
        c1.add_atom(Atom(" C1'", Vector3D(1, 2, 3)));
        c1.add_atom(Atom(" N1 ", Vector3D(2, 3, 4)));
        chain_a_.add_residue(c1);

        Residue g1("  G", 2, 'A');
        g1.add_atom(Atom(" C1'", Vector3D(5, 6, 7)));
        g1.add_atom(Atom(" N9 ", Vector3D(6, 7, 8)));
        chain_a_.add_residue(g1);

        Residue a1("  A", 3, 'A');
        a1.add_atom(Atom(" C1'", Vector3D(9, 10, 11)));
        a1.add_atom(Atom(" N9 ", Vector3D(10, 11, 12)));
        chain_a_.add_residue(a1);
    }

    Chain chain_a_;
};

// Constructor tests
TEST_F(ChainTest, DefaultConstructor) {
    Chain chain;
    EXPECT_EQ(chain.chain_id(), '\0');
    EXPECT_EQ(chain.num_residues(), 0);
    EXPECT_EQ(chain.num_atoms(), 0);
}

TEST_F(ChainTest, ChainIdConstructor) {
    Chain chain('B');
    EXPECT_EQ(chain.chain_id(), 'B');
    EXPECT_EQ(chain.num_residues(), 0);
}

// Residue management tests
TEST_F(ChainTest, AddResidue) {
    Chain chain('A');
    EXPECT_EQ(chain.num_residues(), 0);

    Residue residue("  C", 1, 'A');
    chain.add_residue(residue);
    EXPECT_EQ(chain.num_residues(), 1);

    Residue residue2("  G", 2, 'A');
    chain.add_residue(residue2);
    EXPECT_EQ(chain.num_residues(), 2);
}

TEST_F(ChainTest, NumAtoms) {
    EXPECT_EQ(chain_a_.num_atoms(), 6); // 2 atoms per residue * 3 residues
}

TEST_F(ChainTest, FindResidue) {
    auto residue = chain_a_.find_residue(2);
    ASSERT_TRUE(residue.has_value());
    EXPECT_EQ(residue->name(), "  G");
    EXPECT_EQ(residue->seq_num(), 2);
}

TEST_F(ChainTest, FindResidueNotFound) {
    auto residue = chain_a_.find_residue(99);
    EXPECT_FALSE(residue.has_value());
}

// Sequence tests
TEST_F(ChainTest, Sequence) {
    std::string seq = chain_a_.sequence();
    EXPECT_EQ(seq, "CGA");
}

TEST_F(ChainTest, SequenceWithNonNucleotides) {
    Chain chain('A');
    Residue c1("  C", 1, 'A');
    chain.add_residue(c1);

    Residue unknown("XXX", 2, 'A');
    chain.add_residue(unknown);

    Residue g1("  G", 3, 'A');
    chain.add_residue(g1);

    std::string seq = chain.sequence();
    EXPECT_EQ(seq, "CG"); // Unknown residue skipped
}

TEST_F(ChainTest, Nucleotides) {
    auto nts = chain_a_.nucleotides();
    EXPECT_EQ(nts.size(), 3);
    for (const auto& nt : nts) {
        EXPECT_TRUE(nt.is_nucleotide());
    }
}

TEST_F(ChainTest, NucleotidesWithMixed) {
    Chain chain('A');
    Residue c1("  C", 1, 'A');
    chain.add_residue(c1);

    Residue unknown("XXX", 2, 'A');
    chain.add_residue(unknown);

    Residue g1("  G", 3, 'A');
    chain.add_residue(g1);

    auto nts = chain.nucleotides();
    EXPECT_EQ(nts.size(), 2); // Only C and G
}

// JSON serialization tests - Legacy format
TEST_F(ChainTest, ToJsonLegacy) {
    auto json = chain_a_.to_json_legacy();

    EXPECT_EQ(json["chain_id"], "A");
    EXPECT_EQ(json["num_residues"], 3);
    EXPECT_TRUE(json.contains("residues"));
    EXPECT_TRUE(json["residues"].is_array());
    EXPECT_EQ(json["residues"].size(), 3);
}

TEST_F(ChainTest, FromJsonLegacy) {
    nlohmann::json j = {
        {"chain_id", "B"},
        {"num_residues", 2},
        {"residues",
         {{{"residue_name", "  C"}, {"residue_seq", 1}, {"chain_id", "B"}, {"atoms", {}}},
          {{"residue_name", "  G"}, {"residue_seq", 2}, {"chain_id", "B"}, {"atoms", {}}}}}};

    Chain chain = Chain::from_json_legacy(j);

    EXPECT_EQ(chain.chain_id(), 'B');
    EXPECT_EQ(chain.num_residues(), 2);
}

TEST_F(ChainTest, JsonLegacyRoundTrip) {
    auto json = chain_a_.to_json_legacy();
    Chain chain = Chain::from_json_legacy(json);

    EXPECT_EQ(chain.chain_id(), chain_a_.chain_id());
    EXPECT_EQ(chain.num_residues(), chain_a_.num_residues());
    EXPECT_EQ(chain.sequence(), chain_a_.sequence());
}

// JSON serialization tests - Modern format
TEST_F(ChainTest, ToJsonModern) {
    auto json = chain_a_.to_json();

    EXPECT_EQ(json["chain_id"], "A");
    EXPECT_TRUE(json.contains("residues"));
    EXPECT_TRUE(json["residues"].is_array());
}

TEST_F(ChainTest, FromJsonModern) {
    nlohmann::json j = {
        {"chain_id", "C"},
        {"residues", {{{"name", "  A"}, {"seq_num", 1}, {"chain_id", "C"}, {"atoms", {}}}}}};

    Chain chain = Chain::from_json(j);

    EXPECT_EQ(chain.chain_id(), 'C');
    EXPECT_EQ(chain.num_residues(), 1);
}

TEST_F(ChainTest, JsonModernRoundTrip) {
    auto json = chain_a_.to_json();
    Chain chain = Chain::from_json(json);

    EXPECT_EQ(chain.chain_id(), chain_a_.chain_id());
    EXPECT_EQ(chain.num_residues(), chain_a_.num_residues());
    EXPECT_EQ(chain.sequence(), chain_a_.sequence());
}

// Edge cases
TEST_F(ChainTest, EmptyChain) {
    Chain chain('A');
    EXPECT_EQ(chain.sequence(), "");
    EXPECT_TRUE(chain.nucleotides().empty());
    EXPECT_EQ(chain.num_atoms(), 0);
}

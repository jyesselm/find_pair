/**
 * @file test_residue.cpp
 * @brief Unit tests for Residue class
 */

#include <gtest/gtest.h>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>

using namespace x3dna::core;
using namespace x3dna::geometry;

class ResidueTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a cytosine residue with atoms
        residue_c_ = Residue("  C", 1, 'A');
        residue_c_.add_atom(Atom(" C1'", Vector3D(1.0, 2.0, 3.0), "  C", 'A', 1));
        residue_c_.add_atom(Atom(" N1 ", Vector3D(2.0, 3.0, 4.0), "  C", 'A', 1));
        residue_c_.add_atom(Atom(" C2 ", Vector3D(3.0, 4.0, 5.0), "  C", 'A', 1));
        residue_c_.add_atom(Atom(" N3 ", Vector3D(4.0, 5.0, 6.0), "  C", 'A', 1));

        // Create a guanine residue
        residue_g_ = Residue("  G", 2, 'A');
        residue_g_.add_atom(Atom(" C1'", Vector3D(5.0, 6.0, 7.0), "  G", 'A', 2));
        residue_g_.add_atom(Atom(" N9 ", Vector3D(6.0, 7.0, 8.0), "  G", 'A', 2));
        residue_g_.add_atom(Atom(" C4 ", Vector3D(7.0, 8.0, 9.0), "  G", 'A', 2));
    }

    Residue residue_c_;
    Residue residue_g_;
};

// Constructor tests
TEST_F(ResidueTest, DefaultConstructor) {
    Residue residue;
    EXPECT_EQ(residue.name(), "");
    EXPECT_EQ(residue.seq_num(), 0);
    EXPECT_EQ(residue.chain_id(), '\0');
    EXPECT_EQ(residue.num_atoms(), 0);
}

TEST_F(ResidueTest, NameSeqNumChainConstructor) {
    Residue residue("  A", 10, 'B');
    EXPECT_EQ(residue.name(), "  A");
    EXPECT_EQ(residue.seq_num(), 10);
    EXPECT_EQ(residue.chain_id(), 'B');
}

// Atom management tests
TEST_F(ResidueTest, AddAtom) {
    Residue residue("  C", 1, 'A');
    EXPECT_EQ(residue.num_atoms(), 0);

    residue.add_atom(Atom(" C1'", Vector3D(1, 2, 3)));
    EXPECT_EQ(residue.num_atoms(), 1);

    residue.add_atom(Atom(" N1 ", Vector3D(2, 3, 4)));
    EXPECT_EQ(residue.num_atoms(), 2);
}

TEST_F(ResidueTest, FindAtom) {
    auto atom = residue_c_.find_atom(" N1 ");
    ASSERT_TRUE(atom.has_value());
    EXPECT_EQ(atom->name(), " N1 ");
    EXPECT_EQ(atom->position(), Vector3D(2.0, 3.0, 4.0));
}

TEST_F(ResidueTest, FindAtomNotFound) {
    auto atom = residue_c_.find_atom(" P  ");
    EXPECT_FALSE(atom.has_value());
}

TEST_F(ResidueTest, RingAtoms) {
    auto ring = residue_c_.ring_atoms();
    // N1, C2, N3 are ring atoms
    EXPECT_GT(ring.size(), 0);
    for (const auto& atom : ring) {
        EXPECT_TRUE(atom.is_ring_atom());
    }
}

// Base identification tests
TEST_F(ResidueTest, OneLetterCode) {
    EXPECT_EQ(Residue("  C", 1, 'A').one_letter_code(), 'C');
    EXPECT_EQ(Residue("  G", 1, 'A').one_letter_code(), 'G');
    EXPECT_EQ(Residue("  A", 1, 'A').one_letter_code(), 'A');
    EXPECT_EQ(Residue("  T", 1, 'A').one_letter_code(), 'T');
    EXPECT_EQ(Residue("  U", 1, 'A').one_letter_code(), 'U');
    EXPECT_EQ(Residue("XXX", 1, 'A').one_letter_code(), '?');
}

TEST_F(ResidueTest, IsNucleotide) {
    EXPECT_TRUE(Residue("  C", 1, 'A').is_nucleotide());
    EXPECT_TRUE(Residue("  G", 1, 'A').is_nucleotide());
    EXPECT_TRUE(Residue("  A", 1, 'A').is_nucleotide());
    EXPECT_TRUE(Residue("  T", 1, 'A').is_nucleotide());
    EXPECT_TRUE(Residue("  U", 1, 'A').is_nucleotide());
    EXPECT_FALSE(Residue("XXX", 1, 'A').is_nucleotide());
}

TEST_F(ResidueTest, RYClassification) {
    EXPECT_EQ(Residue("  A", 1, 'A').ry_classification(), 1);  // Purine
    EXPECT_EQ(Residue("  G", 1, 'A').ry_classification(), 1);  // Purine
    EXPECT_EQ(Residue("  C", 1, 'A').ry_classification(), 0);  // Pyrimidine
    EXPECT_EQ(Residue("  T", 1, 'A').ry_classification(), 0);  // Pyrimidine
    EXPECT_EQ(Residue("  U", 1, 'A').ry_classification(), 0);  // Pyrimidine
    EXPECT_EQ(Residue("XXX", 1, 'A').ry_classification(), -1); // Not nucleotide
}

TEST_F(ResidueTest, ResidueType) {
    EXPECT_EQ(Residue("  A", 1, 'A').residue_type(), ResidueType::ADENINE);
    EXPECT_EQ(Residue("  C", 1, 'A').residue_type(), ResidueType::CYTOSINE);
    EXPECT_EQ(Residue("  G", 1, 'A').residue_type(), ResidueType::GUANINE);
    EXPECT_EQ(Residue("  T", 1, 'A').residue_type(), ResidueType::THYMINE);
    EXPECT_EQ(Residue("  U", 1, 'A').residue_type(), ResidueType::URACIL);
}

// Reference frame tests
TEST_F(ResidueTest, ReferenceFrame) {
    Residue residue("  C", 1, 'A');
    EXPECT_FALSE(residue.reference_frame().has_value());

    Matrix3D rotation = Matrix3D::identity();
    Vector3D origin(1.0, 2.0, 3.0);
    ReferenceFrame frame(rotation, origin);

    residue.set_reference_frame(frame);
    ASSERT_TRUE(residue.reference_frame().has_value());

    auto retrieved_frame = residue.reference_frame().value();
    EXPECT_EQ(retrieved_frame.origin(), origin);
}

// JSON serialization tests - Legacy format
TEST_F(ResidueTest, ToJsonLegacy) {
    auto json = residue_c_.to_json_legacy();

    EXPECT_EQ(json["residue_name"], "  C");
    EXPECT_EQ(json["residue_seq"], 1);
    EXPECT_EQ(json["chain_id"], "A");
    EXPECT_TRUE(json.contains("atoms"));
    EXPECT_TRUE(json["atoms"].is_array());
    EXPECT_EQ(json["atoms"].size(), 4);
}

TEST_F(ResidueTest, FromJsonLegacy) {
    nlohmann::json j = {{"residue_name", "  G"},
                        {"residue_seq", 2},
                        {"chain_id", "B"},
                        {"atoms",
                         {{{"atom_name", " C1'"}, {"xyz", {1.0, 2.0, 3.0}}},
                          {{"atom_name", " N9 "}, {"xyz", {2.0, 3.0, 4.0}}}}}};

    Residue residue = Residue::from_json_legacy(j);

    EXPECT_EQ(residue.name(), "  G");
    EXPECT_EQ(residue.seq_num(), 2);
    EXPECT_EQ(residue.chain_id(), 'B');
    EXPECT_EQ(residue.num_atoms(), 2);
}

TEST_F(ResidueTest, JsonLegacyRoundTrip) {
    auto json = residue_c_.to_json_legacy();
    Residue residue = Residue::from_json_legacy(json);

    EXPECT_EQ(residue.name(), residue_c_.name());
    EXPECT_EQ(residue.seq_num(), residue_c_.seq_num());
    EXPECT_EQ(residue.chain_id(), residue_c_.chain_id());
    EXPECT_EQ(residue.num_atoms(), residue_c_.num_atoms());
}

// JSON serialization tests - Modern format
TEST_F(ResidueTest, ToJsonModern) {
    auto json = residue_c_.to_json();

    EXPECT_EQ(json["name"], "  C");
    EXPECT_EQ(json["seq_num"], 1);
    EXPECT_EQ(json["chain_id"], "A");
    EXPECT_TRUE(json.contains("atoms"));
    EXPECT_TRUE(json["atoms"].is_array());
}

TEST_F(ResidueTest, FromJsonModern) {
    nlohmann::json j = {{"name", "  A"},
                        {"seq_num", 3},
                        {"chain_id", "C"},
                        {"atoms",
                         {{{"name", " C1'"}, {"position", {1.0, 2.0, 3.0}}},
                          {{"name", " N9 "}, {"position", {2.0, 3.0, 4.0}}}}}};

    Residue residue = Residue::from_json(j);

    EXPECT_EQ(residue.name(), "  A");
    EXPECT_EQ(residue.seq_num(), 3);
    EXPECT_EQ(residue.chain_id(), 'C');
    EXPECT_EQ(residue.num_atoms(), 2);
}

TEST_F(ResidueTest, JsonModernRoundTrip) {
    auto json = residue_c_.to_json();
    Residue residue = Residue::from_json(json);

    EXPECT_EQ(residue.name(), residue_c_.name());
    EXPECT_EQ(residue.seq_num(), residue_c_.seq_num());
    EXPECT_EQ(residue.chain_id(), residue_c_.chain_id());
    EXPECT_EQ(residue.num_atoms(), residue_c_.num_atoms());
}

// Edge cases
TEST_F(ResidueTest, EmptyResidue) {
    Residue residue("  C", 1, 'A');
    EXPECT_EQ(residue.num_atoms(), 0);
    EXPECT_FALSE(residue.find_atom(" C1'").has_value());
    EXPECT_TRUE(residue.ring_atoms().empty());
}

TEST_F(ResidueTest, ResidueWithReferenceFrame) {
    Residue residue("  C", 1, 'A');
    residue.add_atom(Atom(" C1'", Vector3D(0, 0, 0)));

    Matrix3D rotation = Matrix3D::rotation_z(M_PI / 4.0);
    Vector3D origin(10.0, 20.0, 30.0);
    ReferenceFrame frame(rotation, origin);

    residue.set_reference_frame(frame);

    auto json = residue.to_json_legacy();
    EXPECT_TRUE(json.contains("reference_frame"));

    Residue reconstructed = Residue::from_json_legacy(json);
    ASSERT_TRUE(reconstructed.reference_frame().has_value());
    EXPECT_EQ(reconstructed.reference_frame()->origin(), origin);
}

/**
 * @file test_residue.cpp
 * @brief Unit tests for Residue class
 */

#include <gtest/gtest.h>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/nucleotide_utils.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>
#include <x3dna/io/serializers.hpp>

using namespace x3dna::core;
using namespace x3dna::geometry;
using namespace x3dna::io;

class ResidueTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a cytosine residue with atoms
        residue_c_ = Residue("  C", 1, "A");
        residue_c_.add_atom(Atom(" C1'", Vector3D(1.0, 2.0, 3.0)));
        residue_c_.add_atom(Atom(" N1 ", Vector3D(2.0, 3.0, 4.0)));
        residue_c_.add_atom(Atom(" C2 ", Vector3D(3.0, 4.0, 5.0)));
        residue_c_.add_atom(Atom(" N3 ", Vector3D(4.0, 5.0, 6.0)));

        // Create a guanine residue
        residue_g_ = Residue("  G", 2, "A");
        residue_g_.add_atom(Atom(" C1'", Vector3D(5.0, 6.0, 7.0)));
        residue_g_.add_atom(Atom(" N9 ", Vector3D(6.0, 7.0, 8.0)));
        residue_g_.add_atom(Atom(" C4 ", Vector3D(7.0, 8.0, 9.0)));
    }

    Residue residue_c_;
    Residue residue_g_;
};

// Constructor tests
TEST_F(ResidueTest, DefaultConstructor) {
    Residue residue;
    EXPECT_EQ(residue.name(), "");
    EXPECT_EQ(residue.seq_num(), 0);
    EXPECT_EQ(residue.chain_id(), "");
    EXPECT_EQ(residue.num_atoms(), 0);
}

TEST_F(ResidueTest, NameSeqNumChainConstructor) {
    Residue residue("  A", 10, "B");
    EXPECT_EQ(residue.name(), "A"); // Names are trimmed on construction
    EXPECT_EQ(residue.seq_num(), 10);
    EXPECT_EQ(residue.chain_id(), "B");
}

// Atom management tests
TEST_F(ResidueTest, AddAtom) {
    Residue residue("  C", 1, "A");
    EXPECT_EQ(residue.num_atoms(), 0);

    residue.add_atom(Atom(" C1'", Vector3D(1, 2, 3)));
    EXPECT_EQ(residue.num_atoms(), 1);

    residue.add_atom(Atom(" N1 ", Vector3D(2, 3, 4)));
    EXPECT_EQ(residue.num_atoms(), 2);
}

TEST_F(ResidueTest, FindAtom) {
    // find_atom() accepts both padded and trimmed names
    auto atom = residue_c_.find_atom(" N1 ");
    ASSERT_TRUE(atom.has_value());
    EXPECT_EQ(atom->name(), "N1"); // name() returns trimmed
    EXPECT_EQ(atom->position(), Vector3D(2.0, 3.0, 4.0));
}

TEST_F(ResidueTest, FindAtomNotFound) {
    auto atom = residue_c_.find_atom(" P  ");
    EXPECT_FALSE(atom.has_value());
}

TEST_F(ResidueTest, RingAtoms) {
    auto ring = ring_atoms(residue_c_);
    // N1, C2, N3 are ring atoms
    EXPECT_GT(ring.size(), 0);
    for (const auto& atom : ring) {
        EXPECT_TRUE(atom.is_ring_atom());
    }
}

// Base identification tests (use create_from_atoms for proper initialization)
TEST_F(ResidueTest, OneLetterCode) {
    // one_letter_code() is now a free function that looks up via ModifiedNucleotideRegistry
    EXPECT_EQ(one_letter_code(Residue::create_from_atoms("  C", 1, "A", "", {})), 'C');
    EXPECT_EQ(one_letter_code(Residue::create_from_atoms("  G", 1, "A", "", {})), 'G');
    EXPECT_EQ(one_letter_code(Residue::create_from_atoms("  A", 1, "A", "", {})), 'A');
    EXPECT_EQ(one_letter_code(Residue::create_from_atoms("  T", 1, "A", "", {})), 'T');
    EXPECT_EQ(one_letter_code(Residue::create_from_atoms("  U", 1, "A", "", {})), 'U');
    EXPECT_EQ(one_letter_code(Residue::create_from_atoms("XXX", 1, "A", "", {})), '?');
}

TEST_F(ResidueTest, IsNucleotide) {
    EXPECT_TRUE(Residue::create_from_atoms("  C", 1, "A", "", {}).is_nucleotide());
    EXPECT_TRUE(Residue::create_from_atoms("  G", 1, "A", "", {}).is_nucleotide());
    EXPECT_TRUE(Residue::create_from_atoms("  A", 1, "A", "", {}).is_nucleotide());
    EXPECT_TRUE(Residue::create_from_atoms("  T", 1, "A", "", {}).is_nucleotide());
    EXPECT_TRUE(Residue::create_from_atoms("  U", 1, "A", "", {}).is_nucleotide());
    EXPECT_FALSE(Residue::create_from_atoms("XXX", 1, "A", "", {}).is_nucleotide());
}

TEST_F(ResidueTest, RYClassification) {
    EXPECT_EQ(ry_classification(Residue::create_from_atoms("  A", 1, "A", "", {})), 1);  // Purine
    EXPECT_EQ(ry_classification(Residue::create_from_atoms("  G", 1, "A", "", {})), 1);  // Purine
    EXPECT_EQ(ry_classification(Residue::create_from_atoms("  C", 1, "A", "", {})), 0);  // Pyrimidine
    EXPECT_EQ(ry_classification(Residue::create_from_atoms("  T", 1, "A", "", {})), 0);  // Pyrimidine
    EXPECT_EQ(ry_classification(Residue::create_from_atoms("  U", 1, "A", "", {})), 0);  // Pyrimidine
    EXPECT_EQ(ry_classification(Residue::create_from_atoms("XXX", 1, "A", "", {})), -1); // Not nucleotide
}

TEST_F(ResidueTest, BaseType) {
    EXPECT_EQ(Residue::create_from_atoms("  A", 1, "A", "", {}).base_type(), typing::BaseType::ADENINE);
    EXPECT_EQ(Residue::create_from_atoms("  C", 1, "A", "", {}).base_type(), typing::BaseType::CYTOSINE);
    EXPECT_EQ(Residue::create_from_atoms("  G", 1, "A", "", {}).base_type(), typing::BaseType::GUANINE);
    EXPECT_EQ(Residue::create_from_atoms("  T", 1, "A", "", {}).base_type(), typing::BaseType::THYMINE);
    EXPECT_EQ(Residue::create_from_atoms("  U", 1, "A", "", {}).base_type(), typing::BaseType::URACIL);
}

// Reference frame tests
TEST_F(ResidueTest, ReferenceFrame) {
    Residue residue("  C", 1, "A");
    EXPECT_FALSE(residue.reference_frame().has_value());

    Matrix3D rotation = Matrix3D::identity();
    Vector3D origin(1.0, 2.0, 3.0);
    ReferenceFrame frame(rotation, origin);

    residue.set_reference_frame(frame);
    ASSERT_TRUE(residue.reference_frame().has_value());

    auto retrieved_frame = residue.reference_frame().value();
    EXPECT_EQ(retrieved_frame.origin(), origin);
}

// JSON serialization tests - using ResidueSerializer
TEST_F(ResidueTest, ToJsonLegacy) {
    auto json = ResidueSerializer::to_legacy_json(residue_c_);

    EXPECT_EQ(json["residue_name"], "C"); // Names are trimmed
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

    Residue residue = ResidueSerializer::from_legacy_json(j);

    EXPECT_EQ(residue.name(), "G"); // Names are trimmed on construction
    EXPECT_EQ(residue.seq_num(), 2);
    EXPECT_EQ(residue.chain_id(), "B");
    EXPECT_EQ(residue.num_atoms(), 2);
}

TEST_F(ResidueTest, JsonLegacyRoundTrip) {
    auto json = ResidueSerializer::to_legacy_json(residue_c_);
    Residue residue = ResidueSerializer::from_legacy_json(json);

    EXPECT_EQ(residue.name(), residue_c_.name());
    EXPECT_EQ(residue.seq_num(), residue_c_.seq_num());
    EXPECT_EQ(residue.chain_id(), residue_c_.chain_id());
    EXPECT_EQ(residue.num_atoms(), residue_c_.num_atoms());
}

// JSON serialization tests - Modern format
TEST_F(ResidueTest, ToJsonModern) {
    auto json = ResidueSerializer::to_json(residue_c_);

    EXPECT_EQ(json["name"], "C"); // Names are trimmed
    EXPECT_EQ(json["seq_num"], 1);
    EXPECT_EQ(json["chain_id"], "A");
    EXPECT_TRUE(json.contains("atoms"));
    EXPECT_TRUE(json["atoms"].is_array());
}

TEST_F(ResidueTest, FromJsonModern) {
    // Residue::from_json uses Atom::from_json which expects legacy keys
    nlohmann::json j = {{"name", "  A"},
                        {"seq_num", 3},
                        {"chain_id", "C"},
                        {"atoms",
                         {{{"atom_name", " C1'"}, {"xyz", {1.0, 2.0, 3.0}}},
                          {{"atom_name", " N9 "}, {"xyz", {2.0, 3.0, 4.0}}}}}};

    Residue residue = ResidueSerializer::from_json(j);

    EXPECT_EQ(residue.name(), "A"); // Names are trimmed on construction
    EXPECT_EQ(residue.seq_num(), 3);
    EXPECT_EQ(residue.chain_id(), "C");
    EXPECT_EQ(residue.num_atoms(), 2);
}

TEST_F(ResidueTest, JsonModernRoundTrip) {
    auto json = ResidueSerializer::to_json(residue_c_);
    Residue residue = ResidueSerializer::from_json(json);

    EXPECT_EQ(residue.name(), residue_c_.name());
    EXPECT_EQ(residue.seq_num(), residue_c_.seq_num());
    EXPECT_EQ(residue.chain_id(), residue_c_.chain_id());
    EXPECT_EQ(residue.num_atoms(), residue_c_.num_atoms());
}

// Edge cases
TEST_F(ResidueTest, EmptyResidue) {
    Residue residue("  C", 1, "A");
    EXPECT_EQ(residue.num_atoms(), 0);
    EXPECT_FALSE(residue.find_atom(" C1'").has_value());
    EXPECT_TRUE(ring_atoms(residue).empty());
}

TEST_F(ResidueTest, ResidueWithReferenceFrame) {
    Residue residue("  C", 1, "A");
    residue.add_atom(Atom(" C1'", Vector3D(0, 0, 0)));

    Matrix3D rotation = Matrix3D::rotation_z(M_PI / 4.0);
    Vector3D origin(10.0, 20.0, 30.0);
    ReferenceFrame frame(rotation, origin);

    residue.set_reference_frame(frame);

    auto json = ResidueSerializer::to_legacy_json(residue);
    EXPECT_TRUE(json.contains("reference_frame"));

    Residue reconstructed = ResidueSerializer::from_legacy_json(json);
    ASSERT_TRUE(reconstructed.reference_frame().has_value());
    EXPECT_EQ(reconstructed.reference_frame()->origin(), origin);
}

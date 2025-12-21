/**
 * @file test_atom.cpp
 * @brief Unit tests for Atom class
 */

#include <gtest/gtest.h>
#include <x3dna/core/atom.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/io/serializers.hpp>

using namespace x3dna::core;
using namespace x3dna::geometry;
using namespace x3dna::io;

class AtomTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test atoms - Atom names are trimmed on construction
        atom1_ = Atom(" C1'", Vector3D(1.0, 2.0, 3.0));
        atom2_ = Atom(" N3 ", Vector3D(4.0, 5.0, 6.0));
        atom3_ = Atom(" O2 ", Vector3D(0.0, 0.0, 0.0));

        // Create test residues for serialization tests
        residue1_ = Residue("C", 1, "A");
        residue2_ = Residue("G", 2, "B");
    }

    Atom atom1_;
    Atom atom2_;
    Atom atom3_;
    Residue residue1_;
    Residue residue2_;
};

// Constructor tests
TEST_F(AtomTest, DefaultConstructor) {
    Atom atom;
    EXPECT_EQ(atom.name(), "");
    EXPECT_EQ(atom.position(), Vector3D(0, 0, 0));
}

TEST_F(AtomTest, NamePositionConstructor) {
    Atom atom(" C1'", Vector3D(1.0, 2.0, 3.0));
    // Names are trimmed on construction
    EXPECT_EQ(atom.name(), "C1'");
    EXPECT_EQ(atom.position(), Vector3D(1.0, 2.0, 3.0));
}

TEST_F(AtomTest, FullConstructor) {
    // Names are trimmed on construction
    EXPECT_EQ(atom1_.name(), "C1'");
    EXPECT_EQ(atom1_.position(), Vector3D(1.0, 2.0, 3.0));
}

// Builder pattern tests
TEST_F(AtomTest, BuilderBasic) {
    auto atom = Atom::create(" N1 ", Vector3D(5.0, 6.0, 7.0)).build();

    // Names are trimmed on construction
    EXPECT_EQ(atom.name(), "N1");
    EXPECT_EQ(atom.position(), Vector3D(5.0, 6.0, 7.0));
}

TEST_F(AtomTest, BuilderAllFields) {
    auto atom = Atom::create(" CA ", Vector3D(1.0, 2.0, 3.0))
                    .alt_loc('A')
                    .occupancy(0.75)
                    .atom_serial(100)
                    .model_number(1)
                    .b_factor(25.5)
                    .element("C")
                    .legacy_atom_idx(50)
                    .build();

    EXPECT_EQ(atom.name(), "CA");
    EXPECT_EQ(atom.alt_loc(), 'A');
    EXPECT_NEAR(atom.occupancy(), 0.75, 1e-9);
    EXPECT_EQ(atom.atom_serial(), 100);
    EXPECT_EQ(atom.model_number(), 1);
    EXPECT_NEAR(atom.b_factor(), 25.5, 1e-9);
    EXPECT_EQ(atom.element(), "C");
    EXPECT_EQ(atom.legacy_atom_idx(), 50);
}

// Post-construction setters (for parsing workflow)
TEST_F(AtomTest, PostConstructionSetters) {
    auto atom = Atom::create(" C1'", Vector3D(1.0, 2.0, 3.0)).build();

    // These are the only setters retained for parsing workflow
    atom.set_model_number(2);
    atom.set_legacy_atom_idx(100);

    EXPECT_EQ(atom.model_number(), 2);
    EXPECT_EQ(atom.legacy_atom_idx(), 100);
}

// Distance calculation tests
TEST_F(AtomTest, DistanceTo) {
    double dist = atom1_.distance_to(atom2_);
    double expected = Vector3D(1.0, 2.0, 3.0).distance_to(Vector3D(4.0, 5.0, 6.0));
    EXPECT_NEAR(dist, expected, 1e-9);
    EXPECT_NEAR(dist, std::sqrt(27.0), 1e-9); // sqrt((4-1)^2 + (5-2)^2 + (6-3)^2)
}

TEST_F(AtomTest, DistanceToSelf) {
    double dist = atom1_.distance_to(atom1_);
    EXPECT_NEAR(dist, 0.0, 1e-9);
}

// Ring atom tests - names are now trimmed
TEST_F(AtomTest, IsRingAtom) {
    EXPECT_TRUE(Atom(" N1 ", Vector3D(0, 0, 0)).is_ring_atom());
    EXPECT_TRUE(Atom(" C2 ", Vector3D(0, 0, 0)).is_ring_atom());
    EXPECT_TRUE(Atom(" N3 ", Vector3D(0, 0, 0)).is_ring_atom());
    EXPECT_TRUE(Atom(" C4 ", Vector3D(0, 0, 0)).is_ring_atom());
    EXPECT_TRUE(Atom(" C5 ", Vector3D(0, 0, 0)).is_ring_atom());
    EXPECT_TRUE(Atom(" C6 ", Vector3D(0, 0, 0)).is_ring_atom());
    EXPECT_TRUE(Atom(" N7 ", Vector3D(0, 0, 0)).is_ring_atom());
    EXPECT_TRUE(Atom(" C8 ", Vector3D(0, 0, 0)).is_ring_atom());
    EXPECT_TRUE(Atom(" N9 ", Vector3D(0, 0, 0)).is_ring_atom());

    EXPECT_FALSE(Atom(" C1'", Vector3D(0, 0, 0)).is_ring_atom());
    EXPECT_FALSE(Atom(" O2 ", Vector3D(0, 0, 0)).is_ring_atom());
    EXPECT_FALSE(Atom(" P  ", Vector3D(0, 0, 0)).is_ring_atom());
}

// H-bond donor/acceptor tests
TEST_F(AtomTest, IsHydrogenBondDonor) {
    EXPECT_TRUE(Atom(" N1 ", Vector3D(0, 0, 0)).is_hydrogen_bond_donor());
    EXPECT_TRUE(Atom(" N2 ", Vector3D(0, 0, 0)).is_hydrogen_bond_donor());
    EXPECT_TRUE(Atom(" N3 ", Vector3D(0, 0, 0)).is_hydrogen_bond_donor());
    EXPECT_TRUE(Atom(" N4 ", Vector3D(0, 0, 0)).is_hydrogen_bond_donor());
    EXPECT_TRUE(Atom(" N6 ", Vector3D(0, 0, 0)).is_hydrogen_bond_donor());
    EXPECT_TRUE(Atom(" N7 ", Vector3D(0, 0, 0)).is_hydrogen_bond_donor());
    EXPECT_TRUE(Atom(" N9 ", Vector3D(0, 0, 0)).is_hydrogen_bond_donor());

    EXPECT_FALSE(Atom(" O2 ", Vector3D(0, 0, 0)).is_hydrogen_bond_donor());
    EXPECT_FALSE(Atom(" C1'", Vector3D(0, 0, 0)).is_hydrogen_bond_donor());
}

TEST_F(AtomTest, IsHydrogenBondAcceptor) {
    EXPECT_TRUE(Atom(" O2 ", Vector3D(0, 0, 0)).is_hydrogen_bond_acceptor());
    EXPECT_TRUE(Atom(" O4 ", Vector3D(0, 0, 0)).is_hydrogen_bond_acceptor());
    EXPECT_TRUE(Atom(" O6 ", Vector3D(0, 0, 0)).is_hydrogen_bond_acceptor());
    EXPECT_TRUE(Atom(" N3 ", Vector3D(0, 0, 0)).is_hydrogen_bond_acceptor());
    EXPECT_TRUE(Atom(" N7 ", Vector3D(0, 0, 0)).is_hydrogen_bond_acceptor());

    EXPECT_FALSE(Atom(" N1 ", Vector3D(0, 0, 0)).is_hydrogen_bond_acceptor());
    EXPECT_FALSE(Atom(" C1'", Vector3D(0, 0, 0)).is_hydrogen_bond_acceptor());
}

// JSON serialization tests - using AtomSerializer with Residue context
TEST_F(AtomTest, ToJsonLegacy) {
    // AtomSerializer now requires Residue context for residue-level fields
    auto json = AtomSerializer::to_legacy_json(atom1_, residue1_);

    // Atom names are trimmed, so JSON output contains trimmed names
    EXPECT_EQ(json["atom_name"], "C1'");
    std::vector<double> expected_xyz = {1.0, 2.0, 3.0};
    EXPECT_EQ(json["xyz"], expected_xyz);
    // Residue-level fields come from residue1_
    EXPECT_EQ(json["residue_name"], "C");
    EXPECT_EQ(json["chain_id"], "A");
    EXPECT_EQ(json["residue_seq"], 1);
    EXPECT_EQ(json["record_type"], "A");
}

TEST_F(AtomTest, FromJsonLegacy) {
    nlohmann::json j = {{"atom_name", " N3 "}, {"xyz", {4.0, 5.0, 6.0}}, {"residue_name", "  G"},
                        {"chain_id", "B"},     {"residue_seq", 2},       {"record_type", "A"}};

    Atom atom = AtomSerializer::from_legacy_json(j);

    // Names are trimmed on construction, residue fields are ignored in deserialization
    EXPECT_EQ(atom.name(), "N3");
    EXPECT_EQ(atom.position(), Vector3D(4.0, 5.0, 6.0));
}

TEST_F(AtomTest, JsonLegacyRoundTrip) {
    auto json = AtomSerializer::to_legacy_json(atom1_, residue1_);
    Atom atom = AtomSerializer::from_legacy_json(json);

    // Names are trimmed
    EXPECT_EQ(atom.name(), atom1_.name());
    EXPECT_EQ(atom.position(), atom1_.position());
}

// JSON serialization tests - Modern format
TEST_F(AtomTest, ToJsonModern) {
    auto json = AtomSerializer::to_json(atom1_, residue1_);

    // Atom names are trimmed, so JSON output contains trimmed names
    EXPECT_EQ(json["atom_name"], "C1'");
    std::vector<double> expected_xyz = {1.0, 2.0, 3.0};
    EXPECT_EQ(json["xyz"], expected_xyz);
    // Residue-level fields come from residue1_
    EXPECT_EQ(json["residue_name"], "C");
    EXPECT_EQ(json["chain_id"], "A");
    EXPECT_EQ(json["residue_seq"], 1);
    EXPECT_EQ(json["record_type"], "A");
}

TEST_F(AtomTest, FromJsonModern) {
    nlohmann::json j = {{"atom_name", " N3 "}, {"xyz", {4.0, 5.0, 6.0}}, {"residue_name", "  G"},
                        {"chain_id", "B"},     {"residue_seq", 2},       {"record_type", "A"}};

    Atom atom = AtomSerializer::from_json(j);

    EXPECT_EQ(atom.name(), "N3");
    EXPECT_EQ(atom.position(), Vector3D(4.0, 5.0, 6.0));
}

TEST_F(AtomTest, JsonModernRoundTrip) {
    auto json = AtomSerializer::to_json(atom1_, residue1_);
    Atom atom = AtomSerializer::from_json(json);

    EXPECT_EQ(atom.name(), atom1_.name());
    EXPECT_EQ(atom.position(), atom1_.position());
}

// Edge cases
TEST_F(AtomTest, MinimalJsonLegacy) {
    nlohmann::json j = {{"atom_name", " C1'"}, {"xyz", {1.0, 2.0, 3.0}}};

    Atom atom = AtomSerializer::from_legacy_json(j);
    EXPECT_EQ(atom.name(), "C1'");
    EXPECT_EQ(atom.position(), Vector3D(1.0, 2.0, 3.0));
}

TEST_F(AtomTest, AtomNameWithSpaces) {
    Atom atom("  P  ", Vector3D(0, 0, 0));
    // Name is trimmed
    EXPECT_EQ(atom.name(), "P");
    // Ring atom check uses trimmed name
    EXPECT_FALSE(atom.is_ring_atom());
}

// Test that trimmed names work for comparisons
TEST_F(AtomTest, TrimmedNameComparisons) {
    Atom atom1(" N1 ", Vector3D(0, 0, 0));
    Atom atom2("N1", Vector3D(0, 0, 0));

    // Both should have the same trimmed name
    EXPECT_EQ(atom1.name(), atom2.name());
    EXPECT_EQ(atom1.name(), "N1");
}

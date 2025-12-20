/**
 * @file test_atom.cpp
 * @brief Unit tests for Atom class
 */

#include <gtest/gtest.h>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/io/serializers.hpp>

using namespace x3dna::core;
using namespace x3dna::geometry;
using namespace x3dna::io;

class AtomTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test atoms - names are now trimmed internally
        atom1_ = Atom(" C1'", Vector3D(1.0, 2.0, 3.0), "  C", "A", 1);
        atom2_ = Atom(" N3 ", Vector3D(4.0, 5.0, 6.0), "  G", "A", 2);
        atom3_ = Atom(" O2 ", Vector3D(0.0, 0.0, 0.0));
    }

    Atom atom1_;
    Atom atom2_;
    Atom atom3_;
};

// Constructor tests
TEST_F(AtomTest, DefaultConstructor) {
    Atom atom;
    EXPECT_EQ(atom.name(), "");
    EXPECT_EQ(atom.position(), Vector3D(0, 0, 0));
    EXPECT_EQ(atom.residue_name(), "");
    EXPECT_EQ(atom.chain_id(), "");
    EXPECT_EQ(atom.residue_seq(), 0);
    EXPECT_EQ(atom.record_type(), 'A');
}

TEST_F(AtomTest, NamePositionConstructor) {
    Atom atom(" C1'", Vector3D(1.0, 2.0, 3.0));
    // Names are now trimmed
    EXPECT_EQ(atom.name(), "C1'");
    EXPECT_EQ(atom.original_atom_name(), " C1'");
    EXPECT_EQ(atom.position(), Vector3D(1.0, 2.0, 3.0));
}

TEST_F(AtomTest, FullConstructor) {
    // Names are now trimmed, original names preserved
    EXPECT_EQ(atom1_.name(), "C1'");
    EXPECT_EQ(atom1_.original_atom_name(), " C1'");
    EXPECT_EQ(atom1_.position(), Vector3D(1.0, 2.0, 3.0));
    EXPECT_EQ(atom1_.residue_name(), "C");
    EXPECT_EQ(atom1_.original_residue_name(), "  C");
    EXPECT_EQ(atom1_.chain_id(), "A");
    EXPECT_EQ(atom1_.residue_seq(), 1);
    EXPECT_EQ(atom1_.record_type(), 'A');
}

// Builder pattern tests
TEST_F(AtomTest, BuilderBasic) {
    auto atom = Atom::create(" N1 ", Vector3D(5.0, 6.0, 7.0))
                    .residue_name("  A")
                    .chain_id("B")
                    .residue_seq(10)
                    .record_type('H')
                    .build();

    // Names are trimmed
    EXPECT_EQ(atom.name(), "N1");
    EXPECT_EQ(atom.original_atom_name(), " N1 ");
    EXPECT_EQ(atom.position(), Vector3D(5.0, 6.0, 7.0));
    EXPECT_EQ(atom.residue_name(), "A");
    EXPECT_EQ(atom.chain_id(), "B");
    EXPECT_EQ(atom.residue_seq(), 10);
    EXPECT_EQ(atom.record_type(), 'H');
}

TEST_F(AtomTest, BuilderAllFields) {
    auto atom = Atom::create(" CA ", Vector3D(1.0, 2.0, 3.0))
                    .residue_name("ALA")
                    .chain_id("A")
                    .residue_seq(42)
                    .record_type('A')
                    .alt_loc('A')
                    .insertion("B")
                    .occupancy(0.75)
                    .atom_serial(100)
                    .model_number(1)
                    .b_factor(25.5)
                    .element("C")
                    .original_atom_name(" CA ")
                    .original_residue_name("ALA")
                    .legacy_atom_idx(50)
                    .legacy_residue_idx(10)
                    .build();

    EXPECT_EQ(atom.name(), "CA");
    EXPECT_EQ(atom.residue_name(), "ALA");
    EXPECT_EQ(atom.chain_id(), "A");
    EXPECT_EQ(atom.residue_seq(), 42);
    EXPECT_EQ(atom.alt_loc(), 'A');
    EXPECT_EQ(atom.insertion(), "B");
    EXPECT_NEAR(atom.occupancy(), 0.75, 1e-9);
    EXPECT_EQ(atom.atom_serial(), 100);
    EXPECT_EQ(atom.model_number(), 1);
    EXPECT_NEAR(atom.b_factor(), 25.5, 1e-9);
    EXPECT_EQ(atom.element(), "C");
    EXPECT_EQ(atom.legacy_atom_idx(), 50);
    EXPECT_EQ(atom.legacy_residue_idx(), 10);
}

// Post-construction setters (for parsing workflow)
TEST_F(AtomTest, PostConstructionSetters) {
    auto atom = Atom::create(" C1'", Vector3D(1.0, 2.0, 3.0)).build();

    // These are the only setters retained for parsing workflow
    atom.set_model_number(2);
    atom.set_legacy_atom_idx(100);
    atom.set_legacy_residue_idx(25);

    EXPECT_EQ(atom.model_number(), 2);
    EXPECT_EQ(atom.legacy_atom_idx(), 100);
    EXPECT_EQ(atom.legacy_residue_idx(), 25);
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

// JSON serialization tests - using AtomSerializer
TEST_F(AtomTest, ToJsonLegacy) {
    auto json = AtomSerializer::to_legacy_json(atom1_);

    // Original name is preserved in JSON output
    EXPECT_EQ(json["atom_name"], " C1'");
    std::vector<double> expected_xyz = {1.0, 2.0, 3.0};
    EXPECT_EQ(json["xyz"], expected_xyz);
    EXPECT_EQ(json["residue_name"], "  C");
    EXPECT_EQ(json["chain_id"], "A");
    EXPECT_EQ(json["residue_seq"], 1);
    EXPECT_EQ(json["record_type"], "A");
}

TEST_F(AtomTest, FromJsonLegacy) {
    nlohmann::json j = {{"atom_name", " N3 "}, {"xyz", {4.0, 5.0, 6.0}}, {"residue_name", "  G"},
                        {"chain_id", "B"},     {"residue_seq", 2},       {"record_type", "A"}};

    Atom atom = AtomSerializer::from_legacy_json(j);

    // Names are trimmed on construction
    EXPECT_EQ(atom.name(), "N3");
    EXPECT_EQ(atom.position(), Vector3D(4.0, 5.0, 6.0));
    EXPECT_EQ(atom.residue_name(), "G");
    EXPECT_EQ(atom.chain_id(), "B");
    EXPECT_EQ(atom.residue_seq(), 2);
    EXPECT_EQ(atom.record_type(), 'A');
}

TEST_F(AtomTest, JsonLegacyRoundTrip) {
    auto json = AtomSerializer::to_legacy_json(atom1_);
    Atom atom = AtomSerializer::from_legacy_json(json);

    // Names are trimmed, so compare trimmed versions
    EXPECT_EQ(atom.name(), atom1_.name());
    EXPECT_EQ(atom.position(), atom1_.position());
    EXPECT_EQ(atom.residue_name(), atom1_.residue_name());
    EXPECT_EQ(atom.chain_id(), atom1_.chain_id());
    EXPECT_EQ(atom.residue_seq(), atom1_.residue_seq());
    EXPECT_EQ(atom.record_type(), atom1_.record_type());
}

// JSON serialization tests - Modern format
TEST_F(AtomTest, ToJsonModern) {
    auto json = AtomSerializer::to_json(atom1_);

    // to_json() uses original names for JSON output
    EXPECT_EQ(json["atom_name"], " C1'");
    std::vector<double> expected_xyz = {1.0, 2.0, 3.0};
    EXPECT_EQ(json["xyz"], expected_xyz);
    EXPECT_EQ(json["residue_name"], "  C");
    EXPECT_EQ(json["chain_id"], "A");
    EXPECT_EQ(json["residue_seq"], 1);
    EXPECT_EQ(json["record_type"], "A");
}

TEST_F(AtomTest, FromJsonModern) {
    nlohmann::json j = {{"atom_name", " N3 "}, {"xyz", {4.0, 5.0, 6.0}}, {"residue_name", "  G"}, {"chain_id", "B"},
                        {"residue_seq", 2}, {"record_type", "A"}};

    Atom atom = AtomSerializer::from_json(j);

    EXPECT_EQ(atom.name(), "N3");
    EXPECT_EQ(atom.position(), Vector3D(4.0, 5.0, 6.0));
    EXPECT_EQ(atom.residue_name(), "G");
    EXPECT_EQ(atom.chain_id(), "B");
    EXPECT_EQ(atom.residue_seq(), 2);
    EXPECT_EQ(atom.record_type(), 'A');
}

TEST_F(AtomTest, JsonModernRoundTrip) {
    auto json = AtomSerializer::to_json(atom1_);
    Atom atom = AtomSerializer::from_json(json);

    EXPECT_EQ(atom.name(), atom1_.name());
    EXPECT_EQ(atom.position(), atom1_.position());
    EXPECT_EQ(atom.residue_name(), atom1_.residue_name());
    EXPECT_EQ(atom.chain_id(), atom1_.chain_id());
    EXPECT_EQ(atom.residue_seq(), atom1_.residue_seq());
    EXPECT_EQ(atom.record_type(), atom1_.record_type());
}

// Edge cases
TEST_F(AtomTest, MinimalJsonLegacy) {
    nlohmann::json j = {{"atom_name", " C1'"}, {"xyz", {1.0, 2.0, 3.0}}};

    Atom atom = AtomSerializer::from_legacy_json(j);
    EXPECT_EQ(atom.name(), "C1'");
    EXPECT_EQ(atom.position(), Vector3D(1.0, 2.0, 3.0));
    EXPECT_EQ(atom.residue_name(), "");
    EXPECT_EQ(atom.chain_id(), "");
    EXPECT_EQ(atom.residue_seq(), 0);
}

TEST_F(AtomTest, AtomNameWithSpaces) {
    Atom atom("  P  ", Vector3D(0, 0, 0));
    // Name is trimmed
    EXPECT_EQ(atom.name(), "P");
    EXPECT_EQ(atom.original_atom_name(), "  P  ");
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

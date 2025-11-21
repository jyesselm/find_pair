/**
 * @file test_atom.cpp
 * @brief Unit tests for Atom class
 */

#include <gtest/gtest.h>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>

using namespace x3dna::core;
using namespace x3dna::geometry;

class AtomTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test atoms
        atom1_ = Atom(" C1'", Vector3D(1.0, 2.0, 3.0), "  C", 'A', 1);
        atom2_ = Atom(" N3 ", Vector3D(4.0, 5.0, 6.0), "  G", 'A', 2);
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
    EXPECT_EQ(atom.chain_id(), '\0');
    EXPECT_EQ(atom.residue_seq(), 0);
    EXPECT_EQ(atom.record_type(), 'A');
}

TEST_F(AtomTest, NamePositionConstructor) {
    Atom atom(" C1'", Vector3D(1.0, 2.0, 3.0));
    EXPECT_EQ(atom.name(), " C1'");
    EXPECT_EQ(atom.position(), Vector3D(1.0, 2.0, 3.0));
}

TEST_F(AtomTest, FullConstructor) {
    EXPECT_EQ(atom1_.name(), " C1'");
    EXPECT_EQ(atom1_.position(), Vector3D(1.0, 2.0, 3.0));
    EXPECT_EQ(atom1_.residue_name(), "  C");
    EXPECT_EQ(atom1_.chain_id(), 'A');
    EXPECT_EQ(atom1_.residue_seq(), 1);
    EXPECT_EQ(atom1_.record_type(), 'A');
}

// Getter/Setter tests
TEST_F(AtomTest, Setters) {
    Atom atom;
    atom.set_name(" N1 ");
    atom.set_position(Vector3D(5.0, 6.0, 7.0));
    atom.set_residue_name("  A");
    atom.set_chain_id('B');
    atom.set_residue_seq(10);
    atom.set_record_type('H');
    
    EXPECT_EQ(atom.name(), " N1 ");
    EXPECT_EQ(atom.position(), Vector3D(5.0, 6.0, 7.0));
    EXPECT_EQ(atom.residue_name(), "  A");
    EXPECT_EQ(atom.chain_id(), 'B');
    EXPECT_EQ(atom.residue_seq(), 10);
    EXPECT_EQ(atom.record_type(), 'H');
}

// Distance calculation tests
TEST_F(AtomTest, DistanceTo) {
    double dist = atom1_.distance_to(atom2_);
    double expected = Vector3D(1.0, 2.0, 3.0).distance_to(Vector3D(4.0, 5.0, 6.0));
    EXPECT_NEAR(dist, expected, 1e-9);
    EXPECT_NEAR(dist, std::sqrt(27.0), 1e-9);  // sqrt((4-1)^2 + (5-2)^2 + (6-3)^2)
}

TEST_F(AtomTest, DistanceToSelf) {
    double dist = atom1_.distance_to(atom1_);
    EXPECT_NEAR(dist, 0.0, 1e-9);
}

// Ring atom tests
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

// JSON serialization tests - Legacy format
TEST_F(AtomTest, ToJsonLegacy) {
    auto json = atom1_.to_json_legacy();
    
    EXPECT_EQ(json["atom_name"], " C1'");
    std::vector<double> expected_xyz = {1.0, 2.0, 3.0};
    EXPECT_EQ(json["xyz"], expected_xyz);
    EXPECT_EQ(json["residue_name"], "  C");
    EXPECT_EQ(json["chain_id"], "A");
    EXPECT_EQ(json["residue_seq"], 1);
    EXPECT_EQ(json["record_type"], "A");
}

TEST_F(AtomTest, FromJsonLegacy) {
    nlohmann::json j = {
        {"atom_name", " N3 "},
        {"xyz", {4.0, 5.0, 6.0}},
        {"residue_name", "  G"},
        {"chain_id", "B"},
        {"residue_seq", 2},
        {"record_type", "A"}
    };
    
    Atom atom = Atom::from_json_legacy(j);
    
    EXPECT_EQ(atom.name(), " N3 ");
    EXPECT_EQ(atom.position(), Vector3D(4.0, 5.0, 6.0));
    EXPECT_EQ(atom.residue_name(), "  G");
    EXPECT_EQ(atom.chain_id(), 'B');
    EXPECT_EQ(atom.residue_seq(), 2);
    EXPECT_EQ(atom.record_type(), 'A');
}

TEST_F(AtomTest, JsonLegacyRoundTrip) {
    auto json = atom1_.to_json_legacy();
    Atom atom = Atom::from_json_legacy(json);
    
    EXPECT_EQ(atom.name(), atom1_.name());
    EXPECT_EQ(atom.position(), atom1_.position());
    EXPECT_EQ(atom.residue_name(), atom1_.residue_name());
    EXPECT_EQ(atom.chain_id(), atom1_.chain_id());
    EXPECT_EQ(atom.residue_seq(), atom1_.residue_seq());
    EXPECT_EQ(atom.record_type(), atom1_.record_type());
}

// JSON serialization tests - Modern format
TEST_F(AtomTest, ToJsonModern) {
    auto json = atom1_.to_json();
    
    EXPECT_EQ(json["name"], " C1'");
    std::vector<double> expected_pos = {1.0, 2.0, 3.0};
    EXPECT_EQ(json["position"], expected_pos);
    EXPECT_EQ(json["residue_name"], "  C");
    EXPECT_EQ(json["chain_id"], "A");
    EXPECT_EQ(json["residue_seq"], 1);
    EXPECT_EQ(json["record_type"], "A");
}

TEST_F(AtomTest, FromJsonModern) {
    nlohmann::json j = {
        {"name", " N3 "},
        {"position", {4.0, 5.0, 6.0}},
        {"residue_name", "  G"},
        {"chain_id", "B"},
        {"residue_seq", 2},
        {"record_type", "A"}
    };
    
    Atom atom = Atom::from_json(j);
    
    EXPECT_EQ(atom.name(), " N3 ");
    EXPECT_EQ(atom.position(), Vector3D(4.0, 5.0, 6.0));
    EXPECT_EQ(atom.residue_name(), "  G");
    EXPECT_EQ(atom.chain_id(), 'B');
    EXPECT_EQ(atom.residue_seq(), 2);
    EXPECT_EQ(atom.record_type(), 'A');
}

TEST_F(AtomTest, JsonModernRoundTrip) {
    auto json = atom1_.to_json();
    Atom atom = Atom::from_json(json);
    
    EXPECT_EQ(atom.name(), atom1_.name());
    EXPECT_EQ(atom.position(), atom1_.position());
    EXPECT_EQ(atom.residue_name(), atom1_.residue_name());
    EXPECT_EQ(atom.chain_id(), atom1_.chain_id());
    EXPECT_EQ(atom.residue_seq(), atom1_.residue_seq());
    EXPECT_EQ(atom.record_type(), atom1_.record_type());
}

// Edge cases
TEST_F(AtomTest, MinimalJsonLegacy) {
    nlohmann::json j = {
        {"atom_name", " C1'"},
        {"xyz", {1.0, 2.0, 3.0}}
    };
    
    Atom atom = Atom::from_json_legacy(j);
    EXPECT_EQ(atom.name(), " C1'");
    EXPECT_EQ(atom.position(), Vector3D(1.0, 2.0, 3.0));
    EXPECT_EQ(atom.residue_name(), "");
    EXPECT_EQ(atom.chain_id(), '\0');
    EXPECT_EQ(atom.residue_seq(), 0);
}

TEST_F(AtomTest, AtomNameWithSpaces) {
    Atom atom("  P  ", Vector3D(0, 0, 0));
    EXPECT_EQ(atom.name(), "  P  ");
    // Ring atom check should handle spaces
    EXPECT_FALSE(atom.is_ring_atom());
}


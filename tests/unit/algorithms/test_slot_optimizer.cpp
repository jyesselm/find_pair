/**
 * @file test_slot_optimizer.cpp
 * @brief Unit tests for slot-based H-bond optimizer
 */

#include <gtest/gtest.h>
#include <x3dna/algorithms/hydrogen_bond/slot/slot.hpp>
#include <x3dna/algorithms/hydrogen_bond/slot/atom_capacity.hpp>
#include <x3dna/algorithms/hydrogen_bond/slot/hbond_candidate.hpp>
#include <x3dna/algorithms/hydrogen_bond/slot/slot_optimizer_params.hpp>
#include <cmath>

using namespace x3dna::algorithms::hydrogen_bond::slot;
using namespace x3dna::geometry;

// ============================================================================
// AtomCapacity tests
// ============================================================================

TEST(AtomCapacityTest, DonorCapacityStandardBases) {
    // NH2 amino groups - 2 hydrogens
    EXPECT_EQ(AtomCapacity::get_donor_capacity("A", "N6"), 2);
    EXPECT_EQ(AtomCapacity::get_donor_capacity("C", "N4"), 2);
    EXPECT_EQ(AtomCapacity::get_donor_capacity("G", "N2"), 2);

    // Imino NH - 1 hydrogen
    EXPECT_EQ(AtomCapacity::get_donor_capacity("G", "N1"), 1);
    EXPECT_EQ(AtomCapacity::get_donor_capacity("U", "N3"), 1);
    EXPECT_EQ(AtomCapacity::get_donor_capacity("T", "N3"), 1);

    // Non-donors should return 0
    EXPECT_EQ(AtomCapacity::get_donor_capacity("A", "N1"), 0);
    EXPECT_EQ(AtomCapacity::get_donor_capacity("G", "O6"), 0);
}

TEST(AtomCapacityTest, AcceptorCapacityStandardBases) {
    // sp2 carbonyl oxygens - 2 lone pairs
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("G", "O6"), 2);
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("U", "O2"), 2);
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("U", "O4"), 2);
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("C", "O2"), 2);

    // sp2 ring nitrogens - 1 lone pair
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("A", "N1"), 1);
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("A", "N3"), 1);
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("A", "N7"), 1);
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("G", "N7"), 1);

    // Non-acceptors should return 0
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("A", "N6"), 0);
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("G", "N2"), 0);
}

TEST(AtomCapacityTest, RiboseAtoms) {
    // O2' can both donate (1 H) and accept (2 LP)
    EXPECT_EQ(AtomCapacity::get_donor_capacity("A", "O2'"), 1);
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("A", "O2'"), 2);

    // O4' ring - 1 accessible LP
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("G", "O4'"), 1);
}

TEST(AtomCapacityTest, BackboneAtoms) {
    // Phosphate oxygens - 3 lone pairs (OP1/O1P variants)
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("A", "OP1"), 3);
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("A", "O1P"), 3);
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("G", "OP2"), 3);
    EXPECT_EQ(AtomCapacity::get_acceptor_capacity("G", "O2P"), 3);
}

TEST(AtomCapacityTest, ParentBaseTypeFallback) {
    // DNA variants
    EXPECT_EQ(AtomCapacity::get_donor_capacity("DA", "N6"), 2);
    EXPECT_EQ(AtomCapacity::get_donor_capacity("DG", "N2"), 2);
    EXPECT_EQ(AtomCapacity::get_donor_capacity("DC", "N4"), 2);
    EXPECT_EQ(AtomCapacity::get_donor_capacity("DT", "N3"), 1);

    // Modified bases using parent lookup
    EXPECT_EQ(AtomCapacity::get_donor_capacity("5MC", "N4"), 2); // 5-methylcytosine
}

TEST(AtomCapacityTest, NormalizeAtomName) {
    EXPECT_EQ(AtomCapacity::normalize_atom_name("  N6  "), "N6");
    EXPECT_EQ(AtomCapacity::normalize_atom_name("O2'"), "O2'");
    EXPECT_EQ(AtomCapacity::normalize_atom_name("\tOP1\t"), "OP1");
}

TEST(AtomCapacityTest, IsBackboneAtom) {
    EXPECT_TRUE(AtomCapacity::is_backbone_atom("P"));
    EXPECT_TRUE(AtomCapacity::is_backbone_atom("OP1"));
    EXPECT_TRUE(AtomCapacity::is_backbone_atom("OP2"));
    EXPECT_TRUE(AtomCapacity::is_backbone_atom("O1P"));
    EXPECT_TRUE(AtomCapacity::is_backbone_atom("O2P"));
    EXPECT_TRUE(AtomCapacity::is_backbone_atom("O3'"));
    EXPECT_TRUE(AtomCapacity::is_backbone_atom("O5'"));

    EXPECT_FALSE(AtomCapacity::is_backbone_atom("N1"));
    EXPECT_FALSE(AtomCapacity::is_backbone_atom("O6"));
    EXPECT_FALSE(AtomCapacity::is_backbone_atom("O2'"));
    EXPECT_FALSE(AtomCapacity::is_backbone_atom("C1'"));
}

// ============================================================================
// HSlot tests
// ============================================================================

TEST(HSlotTest, BasicAvailability) {
    Vector3D direction(1.0, 0.0, 0.0);
    HSlot slot(direction, 1);  // Single-bond slot

    EXPECT_TRUE(slot.is_available());
    EXPECT_EQ(slot.bond_count(), 0u);
    EXPECT_EQ(slot.max_bonds(), 1);

    slot.add_bond(direction);
    EXPECT_FALSE(slot.is_available());
    EXPECT_EQ(slot.bond_count(), 1u);
}

TEST(HSlotTest, MultipleBonds) {
    Vector3D direction(0.0, 1.0, 0.0);
    HSlot slot(direction, 2);  // Two-bond slot (like NH2)

    EXPECT_TRUE(slot.is_available());  // Completely unused

    Vector3D first_bond(0.5, 0.866, 0.0);  // 60 degrees off
    slot.add_bond(first_bond);
    EXPECT_FALSE(slot.is_available());  // No longer "available" (has bonds)
    EXPECT_EQ(slot.bond_count(), 1u);   // But only 1 bond so far

    // Can still add a second bond if angularly separated
    Vector3D second_bond(-0.5, 0.866, 0.0);  // ~90 degrees from first
    EXPECT_TRUE(slot.can_add_bond(second_bond, 60.0));  // Far enough apart
    slot.add_bond(second_bond);
    EXPECT_EQ(slot.bond_count(), 2u);  // Now at max capacity
}

TEST(HSlotTest, BifurcationAngleCheck) {
    Vector3D direction(0.0, 0.0, 1.0);
    HSlot slot(direction, 2);

    // First bond pointing mostly along z with slight y
    Vector3D first_bond(0.0, 0.5, 0.866);
    slot.add_bond(first_bond);

    // Too close - angle between (0, 0.5, 0.866) and (0, 0.4, 0.9) is only ~7 degrees
    Vector3D too_close(0.0, 0.4, 0.9165);  // Nearly same direction
    too_close = too_close.normalized();
    EXPECT_FALSE(slot.can_add_bond(too_close, 60.0));

    // Far enough apart - pointing mostly in -y direction (about 90 degrees from first)
    Vector3D far_enough(0.0, -0.866, 0.5);  // dot product = 0.5*(-0.866) + 0.866*0.5 = 0 -> 90 degrees
    EXPECT_TRUE(slot.can_add_bond(far_enough, 60.0));
}

TEST(HSlotTest, Reset) {
    Vector3D direction(1.0, 0.0, 0.0);
    HSlot slot(direction, 1);

    slot.add_bond(direction);
    EXPECT_FALSE(slot.is_available());

    slot.reset();
    EXPECT_TRUE(slot.is_available());
    EXPECT_EQ(slot.bond_count(), 0u);
}

// ============================================================================
// LPSlot tests
// ============================================================================

TEST(LPSlotTest, BasicAvailability) {
    Vector3D direction(0.0, 1.0, 0.0);
    LPSlot slot(direction, 1);  // Single LP

    EXPECT_TRUE(slot.is_available());
    EXPECT_EQ(slot.bond_count(), 0u);

    slot.add_bond(direction);
    EXPECT_FALSE(slot.is_available());
}

TEST(LPSlotTest, MultipleLonePairs) {
    Vector3D direction(-1.0, 0.0, 0.0);
    LPSlot slot(direction, 2);  // sp2 oxygen has 2 LPs

    EXPECT_TRUE(slot.is_available());  // Initially unused

    Vector3D first_bond(-0.866, 0.5, 0.0);
    slot.add_bond(first_bond);
    EXPECT_FALSE(slot.is_available());  // No longer "available" (has bonds)
    EXPECT_EQ(slot.bond_count(), 1u);

    // Can still add second bond if angularly separated (~60 degrees apart)
    Vector3D second_bond(-0.866, -0.5, 0.0);
    EXPECT_TRUE(slot.can_add_bond(second_bond, 43.0));  // These are ~60 degrees apart
    slot.add_bond(second_bond);
    EXPECT_EQ(slot.bond_count(), 2u);  // Now at max
}

// ============================================================================
// HBondCandidate tests
// ============================================================================

TEST(HBondCandidateTest, Direction) {
    HBondCandidate c;
    c.donor_pos = Vector3D(0.0, 0.0, 0.0);
    c.acceptor_pos = Vector3D(3.0, 0.0, 0.0);
    c.distance = 3.0;

    Vector3D dir = c.direction();
    EXPECT_NEAR(dir.x(), 1.0, 1e-6);
    EXPECT_NEAR(dir.y(), 0.0, 1e-6);
    EXPECT_NEAR(dir.z(), 0.0, 1e-6);
}

TEST(HBondCandidateTest, QualityScore) {
    HBondCandidate c1, c2;

    // Shorter distance, lower alignment
    c1.distance = 2.8;
    c1.alignment_score = 1.0;

    // Longer distance, higher alignment
    c2.distance = 3.0;
    c2.alignment_score = 1.5;

    // c1 should score better: -2.8 + 0.4*1.0 = -2.4
    // c2: -3.0 + 0.4*1.5 = -2.4
    // They're equal in this case
    EXPECT_NEAR(c1.quality_score(), c2.quality_score(), 1e-6);

    // Now make c2's alignment even better
    c2.alignment_score = 2.0;
    // c2: -3.0 + 0.4*2.0 = -2.2 (better than -2.4)
    EXPECT_GT(c2.quality_score(), c1.quality_score());
}

// ============================================================================
// SlotOptimizerParams tests
// ============================================================================

TEST(SlotOptimizerParamsTest, DefaultParams) {
    auto params = SlotOptimizerParams::optimized();

    EXPECT_EQ(params.max_distance, 4.0);
    EXPECT_EQ(params.short_distance_threshold, 3.5);
    EXPECT_EQ(params.min_alignment, 0.3);
    EXPECT_EQ(params.min_bifurcation_alignment, 0.5);
    EXPECT_EQ(params.min_bifurcation_angle, 43.0);
    EXPECT_FALSE(params.baseline_mode);
}

TEST(SlotOptimizerParamsTest, BaselineParams) {
    auto params = SlotOptimizerParams::baseline();

    EXPECT_TRUE(params.baseline_mode);
    EXPECT_EQ(params.baseline_min_distance, 2.5);
    EXPECT_EQ(params.baseline_max_distance, 3.5);
}

TEST(SlotOptimizerParamsTest, StrictParams) {
    auto params = SlotOptimizerParams::strict();

    EXPECT_EQ(params.min_alignment, 0.5);
    EXPECT_EQ(params.min_bifurcation_alignment, 0.7);
    EXPECT_FALSE(params.baseline_mode);
}

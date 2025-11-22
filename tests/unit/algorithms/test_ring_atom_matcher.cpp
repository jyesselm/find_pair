/**
 * @file test_ring_atom_matcher.cpp
 * @brief Unit tests for RingAtomMatcher class
 */

#include <gtest/gtest.h>
#include <x3dna/algorithms/ring_atom_matcher.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>

using namespace x3dna::algorithms;
using namespace x3dna::core;
using namespace x3dna::geometry;

class RingAtomMatcherTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a test experimental residue (adenine)
        experimental_residue_ = Residue("  A", 1, 'A');
        
        // Add ring atoms for adenine (purine): 9 atoms
        experimental_residue_.add_atom(Atom(" C4 ", Vector3D(-1.267, 3.124, 0.000), "  A", 'A', 1));
        experimental_residue_.add_atom(Atom(" N3 ", Vector3D(-2.320, 2.290, 0.000), "  A", 'A', 1));
        experimental_residue_.add_atom(Atom(" C2 ", Vector3D(-1.912, 1.023, 0.000), "  A", 'A', 1));
        experimental_residue_.add_atom(Atom(" N1 ", Vector3D(-0.668, 0.532, 0.000), "  A", 'A', 1));
        experimental_residue_.add_atom(Atom(" C6 ", Vector3D(0.369, 1.398, 0.000), "  A", 'A', 1));
        experimental_residue_.add_atom(Atom(" C5 ", Vector3D(0.071, 2.771, 0.000), "  A", 'A', 1));
        experimental_residue_.add_atom(Atom(" N7 ", Vector3D(0.877, 3.902, 0.000), "  A", 'A', 1));
        experimental_residue_.add_atom(Atom(" C8 ", Vector3D(0.024, 4.897, 0.000), "  A", 'A', 1));
        experimental_residue_.add_atom(Atom(" N9 ", Vector3D(-1.291, 4.498, 0.000), "  A", 'A', 1));
        
        // Create standard template structure (same atoms, different coordinates)
        standard_template_ = Structure("ATOMIC_A");
        Chain chain('A');
        Residue template_residue("  A", 1, 'A');
        
        // Standard template coordinates (from Atomic_A.pdb format)
        template_residue.add_atom(Atom(" C1'", Vector3D(-2.479, 5.346, 0.000), "  A", 'A', 1));
        template_residue.add_atom(Atom(" N9 ", Vector3D(-1.291, 4.498, 0.000), "  A", 'A', 1));
        template_residue.add_atom(Atom(" C8 ", Vector3D(0.024, 4.897, 0.000), "  A", 'A', 1));
        template_residue.add_atom(Atom(" N7 ", Vector3D(0.877, 3.902, 0.000), "  A", 'A', 1));
        template_residue.add_atom(Atom(" C5 ", Vector3D(0.071, 2.771, 0.000), "  A", 'A', 1));
        template_residue.add_atom(Atom(" C6 ", Vector3D(0.369, 1.398, 0.000), "  A", 'A', 1));
        template_residue.add_atom(Atom(" N6 ", Vector3D(1.611, 0.909, 0.000), "  A", 'A', 1));
        template_residue.add_atom(Atom(" N1 ", Vector3D(-0.668, 0.532, 0.000), "  A", 'A', 1));
        template_residue.add_atom(Atom(" C2 ", Vector3D(-1.912, 1.023, 0.000), "  A", 'A', 1));
        template_residue.add_atom(Atom(" N3 ", Vector3D(-2.320, 2.290, 0.000), "  A", 'A', 1));
        template_residue.add_atom(Atom(" C4 ", Vector3D(-1.267, 3.124, 0.000), "  A", 'A', 1));
        
        chain.add_residue(template_residue);
        standard_template_.add_chain(chain);
    }

    Residue experimental_residue_;
    Structure standard_template_;
};

// Test matching for purine (adenine)
TEST_F(RingAtomMatcherTest, MatchPurineAtoms) {
    MatchedAtoms matched = RingAtomMatcher::match(experimental_residue_, standard_template_, false);
    
    EXPECT_GE(matched.num_matched, 9);  // Should match all 9 purine ring atoms
    EXPECT_TRUE(matched.is_valid());
    EXPECT_EQ(matched.experimental.size(), matched.standard.size());
    EXPECT_EQ(matched.experimental.size(), matched.atom_names.size());
}

// Test matching for pyrimidine (cytosine)
TEST_F(RingAtomMatcherTest, MatchPyrimidineAtoms) {
    Residue cytosine("  C", 1, 'A');
    
    // Add pyrimidine ring atoms: 6 atoms
    cytosine.add_atom(Atom(" C4 ", Vector3D(0.0, 0.0, 0.0), "  C", 'A', 1));
    cytosine.add_atom(Atom(" N3 ", Vector3D(1.0, 0.0, 0.0), "  C", 'A', 1));
    cytosine.add_atom(Atom(" C2 ", Vector3D(2.0, 0.0, 0.0), "  C", 'A', 1));
    cytosine.add_atom(Atom(" N1 ", Vector3D(3.0, 0.0, 0.0), "  C", 'A', 1));
    cytosine.add_atom(Atom(" C6 ", Vector3D(4.0, 0.0, 0.0), "  C", 'A', 1));
    cytosine.add_atom(Atom(" C5 ", Vector3D(5.0, 0.0, 0.0), "  C", 'A', 1));
    
    Structure template_c("ATOMIC_C");
    Chain chain('A');
    Residue template_residue("  C", 1, 'A');
    template_residue.add_atom(Atom(" C4 ", Vector3D(0.0, 0.0, 0.0), "  C", 'A', 1));
    template_residue.add_atom(Atom(" N3 ", Vector3D(1.0, 0.0, 0.0), "  C", 'A', 1));
    template_residue.add_atom(Atom(" C2 ", Vector3D(2.0, 0.0, 0.0), "  C", 'A', 1));
    template_residue.add_atom(Atom(" N1 ", Vector3D(3.0, 0.0, 0.0), "  C", 'A', 1));
    template_residue.add_atom(Atom(" C6 ", Vector3D(4.0, 0.0, 0.0), "  C", 'A', 1));
    template_residue.add_atom(Atom(" C5 ", Vector3D(5.0, 0.0, 0.0), "  C", 'A', 1));
    chain.add_residue(template_residue);
    template_c.add_chain(chain);
    
    MatchedAtoms matched = RingAtomMatcher::match(cytosine, template_c, false);
    
    EXPECT_GE(matched.num_matched, 6);  // Should match all 6 pyrimidine ring atoms
    EXPECT_TRUE(matched.is_valid());
}

// Test RNA matching (includes C1')
TEST_F(RingAtomMatcherTest, MatchRNAAtoms) {
    Residue rna_residue("  A", 1, 'A');
    rna_residue.add_atom(Atom(" C1'", Vector3D(0.0, 0.0, 0.0), "  A", 'A', 1));
    rna_residue.add_atom(Atom(" C4 ", Vector3D(1.0, 0.0, 0.0), "  A", 'A', 1));
    rna_residue.add_atom(Atom(" N3 ", Vector3D(2.0, 0.0, 0.0), "  A", 'A', 1));
    
    Structure template_rna("ATOMIC_A_RNA");
    Chain chain('A');
    Residue template_residue("  A", 1, 'A');
    template_residue.add_atom(Atom(" C1'", Vector3D(0.0, 0.0, 0.0), "  A", 'A', 1));
    template_residue.add_atom(Atom(" C4 ", Vector3D(1.0, 0.0, 0.0), "  A", 'A', 1));
    template_residue.add_atom(Atom(" N3 ", Vector3D(2.0, 0.0, 0.0), "  A", 'A', 1));
    chain.add_residue(template_residue);
    template_rna.add_chain(chain);
    
    MatchedAtoms matched = RingAtomMatcher::match(rna_residue, template_rna, true);
    
    // Should include C1' in the match
    EXPECT_GE(matched.num_matched, 3);
    bool has_c1_prime = false;
    for (const auto& name : matched.atom_names) {
        if (name == " C1'") {
            has_c1_prime = true;
            break;
        }
    }
    EXPECT_TRUE(has_c1_prime);
}

// Test with missing atoms
TEST_F(RingAtomMatcherTest, MatchWithMissingAtoms) {
    Residue incomplete_residue("  A", 1, 'A');
    // Only add 5 atoms (missing 4)
    incomplete_residue.add_atom(Atom(" C4 ", Vector3D(0.0, 0.0, 0.0), "  A", 'A', 1));
    incomplete_residue.add_atom(Atom(" N3 ", Vector3D(1.0, 0.0, 0.0), "  A", 'A', 1));
    incomplete_residue.add_atom(Atom(" C2 ", Vector3D(2.0, 0.0, 0.0), "  A", 'A', 1));
    incomplete_residue.add_atom(Atom(" N1 ", Vector3D(3.0, 0.0, 0.0), "  A", 'A', 1));
    incomplete_residue.add_atom(Atom(" C6 ", Vector3D(4.0, 0.0, 0.0), "  A", 'A', 1));
    
    MatchedAtoms matched = RingAtomMatcher::match(incomplete_residue, standard_template_, false);
    
    // Should still match what's available (at least 5)
    EXPECT_GE(matched.num_matched, 5);
    // May or may not be valid depending on threshold
}

// Test ring atom names retrieval
TEST_F(RingAtomMatcherTest, GetRingAtomNames) {
    // Purine
    auto purine_names = RingAtomMatcher::get_ring_atom_names(ResidueType::ADENINE, false);
    EXPECT_EQ(purine_names.size(), 9);
    
    // Pyrimidine
    auto pyrimidine_names = RingAtomMatcher::get_ring_atom_names(ResidueType::CYTOSINE, false);
    EXPECT_EQ(pyrimidine_names.size(), 6);
    
    // RNA purine (includes C1')
    auto rna_purine_names = RingAtomMatcher::get_ring_atom_names(ResidueType::ADENINE, true);
    EXPECT_EQ(rna_purine_names.size(), 10);  // 9 ring atoms + C1'
    EXPECT_EQ(rna_purine_names[0], " C1'");
}


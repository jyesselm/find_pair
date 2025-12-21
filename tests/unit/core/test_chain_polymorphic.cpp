/**
 * @file test_chain_polymorphic.cpp
 * @brief Unit tests for the polymorphic Chain class
 */

#include <gtest/gtest.h>
#include <x3dna/core/structure/residue.hpp>
#include <x3dna/geometry/vector3d.hpp>

using namespace x3dna::core::structure;
using namespace x3dna::core;  // For Atom
using namespace x3dna::geometry;

class ChainPolymorphicTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

// === Basic Construction ===

TEST_F(ChainPolymorphicTest, DefaultConstruction) {
    Chain chain;
    EXPECT_TRUE(chain.empty());
    EXPECT_EQ(chain.size(), 0);
    EXPECT_EQ(chain.chain_id(), "");
}

TEST_F(ChainPolymorphicTest, ConstructWithId) {
    Chain chain("A");
    EXPECT_EQ(chain.chain_id(), "A");
    EXPECT_TRUE(chain.empty());
}

// === Adding Residues ===

TEST_F(ChainPolymorphicTest, AddRNAResidue) {
    Chain chain("A");

    chain.add_residue(ResidueFactory::create("A", 1, "A", ""));
    chain.add_residue(ResidueFactory::create("G", 2, "A", ""));
    chain.add_residue(ResidueFactory::create("C", 3, "A", ""));

    EXPECT_EQ(chain.size(), 3);
    EXPECT_TRUE(chain[0].is_rna());
    EXPECT_TRUE(chain[1].is_rna());
    EXPECT_TRUE(chain[2].is_rna());
}

TEST_F(ChainPolymorphicTest, AddMixedResidues) {
    Chain chain("A");

    chain.add_residue(ResidueFactory::create("A", 1, "A", ""));      // RNA
    chain.add_residue(ResidueFactory::create("DA", 2, "A", ""));     // DNA
    chain.add_residue(ResidueFactory::create("ALA", 3, "A", ""));    // Protein
    chain.add_residue(ResidueFactory::create("HOH", 4, "A", ""));    // Water

    EXPECT_EQ(chain.size(), 4);
    EXPECT_TRUE(chain[0].is_rna());
    EXPECT_TRUE(chain[1].is_dna());
    EXPECT_TRUE(chain[2].is_protein());
    EXPECT_TRUE(chain[3].is_ligand());
}

// === Iteration ===

TEST_F(ChainPolymorphicTest, IterateResidues) {
    Chain chain("A");
    chain.add_residue(ResidueFactory::create("A", 1, "A", ""));
    chain.add_residue(ResidueFactory::create("G", 2, "A", ""));
    chain.add_residue(ResidueFactory::create("C", 3, "A", ""));

    int count = 0;
    for (const auto& res : chain) {
        EXPECT_TRUE(res.is_nucleotide());
        count++;
    }
    EXPECT_EQ(count, 3);
}

TEST_F(ChainPolymorphicTest, ConstIteration) {
    Chain chain("A");
    chain.add_residue(ResidueFactory::create("A", 1, "A", ""));
    chain.add_residue(ResidueFactory::create("G", 2, "A", ""));

    const Chain& const_chain = chain;
    int count = 0;
    for (const auto& res : const_chain) {
        EXPECT_TRUE(res.is_nucleotide());
        count++;
    }
    EXPECT_EQ(count, 2);
}

// === Nucleotide Access ===

TEST_F(ChainPolymorphicTest, GetNucleotides) {
    Chain chain("A");
    chain.add_residue(ResidueFactory::create("A", 1, "A", ""));      // RNA
    chain.add_residue(ResidueFactory::create("ALA", 2, "A", ""));    // Protein
    chain.add_residue(ResidueFactory::create("G", 3, "A", ""));      // RNA
    chain.add_residue(ResidueFactory::create("HOH", 4, "A", ""));    // Water

    auto nucs = chain.nucleotides();
    EXPECT_EQ(nucs.size(), 2);
    EXPECT_EQ(nucs[0]->one_letter_code(), 'A');
    EXPECT_EQ(nucs[1]->one_letter_code(), 'G');
}

TEST_F(ChainPolymorphicTest, GetNucleotidesConst) {
    Chain chain("A");
    chain.add_residue(ResidueFactory::create("C", 1, "A", ""));
    chain.add_residue(ResidueFactory::create("U", 2, "A", ""));

    const Chain& const_chain = chain;
    auto nucs = const_chain.nucleotides();
    EXPECT_EQ(nucs.size(), 2);
}

// === Sequence ===

TEST_F(ChainPolymorphicTest, GetSequence) {
    Chain chain("A");
    chain.add_residue(ResidueFactory::create("A", 1, "A", ""));
    chain.add_residue(ResidueFactory::create("G", 2, "A", ""));
    chain.add_residue(ResidueFactory::create("C", 3, "A", ""));
    chain.add_residue(ResidueFactory::create("U", 4, "A", ""));

    EXPECT_EQ(chain.sequence(), "AGCU");
}

TEST_F(ChainPolymorphicTest, SequenceSkipsNonNucleotides) {
    Chain chain("A");
    chain.add_residue(ResidueFactory::create("A", 1, "A", ""));
    chain.add_residue(ResidueFactory::create("ALA", 2, "A", ""));   // Protein - skipped
    chain.add_residue(ResidueFactory::create("G", 3, "A", ""));

    EXPECT_EQ(chain.sequence(), "AG");
}

// === Find Residue ===

TEST_F(ChainPolymorphicTest, FindResidueBySeqNum) {
    Chain chain("A");
    chain.add_residue(ResidueFactory::create("A", 10, "A", ""));
    chain.add_residue(ResidueFactory::create("G", 20, "A", ""));
    chain.add_residue(ResidueFactory::create("C", 30, "A", ""));

    auto* found = chain.find_residue(20);
    ASSERT_NE(found, nullptr);
    EXPECT_EQ(found->name(), "G");

    auto* not_found = chain.find_residue(99);
    EXPECT_EQ(not_found, nullptr);
}

// === Clone ===

TEST_F(ChainPolymorphicTest, CloneChain) {
    Chain original("A");
    original.add_residue(ResidueFactory::create("A", 1, "A", ""));
    original.add_residue(ResidueFactory::create("G", 2, "A", ""));

    Chain cloned = original.clone();

    EXPECT_EQ(cloned.chain_id(), "A");
    EXPECT_EQ(cloned.size(), 2);
    EXPECT_EQ(cloned[0].name(), "A");
    EXPECT_EQ(cloned[1].name(), "G");

    // Verify it's a deep copy (different addresses)
    EXPECT_NE(&cloned[0], &original[0]);
}

// === Move Semantics ===

TEST_F(ChainPolymorphicTest, MoveConstruction) {
    Chain original("A");
    original.add_residue(ResidueFactory::create("A", 1, "A", ""));

    Chain moved(std::move(original));

    EXPECT_EQ(moved.chain_id(), "A");
    EXPECT_EQ(moved.size(), 1);
}

TEST_F(ChainPolymorphicTest, MoveAssignment) {
    Chain original("A");
    original.add_residue(ResidueFactory::create("A", 1, "A", ""));

    Chain target("B");
    target = std::move(original);

    EXPECT_EQ(target.chain_id(), "A");
    EXPECT_EQ(target.size(), 1);
}

// === Dynamic Cast Access ===

TEST_F(ChainPolymorphicTest, DynamicCastToINucleotide) {
    Chain chain("A");
    chain.add_residue(ResidueFactory::create("A", 1, "A", ""));
    chain.add_residue(ResidueFactory::create("ALA", 2, "A", ""));

    // First residue should cast to INucleotide
    auto* nuc = dynamic_cast<INucleotide*>(&chain[0]);
    ASSERT_NE(nuc, nullptr);
    EXPECT_EQ(nuc->one_letter_code(), 'A');
    EXPECT_TRUE(nuc->is_purine());

    // Second residue should NOT cast to INucleotide
    auto* not_nuc = dynamic_cast<INucleotide*>(&chain[1]);
    EXPECT_EQ(not_nuc, nullptr);
}

TEST_F(ChainPolymorphicTest, DynamicCastToRNA) {
    Chain chain("A");
    chain.add_residue(ResidueFactory::create("A", 1, "A", ""));
    chain.add_residue(ResidueFactory::create("DA", 2, "A", ""));

    // First is RNA
    auto* rna = dynamic_cast<RNA*>(&chain[0]);
    ASSERT_NE(rna, nullptr);

    // Second is DNA, not RNA
    auto* not_rna = dynamic_cast<RNA*>(&chain[1]);
    EXPECT_EQ(not_rna, nullptr);

    auto* dna = dynamic_cast<DNA*>(&chain[1]);
    ASSERT_NE(dna, nullptr);
}

// === Atom Count ===

TEST_F(ChainPolymorphicTest, NumAtoms) {
    Chain chain("A");

    auto res1 = ResidueFactory::create("A", 1, "A", "");
    res1->add_atom(Atom("N9", Vector3D(0, 0, 0)));
    res1->add_atom(Atom("C8", Vector3D(1, 0, 0)));

    auto res2 = ResidueFactory::create("G", 2, "A", "");
    res2->add_atom(Atom("N9", Vector3D(2, 0, 0)));

    chain.add_residue(std::move(res1));
    chain.add_residue(std::move(res2));

    EXPECT_EQ(chain.num_atoms(), 3);
}

/**
 * @file test_residue_hierarchy.cpp
 * @brief Unit tests for the polymorphic residue hierarchy
 */

#include <gtest/gtest.h>
#include <x3dna/core/residue/residue.hpp>
#include <x3dna/geometry/matrix3d.hpp>

using namespace x3dna::core::poly;
using namespace x3dna::core;  // For ReferenceFrame, Atom
using namespace x3dna::geometry;

class ResidueHierarchyTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

// === Factory Tests ===

TEST_F(ResidueHierarchyTest, FactoryCreatesRNAForAdenine) {
    auto residue = ResidueFactory::create("A", 1, "A", "");

    ASSERT_NE(residue, nullptr);
    EXPECT_TRUE(residue->is_nucleotide());
    EXPECT_TRUE(residue->is_rna());
    EXPECT_FALSE(residue->is_dna());
    EXPECT_FALSE(residue->is_protein());
    EXPECT_FALSE(residue->is_ligand());
    EXPECT_EQ(residue->name(), "A");
}

TEST_F(ResidueHierarchyTest, FactoryCreatesRNAForGuanine) {
    auto residue = ResidueFactory::create("G", 2, "A", "");

    ASSERT_NE(residue, nullptr);
    EXPECT_TRUE(residue->is_rna());
    EXPECT_EQ(residue->name(), "G");
}

TEST_F(ResidueHierarchyTest, FactoryCreatesDNAForDA) {
    auto residue = ResidueFactory::create("DA", 1, "A", "");

    ASSERT_NE(residue, nullptr);
    EXPECT_TRUE(residue->is_nucleotide());
    EXPECT_TRUE(residue->is_dna());
    EXPECT_FALSE(residue->is_rna());
    EXPECT_EQ(residue->name(), "DA");
}

TEST_F(ResidueHierarchyTest, FactoryCreatesDNAForDT) {
    auto residue = ResidueFactory::create("DT", 1, "B", "");

    ASSERT_NE(residue, nullptr);
    EXPECT_TRUE(residue->is_dna());
    EXPECT_FALSE(residue->is_rna());
}

TEST_F(ResidueHierarchyTest, FactoryCreatesProteinForALA) {
    auto residue = ResidueFactory::create("ALA", 1, "A", "");

    ASSERT_NE(residue, nullptr);
    EXPECT_TRUE(residue->is_protein());
    EXPECT_FALSE(residue->is_nucleotide());
    EXPECT_FALSE(residue->is_ligand());
    EXPECT_EQ(residue->name(), "ALA");
}

TEST_F(ResidueHierarchyTest, FactoryCreatesProteinForGLY) {
    auto residue = ResidueFactory::create("GLY", 5, "A", "");

    ASSERT_NE(residue, nullptr);
    EXPECT_TRUE(residue->is_protein());
}

TEST_F(ResidueHierarchyTest, FactoryCreatesLigandForWater) {
    auto residue = ResidueFactory::create("HOH", 1, "W", "");

    ASSERT_NE(residue, nullptr);
    EXPECT_TRUE(residue->is_ligand());
    EXPECT_FALSE(residue->is_nucleotide());
    EXPECT_FALSE(residue->is_protein());

    // Cast to Ligand to check water-specific method
    auto* ligand = dynamic_cast<Ligand*>(residue.get());
    ASSERT_NE(ligand, nullptr);
    EXPECT_TRUE(ligand->is_water());
}

TEST_F(ResidueHierarchyTest, FactoryCreatesLigandForMagnesiumIon) {
    auto residue = ResidueFactory::create("MG", 1, "A", "");

    ASSERT_NE(residue, nullptr);
    EXPECT_TRUE(residue->is_ligand());

    auto* ligand = dynamic_cast<Ligand*>(residue.get());
    ASSERT_NE(ligand, nullptr);
    EXPECT_TRUE(ligand->is_ion());
}

// === INucleotide Interface Tests ===

TEST_F(ResidueHierarchyTest, RNACastToINucleotide) {
    auto residue = ResidueFactory::create("A", 1, "A", "");

    auto* nucleotide = dynamic_cast<INucleotide*>(residue.get());
    ASSERT_NE(nucleotide, nullptr);
    EXPECT_EQ(nucleotide->one_letter_code(), 'A');
    EXPECT_TRUE(nucleotide->is_purine());
    EXPECT_FALSE(nucleotide->is_pyrimidine());
}

TEST_F(ResidueHierarchyTest, DNACastToINucleotide) {
    auto residue = ResidueFactory::create("DC", 1, "A", "");

    auto* nucleotide = dynamic_cast<INucleotide*>(residue.get());
    ASSERT_NE(nucleotide, nullptr);
    EXPECT_EQ(nucleotide->one_letter_code(), 'C');
    EXPECT_FALSE(nucleotide->is_purine());
    EXPECT_TRUE(nucleotide->is_pyrimidine());
}

TEST_F(ResidueHierarchyTest, RNAPurineClassification) {
    auto a = ResidueFactory::create_rna("A", 1, "A");
    auto g = ResidueFactory::create_rna("G", 2, "A");

    EXPECT_TRUE(a->is_purine());
    EXPECT_TRUE(g->is_purine());
    EXPECT_EQ(a->ry_classification(), 1);
    EXPECT_EQ(g->ry_classification(), 1);
}

TEST_F(ResidueHierarchyTest, RNAPyrimidineClassification) {
    auto c = ResidueFactory::create_rna("C", 1, "A");
    auto u = ResidueFactory::create_rna("U", 2, "A");

    EXPECT_TRUE(c->is_pyrimidine());
    EXPECT_TRUE(u->is_pyrimidine());
    EXPECT_EQ(c->ry_classification(), 0);
    EXPECT_EQ(u->ry_classification(), 0);
}

// === Clone Tests ===

TEST_F(ResidueHierarchyTest, RNAClone) {
    auto original = ResidueFactory::create_rna("A", 1, "A");
    original->set_legacy_residue_idx(42);

    auto cloned = original->clone();

    ASSERT_NE(cloned, nullptr);
    EXPECT_EQ(cloned->name(), "A");
    EXPECT_EQ(cloned->seq_num(), 1);
    EXPECT_EQ(cloned->legacy_residue_idx(), 42);
    EXPECT_TRUE(cloned->is_rna());
}

TEST_F(ResidueHierarchyTest, DNAClone) {
    auto original = ResidueFactory::create_dna("DG", 5, "B");

    auto cloned = original->clone();

    ASSERT_NE(cloned, nullptr);
    EXPECT_EQ(cloned->name(), "DG");
    EXPECT_TRUE(cloned->is_dna());
}

// === Atom Management Tests ===

TEST_F(ResidueHierarchyTest, AddAtomToResidue) {
    auto residue = ResidueFactory::create("A", 1, "A", "");

    EXPECT_EQ(residue->num_atoms(), 0);

    Atom atom("N9", Vector3D(0.0, 0.0, 0.0));
    residue->add_atom(atom);

    EXPECT_EQ(residue->num_atoms(), 1);
    EXPECT_EQ(residue->atoms()[0].name(), "N9");
}

TEST_F(ResidueHierarchyTest, FindAtomInResidue) {
    auto residue = ResidueFactory::create("A", 1, "A", "");

    Atom n9("N9", Vector3D(1.0, 2.0, 3.0));
    Atom c8("C8", Vector3D(4.0, 5.0, 6.0));
    residue->add_atom(n9);
    residue->add_atom(c8);

    auto found = residue->find_atom("N9");
    ASSERT_TRUE(found.has_value());
    EXPECT_EQ(found->name(), "N9");

    auto not_found = residue->find_atom("XYZ");
    EXPECT_FALSE(not_found.has_value());
}

// === Reference Frame Tests (Nucleotides only) ===

TEST_F(ResidueHierarchyTest, NucleotideReferenceFrame) {
    auto residue = ResidueFactory::create_rna("G", 1, "A");

    EXPECT_FALSE(residue->reference_frame().has_value());

    // Create a simple reference frame (identity rotation, origin at (1,0,0))
    ReferenceFrame frame(Matrix3D::identity(), Vector3D(1.0, 0.0, 0.0));

    residue->set_reference_frame(frame);

    EXPECT_TRUE(residue->reference_frame().has_value());
}

// === Modified Nucleotide Tests ===

TEST_F(ResidueHierarchyTest, ModifiedNucleotideOneLetterCode) {
    // PSU (pseudouridine) should map to lowercase 'P' or similar
    auto psu = ResidueFactory::create("PSU", 1, "A", "");

    ASSERT_NE(psu, nullptr);
    EXPECT_TRUE(psu->is_nucleotide());

    auto* nuc = dynamic_cast<INucleotide*>(psu.get());
    ASSERT_NE(nuc, nullptr);
    // Modified nucleotides typically have lowercase codes
    char code = nuc->one_letter_code();
    EXPECT_TRUE(code == 'P' || code == 'p' || code == 'U' || code == 'u');
}

// === Direct Type Construction Tests ===

TEST_F(ResidueHierarchyTest, DirectRNAConstruction) {
    RNA rna("A", 1, "A", "");

    EXPECT_TRUE(rna.is_rna());
    EXPECT_FALSE(rna.is_dna());
    EXPECT_EQ(rna.name(), "A");
}

TEST_F(ResidueHierarchyTest, DirectDNAConstruction) {
    DNA dna("DT", 1, "B", "");

    EXPECT_TRUE(dna.is_dna());
    EXPECT_FALSE(dna.is_rna());
    EXPECT_EQ(dna.name(), "DT");
}

TEST_F(ResidueHierarchyTest, DirectProteinConstruction) {
    Protein protein("ALA", 1, "A", "");

    EXPECT_TRUE(protein.is_protein());
    EXPECT_FALSE(protein.is_nucleotide());
    EXPECT_EQ(protein.name(), "ALA");
}

TEST_F(ResidueHierarchyTest, DirectLigandConstruction) {
    Ligand ligand("HOH", 1, "W", "");

    EXPECT_TRUE(ligand.is_ligand());
    EXPECT_FALSE(ligand.is_nucleotide());
    EXPECT_FALSE(ligand.is_protein());
}

/**
 * @file test_role_classifier.cpp
 * @brief Unit tests for H-bond role classifier
 */

#include <gtest/gtest.h>
#include <x3dna/algorithms/hydrogen_bond/role_classifier.hpp>

using namespace x3dna::algorithms;
using namespace x3dna::core;

// ============================================================================
// Nucleotide atom role tests
// ============================================================================

TEST(RoleClassifierTest, AdenineAtomRoles) {
    // Backbone
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('A', " O1P"), HBondAtomRole::ACCEPTOR);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('A', " O2P"), HBondAtomRole::ACCEPTOR);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('A', " O2'"), HBondAtomRole::EITHER);

    // Base atoms
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('A', " N9 "), HBondAtomRole::EITHER); // Glycosidic
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('A', " N7 "), HBondAtomRole::ACCEPTOR);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('A', " N6 "), HBondAtomRole::DONOR);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('A', " N1 "), HBondAtomRole::ACCEPTOR);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('A', " N3 "), HBondAtomRole::ACCEPTOR);
}

TEST(RoleClassifierTest, GuanineAtomRoles) {
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('G', " N9 "), HBondAtomRole::EITHER);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('G', " O6 "), HBondAtomRole::ACCEPTOR);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('G', " N1 "), HBondAtomRole::DONOR);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('G', " N2 "), HBondAtomRole::DONOR);
}

TEST(RoleClassifierTest, CytosineAtomRoles) {
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('C', " N1 "), HBondAtomRole::EITHER);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('C', " O2 "), HBondAtomRole::ACCEPTOR);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('C', " N3 "), HBondAtomRole::ACCEPTOR);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('C', " N4 "), HBondAtomRole::DONOR);
}

TEST(RoleClassifierTest, UracilAtomRoles) {
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('U', " N1 "), HBondAtomRole::EITHER);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('U', " O2 "), HBondAtomRole::ACCEPTOR);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('U', " N3 "), HBondAtomRole::DONOR);
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('U', " O4 "), HBondAtomRole::ACCEPTOR);
}

TEST(RoleClassifierTest, UnknownBaseUsesElementFallback) {
    // Unknown base, N atoms fall back to element-based EITHER
    // (enables H-bond detection for modified nucleotides)
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('X', " N1 "), HBondAtomRole::EITHER);

    // Unknown base (PSU='P'), backbone atoms still return their role
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('P', " O2'"), HBondAtomRole::EITHER);

    // Unknown base, N atoms fall back to EITHER
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('P', " N1 "), HBondAtomRole::EITHER);

    // Carbon atoms return UNKNOWN (not H-bond capable)
    EXPECT_EQ(HBondRoleClassifier::get_nucleotide_atom_role('X', " C1 "), HBondAtomRole::UNKNOWN);
}

// ============================================================================
// Nucleotide bond classification tests
// ============================================================================

TEST(RoleClassifierTest, StandardBondClassification) {
    // A-D (Acceptor-Donor)
    auto result = HBondRoleClassifier::classify_nucleotide_bond('A', 'G', " N1 ", " N2 ");
    EXPECT_EQ(result, HBondClassification::STANDARD);

    // D-A (Donor-Acceptor)
    result = HBondRoleClassifier::classify_nucleotide_bond('G', 'A', " N2 ", " N1 ");
    EXPECT_EQ(result, HBondClassification::STANDARD);

    // X-D (Either-Donor)
    result = HBondRoleClassifier::classify_nucleotide_bond('A', 'G', " O2'", " N2 ");
    EXPECT_EQ(result, HBondClassification::STANDARD);
}

TEST(RoleClassifierTest, NonStandardBondClassification) {
    // A-A (Acceptor-Acceptor) - invalid
    auto result = HBondRoleClassifier::classify_nucleotide_bond('A', 'A', " N1 ", " N3 ");
    EXPECT_EQ(result, HBondClassification::NON_STANDARD);

    // D-D (Donor-Donor) - invalid
    result = HBondRoleClassifier::classify_nucleotide_bond('G', 'G', " N2 ", " N1 ");
    EXPECT_EQ(result, HBondClassification::NON_STANDARD);

    // Unknown base
    result = HBondRoleClassifier::classify_nucleotide_bond('X', 'A', " N1 ", " N1 ");
    EXPECT_EQ(result, HBondClassification::NON_STANDARD);
}

// ============================================================================
// Protein atom role tests
// ============================================================================

TEST(RoleClassifierTest, ProteinMainchainRoles) {
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("ALA", " N  "), HBondAtomRole::DONOR);
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("ALA", " O  "), HBondAtomRole::ACCEPTOR);
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("VAL", " OXT"), HBondAtomRole::ACCEPTOR);
}

TEST(RoleClassifierTest, ProteinSidechainRoles) {
    // Serine - hydroxyl
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("SER", " OG "), HBondAtomRole::EITHER);

    // Asparagine - amide
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("ASN", " OD1"), HBondAtomRole::ACCEPTOR);
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("ASN", " ND2"), HBondAtomRole::DONOR);

    // Aspartate - carboxyl
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("ASP", " OD1"), HBondAtomRole::ACCEPTOR);
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("ASP", " OD2"), HBondAtomRole::ACCEPTOR);

    // Lysine - amino
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("LYS", " NZ "), HBondAtomRole::DONOR);

    // Arginine - guanidinium
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("ARG", " NH1"), HBondAtomRole::DONOR);
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("ARG", " NH2"), HBondAtomRole::DONOR);

    // Histidine - imidazole
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("HIS", " ND1"), HBondAtomRole::EITHER);
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("HIS", " NE2"), HBondAtomRole::EITHER);
}

TEST(RoleClassifierTest, ProteinCaseInsensitive) {
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("ser", " OG "), HBondAtomRole::EITHER);
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("Ser", " OG "), HBondAtomRole::EITHER);
    EXPECT_EQ(HBondRoleClassifier::get_protein_atom_role("SER", " OG "), HBondAtomRole::EITHER);
}

TEST(RoleClassifierTest, IsMainchainAtom) {
    EXPECT_TRUE(HBondRoleClassifier::is_mainchain_atom(" N  "));
    EXPECT_TRUE(HBondRoleClassifier::is_mainchain_atom(" O  "));
    EXPECT_TRUE(HBondRoleClassifier::is_mainchain_atom(" OXT"));
    EXPECT_FALSE(HBondRoleClassifier::is_mainchain_atom(" CA "));
    EXPECT_FALSE(HBondRoleClassifier::is_mainchain_atom(" OG "));
}

// ============================================================================
// Ligand atom role tests
// ============================================================================

TEST(RoleClassifierTest, LigandElementBasedRoles) {
    // Nitrogen - could be donor or acceptor
    EXPECT_EQ(HBondRoleClassifier::get_ligand_atom_role(" N1 "), HBondAtomRole::EITHER);

    // Oxygen - could be donor or acceptor
    EXPECT_EQ(HBondRoleClassifier::get_ligand_atom_role(" O2 "), HBondAtomRole::EITHER);

    // Sulfur - thiol
    EXPECT_EQ(HBondRoleClassifier::get_ligand_atom_role(" SG "), HBondAtomRole::EITHER);

    // Carbon - not typically H-bond participant
    EXPECT_EQ(HBondRoleClassifier::get_ligand_atom_role(" C1 "), HBondAtomRole::UNKNOWN);
}

// ============================================================================
// General classification tests
// ============================================================================

TEST(RoleClassifierTest, GetAtomRoleByMoleculeType) {
    // Nucleic acid
    auto role = HBondRoleClassifier::get_atom_role(MoleculeType::NUCLEIC_ACID, "A", " N1 ");
    EXPECT_EQ(role, HBondAtomRole::ACCEPTOR);

    // Protein
    role = HBondRoleClassifier::get_atom_role(MoleculeType::PROTEIN, "SER", " OG ");
    EXPECT_EQ(role, HBondAtomRole::EITHER);

    // Ligand
    role = HBondRoleClassifier::get_atom_role(MoleculeType::LIGAND, "UNK", " N1 ");
    EXPECT_EQ(role, HBondAtomRole::EITHER);

    // Unknown type
    role = HBondRoleClassifier::get_atom_role(MoleculeType::UNKNOWN, "UNK", " N1 ");
    EXPECT_EQ(role, HBondAtomRole::UNKNOWN);
}

TEST(RoleClassifierTest, ClassifyByRoles) {
    // Valid combinations
    EXPECT_EQ(HBondRoleClassifier::classify_by_roles(HBondAtomRole::ACCEPTOR, HBondAtomRole::DONOR),
              HBondClassification::STANDARD);
    EXPECT_EQ(HBondRoleClassifier::classify_by_roles(HBondAtomRole::DONOR, HBondAtomRole::ACCEPTOR),
              HBondClassification::STANDARD);
    EXPECT_EQ(HBondRoleClassifier::classify_by_roles(HBondAtomRole::EITHER, HBondAtomRole::DONOR),
              HBondClassification::STANDARD);
    EXPECT_EQ(HBondRoleClassifier::classify_by_roles(HBondAtomRole::EITHER, HBondAtomRole::EITHER),
              HBondClassification::STANDARD);

    // Invalid combinations
    EXPECT_EQ(HBondRoleClassifier::classify_by_roles(HBondAtomRole::ACCEPTOR, HBondAtomRole::ACCEPTOR),
              HBondClassification::NON_STANDARD);
    EXPECT_EQ(HBondRoleClassifier::classify_by_roles(HBondAtomRole::DONOR, HBondAtomRole::DONOR),
              HBondClassification::NON_STANDARD);
    EXPECT_EQ(HBondRoleClassifier::classify_by_roles(HBondAtomRole::UNKNOWN, HBondAtomRole::DONOR),
              HBondClassification::NON_STANDARD);
}

// ============================================================================
// Utility tests
// ============================================================================

TEST(RoleClassifierTest, IsGoodHBondDistance) {
    EXPECT_TRUE(HBondRoleClassifier::is_good_hbond_distance(2.8));
    EXPECT_TRUE(HBondRoleClassifier::is_good_hbond_distance(2.5));
    EXPECT_TRUE(HBondRoleClassifier::is_good_hbond_distance(3.5));
    EXPECT_FALSE(HBondRoleClassifier::is_good_hbond_distance(2.4));
    EXPECT_FALSE(HBondRoleClassifier::is_good_hbond_distance(3.6));
}

TEST(RoleClassifierTest, CountGoodHBonds) {
    std::vector<HBond> bonds;

    // Good bond
    HBond bond1;
    bond1.distance = 2.8;
    bond1.classification = HBondClassification::STANDARD;
    bonds.push_back(bond1);

    // Too short
    HBond bond2;
    bond2.distance = 2.3;
    bond2.classification = HBondClassification::STANDARD;
    bonds.push_back(bond2);

    // Non-standard
    HBond bond3;
    bond3.distance = 2.9;
    bond3.classification = HBondClassification::NON_STANDARD;
    bonds.push_back(bond3);

    // Good bond
    HBond bond4;
    bond4.distance = 3.2;
    bond4.classification = HBondClassification::STANDARD;
    bonds.push_back(bond4);

    EXPECT_EQ(HBondRoleClassifier::count_good_hbonds(bonds), 2);
}

// ============================================================================
// Legacy compatibility tests
// ============================================================================

TEST(RoleClassifierTest, LegacyCompatibility) {
    // Legacy get_atom_role should call get_nucleotide_atom_role
    auto role1 = HBondRoleClassifier::get_atom_role('A', " N1 ");
    auto role2 = HBondRoleClassifier::get_nucleotide_atom_role('A', " N1 ");
    EXPECT_EQ(role1, role2);

    // Legacy classify_bond should call classify_nucleotide_bond
    auto class1 = HBondRoleClassifier::classify_bond('A', 'G', " N1 ", " N2 ");
    auto class2 = HBondRoleClassifier::classify_nucleotide_bond('A', 'G', " N1 ", " N2 ");
    EXPECT_EQ(class1, class2);
}

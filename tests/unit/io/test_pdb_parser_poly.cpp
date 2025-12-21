/**
 * @file test_pdb_parser_poly.cpp
 * @brief Unit tests for PDB parser with polymorphic types
 */

#include <gtest/gtest.h>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/residue/residue.hpp>
#include <sstream>

using namespace x3dna::io;
using namespace x3dna::core::poly;  // Polymorphic types

class PdbParserPolyTest : public ::testing::Test {
protected:
    PdbParser parser;

    // Simple PDB content with RNA nucleotides (properly formatted)
    const std::string rna_pdb = R"(HEADER    RNA STRUCTURE
ATOM      1  C1'   A A   1       1.000   2.000   3.000  1.00 20.00           C
ATOM      2  N9    A A   1       1.100   2.100   3.100  1.00 20.00           N
ATOM      3  C1'   G A   2       2.000   3.000   4.000  1.00 20.00           C
ATOM      4  N9    G A   2       2.100   3.100   4.100  1.00 20.00           N
)";

    // DNA content
    const std::string dna_pdb = R"(HEADER    DNA STRUCTURE
ATOM      1  C1'  DA B   1       1.000   2.000   3.000  1.00 20.00           C
ATOM      2  N9   DA B   1       1.100   2.100   3.100  1.00 20.00           N
ATOM      3  C1'  DT B   2       2.000   3.000   4.000  1.00 20.00           C
ATOM      4  N1   DT B   2       2.100   3.100   4.100  1.00 20.00           N
)";

    // Protein content
    const std::string protein_pdb = R"(HEADER    PROTEIN STRUCTURE
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.100   2.100   3.100  1.00 20.00           C
ATOM      3  N   GLY A   2       2.000   3.000   4.000  1.00 20.00           N
ATOM      4  CA  GLY A   2       2.100   3.100   4.100  1.00 20.00           C
)";

    // Mixed content (RNA + protein)
    const std::string mixed_pdb = R"(HEADER    MIXED STRUCTURE
ATOM      1  C1'   A A   1       1.000   2.000   3.000  1.00 20.00           C
ATOM      2  N9    A A   1       1.100   2.100   3.100  1.00 20.00           N
ATOM      3  N   ALA A   2       2.000   3.000   4.000  1.00 20.00           N
ATOM      4  CA  ALA A   2       2.100   3.100   4.100  1.00 20.00           C
)";
};

// === Basic Parsing ===

TEST_F(PdbParserPolyTest, ParseStringReturnsPolyStructure) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    EXPECT_FALSE(structure.empty());
    EXPECT_EQ(structure.size(), 1);  // One chain
}

TEST_F(PdbParserPolyTest, ParseStringCreatesCorrectChains) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    EXPECT_EQ(structure.size(), 1);
    EXPECT_EQ(structure[0].chain_id(), "A");
}

TEST_F(PdbParserPolyTest, ParseStringCreatesCorrectResidueCount) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    EXPECT_EQ(structure.num_residues(), 2);  // A, G
}

// === RNA Type Detection ===

TEST_F(PdbParserPolyTest, ParseRNACreatesRNAType) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    // Both residues should be RNA
    EXPECT_TRUE(structure[0][0].is_rna());
    EXPECT_TRUE(structure[0][0].is_nucleotide());
    EXPECT_FALSE(structure[0][0].is_dna());
    EXPECT_FALSE(structure[0][0].is_protein());

    EXPECT_TRUE(structure[0][1].is_rna());
    EXPECT_TRUE(structure[0][1].is_nucleotide());
}

TEST_F(PdbParserPolyTest, RNAResidueNames) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    EXPECT_EQ(structure[0][0].name(), "A");
    EXPECT_EQ(structure[0][1].name(), "G");
}

// === DNA Type Detection ===

TEST_F(PdbParserPolyTest, ParseDNACreatesDNAType) {
    Structure structure = parser.parse_string_poly(dna_pdb);

    EXPECT_EQ(structure.num_residues(), 2);

    // First residue should be DNA Adenine
    EXPECT_TRUE(structure[0][0].is_dna());
    EXPECT_TRUE(structure[0][0].is_nucleotide());
    EXPECT_FALSE(structure[0][0].is_rna());
    EXPECT_EQ(structure[0][0].name(), "DA");

    // Second residue should be DNA Thymine
    EXPECT_TRUE(structure[0][1].is_dna());
    EXPECT_EQ(structure[0][1].name(), "DT");
}

// === Protein Type Detection ===

TEST_F(PdbParserPolyTest, ParseProteinCreatesProteinType) {
    Structure structure = parser.parse_string_poly(protein_pdb);

    EXPECT_EQ(structure.num_residues(), 2);

    // Both should be protein
    EXPECT_TRUE(structure[0][0].is_protein());
    EXPECT_FALSE(structure[0][0].is_nucleotide());
    EXPECT_EQ(structure[0][0].name(), "ALA");

    EXPECT_TRUE(structure[0][1].is_protein());
    EXPECT_EQ(structure[0][1].name(), "GLY");
}

// === Mixed Types ===

TEST_F(PdbParserPolyTest, ParseMixedTypesCorrectly) {
    Structure structure = parser.parse_string_poly(mixed_pdb);

    EXPECT_EQ(structure.num_residues(), 2);

    // First is RNA
    EXPECT_TRUE(structure[0][0].is_rna());
    EXPECT_EQ(structure[0][0].name(), "A");

    // Second is protein
    EXPECT_TRUE(structure[0][1].is_protein());
    EXPECT_EQ(structure[0][1].name(), "ALA");
}

// === Nucleotide Interface ===

TEST_F(PdbParserPolyTest, RNANucleotideOneLetterCode) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    auto nucs = structure.nucleotides();
    ASSERT_EQ(nucs.size(), 2);

    EXPECT_EQ(nucs[0]->one_letter_code(), 'A');
    EXPECT_EQ(nucs[1]->one_letter_code(), 'G');
}

TEST_F(PdbParserPolyTest, RNAPurinePyrimidine) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    auto nucs = structure.nucleotides();
    ASSERT_EQ(nucs.size(), 2);

    // A and G are purines
    EXPECT_TRUE(nucs[0]->is_purine());
    EXPECT_TRUE(nucs[1]->is_purine());
    EXPECT_FALSE(nucs[0]->is_pyrimidine());
}

TEST_F(PdbParserPolyTest, DNANucleotideOneLetterCode) {
    Structure structure = parser.parse_string_poly(dna_pdb);

    auto nucs = structure.nucleotides();
    ASSERT_EQ(nucs.size(), 2);

    EXPECT_EQ(nucs[0]->one_letter_code(), 'A');  // DA -> A
    EXPECT_EQ(nucs[1]->one_letter_code(), 'T');  // DT -> T
}

// === Atom Preservation ===

TEST_F(PdbParserPolyTest, AtomsArePreserved) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    // Each RNA residue has 2 atoms
    EXPECT_EQ(structure[0][0].num_atoms(), 2);
    EXPECT_EQ(structure[0][1].num_atoms(), 2);
}

TEST_F(PdbParserPolyTest, AtomCoordinatesAreCorrect) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    auto atom_opt = structure[0][0].find_atom("N9");
    ASSERT_TRUE(atom_opt.has_value());
    EXPECT_DOUBLE_EQ(atom_opt->position().x(), 1.1);
    EXPECT_DOUBLE_EQ(atom_opt->position().y(), 2.1);
    EXPECT_DOUBLE_EQ(atom_opt->position().z(), 3.1);
}

// === Legacy Index Support ===

TEST_F(PdbParserPolyTest, LegacyResidueIdxIsSet) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    // Legacy indices should be 1-based and in order
    EXPECT_EQ(structure[0][0].legacy_residue_idx(), 1);
    EXPECT_EQ(structure[0][1].legacy_residue_idx(), 2);
}

TEST_F(PdbParserPolyTest, LegacyAtomIdxIsSet) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    // First atom should have legacy idx 1
    EXPECT_EQ(structure[0][0].atoms()[0].legacy_atom_idx(), 1);
    EXPECT_EQ(structure[0][0].atoms()[1].legacy_atom_idx(), 2);
}

// === Dynamic Cast Access ===

TEST_F(PdbParserPolyTest, DynamicCastToRNA) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    auto* rna = dynamic_cast<RNA*>(&structure[0][0]);
    ASSERT_NE(rna, nullptr);
    EXPECT_EQ(rna->one_letter_code(), 'A');
}

TEST_F(PdbParserPolyTest, DynamicCastToDNA) {
    Structure structure = parser.parse_string_poly(dna_pdb);

    auto* dna = dynamic_cast<DNA*>(&structure[0][0]);
    ASSERT_NE(dna, nullptr);
    EXPECT_EQ(dna->one_letter_code(), 'A');
}

TEST_F(PdbParserPolyTest, DynamicCastToProtein) {
    Structure structure = parser.parse_string_poly(protein_pdb);

    auto* protein = dynamic_cast<Protein*>(&structure[0][0]);
    ASSERT_NE(protein, nullptr);
    EXPECT_EQ(protein->name(), "ALA");
}

// === Stream Parsing ===

TEST_F(PdbParserPolyTest, ParseStreamPoly) {
    std::istringstream stream(rna_pdb);
    Structure structure = parser.parse_stream_poly(stream);

    EXPECT_EQ(structure.num_residues(), 2);
    EXPECT_TRUE(structure[0][0].is_rna());
}

// === Sequence from Structure ===

TEST_F(PdbParserPolyTest, ChainSequence) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    // Sequence should only include nucleotides
    EXPECT_EQ(structure[0].sequence(), "AG");
}

TEST_F(PdbParserPolyTest, MixedChainSequenceExcludesProtein) {
    Structure structure = parser.parse_string_poly(mixed_pdb);

    // Sequence should only include the RNA nucleotide
    EXPECT_EQ(structure[0].sequence(), "A");
}

// === Nucleotides filtering ===

TEST_F(PdbParserPolyTest, NucleotidesReturnsOnlyNucleotides) {
    Structure structure = parser.parse_string_poly(mixed_pdb);

    auto nucs = structure.nucleotides();
    EXPECT_EQ(nucs.size(), 1);  // Only A, not ALA
    EXPECT_EQ(nucs[0]->name(), "A");
}

TEST_F(PdbParserPolyTest, NucleotidesFromRNAOnly) {
    Structure structure = parser.parse_string_poly(rna_pdb);

    auto nucs = structure.nucleotides();
    EXPECT_EQ(nucs.size(), 2);
}

TEST_F(PdbParserPolyTest, NucleotidesFromProteinReturnsEmpty) {
    Structure structure = parser.parse_string_poly(protein_pdb);

    auto nucs = structure.nucleotides();
    EXPECT_EQ(nucs.size(), 0);
}


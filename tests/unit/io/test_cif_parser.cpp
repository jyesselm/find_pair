/**
 * @file test_cif_parser.cpp
 * @brief Unit tests for CifParser
 */

#include <gtest/gtest.h>
#include <x3dna/io/cif_parser.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/chain.hpp>
#include <filesystem>
#include <set>

namespace x3dna::test {
namespace io {

using namespace x3dna::io;
using namespace x3dna::core;

/**
 * @brief Test parsing a simple CIF string
 */
TEST(CifParserTest, ParseSimpleString) {
    // Simple mmCIF content with one residue
    std::string cif_content = R"(data_TEST
_entry.id TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C "C1'" C A 1 A 1 . . 1.000 2.000 3.000 1.00 20.00
ATOM 2 N N1 C A 1 A 1 . . 1.100 2.100 3.100 1.00 20.00
ATOM 3 C "C1'" G A 2 A 2 . . 2.000 3.000 4.000 1.00 20.00
ATOM 4 N N1 G A 2 A 2 . . 2.100 3.100 4.100 1.00 20.00
)";

    CifParser parser;
    Structure structure = parser.parse_string(cif_content);

    EXPECT_EQ(structure.num_atoms(), 4);
    EXPECT_EQ(structure.num_residues(), 2);
    EXPECT_EQ(structure.num_chains(), 1);
}

/**
 * @brief Test parsing ATOM records
 */
TEST(CifParserTest, ParseAtomRecords) {
    std::string cif_content = R"(data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C "C1'" C A 1 A 1 . . 1.000 2.000 3.000 1.00 20.00
ATOM 2 N N1 C A 1 A 1 . . 1.100 2.100 3.100 1.00 20.00
)";

    CifParser parser;
    Structure structure = parser.parse_string(cif_content);

    EXPECT_EQ(structure.num_atoms(), 2);

    // Check first atom
    auto chain = structure.find_chain("A");
    ASSERT_TRUE(chain.has_value());

    auto residues = chain.value().residues();
    ASSERT_EQ(residues.size(), 1);

    auto atoms = residues[0].atoms();
    ASSERT_EQ(atoms.size(), 2);

    // Check residue-level fields from residue (not atom)
    EXPECT_EQ(residues[0].chain_id(), "A");
    EXPECT_EQ(residues[0].seq_num(), 1);

    // Check atom-level fields
    EXPECT_DOUBLE_EQ(atoms[0].position().x(), 1.0);
    EXPECT_DOUBLE_EQ(atoms[0].position().y(), 2.0);
    EXPECT_DOUBLE_EQ(atoms[0].position().z(), 3.0);
}

/**
 * @brief Test parsing HETATM records (when enabled)
 */
TEST(CifParserTest, ParseHetatmRecords) {
    std::string cif_content = R"(data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C "C1'" C A 1 A 1 . . 1.000 2.000 3.000 1.00 20.00
HETATM 2 N N1 SPM A 21 A 21 . . 10.683 -8.783 22.839 1.00 40.13
)";

    CifParser parser;
    parser.set_include_hetatm(true);
    Structure structure = parser.parse_string(cif_content);

    EXPECT_EQ(structure.num_atoms(), 2);

    // Check HETATM atom
    auto chain = structure.find_chain("A");
    ASSERT_TRUE(chain.has_value());

    // Should have residues 1 and 21
    auto residues = chain.value().residues();
    ASSERT_GE(residues.size(), 2);

    // Find residue 21
    auto it = std::find_if(residues.begin(), residues.end(), [](const Residue& r) {
        return r.seq_num() == 21;
    });
    ASSERT_NE(it, residues.end());

    auto atoms = it->atoms();
    ASSERT_GE(atoms.size(), 1);
    }

/**
 * @brief Test HETATM exclusion (default)
 */
TEST(CifParserTest, ExcludeHetatmByDefault) {
    std::string cif_content = R"(data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C "C1'" C A 1 A 1 . . 1.000 2.000 3.000 1.00 20.00
HETATM 2 N N1 SPM A 21 A 21 . . 10.683 -8.783 22.839 1.00 40.13
)";

    CifParser parser;
    // Default: include_hetatm = false
    Structure structure = parser.parse_string(cif_content);

    EXPECT_EQ(structure.num_atoms(), 1); // Only ATOM, not HETATM
}

/**
 * @brief Test water exclusion
 */
TEST(CifParserTest, ExcludeWaters) {
    std::string cif_content = R"(data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C "C1'" C A 1 A 1 . . 1.000 2.000 3.000 1.00 20.00
HETATM 2 O O HOH A 22 A 22 . . 5.000 6.000 7.000 1.00 30.00
)";

    CifParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(false); // Exclude waters
    Structure structure = parser.parse_string(cif_content);

    EXPECT_EQ(structure.num_atoms(), 1); // Only ATOM, not HOH
}

/**
 * @brief Test chain identification
 */
TEST(CifParserTest, ParseMultipleChains) {
    std::string cif_content = R"(data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C "C1'" C A 1 A 1 . . 1.000 2.000 3.000 1.00 20.00
ATOM 2 C "C1'" G B 1 B 1 . . 2.000 3.000 4.000 1.00 20.00
)";

    CifParser parser;
    Structure structure = parser.parse_string(cif_content);

    EXPECT_EQ(structure.num_chains(), 2);
    EXPECT_TRUE(structure.find_chain("A").has_value());
    EXPECT_TRUE(structure.find_chain("B").has_value());
}

/**
 * @brief Test residue numbering
 * Note: Residue order in chain may not match CIF file order due to internal sorting
 */
TEST(CifParserTest, ParseResidueNumbering) {
    std::string cif_content = R"(data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C "C1'" C A 1 A 1 . . 1.000 2.000 3.000 1.00 20.00
ATOM 2 C "C1'" G A 2 A 2 . . 2.000 3.000 4.000 1.00 20.00
ATOM 3 C "C1'" A A 3 A 3 . . 3.000 4.000 5.000 1.00 20.00
)";

    CifParser parser;
    Structure structure = parser.parse_string(cif_content);

    auto chain = structure.find_chain("A");
    ASSERT_TRUE(chain.has_value());

    auto residues = chain.value().residues();
    EXPECT_EQ(residues.size(), 3);

    // Verify all three sequence numbers are present (order may vary)
    std::set<int> seq_nums;
    for (const auto& r : residues) {
        seq_nums.insert(r.seq_num());
    }
    EXPECT_TRUE(seq_nums.count(1) > 0);
    EXPECT_TRUE(seq_nums.count(2) > 0);
    EXPECT_TRUE(seq_nums.count(3) > 0);
}

/**
 * @brief Test parsing real CIF file
 */
TEST(CifParserTest, ParseRealCifFile) {
    std::filesystem::path cif_file = "data/cif/100D.cif";

    if (!std::filesystem::exists(cif_file)) {
        GTEST_SKIP() << "CIF file not found: " << cif_file;
    }

    CifParser parser;
    Structure structure = parser.parse_file(cif_file);

    EXPECT_GT(structure.num_atoms(), 0);
    EXPECT_GT(structure.num_residues(), 0);
    EXPECT_GT(structure.num_chains(), 0);

    // Verify we can find specific atoms
    auto chain = structure.find_chain("A");
    ASSERT_TRUE(chain.has_value());

    auto residues = chain.value().residues();
    ASSERT_GT(residues.size(), 0);

    // Check first residue has atoms
    auto atoms = residues[0].atoms();
    EXPECT_GT(atoms.size(), 0);
}

/**
 * @brief Test error handling for missing file
 */
TEST(CifParserTest, ErrorOnMissingFile) {
    CifParser parser;
    std::filesystem::path missing_file = "data/cif/nonexistent.cif";

    EXPECT_THROW(parser.parse_file(missing_file), CifParser::ParseError);
}

/**
 * @brief Test atom name normalization
 * CIF atom names should be converted to 4-character PDB format
 */
TEST(CifParserTest, AtomNameNormalization) {
    std::string cif_content = R"(data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C "C1'" C A 1 A 1 . . 1.000 2.000 3.000 1.00 20.00
ATOM 2 N N1 C A 1 A 1 . . 1.100 2.100 3.100 1.00 20.00
)";

    CifParser parser;
    Structure structure = parser.parse_string(cif_content);

    auto chain = structure.find_chain("A");
    ASSERT_TRUE(chain.has_value());

    auto residues = chain.value().residues();
    ASSERT_EQ(residues.size(), 1);

    auto atoms = residues[0].atoms();
    ASSERT_GE(atoms.size(), 2);

    // Atom names are now trimmed; original names preserved for PDB output
    EXPECT_EQ(atoms[0].name(), "C1'");
    EXPECT_EQ(atoms[1].name(), "N1");
    // Original names are preserved in PDB format for round-trip
        }

/**
 * @brief Test phosphate atom name conversion
 * OP1 -> O1P, OP2 -> O2P (CIF to PDB convention)
 */
TEST(CifParserTest, PhosphateAtomNameConversion) {
    std::string cif_content = R"(data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 P P G A 1 A 1 . . 1.000 2.000 3.000 1.00 20.00
ATOM 2 O OP1 G A 1 A 1 . . 1.100 2.100 3.100 1.00 20.00
ATOM 3 O OP2 G A 1 A 1 . . 1.200 2.200 3.200 1.00 20.00
)";

    CifParser parser;
    Structure structure = parser.parse_string(cif_content);

    auto chain = structure.find_chain("A");
    ASSERT_TRUE(chain.has_value());

    auto residues = chain.value().residues();
    ASSERT_EQ(residues.size(), 1);

    auto atoms = residues[0].atoms();
    ASSERT_EQ(atoms.size(), 3);

    // P (trimmed name)
    EXPECT_EQ(atoms[0].name(), "P"); // name() returns trimmed
    // OP1 becomes O1P (trimmed)
    EXPECT_EQ(atoms[1].name(), "O1P");
    // OP2 becomes O2P (trimmed)
    EXPECT_EQ(atoms[2].name(), "O2P");
}

/**
 * @brief Test alternate conformation handling
 * Should keep atoms with alt_loc ' ', 'A', or '1' and skip others
 */
TEST(CifParserTest, AlternateConformationFilter) {
    std::string cif_content = R"(data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C "C1'" C A 1 A 1 . . 1.000 2.000 3.000 1.00 20.00
ATOM 2 C "C1'" C A 1 A 1 . A 1.100 2.100 3.100 0.50 20.00
ATOM 3 C "C1'" C A 1 A 1 . B 1.200 2.200 3.200 0.50 20.00
)";

    CifParser parser;
    Structure structure = parser.parse_string(cif_content);

    // Should have 2 atoms: one with no alt_loc and one with 'A'
    // The 'B' alt_loc atom should be skipped
    EXPECT_EQ(structure.num_atoms(), 2);
}

/**
 * @brief Test modified nucleotide auto-inclusion
 * Modified nucleotides should be included even without include_hetatm flag
 */
TEST(CifParserTest, ModifiedNucleotideAutoInclude) {
    std::string cif_content = R"(data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C "C1'" C A 1 A 1 . . 1.000 2.000 3.000 1.00 20.00
HETATM 2 N N1 PSU A 2 A 2 . . 2.000 3.000 4.000 1.00 20.00
)";

    CifParser parser;
    // include_hetatm is false by default
    Structure structure = parser.parse_string(cif_content);

    // Should have 2 atoms: regular C and PSU (modified nucleotide)
    EXPECT_EQ(structure.num_atoms(), 2);
}

/**
 * @brief Test empty content handling
 */
TEST(CifParserTest, ErrorOnEmptyContent) {
    std::string cif_content = "";

    CifParser parser;
    EXPECT_THROW(parser.parse_string(cif_content), CifParser::ParseError);
}

/**
 * @brief Test legacy index assignment
 * Legacy indices should be assigned sequentially
 */
TEST(CifParserTest, LegacyIndexAssignment) {
    std::string cif_content = R"(data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C "C1'" C A 1 A 1 . . 1.000 2.000 3.000 1.00 20.00
ATOM 2 N N1 C A 1 A 1 . . 1.100 2.100 3.100 1.00 20.00
ATOM 3 C "C1'" G A 2 A 2 . . 2.000 3.000 4.000 1.00 20.00
)";

    CifParser parser;
    Structure structure = parser.parse_string(cif_content);

    auto chain = structure.find_chain("A");
    ASSERT_TRUE(chain.has_value());

    auto residues = chain.value().residues();
    ASSERT_EQ(residues.size(), 2);

    // First residue atoms
    auto atoms1 = residues[0].atoms();
    EXPECT_EQ(atoms1[0].legacy_atom_idx(), 1);
    EXPECT_EQ(atoms1[1].legacy_atom_idx(), 2);
    EXPECT_EQ(residues[0].legacy_residue_idx(), 1);

    // Second residue atoms
    auto atoms2 = residues[1].atoms();
    EXPECT_EQ(atoms2[0].legacy_atom_idx(), 3);
    EXPECT_EQ(residues[1].legacy_residue_idx(), 2);
}

} // namespace io
} // namespace x3dna::test

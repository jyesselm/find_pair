/**
 * @file test_pdb_parser.cpp
 * @brief Unit tests for PdbParser
 */

#include <gtest/gtest.h>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/chain.hpp>
#include <sstream>
#include <filesystem>
#include <set>

namespace x3dna::test {
namespace io {

using namespace x3dna::io;
using namespace x3dna::core;

/**
 * @brief Test parsing a simple PDB string
 */
TEST(PdbParserTest, ParseSimpleString) {
    std::string pdb_content = R"(HEADER    TEST STRUCTURE
ATOM      1  C1'   C A   1       1.000   2.000   3.000  1.00 20.00           C  
ATOM      2  N1    C A   1       1.100   2.100   3.100  1.00 20.00           N  
ATOM      3  C1'   G A   2       2.000   3.000   4.000  1.00 20.00           C  
ATOM      4  N1    G A   2       2.100   3.100   4.100  1.00 20.00           N  
)";

    PdbParser parser;
    Structure structure = parser.parse_string(pdb_content);

    EXPECT_EQ(structure.num_atoms(), 4);
    EXPECT_EQ(structure.num_residues(), 2);
    EXPECT_EQ(structure.num_chains(), 1);
}

/**
 * @brief Test parsing ATOM records
 */
TEST(PdbParserTest, ParseAtomRecords) {
    std::string pdb_content =
        R"(ATOM      1  C1'   C A   1       1.000   2.000   3.000  1.00 20.00           C  
ATOM      2  N1    C A   1       1.100   2.100   3.100  1.00 20.00           N  
)";

    PdbParser parser;
    Structure structure = parser.parse_string(pdb_content);

    EXPECT_EQ(structure.num_atoms(), 2);

    // Check first atom
    auto chain = structure.find_chain("A");
    ASSERT_TRUE(chain.has_value());

    auto residues = chain.value().residues();
    ASSERT_EQ(residues.size(), 1);

    auto atoms = residues[0].atoms();
    ASSERT_EQ(atoms.size(), 2);

    // Check atom properties
    EXPECT_EQ(atoms[0].name(), "C1'"); // name() returns trimmed
    EXPECT_DOUBLE_EQ(atoms[0].position().x(), 1.0);
    EXPECT_DOUBLE_EQ(atoms[0].position().y(), 2.0);
    EXPECT_DOUBLE_EQ(atoms[0].position().z(), 3.0);

    // Check residue properties (residue-level fields now in Residue, not Atom)
    EXPECT_EQ(residues[0].name(), "C"); // GEMMI trims residue names
    EXPECT_EQ(residues[0].chain_id(), "A");
    EXPECT_EQ(residues[0].seq_num(), 1);
}

/**
 * @brief Test parsing HETATM records (when enabled)
 */
TEST(PdbParserTest, ParseHetatmRecords) {
    std::string pdb_content =
        R"(ATOM      1  C1'   C A   1       1.000   2.000   3.000  1.00 20.00           C  
HETATM    2  N1  SPM A  21      10.683  -8.783  22.839  1.00 40.13           N  
)";

    PdbParser parser;
    parser.set_include_hetatm(true);
    Structure structure = parser.parse_string(pdb_content);

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
TEST(PdbParserTest, ExcludeHetatmByDefault) {
    std::string pdb_content =
        R"(ATOM      1  C1'   C A   1       1.000   2.000   3.000  1.00 20.00           C  
HETATM    2  N1  SPM A  21      10.683  -8.783  22.839  1.00 40.13           N  
)";

    PdbParser parser;
    // Default: include_hetatm = false
    Structure structure = parser.parse_string(pdb_content);

    EXPECT_EQ(structure.num_atoms(), 1); // Only ATOM, not HETATM
}

/**
 * @brief Test water exclusion
 */
TEST(PdbParserTest, ExcludeWaters) {
    std::string pdb_content =
        R"(ATOM      1  C1'   C A   1       1.000   2.000   3.000  1.00 20.00           C  
HETATM    2  O   HOH A  22       5.000   6.000   7.000  1.00 30.00           O  
)";

    PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(false); // Exclude waters
    Structure structure = parser.parse_string(pdb_content);

    EXPECT_EQ(structure.num_atoms(), 1); // Only ATOM, not HOH
}

/**
 * @brief Test chain identification
 */
TEST(PdbParserTest, ParseMultipleChains) {
    std::string pdb_content =
        R"(ATOM      1  C1'   C A   1       1.000   2.000   3.000  1.00 20.00           C  
ATOM      2  C1'   G B   1       2.000   3.000   4.000  1.00 20.00           C  
)";

    PdbParser parser;
    Structure structure = parser.parse_string(pdb_content);

    EXPECT_EQ(structure.num_chains(), 2);
    EXPECT_TRUE(structure.find_chain("A").has_value());
    EXPECT_TRUE(structure.find_chain("B").has_value());
}

/**
 * @brief Test residue numbering
 */
TEST(PdbParserTest, ParseResidueNumbering) {
    std::string pdb_content =
        R"(ATOM      1  C1'   C A   1       1.000   2.000   3.000  1.00 20.00           C
ATOM      2  C1'   G A   2       2.000   3.000   4.000  1.00 20.00           C
ATOM      3  C1'   A A   3       3.000   4.000   5.000  1.00 20.00           C
)";

    PdbParser parser;
    Structure structure = parser.parse_string(pdb_content);

    auto chain = structure.find_chain("A");
    ASSERT_TRUE(chain.has_value());

    auto residues = chain.value().residues();
    EXPECT_EQ(residues.size(), 3);

    // Check that all expected sequence numbers are present
    // (order may vary due to map key sorting)
    std::set<int> seq_nums;
    for (const auto& res : residues) {
        seq_nums.insert(res.seq_num());
    }
    EXPECT_TRUE(seq_nums.count(1) > 0);
    EXPECT_TRUE(seq_nums.count(2) > 0);
    EXPECT_TRUE(seq_nums.count(3) > 0);
}

/**
 * @brief Test parsing real PDB file
 */
TEST(PdbParserTest, ParseRealPdbFile) {
    std::filesystem::path pdb_file = "data/pdb/100D.pdb";

    if (!std::filesystem::exists(pdb_file)) {
        GTEST_SKIP() << "PDB file not found: " << pdb_file;
    }

    PdbParser parser;
    Structure structure = parser.parse_file(pdb_file);

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
TEST(PdbParserTest, ErrorOnMissingFile) {
    PdbParser parser;
    std::filesystem::path missing_file = "data/pdb/nonexistent.pdb";

    EXPECT_THROW(parser.parse_file(missing_file), PdbParser::ParseError);
}

/**
 * @brief Test error handling for malformed coordinates
 */
TEST(PdbParserTest, ErrorOnMalformedCoordinates) {
    std::string pdb_content =
        R"(ATOM      1  C1'   C A   1       invalid   2.000   3.000  1.00 20.00           C  
)";

    PdbParser parser;
    // Should skip malformed lines (current implementation)
    // In the future, we might want to throw or collect errors
    EXPECT_NO_THROW(parser.parse_string(pdb_content));
}

/**
 * @brief Test atom name normalization
 */
TEST(PdbParserTest, AtomNameNormalization) {
    std::string pdb_content =
        R"(ATOM      1  C1'   C A   1       1.000   2.000   3.000  1.00 20.00           C
ATOM      2  N1    C A   1       1.100   2.100   3.100  1.00 20.00           N
)";

    PdbParser parser;
    Structure structure = parser.parse_string(pdb_content);

    auto chain = structure.find_chain("A");
    ASSERT_TRUE(chain.has_value());

    auto residues = chain.value().residues();
    ASSERT_EQ(residues.size(), 1);

    auto atoms = residues[0].atoms();
    ASSERT_GE(atoms.size(), 2);

    // Atom names are now trimmed; original names preserved for PDB output
    EXPECT_EQ(atoms[0].name(), "C1'");
    EXPECT_EQ(atoms[1].name(), "N1");
    // Original names are preserved for round-trip
        }

/**
 * @brief Test residue name normalization
 */
TEST(PdbParserTest, ResidueNameNormalization) {
    std::string pdb_content =
        R"(ATOM      1  C1'   C A   1       1.000   2.000   3.000  1.00 20.00           C
ATOM      2  C1'   G A   2       2.000   3.000   4.000  1.00 20.00           C
)";

    PdbParser parser;
    Structure structure = parser.parse_string(pdb_content);

    auto chain = structure.find_chain("A");
    ASSERT_TRUE(chain.has_value());

    auto residues = chain.value().residues();
    ASSERT_EQ(residues.size(), 2);

    // GEMMI trims residue names, so they are no longer padded to 3 characters
    // Just verify names are non-empty and correctly identified
    std::set<std::string> names;
    for (const auto& res : residues) {
        EXPECT_FALSE(res.name().empty());
        names.insert(res.name());
    }
    EXPECT_TRUE(names.count("C") > 0 || names.count("  C") > 0);
    EXPECT_TRUE(names.count("G") > 0 || names.count("  G") > 0);
}

/**
 * @brief Test parsing from stream
 */
TEST(PdbParserTest, ParseFromStream) {
    std::string pdb_content =
        R"(ATOM      1  C1'   C A   1       1.000   2.000   3.000  1.00 20.00           C  
ATOM      2  N1    C A   1       1.100   2.100   3.100  1.00 20.00           N  
)";

    std::istringstream stream(pdb_content);
    PdbParser parser;
    Structure structure = parser.parse_stream(stream);

    EXPECT_EQ(structure.num_atoms(), 2);
}

} // namespace io
} // namespace x3dna::test

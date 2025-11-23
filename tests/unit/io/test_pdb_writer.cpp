/**
 * @file test_pdb_writer.cpp
 * @brief Unit tests for PdbWriter class
 */

#include <gtest/gtest.h>
#include <x3dna/io/pdb_writer.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/structure.hpp>
#include <filesystem>
#include <fstream>
#include <sstream>

using namespace x3dna::io;
using namespace x3dna::core;
using namespace x3dna::geometry;

class PdbWriterTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a test structure
        structure_ = Structure("TEST");
        Chain chain('A');
        Residue residue("  C", 1, 'A');

        Atom atom(" C1'", Vector3D(1.234, 5.678, 9.012), "  C", 'A', 1);
        atom.set_occupancy(1.0);
        atom.set_b_factor(20.5);
        residue.add_atom(atom);

        chain.add_residue(residue);
        structure_.add_chain(chain);

        output_file_ = "test_output.pdb";
    }

    void TearDown() override {
        if (std::filesystem::exists(output_file_)) {
            std::filesystem::remove(output_file_);
        }
    }

    Structure structure_;
    std::string output_file_;
};

// File writing tests
TEST_F(PdbWriterTest, WriteFile) {
    PdbWriter writer;
    writer.write_file(structure_, output_file_);

    EXPECT_TRUE(std::filesystem::exists(output_file_));

    // Read back and verify format
    std::ifstream file(output_file_);
    std::string line;
    bool found_atom = false;

    while (std::getline(file, line)) {
        if (line.substr(0, 4) == "ATOM") {
            found_atom = true;
            // Verify format
            EXPECT_GE(line.length(), 54); // Minimum for coordinates
            break;
        }
        if (line.substr(0, 3) == "END") {
            break;
        }
    }

    EXPECT_TRUE(found_atom);
}

// Stream writing tests
TEST_F(PdbWriterTest, WriteStream) {
    PdbWriter writer;
    std::ostringstream stream;
    writer.write_stream(structure_, stream);

    std::string output = stream.str();
    EXPECT_FALSE(output.empty());
    EXPECT_NE(output.find("ATOM"), std::string::npos);
    EXPECT_NE(output.find("END"), std::string::npos);
}

// String conversion tests
TEST_F(PdbWriterTest, ToString) {
    PdbWriter writer;
    std::string pdb_string = writer.to_string(structure_);

    EXPECT_FALSE(pdb_string.empty());
    EXPECT_NE(pdb_string.find("ATOM"), std::string::npos);
}

// Round-trip test
TEST_F(PdbWriterTest, RoundTrip) {
    PdbWriter writer;
    writer.write_file(structure_, output_file_);

    // Parse back
    PdbParser parser;
    Structure parsed = parser.parse_file(output_file_);

    EXPECT_EQ(parsed.num_atoms(), structure_.num_atoms());
    EXPECT_EQ(parsed.num_residues(), structure_.num_residues());
}

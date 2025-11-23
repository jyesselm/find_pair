/**
 * @file test_json_reader.cpp
 * @brief Unit tests for JsonReader class
 */

#include <gtest/gtest.h>
#include <x3dna/io/json_reader.hpp>
#include <x3dna/core/structure.hpp>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>

using namespace x3dna::io;
using namespace x3dna::core;

class JsonReaderTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a test JSON file
        test_json_file_ = "test_reader.json";
        create_test_json();
    }

    void TearDown() override {
        if (std::filesystem::exists(test_json_file_)) {
            std::filesystem::remove(test_json_file_);
        }
    }

    void create_test_json() {
        nlohmann::json json;
        json["pdb_file"] = "test.pdb";
        json["pdb_name"] = "TEST";
        json["calculations"] = nlohmann::json::array();

        // Add pdb_atoms record
        nlohmann::json atoms_record;
        atoms_record["type"] = "pdb_atoms";
        atoms_record["pdb_id"] = "TEST";
        atoms_record["num_atoms"] = 2;

        nlohmann::json atoms_array = nlohmann::json::array();
        nlohmann::json atom1;
        atom1["atom_name"] = " C1'";
        atom1["xyz"] = {1.0, 2.0, 3.0};
        atom1["residue_name"] = "  C";
        atom1["chain_id"] = "A";
        atom1["residue_seq"] = 1;
        atom1["record_type"] = "A";
        atoms_array.push_back(atom1);

        nlohmann::json atom2;
        atom2["atom_name"] = " N3 ";
        atom2["xyz"] = {4.0, 5.0, 6.0};
        atom2["residue_name"] = "  G";
        atom2["chain_id"] = "A";
        atom2["residue_seq"] = 2;
        atom2["record_type"] = "A";
        atoms_array.push_back(atom2);

        atoms_record["atoms"] = atoms_array;
        json["calculations"].push_back(atoms_record);

        // Write to file
        std::ofstream file(test_json_file_);
        file << json.dump(2);
    }

    std::string test_json_file_;
};

// File reading tests
TEST_F(JsonReaderTest, ReadStructureLegacyFromFile) {
    std::filesystem::path file_path(test_json_file_);
    Structure structure = JsonReader::read_structure_legacy(file_path);
    EXPECT_EQ(structure.num_atoms(), 2);
    EXPECT_GT(structure.num_residues(), 0);
}

// Record finding tests
TEST_F(JsonReaderTest, FindRecordsByType) {
    std::ifstream file(test_json_file_);
    nlohmann::json json;
    file >> json;

    auto records = JsonReader::find_records_by_type(json, "pdb_atoms");
    EXPECT_EQ(records.size(), 1);
    EXPECT_EQ(records[0]["type"], "pdb_atoms");
}

// Multiple record types
TEST_F(JsonReaderTest, FindMultipleRecordTypes) {
    nlohmann::json json;
    json["calculations"] = nlohmann::json::array();

    nlohmann::json record1;
    record1["type"] = "base_frame_calc";
    json["calculations"].push_back(record1);

    nlohmann::json record2;
    record2["type"] = "ls_fitting";
    json["calculations"].push_back(record2);

    nlohmann::json record3;
    record3["type"] = "base_frame_calc";
    json["calculations"].push_back(record3);

    auto records = JsonReader::find_records_by_type(json, "base_frame_calc");
    EXPECT_EQ(records.size(), 2);

    auto ls_records = JsonReader::find_records_by_type(json, "ls_fitting");
    EXPECT_EQ(ls_records.size(), 1);
}

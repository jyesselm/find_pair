/**
 * @file test_json_writer.cpp
 * @brief Unit tests for JsonWriter class
 */

#include <gtest/gtest.h>
#include <x3dna/io/json_writer.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/core/parameters.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>
#include <filesystem>
#include <fstream>

using namespace x3dna::io;
using namespace x3dna::core;
using namespace x3dna::geometry;

class JsonWriterTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_pdb_path_ = std::filesystem::path("test.pdb");
        writer_ = std::make_unique<JsonWriter>(test_pdb_path_);
    }

    void TearDown() override {
        // Clean up any test files
        if (std::filesystem::exists("test_output.json")) {
            std::filesystem::remove("test_output.json");
        }
    }

    std::filesystem::path test_pdb_path_;
    std::unique_ptr<JsonWriter> writer_;
};

// Constructor tests
TEST_F(JsonWriterTest, Constructor) {
    EXPECT_EQ(writer_->json()["pdb_file"], test_pdb_path_.string());
    EXPECT_EQ(writer_->json()["pdb_name"], "test");
    EXPECT_TRUE(writer_->json()["calculations"].is_array());
}

// PDB atoms recording
TEST_F(JsonWriterTest, RecordPdbAtoms) {
    Structure structure("TEST");
    Chain chain('A');
    Residue residue("  C", 1, 'A');
    residue.add_atom(Atom(" C1'", Vector3D(1.0, 2.0, 3.0), "  C", 'A', 1));
    chain.add_residue(residue);
    structure.add_chain(chain);

    writer_->record_pdb_atoms(structure);

    auto& json = writer_->json();
    EXPECT_TRUE(json.contains("calculations"));

    // Find pdb_atoms record
    bool found = false;
    for (const auto& record : json["calculations"]) {
        if (record.contains("type") && record["type"] == "pdb_atoms") {
            found = true;
            EXPECT_EQ(record["num_atoms"], 1);
            EXPECT_TRUE(record.contains("atoms"));
            EXPECT_EQ(record["atoms"].size(), 1);

            const auto& atom = record["atoms"][0];
            EXPECT_EQ(atom["atom_name"], " C1'");
            EXPECT_EQ(atom["residue_name"], "  C");
            break;
        }
    }
    EXPECT_TRUE(found);
}

// Base frame calc recording
TEST_F(JsonWriterTest, RecordBaseFrameCalc) {
    std::vector<std::string> matched_atoms = {" N3 ", " C2 ", " N1 "};
    writer_->record_base_frame_calc(0, 'A', "data/templates/Atomic_A.pdb", 0.001234, matched_atoms);

    auto& json = writer_->json();
    bool found = false;
    for (const auto& record : json["calculations"]) {
        if (record.contains("type") && record["type"] == "base_frame_calc") {
            found = true;
            EXPECT_EQ(record["residue_idx"], 0);
            EXPECT_EQ(record["base_type"], "A");
            EXPECT_EQ(record["num_matched_atoms"], 3);
            EXPECT_TRUE(record.contains("matched_atoms"));
            break;
        }
    }
    EXPECT_TRUE(found);
}

// LS fitting recording
TEST_F(JsonWriterTest, RecordLsFitting) {
    Matrix3D rotation = Matrix3D::identity();
    Vector3D translation(1.0, 2.0, 3.0);
    writer_->record_ls_fitting(0, 9, 0.001234, rotation, translation);

    auto& json = writer_->json();
    bool found = false;
    for (const auto& record : json["calculations"]) {
        if (record.contains("type") && record["type"] == "ls_fitting") {
            found = true;
            EXPECT_EQ(record["residue_idx"], 0);
            EXPECT_EQ(record["num_points"], 9);
            EXPECT_TRUE(record.contains("rotation_matrix"));
            EXPECT_TRUE(record.contains("translation"));
            break;
        }
    }
    EXPECT_TRUE(found);
}

// Base pair recording
TEST_F(JsonWriterTest, RecordBasePair) {
    BasePair bp(0, 1, BasePairType::WATSON_CRICK);
    bp.set_bp_type("CG");

    Matrix3D rot = Matrix3D::identity();
    Vector3D org1(0, 0, 0);
    Vector3D org2(10, 0, 0);
    ReferenceFrame frame1(rot, org1);
    ReferenceFrame frame2(rot, org2);
    bp.set_frame1(frame1);
    bp.set_frame2(frame2);

    writer_->record_base_pair(bp);

    auto& json = writer_->json();
    bool found = false;
    for (const auto& record : json["calculations"]) {
        if (record.contains("type") && record["type"] == "base_pair") {
            found = true;
            EXPECT_EQ(record["residue_idx1"], 0);
            EXPECT_EQ(record["residue_idx2"], 1);
            EXPECT_EQ(record["bp_type"], "CG");
            break;
        }
    }
    EXPECT_TRUE(found);
}

// Removed atom recording
TEST_F(JsonWriterTest, RecordRemovedAtom) {
    Vector3D xyz(1.0, 2.0, 3.0);
    writer_->record_removed_atom("ATOM   1234  C1'  C   A   1 ", "line_too_short", 1234, " C1'",
                                 "  C", 'A', 1, &xyz, 0);

    auto& json = writer_->json();
    bool found = false;
    for (const auto& record : json["calculations"]) {
        if (record.contains("type") && record["type"] == "removed_atom") {
            found = true;
            EXPECT_EQ(record["reason"], "line_too_short");
            EXPECT_EQ(record["atom_serial"], 1234);
            break;
        }
    }
    EXPECT_TRUE(found);
}

// File writing
TEST_F(JsonWriterTest, WriteToFile) {
    writer_->record_base_frame_calc(0, 'A', "Atomic_A.pdb", 0.001, {" N3 "});

    std::filesystem::path output_file = "test_output.json";
    writer_->write_to_file(output_file);

    EXPECT_TRUE(std::filesystem::exists(output_file));

    // Read back and verify
    std::ifstream file(output_file);
    nlohmann::json read_json;
    file >> read_json;

    EXPECT_EQ(read_json["pdb_name"], "test");
    EXPECT_TRUE(read_json.contains("calculations"));
}

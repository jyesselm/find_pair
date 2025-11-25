/**
 * @file test_residue_ordering_json_comparison.cpp
 * @brief Test residue ordering by comparing JSON outputs
 * 
 * This test generates JSON files for residue ordering and compares them
 * to ensure ordering matches legacy.
 */

#include <gtest/gtest.h>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/structure_legacy_order.hpp>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <vector>
#include <map>

using namespace x3dna::core;
using namespace x3dna::io;
using json = nlohmann::json;

class ResidueOrderingJsonTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_pdb_ = "data/pdb/3G8T.pdb";
        output_dir_ = "data/residue_ordering";
        
        if (!std::filesystem::exists(test_pdb_)) {
            GTEST_SKIP() << "Test PDB file not found: " << test_pdb_;
        }
        
        // Create output directory
        std::filesystem::create_directories(output_dir_);
    }
    
    std::filesystem::path test_pdb_;
    std::filesystem::path output_dir_;
    
    struct ResidueInfo {
        int legacy_index;
        std::string residue_name;
        std::string chain_id;
        int residue_seq;
        std::string insertion_code;
        
        bool operator==(const ResidueInfo& other) const {
            return residue_name == other.residue_name &&
                   chain_id == other.chain_id &&
                   residue_seq == other.residue_seq &&
                   insertion_code == other.insertion_code;
        }
    };
    
    std::vector<ResidueInfo> generate_residue_ordering_json(const std::filesystem::path& pdb_file,
                                                             const std::filesystem::path& json_file) {
        // Parse PDB
        PdbParser parser;
        parser.set_include_hetatm(true);
        parser.set_include_waters(true);
        Structure structure = parser.parse_file(pdb_file);
        
        // Get residues in legacy order (using Structure's built-in method)
        auto residues = structure.residues_in_legacy_order();
        
        // Build JSON
        json output;
        output["pdb_id"] = pdb_file.stem().string();
        output["total_residues"] = residues.size();
        output["residues"] = json::array();
        
        std::vector<ResidueInfo> residue_infos;
        
        for (size_t i = 0; i < residues.size(); i++) {
            const Residue* res = residues[i];
            ResidueInfo info;
            info.legacy_index = static_cast<int>(i + 1);
            info.residue_name = res->name();
            info.chain_id = std::string(1, res->chain_id());
            info.residue_seq = res->seq_num();
            info.insertion_code = std::string(1, res->insertion());
            residue_infos.push_back(info);
            
            json residue_json;
            residue_json["legacy_index"] = info.legacy_index;
            residue_json["residue_name"] = info.residue_name;
            residue_json["chain_id"] = info.chain_id;
            residue_json["residue_seq"] = info.residue_seq;
            residue_json["insertion_code"] = info.insertion_code;
            residue_json["num_atoms"] = res->num_atoms();
            
            output["residues"].push_back(residue_json);
        }
        
        // Write JSON
        std::ofstream out_file(json_file);
        if (out_file.is_open()) {
            out_file << std::setw(2) << output << "\n";
        }
        
        return residue_infos;
    }
    
    std::vector<ResidueInfo> load_residue_ordering_json(const std::filesystem::path& json_file) {
        std::ifstream in_file(json_file);
        if (!in_file.is_open()) {
            return {};
        }
        
        json j;
        in_file >> j;
        
        std::vector<ResidueInfo> residue_infos;
        
        if (j.contains("residues") && j["residues"].is_array()) {
            for (const auto& res_json : j["residues"]) {
                ResidueInfo info;
                info.legacy_index = res_json.value("legacy_index", 0);
                info.residue_name = res_json.value("residue_name", "");
                std::string chain_str = res_json.value("chain_id", " ");
                info.chain_id = chain_str.empty() ? " " : std::string(1, chain_str[0]);
                info.residue_seq = res_json.value("residue_seq", 0);
                std::string ins_str = res_json.value("insertion_code", " ");
                info.insertion_code = ins_str.empty() ? " " : std::string(1, ins_str[0]);
                residue_infos.push_back(info);
            }
        }
        
        return residue_infos;
    }
};

/**
 * @test Generate and verify residue ordering JSON
 */
TEST_F(ResidueOrderingJsonTest, GenerateAndVerifyJson) {
    std::filesystem::path json_file = output_dir_ / "3G8T.json";
    
    // Generate JSON
    auto residues = generate_residue_ordering_json(test_pdb_, json_file);
    
    EXPECT_GT(residues.size(), 0u) << "Should have at least one residue";
    EXPECT_EQ(residues.size(), 1070u) << "3G8T should have 1070 residues";
    
    // Verify JSON file was created
    EXPECT_TRUE(std::filesystem::exists(json_file)) << "JSON file should be created";
    
    // Load and verify JSON
    auto loaded_residues = load_residue_ordering_json(json_file);
    EXPECT_EQ(loaded_residues.size(), residues.size()) << "Loaded residues should match generated";
    
    // Verify specific known residues
    if (residues.size() >= 946) {
        EXPECT_EQ(residues[945].residue_name, "  C") << "Residue 946 should be C";
        EXPECT_EQ(residues[945].chain_id, "S") << "Residue 946 should be in chain S";
        EXPECT_EQ(residues[945].residue_seq, 113) << "Residue 946 should have seq 113";
    }
    
    if (residues.size() >= 947) {
        EXPECT_EQ(residues[946].residue_name, "  U") << "Residue 947 should be U";
        EXPECT_EQ(residues[946].chain_id, "S") << "Residue 947 should be in chain S";
        EXPECT_EQ(residues[946].residue_seq, 114) << "Residue 947 should have seq 114";
    }
}

/**
 * @test Verify JSON consistency across multiple generations
 */
TEST_F(ResidueOrderingJsonTest, JsonConsistency) {
    std::filesystem::path json_file1 = output_dir_ / "3G8T_test1.json";
    std::filesystem::path json_file2 = output_dir_ / "3G8T_test2.json";
    
    // Generate JSON twice
    auto residues1 = generate_residue_ordering_json(test_pdb_, json_file1);
    auto residues2 = generate_residue_ordering_json(test_pdb_, json_file2);
    
    EXPECT_EQ(residues1.size(), residues2.size()) << "Both generations should have same count";
    
    // Verify all residues match
    for (size_t i = 0; i < residues1.size(); i++) {
        EXPECT_EQ(residues1[i], residues2[i]) 
            << "Residue at index " << i << " should match between generations";
    }
}

/**
 * @test Verify JSON structure is correct
 */
TEST_F(ResidueOrderingJsonTest, JsonStructure) {
    std::filesystem::path json_file = output_dir_ / "3G8T_structure_test.json";
    
    generate_residue_ordering_json(test_pdb_, json_file);
    
    std::ifstream in_file(json_file);
    ASSERT_TRUE(in_file.is_open()) << "JSON file should be readable";
    
    json j;
    in_file >> j;
    
    // Verify structure
    EXPECT_TRUE(j.contains("pdb_id")) << "JSON should have pdb_id";
    EXPECT_TRUE(j.contains("total_residues")) << "JSON should have total_residues";
    EXPECT_TRUE(j.contains("residues")) << "JSON should have residues array";
    
    EXPECT_EQ(j["pdb_id"], "3G8T") << "pdb_id should be 3G8T";
    EXPECT_EQ(j["total_residues"], 1070) << "total_residues should be 1070";
    EXPECT_TRUE(j["residues"].is_array()) << "residues should be an array";
    EXPECT_EQ(j["residues"].size(), 1070u) << "residues array should have 1070 elements";
    
    // Verify first residue structure
    if (j["residues"].size() > 0) {
        const auto& first = j["residues"][0];
        EXPECT_TRUE(first.contains("legacy_index")) << "Residue should have legacy_index";
        EXPECT_TRUE(first.contains("residue_name")) << "Residue should have residue_name";
        EXPECT_TRUE(first.contains("chain_id")) << "Residue should have chain_id";
        EXPECT_TRUE(first.contains("residue_seq")) << "Residue should have residue_seq";
        EXPECT_TRUE(first.contains("insertion_code")) << "Residue should have insertion_code";
        EXPECT_TRUE(first.contains("num_atoms")) << "Residue should have num_atoms";
    }
}


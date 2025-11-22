/**
 * @file test_problem_pdbs.cpp
 * @brief Test frame calculation on known problematic PDBs with insertion codes
 */

#include <gtest/gtest.h>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <nlohmann/json.hpp>
#include <fstream>
#include <vector>
#include <tuple>
#include <set>
#include <filesystem>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;

class ProblemPDBsTest : public ::testing::Test {
protected:
    void SetUp() override {
        calculator_ = std::make_unique<BaseFrameCalculator>("data/templates");
        
        // List of problematic PDBs to test
        problem_pdbs_ = {"8ZYD", "4PWD", "6CAR", "3G96", "4B3M"};
    }

    nlohmann::json load_legacy_json(const std::string& pdb_name) {
        std::filesystem::path json_file = std::filesystem::path("data/json_legacy") / (pdb_name + ".json");
        std::ifstream file(json_file);
        if (!file.is_open()) {
            return nlohmann::json();
        }
        nlohmann::json json;
        file >> json;
        return json;
    }

    std::vector<nlohmann::json> find_records_by_type(const nlohmann::json& json, const std::string& type) {
        std::vector<nlohmann::json> results;
        if (json.contains("calculations")) {
            for (const auto& calc : json["calculations"]) {
                if (calc.contains("type") && calc["type"] == type) {
                    results.push_back(calc);
                }
            }
        }
        return results;
    }

    std::vector<std::tuple<char, int, char, std::string>> build_ordered_residue_list(const nlohmann::json& legacy_json) {
        std::vector<std::tuple<char, int, char, std::string>> ordered_residues;
        if (legacy_json.contains("calculations")) {
            for (const auto& calc : legacy_json["calculations"]) {
                if (calc.contains("type") && calc["type"] == "pdb_atoms" &&
                    calc.contains("atoms") && calc["atoms"].is_array()) {
                    std::set<std::tuple<char, int, char, std::string>> seen;
                    for (const auto& atom : calc["atoms"]) {
                        std::string chain_str = atom.value("chain_id", "");
                        char chain_id = chain_str.empty() ? '\0' : chain_str[0];
                        int seq_num = atom.value("residue_seq", 0);
                        std::string insertion_str = atom.value("insertion", " ");
                        char insertion = insertion_str.empty() ? ' ' : insertion_str[0];
                        std::string res_name = atom.value("residue_name", "");
                        auto key = std::make_tuple(chain_id, seq_num, insertion, res_name);
                        if (seen.find(key) == seen.end()) {
                            ordered_residues.push_back(key);
                            seen.insert(key);
                        }
                    }
                    break;
                }
            }
        }
        return ordered_residues;
    }

    std::optional<const Residue*> find_residue_by_legacy_idx(
        const Structure& structure, 
        size_t legacy_residue_idx,
        const std::vector<std::tuple<char, int, char, std::string>>& ordered_residues) {
        
        if (legacy_residue_idx == 0 || legacy_residue_idx > ordered_residues.size()) {
            return std::nullopt;
        }
        
        auto [legacy_chain, legacy_seq, legacy_insertion, legacy_name] = ordered_residues[legacy_residue_idx - 1];
        
        for (const auto& chain : structure.chains()) {
            if (chain.chain_id() != legacy_chain) continue;
            
            for (const auto& residue : chain.residues()) {
                if (residue.seq_num() == legacy_seq && residue.insertion() == legacy_insertion) {
                    return &residue;
                }
            }
        }
        
        return std::nullopt;
    }

    std::unique_ptr<BaseFrameCalculator> calculator_;
    std::vector<std::string> problem_pdbs_;
};

TEST_F(ProblemPDBsTest, TestResidueMatchingWithInsertionCodes) {
    for (const auto& pdb_name : problem_pdbs_) {
        std::filesystem::path pdb_file = std::filesystem::path("data/pdb") / (pdb_name + ".pdb");
        if (!std::filesystem::exists(pdb_file)) {
            GTEST_SKIP() << "PDB file not found: " << pdb_file;
        }

        // Load PDB
        PdbParser parser;
        Structure structure = parser.parse_file(pdb_file);
        
        // Load legacy JSON
        nlohmann::json legacy_json = load_legacy_json(pdb_name);
        if (legacy_json.empty()) {
            GTEST_SKIP() << "Legacy JSON not found for " << pdb_name;
        }

        // Build ordered residue list
        auto ordered_residues = build_ordered_residue_list(legacy_json);
        
        // Calculate frames
        calculator_->calculate_all_frames(structure);

        // Get ls_fitting records
        auto ls_records = find_records_by_type(legacy_json, "ls_fitting");
        
        size_t matched_residues = 0;
        size_t unmatched_residues = 0;
        size_t residues_with_insertion_codes = 0;

        std::cout << "\n=== Testing " << pdb_name << " ===" << std::endl;
        std::cout << "Total residues in legacy JSON: " << ordered_residues.size() << std::endl;

        for (const auto& ls_record : ls_records) {
            if (!ls_record.contains("residue_idx")) continue;
            
            size_t legacy_residue_idx = ls_record["residue_idx"].get<size_t>();
            if (legacy_residue_idx == 0 || legacy_residue_idx > ordered_residues.size()) continue;
            
            auto [legacy_chain, legacy_seq, legacy_insertion, legacy_name] = ordered_residues[legacy_residue_idx - 1];
            
            // Check for insertion codes
            if (legacy_insertion != ' ') {
                residues_with_insertion_codes++;
            }
            
            auto residue_opt = find_residue_by_legacy_idx(structure, legacy_residue_idx, ordered_residues);
            
            if (residue_opt.has_value()) {
                const Residue* residue_ptr = residue_opt.value();
                
                // Verify insertion code matches
                if (residue_ptr->insertion() == legacy_insertion && 
                    residue_ptr->seq_num() == legacy_seq &&
                    residue_ptr->chain_id() == legacy_chain) {
                    matched_residues++;
                } else {
                    unmatched_residues++;
                    std::cout << "  Mismatch at residue_idx " << legacy_residue_idx 
                              << " (" << legacy_chain << ":" << legacy_seq 
                              << (legacy_insertion != ' ' ? std::string(1, legacy_insertion) : "")
                              << " " << legacy_name << ")" << std::endl;
                    std::cout << "    Legacy: " << legacy_chain << ":" << legacy_seq 
                              << (legacy_insertion != ' ' ? std::string(1, legacy_insertion) : "") << std::endl;
                    std::cout << "    Our:    " << residue_ptr->chain_id() << ":" << residue_ptr->seq_num()
                              << (residue_ptr->insertion() != ' ' ? std::string(1, residue_ptr->insertion()) : "") << std::endl;
                }
            } else {
                unmatched_residues++;
                if (legacy_insertion != ' ') {
                    std::cout << "  Not found: residue_idx " << legacy_residue_idx 
                              << " (" << legacy_chain << ":" << legacy_seq << legacy_insertion
                              << " " << legacy_name << ") [has insertion code]" << std::endl;
                }
            }
        }

        std::cout << "Matched residues: " << matched_residues << std::endl;
        std::cout << "Unmatched residues: " << unmatched_residues << std::endl;
        std::cout << "Residues with insertion codes: " << residues_with_insertion_codes << std::endl;
        
        // Check specifically for 8ZYD C:21 issue
        if (pdb_name == "8ZYD") {
            bool found_c21_blank = false;
            bool found_c21a = false;
            
            for (const auto& chain : structure.chains()) {
                if (chain.chain_id() != 'C') continue;
                for (const auto& residue : chain.residues()) {
                    if (residue.seq_num() == 21) {
                        if (residue.insertion() == ' ') {
                            found_c21_blank = true;
                            std::cout << "  Found C:21 (blank): " << residue.name() 
                                      << " with " << residue.num_atoms() << " atoms" << std::endl;
                        } else if (residue.insertion() == 'A') {
                            found_c21a = true;
                            std::cout << "  Found C:21A (insertion A): " << residue.name() 
                                      << " with " << residue.num_atoms() << " atoms" << std::endl;
                        }
                    }
                }
            }
            
            std::cout << "C:21 (blank) found: " << (found_c21_blank ? "YES" : "NO") << std::endl;
            std::cout << "C:21A (insertion A) found: " << (found_c21a ? "YES" : "NO") << std::endl;
            
            EXPECT_TRUE(found_c21_blank || found_c21a) << "Should find at least one C:21 residue";
        }
    }
}


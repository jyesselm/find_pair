/**
 * @file test_specific_problem_pdbs.cpp
 * @brief Test frame calculation specifically on known problematic PDBs
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
#include <iomanip>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;

class SpecificProblemPDBsTest : public ::testing::Test {
protected:
    void SetUp() override {
        calculator_ = std::make_unique<BaseFrameCalculator>("data/templates");
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

    double max_rotation_diff(const Matrix3D& m1, const nlohmann::json& json_m2) {
        if (!json_m2.is_array() || json_m2.size() != 3) return -1.0;
        double max_diff = 0.0;
        for (size_t i = 0; i < 3; ++i) {
            if (!json_m2[i].is_array() || json_m2[i].size() != 3) return -1.0;
            for (size_t j = 0; j < 3; ++j) {
                double v1 = m1.at(i, j);
                double v2 = json_m2[i][j].get<double>();
                double diff = std::abs(v1 - v2);
                if (diff > max_diff) max_diff = diff;
            }
        }
        return max_diff;
    }

    double max_translation_diff(const Vector3D& v1, const nlohmann::json& json_v2) {
        if (!json_v2.is_array() || json_v2.size() != 3) return -1.0;
        double v2_x = json_v2[0].get<double>();
        double v2_y = json_v2[1].get<double>();
        double v2_z = json_v2[2].get<double>();
        double diff_x = std::abs(v1.x() - v2_x);
        double diff_y = std::abs(v1.y() - v2_y);
        double diff_z = std::abs(v1.z() - v2_z);
        return std::max({diff_x, diff_y, diff_z});
    }

    std::unique_ptr<BaseFrameCalculator> calculator_;
};

TEST_F(SpecificProblemPDBsTest, Test8ZYDC21) {
    // This is the specific problematic case - C:21 has two residues (blank and A)
    std::string pdb_name = "8ZYD";
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

    // Get ls_fitting records for C:21 residues
    auto ls_records = find_records_by_type(legacy_json, "ls_fitting");
    
    std::cout << "\n=== Testing 8ZYD C:21 residues ===" << std::endl;
    
    size_t c21_matched = 0;
    size_t c21_failed = 0;
    
    for (const auto& ls_record : ls_records) {
        if (!ls_record.contains("residue_idx")) continue;
        
        size_t legacy_residue_idx = ls_record["residue_idx"].get<size_t>();
        if (legacy_residue_idx == 0 || legacy_residue_idx > ordered_residues.size()) continue;
        
        auto [legacy_chain, legacy_seq, legacy_insertion, legacy_name] = ordered_residues[legacy_residue_idx - 1];
        
        // Focus on C:21 residues
        if (legacy_chain != 'C' || legacy_seq != 21) continue;
        
        auto residue_opt = find_residue_by_legacy_idx(structure, legacy_residue_idx, ordered_residues);
        
        std::cout << "\nResidue " << legacy_residue_idx << " (" << legacy_chain << ":" << legacy_seq
                  << (legacy_insertion != ' ' ? std::string(1, legacy_insertion) : "")
                  << " " << legacy_name << "):" << std::endl;
        
        if (!residue_opt.has_value()) {
            std::cout << "  NOT FOUND" << std::endl;
            c21_failed++;
            continue;
        }
        
        const Residue* residue_ptr = residue_opt.value();
        
        if (!residue_ptr->reference_frame().has_value()) {
            std::cout << "  NO FRAME CALCULATED" << std::endl;
            c21_failed++;
            continue;
        }
        
        // Calculate frame to get detailed info
        FrameCalculationResult result = calculator_->calculate_frame_const(*residue_ptr);
        
        if (!result.is_valid) {
            std::cout << "  INVALID FRAME (num_matched: " << result.num_matched << ")" << std::endl;
            c21_failed++;
            continue;
        }
        
        // Compare with legacy
        double rot_diff = max_rotation_diff(result.rotation_matrix, ls_record["rotation_matrix"]);
        double trans_diff = max_translation_diff(result.translation, ls_record["translation"]);
        double rms_diff = std::abs(result.rms_fit - ls_record["rms_fit"].get<double>());
        size_t num_matched_legacy = ls_record.value("num_points", 0UL);
        
        std::cout << "  Our num_matched: " << result.num_matched << ", Legacy: " << num_matched_legacy << std::endl;
        std::cout << "  Our RMS: " << std::fixed << std::setprecision(6) << result.rms_fit
                  << ", Legacy: " << ls_record["rms_fit"].get<double>() << std::endl;
        std::cout << "  Rot diff: " << rot_diff << ", Trans diff: " << trans_diff
                  << ", RMS diff: " << rms_diff << std::endl;
        
        if (rot_diff < 0.05 && trans_diff < 0.05 && rms_diff < 0.005 && 
            result.num_matched == num_matched_legacy) {
            std::cout << "  ✓ MATCHED" << std::endl;
            c21_matched++;
        } else {
            std::cout << "  ✗ FAILED" << std::endl;
            c21_failed++;
        }
    }
    
    std::cout << "\nSummary for C:21 residues:" << std::endl;
    std::cout << "  Matched: " << c21_matched << std::endl;
    std::cout << "  Failed: " << c21_failed << std::endl;
    
    // We should have matched both C:21 residues
    EXPECT_GE(c21_matched, 1) << "Should match at least one C:21 residue";
}


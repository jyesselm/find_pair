/**
 * @file check_failing_residues.cpp
 * @brief Check actual differences for failing residues
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <nlohmann/json.hpp>
#include "test_data_discovery.hpp"
#include <cmath>
#include <algorithm>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::test;
using namespace x3dna::geometry;

// Copy comparison functions from test
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

std::vector<std::tuple<char, int, std::string>> build_ordered_residue_list(const nlohmann::json& legacy_json) {
    std::vector<std::tuple<char, int, std::string>> ordered_residues;
    if (legacy_json.contains("calculations")) {
        for (const auto& calc : legacy_json["calculations"]) {
            if (calc.contains("type") && calc["type"] == "pdb_atoms" &&
                calc.contains("atoms") && calc["atoms"].is_array()) {
                std::vector<std::tuple<char, int, std::string>> seen;
                for (const auto& atom : calc["atoms"]) {
                    std::string chain_str = atom.value("chain_id", "");
                    char chain_id = chain_str.empty() ? '\0' : chain_str[0];
                    int seq_num = atom.value("residue_seq", 0);
                    std::string res_name = atom.value("residue_name", "");
                    auto key = std::make_tuple(chain_id, seq_num, res_name);
                    bool found = false;
                    for (const auto& s : seen) {
                        if (s == key) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        ordered_residues.push_back(key);
                        seen.push_back(key);
                    }
                }
                break;
            }
        }
    }
    return ordered_residues;
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

std::optional<const Residue*> find_residue_by_legacy_idx(
    const Structure& structure, 
    size_t legacy_residue_idx,
    const std::vector<std::tuple<char, int, std::string>>& ordered_residues) {
    
    if (legacy_residue_idx == 0 || legacy_residue_idx > ordered_residues.size()) {
        return std::nullopt;
    }
    
    auto [legacy_chain, legacy_seq, legacy_name] = ordered_residues[legacy_residue_idx - 1];
    
    for (const auto& chain : structure.chains()) {
        if (chain.chain_id() != legacy_chain) continue;
        for (const auto& residue : chain.residues()) {
            if (residue.seq_num() == legacy_seq) {
                return &residue;
            }
        }
    }
    
    return std::nullopt;
}

struct DiffStats {
    double max_rot_diff = 0.0;
    double max_trans_diff = 0.0;
    double max_rms_diff = 0.0;
    size_t count = 0;
};

int main() {
    auto pairs = test_data_discovery::discover_pairs();
    if (pairs.empty()) {
        std::cerr << "No PDB/JSON pairs found" << std::endl;
        return 1;
    }
    
    BaseFrameCalculator calculator("data/templates");
    DiffStats stats;
    std::vector<std::tuple<std::string, size_t, double, double, double>> failing;
    
    size_t num_to_test = std::min(pairs.size(), size_t(5));
    
    for (size_t i = 0; i < num_to_test; ++i) {
        const auto& pair = pairs[i];
        
        // Load PDB
        PdbParser parser;
        Structure structure = parser.parse_file(pair.pdb_file);
        
        // Load legacy JSON
        std::ifstream json_file(pair.json_file);
        nlohmann::json legacy_json;
        json_file >> legacy_json;
        
        // Get records
        auto ls_records = find_records_by_type(legacy_json, "ls_fitting");
        auto ordered_residues = build_ordered_residue_list(legacy_json);
        
        // Calculate frames
        calculator.calculate_all_frames(structure);
        
        // Check each residue
        for (const auto& ls_record : ls_records) {
            if (!ls_record.contains("residue_idx")) continue;
            
            size_t legacy_residue_idx = ls_record["residue_idx"].get<size_t>();
            auto residue_opt = find_residue_by_legacy_idx(structure, legacy_residue_idx, ordered_residues);
            if (!residue_opt.has_value()) continue;
            
            const Residue* residue_ptr = residue_opt.value();
            if (!residue_ptr->reference_frame().has_value()) continue;
            
            FrameCalculationResult result = calculator.calculate_frame_const(*residue_ptr);
            if (!result.is_valid) continue;
            
            // Calculate differences
            double rot_diff = 0.0;
            double trans_diff = 0.0;
            double rms_diff = 0.0;
            bool failed = false;
            
            if (ls_record.contains("rotation_matrix")) {
                rot_diff = max_rotation_diff(result.rotation_matrix, ls_record["rotation_matrix"]);
                if (rot_diff > 0.05) failed = true;
            }
            
            if (ls_record.contains("translation")) {
                trans_diff = max_translation_diff(result.translation, ls_record["translation"]);
                if (trans_diff > 0.05) failed = true;
            }
            
            if (ls_record.contains("rms_fit") && !ls_record["rms_fit"].is_null()) {
                double legacy_rms = ls_record["rms_fit"].get<double>();
                rms_diff = std::abs(result.rms_fit - legacy_rms);
                if (rms_diff > 0.005) failed = true;
            }
            
            if (failed) {
                stats.max_rot_diff = std::max(stats.max_rot_diff, rot_diff);
                stats.max_trans_diff = std::max(stats.max_trans_diff, trans_diff);
                stats.max_rms_diff = std::max(stats.max_rms_diff, rms_diff);
                stats.count++;
                
                if (failing.size() < 20) {  // Show first 20
                    failing.push_back({pair.pdb_name, legacy_residue_idx, rot_diff, trans_diff, rms_diff});
                }
            }
        }
    }
    
    std::cout << "\n=== Failing Residue Analysis ===" << std::endl;
    std::cout << "Total failing residues found: " << stats.count << std::endl;
    std::cout << "\nMaximum differences:" << std::endl;
    std::cout << "  Rotation matrix: " << std::fixed << std::setprecision(6) << stats.max_rot_diff << std::endl;
    std::cout << "  Translation: " << stats.max_trans_diff << std::endl;
    std::cout << "  RMS: " << stats.max_rms_diff << std::endl;
    
    std::cout << "\nFirst 20 failing residues:" << std::endl;
    std::cout << std::setw(10) << "PDB" << std::setw(8) << "ResIDX" 
              << std::setw(12) << "Rot Diff" 
              << std::setw(12) << "Trans Diff"
              << std::setw(12) << "RMS Diff" << std::endl;
    std::cout << std::string(54, '-') << std::endl;
    
    for (const auto& [pdb, idx, rot, trans, rms] : failing) {
        std::cout << std::setw(10) << pdb 
                  << std::setw(8) << idx
                  << std::setw(12) << std::fixed << std::setprecision(6) << rot
                  << std::setw(12) << trans
                  << std::setw(12) << rms << std::endl;
    }
    
    return 0;
}


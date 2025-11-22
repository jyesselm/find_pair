/**
 * @file analyze_frame_mismatches.cpp
 * @brief Analyze frame calculation mismatches in detail
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/ring_atom_matcher.hpp>
#include <x3dna/algorithms/standard_base_templates.hpp>
#include <nlohmann/json.hpp>
#include "test_data_discovery.hpp"
#include <filesystem>
#include <cmath>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::test;
using namespace x3dna::geometry;

struct MismatchInfo {
    std::string pdb_name;
    size_t legacy_residue_idx;
    char chain_id;
    int seq_num;
    std::string residue_name;
    std::string base_type;
    size_t our_num_matched;
    size_t legacy_num_matched;
    double our_rms;
    double legacy_rms;
    double rms_diff;
    double max_rot_diff;
    double max_trans_diff;
    std::vector<std::string> our_atoms;
    std::vector<std::string> legacy_atoms;
};

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
                    // Check if already seen
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

int main() {
    auto pairs = test_data_discovery::discover_pairs();
    if (pairs.empty()) {
        std::cerr << "No PDB/JSON pairs found" << std::endl;
        return 1;
    }
    
    BaseFrameCalculator calculator("data/templates");
    std::vector<MismatchInfo> mismatches;
    
    // Analyze first 5 PDBs
    size_t num_to_analyze = std::min(pairs.size(), size_t(5));
    
    for (size_t i = 0; i < num_to_analyze; ++i) {
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
        auto base_frame_records = find_records_by_type(legacy_json, "base_frame_calc");
        
        // Build ordered residue list
        auto ordered_residues = build_ordered_residue_list(legacy_json);
        
        // Calculate frames
        calculator.calculate_all_frames(structure);
        
        // Analyze each residue
        for (const auto& ls_record : ls_records) {
            if (!ls_record.contains("residue_idx")) continue;
            
            size_t legacy_residue_idx = ls_record["residue_idx"].get<size_t>();
            if (legacy_residue_idx == 0 || legacy_residue_idx > ordered_residues.size()) continue;
            
            auto [legacy_chain, legacy_seq, legacy_name] = ordered_residues[legacy_residue_idx - 1];
            
            // Find residue in structure
            const Residue* residue_ptr = nullptr;
            for (const auto& chain : structure.chains()) {
                if (chain.chain_id() != legacy_chain) continue;
                for (const auto& residue : chain.residues()) {
                    if (residue.seq_num() == legacy_seq) {
                        residue_ptr = &residue;
                        break;
                    }
                }
                if (residue_ptr) break;
            }
            
            if (!residue_ptr || !residue_ptr->reference_frame().has_value()) continue;
            
            // Calculate frame to get metrics
            FrameCalculationResult result = calculator.calculate_frame_const(*residue_ptr);
            if (!result.is_valid) continue;
            
            // Find corresponding base_frame_calc record
            const nlohmann::json* base_frame_record = nullptr;
            for (const auto& bf : base_frame_records) {
                if (bf.contains("residue_idx") && 
                    bf["residue_idx"].get<size_t>() == legacy_residue_idx) {
                    base_frame_record = &bf;
                    break;
                }
            }
            
            // Compare values (use same tolerances as test)
            bool mismatch = false;
            MismatchInfo info;
            info.pdb_name = pair.pdb_name;
            info.legacy_residue_idx = legacy_residue_idx;
            info.chain_id = legacy_chain;
            info.seq_num = legacy_seq;
            info.residue_name = legacy_name;
            info.our_num_matched = result.num_matched;
            info.our_atoms = result.matched_atoms;
            info.our_rms = result.rms_fit;
            
            // Compare with ls_fitting record
            if (ls_record.contains("rms_fit") && !ls_record["rms_fit"].is_null()) {
                info.legacy_rms = ls_record["rms_fit"].get<double>();
                info.rms_diff = std::abs(result.rms_fit - info.legacy_rms);
                if (info.rms_diff > 0.001) mismatch = true;  // RMS tolerance: 0.001
            }
            
            if (ls_record.contains("num_points")) {
                info.legacy_num_matched = ls_record["num_points"].get<size_t>();
                if (info.our_num_matched != info.legacy_num_matched) mismatch = true;
            }
            
            if (ls_record.contains("rotation_matrix")) {
                info.max_rot_diff = max_rotation_diff(result.rotation_matrix, ls_record["rotation_matrix"]);
                if (info.max_rot_diff > 0.01) mismatch = true;  // Rotation tolerance: 0.01
            }
            
            if (ls_record.contains("translation")) {
                info.max_trans_diff = max_translation_diff(result.translation, ls_record["translation"]);
                if (info.max_trans_diff > 0.01) mismatch = true;  // Translation tolerance: 0.01
            }
            
            // Always record info if there are ANY differences (even small ones)
            if (info.rms_diff > 0.0 || info.max_rot_diff > 0.0 || info.max_trans_diff > 0.0 ||
                info.our_num_matched != info.legacy_num_matched) {
                mismatch = true;
            }
            
            if (base_frame_record) {
                if (base_frame_record->contains("base_type")) {
                    info.base_type = base_frame_record->value("base_type", "");
                }
                if (base_frame_record->contains("matched_atoms")) {
                    info.legacy_atoms = base_frame_record->value("matched_atoms", std::vector<std::string>());
                }
                if (base_frame_record->contains("num_matched_atoms")) {
                    info.legacy_num_matched = base_frame_record->value("num_matched_atoms", 0UL);
                }
            }
            
            if (mismatch) {
                mismatches.push_back(info);
            }
        }
    }
    
    // Report mismatches
    std::cout << "\n=== Frame Calculation Mismatch Analysis ===" << std::endl;
    std::cout << "Total mismatches found: " << mismatches.size() << std::endl;
    std::cout << "\nMismatch Categories:\n" << std::endl;
    
    // Group by mismatch type
    size_t num_matched_atom_diff = 0;
    size_t num_rms_diff = 0;
    size_t num_rot_diff = 0;
    size_t num_trans_diff = 0;
    
    std::map<std::string, size_t> base_type_mismatches;
    std::map<std::pair<size_t, size_t>, size_t> atom_count_mismatches;
    
    for (const auto& mismatch : mismatches) {
        if (mismatch.our_num_matched != mismatch.legacy_num_matched) {
            num_matched_atom_diff++;
            atom_count_mismatches[{mismatch.our_num_matched, mismatch.legacy_num_matched}]++;
        }
        if (mismatch.rms_diff > 0.001) num_rms_diff++;
        if (mismatch.max_rot_diff > 0.01) num_rot_diff++;
        if (mismatch.max_trans_diff > 0.01) num_trans_diff++;
        
        base_type_mismatches[mismatch.base_type]++;
    }
    
    std::cout << "  - Number of matched atoms differs: " << num_matched_atom_diff << std::endl;
    std::cout << "  - RMS differs (>0.001): " << num_rms_diff << std::endl;
    std::cout << "  - Rotation matrix differs (>0.01): " << num_rot_diff << std::endl;
    std::cout << "  - Translation differs (>0.01): " << num_trans_diff << std::endl;
    
    std::cout << "\nMismatches by base type:\n" << std::endl;
    for (const auto& [base_type, count] : base_type_mismatches) {
        std::cout << "  " << base_type << ": " << count << std::endl;
    }
    
    std::cout << "\nAtom count differences:\n" << std::endl;
    for (const auto& [pair, count] : atom_count_mismatches) {
        std::cout << "  Our: " << pair.first << " vs Legacy: " << pair.second << " -> " << count << " residues" << std::endl;
    }
    
    // Show first 10 detailed mismatches
    std::cout << "\n=== First 10 Detailed Mismatches ===\n" << std::endl;
    size_t num_to_show = std::min(mismatches.size(), size_t(10));
    for (size_t i = 0; i < num_to_show; ++i) {
        const auto& m = mismatches[i];
        std::cout << i+1 << ". " << m.pdb_name 
                  << " residue_idx " << m.legacy_residue_idx
                  << " (" << m.chain_id << ":" << m.seq_num << " " << m.residue_name << ")"
                  << " base_type: " << m.base_type << std::endl;
        
        if (m.our_num_matched != m.legacy_num_matched) {
            std::cout << "   Matched atoms: Our=" << m.our_num_matched 
                      << ", Legacy=" << m.legacy_num_matched << std::endl;
        }
        
        if (m.rms_diff > 0.001) {
            std::cout << "   RMS: Our=" << std::fixed << std::setprecision(6) << m.our_rms
                      << ", Legacy=" << m.legacy_rms
                      << ", Diff=" << m.rms_diff << std::endl;
        }
        
        if (m.max_rot_diff > 0.01) {
            std::cout << "   Max rotation diff: " << std::fixed << std::setprecision(6) 
                      << m.max_rot_diff << std::endl;
        }
        
        if (m.max_trans_diff > 0.01) {
            std::cout << "   Max translation diff: " << std::fixed << std::setprecision(6) 
                      << m.max_trans_diff << std::endl;
        }
        
        if (m.our_atoms != m.legacy_atoms) {
            std::cout << "   Our atoms: [";
            for (const auto& atom : m.our_atoms) std::cout << atom << " ";
            std::cout << "]" << std::endl;
            std::cout << "   Legacy atoms: [";
            for (const auto& atom : m.legacy_atoms) std::cout << atom << " ";
            std::cout << "]" << std::endl;
        }
        std::cout << std::endl;
    }
    
    return 0;
}


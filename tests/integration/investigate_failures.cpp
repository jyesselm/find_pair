/**
 * @file investigate_failures.cpp
 * @brief Investigate why test reports failures when calculations match
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <nlohmann/json.hpp>
#include "test_data_discovery.hpp"
#include <cmath>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::test;

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

enum class FailureReason {
    RESIDUE_NOT_FOUND,
    NO_FRAME,
    INVALID_CALCULATION,
    COMPARISON_FAILED,
    MATCHED
};

struct FailureAnalysis {
    size_t legacy_residue_idx;
    char chain_id;
    int seq_num;
    std::string residue_name;
    FailureReason reason;
    std::string details;
};

int main() {
    auto pairs = test_data_discovery::discover_pairs();
    if (pairs.empty()) {
        std::cerr << "No PDB/JSON pairs found" << std::endl;
        return 1;
    }
    
    BaseFrameCalculator calculator("data/templates");
    
    size_t num_to_test = pairs.size();  // Test all PDBs
    std::vector<FailureAnalysis> failures;
    size_t total_checked = 0;
    size_t matched = 0;
    
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
        
        // Analyze each residue
        for (const auto& ls_record : ls_records) {
            if (!ls_record.contains("residue_idx")) continue;
            
            total_checked++;
            size_t legacy_residue_idx = ls_record["residue_idx"].get<size_t>();
            
            auto [legacy_chain, legacy_seq, legacy_name] = ordered_residues[legacy_residue_idx - 1];
            
            FailureAnalysis analysis;
            analysis.legacy_residue_idx = legacy_residue_idx;
            analysis.chain_id = legacy_chain;
            analysis.seq_num = legacy_seq;
            analysis.residue_name = legacy_name;
            
            auto residue_opt = find_residue_by_legacy_idx(structure, legacy_residue_idx, ordered_residues);
            
            if (!residue_opt.has_value()) {
                analysis.reason = FailureReason::RESIDUE_NOT_FOUND;
                analysis.details = "Residue not found in structure";
                failures.push_back(analysis);
                continue;
            }
            
            const Residue* residue_ptr = residue_opt.value();
            
            if (!residue_ptr->reference_frame().has_value()) {
                analysis.reason = FailureReason::NO_FRAME;
                
                // Investigate why no frame
                ResidueType res_type = residue_ptr->residue_type();
                std::string details = "Residue type: ";
                if (res_type == ResidueType::UNKNOWN) details += "UNKNOWN";
                else if (res_type == ResidueType::AMINO_ACID) details += "AMINO_ACID";
                else if (res_type == ResidueType::ADENINE) details += "ADENINE";
                else if (res_type == ResidueType::CYTOSINE) details += "CYTOSINE";
                else if (res_type == ResidueType::GUANINE) details += "GUANINE";
                else if (res_type == ResidueType::THYMINE) details += "THYMINE";
                else if (res_type == ResidueType::URACIL) details += "URACIL";
                else details += "OTHER";
                
                // Try to calculate frame to see why it fails
                FrameCalculationResult test_result = calculator.calculate_frame_const(*residue_ptr);
                if (!test_result.is_valid) {
                    details += "; Invalid calculation (num_matched=" + std::to_string(test_result.num_matched) + ")";
                } else {
                    details += "; Frame calculation succeeded but not stored in residue";
                }
                
                analysis.details = details;
                failures.push_back(analysis);
                continue;
            }
            
            FrameCalculationResult result = calculator.calculate_frame_const(*residue_ptr);
            
            if (!result.is_valid) {
                analysis.reason = FailureReason::INVALID_CALCULATION;
                analysis.details = "Frame calculation is invalid";
                failures.push_back(analysis);
                continue;
            }
            
            // Check if comparison would fail
            bool failed = false;
            std::string failure_details;
            
            if (ls_record.contains("rotation_matrix")) {
                auto leg_rot = ls_record["rotation_matrix"];
                double max_diff = 0.0;
                for (size_t i = 0; i < 3; ++i) {
                    for (size_t j = 0; j < 3; ++j) {
                        double our_val = result.rotation_matrix.at(i, j);
                        double leg_val = leg_rot[i][j].get<double>();
                        double diff = std::abs(our_val - leg_val);
                        if (diff > max_diff) max_diff = diff;
                    }
                }
                if (max_diff > 0.05) {
                    failed = true;
                    failure_details += "Rotation diff: " + std::to_string(max_diff) + "; ";
                }
            }
            
            if (ls_record.contains("translation")) {
                auto leg_trans = ls_record["translation"];
                double max_diff = 0.0;
                for (size_t i = 0; i < 3; ++i) {
                    double our_val = (i == 0) ? result.translation.x() : 
                                    (i == 1) ? result.translation.y() : result.translation.z();
                    double leg_val = leg_trans[i].get<double>();
                    double diff = std::abs(our_val - leg_val);
                    if (diff > max_diff) max_diff = diff;
                }
                if (max_diff > 0.05) {
                    failed = true;
                    failure_details += "Translation diff: " + std::to_string(max_diff) + "; ";
                }
            }
            
            if (ls_record.contains("rms_fit") && !ls_record["rms_fit"].is_null()) {
                double legacy_rms = ls_record["rms_fit"].get<double>();
                double diff = std::abs(result.rms_fit - legacy_rms);
                if (diff > 0.005) {
                    failed = true;
                    failure_details += "RMS diff: " + std::to_string(diff) + "; ";
                }
            }
            
            if (failed) {
                analysis.reason = FailureReason::COMPARISON_FAILED;
                analysis.details = failure_details;
                failures.push_back(analysis);
            } else {
                matched++;
            }
        }
    }
    
    std::cout << "\n=== Failure Analysis ===" << std::endl;
    std::cout << "Total residues checked: " << total_checked << std::endl;
    std::cout << "Matched: " << matched << std::endl;
    std::cout << "Failed: " << failures.size() << std::endl;
    
    // Group failures by reason
    size_t not_found = 0, no_frame = 0, invalid = 0, comparison_failed = 0;
    for (const auto& f : failures) {
        switch (f.reason) {
            case FailureReason::RESIDUE_NOT_FOUND: not_found++; break;
            case FailureReason::NO_FRAME: no_frame++; break;
            case FailureReason::INVALID_CALCULATION: invalid++; break;
            case FailureReason::COMPARISON_FAILED: comparison_failed++; break;
            default: break;
        }
    }
    
    std::cout << "\nFailure breakdown:" << std::endl;
    std::cout << "  Residue not found: " << not_found << std::endl;
    std::cout << "  No frame calculated: " << no_frame << std::endl;
    std::cout << "  Invalid calculation: " << invalid << std::endl;
    std::cout << "  Comparison failed: " << comparison_failed << std::endl;
    
    // Show first 20 failures with details
    std::cout << "\nFirst 20 failures with details:" << std::endl;
    size_t shown = 0;
    for (const auto& f : failures) {
        if (shown < 20) {
            std::string reason_str;
            switch (f.reason) {
                case FailureReason::RESIDUE_NOT_FOUND: reason_str = "RESIDUE_NOT_FOUND"; break;
                case FailureReason::NO_FRAME: reason_str = "NO_FRAME"; break;
                case FailureReason::INVALID_CALCULATION: reason_str = "INVALID_CALCULATION"; break;
                case FailureReason::COMPARISON_FAILED: reason_str = "COMPARISON_FAILED"; break;
                default: reason_str = "UNKNOWN"; break;
            }
            std::cout << "  Residue " << f.legacy_residue_idx 
                      << " (" << f.chain_id << ":" << f.seq_num << " " << f.residue_name << "): "
                      << reason_str << " - " << f.details << std::endl;
            shown++;
        }
    }
    
    return 0;
}


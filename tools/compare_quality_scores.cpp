/**
 * @file compare_quality_scores.cpp
 * @brief Compare quality scores between legacy and modern for specific pairs
 * 
 * This tool helps debug quality score differences that cause pair selection mismatches.
 */

#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <iomanip>
#include <vector>
#include <map>

using json = nlohmann::json;

struct QualityScoreInfo {
    int residue1_idx;
    int residue2_idx;
    double base_score;          // dorg + 2.0*d_v + plane_angle/20.0
    double adjust_pairQuality;   // Adjustment from H-bonds
    int bp_type_id;              // bp_type_id value
    double final_score;          // Final adjusted quality_score
    int num_good_hb;             // Number of good H-bonds (distance in [2.5, 3.5])
    int num_total_hb;            // Total number of H-bonds
    bool is_valid;
    bool is_selected;
};

QualityScoreInfo extract_modern_quality(const json& validation_records, int idx1, int idx2) {
    QualityScoreInfo info;
    info.residue1_idx = idx1;
    info.residue2_idx = idx2;
    info.base_score = 0.0;
    info.adjust_pairQuality = 0.0;
    info.bp_type_id = 0;
    info.final_score = 0.0;
    info.num_good_hb = 0;
    info.num_total_hb = 0;
    info.is_valid = false;
    info.is_selected = false;
    
    // Find matching validation record
    for (const auto& record : validation_records) {
        int r1 = -1, r2 = -1;
        // Try both field name formats
        if (record.contains("residue1_idx") && record.contains("residue2_idx")) {
            r1 = record["residue1_idx"];
            r2 = record["residue2_idx"];
        } else if (record.contains("base_i") && record.contains("base_j")) {
            r1 = record["base_i"];
            r2 = record["base_j"];
        }
        
        if (r1 > 0 && r2 > 0 && ((r1 == idx1 && r2 == idx2) || (r1 == idx2 && r2 == idx1))) {
            // Found matching record
            // Extract rtn_val from calculated_values or directly
            std::array<double, 5> rtn_val = {0.0, 0.0, 0.0, 0.0, 0.0};
            if (record.contains("calculated_values") && record["calculated_values"].is_object()) {
                auto calc = record["calculated_values"];
                rtn_val[0] = calc.value("dorg", 0.0);
                rtn_val[1] = calc.value("d_v", 0.0);
                rtn_val[2] = calc.value("plane_angle", 0.0);
                rtn_val[3] = calc.value("dNN", 0.0);
                rtn_val[4] = calc.value("quality_score", 0.0);
            } else if (record.contains("rtn_val") && record["rtn_val"].is_array() && record["rtn_val"].size() >= 5) {
                auto rtn_val_array = record["rtn_val"];
                for (size_t i = 0; i < 5 && i < rtn_val_array.size(); i++) {
                    rtn_val[i] = rtn_val_array[i].get<double>();
                }
            }
            
            info.base_score = rtn_val[0] + 2.0 * rtn_val[1] + rtn_val[2] / 20.0;
            info.final_score = rtn_val[4]; // rtn_val[5] is index 4 (0-based)
            info.bp_type_id = record.value("bp_type_id", 0);
            
            // Calculate adjust_pairQuality: final_score = base_score + adjust_pairQuality - (bp_type_id == 2 ? 2.0 : 0.0)
            // So: adjust_pairQuality = final_score - base_score + (bp_type_id == 2 ? 2.0 : 0.0)
            info.adjust_pairQuality = info.final_score - info.base_score;
            if (info.bp_type_id == 2) {
                info.adjust_pairQuality += 2.0; // Adjust back to get true adjust_pairQuality
            }
            
            // Handle is_valid as either boolean or number
            if (record.contains("is_valid")) {
                if (record["is_valid"].is_boolean()) {
                    info.is_valid = record["is_valid"].get<bool>();
                } else if (record["is_valid"].is_number()) {
                    info.is_valid = (record["is_valid"].get<int>() != 0);
                }
            }
            break;
        }
    }
    
    return info;
}

QualityScoreInfo extract_legacy_quality(const json& validation_records, int idx1, int idx2) {
    QualityScoreInfo info;
    info.residue1_idx = idx1;
    info.residue2_idx = idx2;
    info.base_score = 0.0;
    info.adjust_pairQuality = 0.0;
    info.bp_type_id = 0;
    info.final_score = 0.0;
    info.num_good_hb = 0;
    info.num_total_hb = 0;
    info.is_valid = false;
    info.is_selected = false;
    
    // Find matching validation record
    for (const auto& record : validation_records) {
        int r1 = -1, r2 = -1;
        // Try both field name formats
        if (record.contains("residue1_idx") && record.contains("residue2_idx")) {
            r1 = record["residue1_idx"];
            r2 = record["residue2_idx"];
        } else if (record.contains("base_i") && record.contains("base_j")) {
            r1 = record["base_i"];
            r2 = record["base_j"];
        }
        
        if (r1 > 0 && r2 > 0 && ((r1 == idx1 && r2 == idx2) || (r1 == idx2 && r2 == idx1))) {
            // Found matching record
            // Extract rtn_val from calculated_values or directly
            std::array<double, 5> rtn_val = {0.0, 0.0, 0.0, 0.0, 0.0};
            if (record.contains("calculated_values") && record["calculated_values"].is_object()) {
                auto calc = record["calculated_values"];
                rtn_val[0] = calc.value("dorg", 0.0);
                rtn_val[1] = calc.value("d_v", 0.0);
                rtn_val[2] = calc.value("plane_angle", 0.0);
                rtn_val[3] = calc.value("dNN", 0.0);
                rtn_val[4] = calc.value("quality_score", 0.0);
            } else if (record.contains("rtn_val") && record["rtn_val"].is_array() && record["rtn_val"].size() >= 5) {
                auto rtn_val_array = record["rtn_val"];
                for (size_t i = 0; i < 5 && i < rtn_val_array.size(); i++) {
                    rtn_val[i] = rtn_val_array[i].get<double>();
                }
            }
            
            info.base_score = rtn_val[0] + 2.0 * rtn_val[1] + rtn_val[2] / 20.0;
            info.final_score = rtn_val[4]; // rtn_val[5] is index 4 (0-based)
            info.bp_type_id = record.value("bp_type_id", 0);
            
            // Calculate adjust_pairQuality: final_score = base_score + adjust_pairQuality - (bp_type_id == 2 ? 2.0 : 0.0)
            // So: adjust_pairQuality = final_score - base_score + (bp_type_id == 2 ? 2.0 : 0.0)
            info.adjust_pairQuality = info.final_score - info.base_score;
            if (info.bp_type_id == 2) {
                info.adjust_pairQuality += 2.0; // Adjust back to get true adjust_pairQuality
            }
            
            // Handle is_valid as either boolean or number
            if (record.contains("is_valid")) {
                if (record["is_valid"].is_boolean()) {
                    info.is_valid = record["is_valid"].get<bool>();
                } else if (record["is_valid"].is_number()) {
                    info.is_valid = (record["is_valid"].get<int>() != 0);
                }
            }
            break;
        }
    }
    
    return info;
}

void print_quality_comparison(const QualityScoreInfo& modern, const QualityScoreInfo& legacy) {
    std::cout << "\n========================================\n";
    std::cout << "Quality Score Comparison\n";
    std::cout << "========================================\n";
    std::cout << "Pair: (" << modern.residue1_idx << ", " << modern.residue2_idx << ")\n\n";
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Base Score (dorg + 2.0*d_v + plane_angle/20.0):\n";
    std::cout << "  Modern: " << modern.base_score << "\n";
    std::cout << "  Legacy: " << legacy.base_score << "\n";
    std::cout << "  Difference: " << (modern.base_score - legacy.base_score) << "\n\n";
    
    std::cout << "adjust_pairQuality:\n";
    std::cout << "  Modern: " << modern.adjust_pairQuality << "\n";
    std::cout << "  Legacy: " << legacy.adjust_pairQuality << "\n";
    std::cout << "  Difference: " << (modern.adjust_pairQuality - legacy.adjust_pairQuality) << "\n\n";
    
    std::cout << "bp_type_id:\n";
    std::cout << "  Modern: " << modern.bp_type_id << "\n";
    std::cout << "  Legacy: " << legacy.bp_type_id << "\n";
    std::cout << "  Match: " << (modern.bp_type_id == legacy.bp_type_id ? "✓" : "✗") << "\n\n";
    
    std::cout << "Final Quality Score (rtn_val[5]):\n";
    std::cout << "  Modern: " << modern.final_score << "\n";
    std::cout << "  Legacy: " << legacy.final_score << "\n";
    std::cout << "  Difference: " << (modern.final_score - legacy.final_score) << "\n\n";
    
    std::cout << "Is Valid:\n";
    std::cout << "  Modern: " << (modern.is_valid ? "yes" : "no") << "\n";
    std::cout << "  Legacy: " << (legacy.is_valid ? "yes" : "no") << "\n";
    std::cout << "  Match: " << (modern.is_valid == legacy.is_valid ? "✓" : "✗") << "\n\n";
    
    if (std::abs(modern.final_score - legacy.final_score) > 0.001) {
        std::cout << "⚠️  QUALITY SCORE MISMATCH!\n";
        if (std::abs(modern.adjust_pairQuality - legacy.adjust_pairQuality) > 0.001) {
            std::cout << "   Root cause: adjust_pairQuality difference\n";
            std::cout << "   This suggests H-bond detection differences for quality adjustment\n";
        }
    } else {
        std::cout << "✓ Quality scores match\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_id> <residue1_idx> <residue2_idx>\n";
        std::cerr << "Example: " << argv[0] << " 3G8T 92 160\n";
        std::cerr << "Example: " << argv[0] << " 6CAQ 75 78\n";
        return 1;
    }
    
    std::string pdb_id = argv[1];
    int idx1 = std::stoi(argv[2]);
    int idx2 = std::stoi(argv[3]);
    
    std::filesystem::path modern_file = "data/json/" + pdb_id + "_pair_validation.json";
    std::filesystem::path legacy_file = "data/json_legacy/" + pdb_id + "_pair_validation.json";
    
    if (!std::filesystem::exists(modern_file)) {
        std::cerr << "Error: Modern JSON not found: " << modern_file << "\n";
        return 1;
    }
    
    if (!std::filesystem::exists(legacy_file)) {
        std::cerr << "Error: Legacy JSON not found: " << legacy_file << "\n";
        return 1;
    }
    
    try {
        // Load modern validation records
        std::ifstream modern_stream(modern_file);
        json modern_data = json::parse(modern_stream);
        
        // Load legacy validation records
        std::ifstream legacy_stream(legacy_file);
        json legacy_data = json::parse(legacy_stream);
        
        // Extract quality scores
        QualityScoreInfo modern_info = extract_modern_quality(modern_data, idx1, idx2);
        QualityScoreInfo legacy_info = extract_legacy_quality(legacy_data, idx1, idx2);
        
        // Print comparison
        print_quality_comparison(modern_info, legacy_info);
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}


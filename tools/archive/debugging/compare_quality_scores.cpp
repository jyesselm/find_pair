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
    double base_score;         // dorg + 2.0*d_v + plane_angle/20.0
    double adjust_pairQuality; // Adjustment from H-bonds
    int bp_type_id;            // bp_type_id value
    double final_score;        // Final adjusted quality_score
    int num_good_hb;           // Number of good H-bonds (distance in [2.5, 3.5])
    int num_total_hb;          // Total number of H-bonds
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
    // CRITICAL: base_i/base_j are legacy indices (1-based) - prefer these
    // residue1_idx/residue2_idx may also be legacy indices, but base_i/base_j are guaranteed
    for (const auto& record : validation_records) {
        int r1 = -1, r2 = -1;
        // Prefer base_i/base_j (these are always legacy indices)
        if (record.contains("base_i") && record.contains("base_j")) {
            r1 = record["base_i"];
            r2 = record["base_j"];
        } else if (record.contains("residue1_idx") && record.contains("residue2_idx")) {
            r1 = record["residue1_idx"];
            r2 = record["residue2_idx"];
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

            // Note: quality_score in calculated_values is the BASE score, not the final adjusted
            // score
            info.base_score = rtn_val[4]; // This is dorg + 2.0*d_v + plane_angle/20.0 (BASE score)
            info.bp_type_id = record.value("bp_type_id", 0);

            // The pair_validation record only contains the BASE quality_score
            // The FINAL adjusted score = base_score + adjust_pairQuality - (bp_type_id == 2 ? 2.0 :
            // 0.0) We don't have adjust_pairQuality or final_score in pair_validation records These
            // would need to come from find_bestpair_selection or be calculated from H-bonds
            info.final_score = info.base_score; // Placeholder - actual final score not in pair_validation
            info.adjust_pairQuality = 0.0;      // Would need to be calculated from H-bonds or found elsewhere

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
    // CRITICAL: base_i/base_j are legacy indices (1-based) - prefer these
    // residue1_idx/residue2_idx may also be legacy indices, but base_i/base_j are guaranteed
    for (const auto& record : validation_records) {
        int r1 = -1, r2 = -1;
        // Prefer base_i/base_j (these are always legacy indices)
        if (record.contains("base_i") && record.contains("base_j")) {
            r1 = record["base_i"];
            r2 = record["base_j"];
        } else if (record.contains("residue1_idx") && record.contains("residue2_idx")) {
            r1 = record["residue1_idx"];
            r2 = record["residue2_idx"];
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

            // Note: quality_score in calculated_values is the BASE score, not the final adjusted
            // score
            info.base_score = rtn_val[4]; // This is dorg + 2.0*d_v + plane_angle/20.0 (BASE score)
            info.bp_type_id = record.value("bp_type_id", 0);

            // The pair_validation record only contains the BASE quality_score
            // The FINAL adjusted score = base_score + adjust_pairQuality - (bp_type_id == 2 ? 2.0 :
            // 0.0) We don't have adjust_pairQuality or final_score in pair_validation records These
            // would need to come from find_bestpair_selection or be calculated from H-bonds
            info.final_score = info.base_score; // Placeholder - actual final score not in pair_validation
            info.adjust_pairQuality = 0.0;      // Would need to be calculated from H-bonds or found elsewhere

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

    std::cout << "Final Quality Score (BASE score from pair_validation):\n";
    std::cout << "  Modern: " << modern.final_score << "\n";
    std::cout << "  Legacy: " << legacy.final_score << "\n";
    std::cout << "  Difference: " << (modern.final_score - legacy.final_score) << "\n";
    std::cout << "  Note: This is the BASE score, not the final adjusted score used for pair "
                 "selection\n\n";

    std::cout << "Is Valid:\n";
    std::cout << "  Modern: " << (modern.is_valid ? "yes" : "no") << "\n";
    std::cout << "  Legacy: " << (legacy.is_valid ? "yes" : "no") << "\n";
    std::cout << "  Match: " << (modern.is_valid == legacy.is_valid ? "✓" : "✗") << "\n\n";

    if (std::abs(modern.base_score - legacy.base_score) > 0.001) {
        std::cout << "⚠️  BASE QUALITY SCORE MISMATCH!\n";
        std::cout << "   This suggests differences in geometric calculations (dorg, d_v, plane_angle)\n";
    } else {
        std::cout << "✓ Base quality scores match\n";
    }

    std::cout << "\n";
    std::cout << "NOTE: To get the FINAL adjusted quality score used for pair selection,\n";
    std::cout << "      we need to look at find_bestpair_selection records or calculate\n";
    std::cout << "      adjust_pairQuality from H-bonds (good H-bonds in [2.5, 3.5] range).\n";
}

std::filesystem::path find_json_file(const std::string& pdb_id, bool is_legacy) {
    // Try segmented directory structure first (new format)
    std::filesystem::path base_dir = is_legacy ? "data/json_legacy" : "data/json";
    std::filesystem::path segmented_file = base_dir / "pair_validation" / (pdb_id + ".json");

    if (std::filesystem::exists(segmented_file)) {
        return segmented_file;
    }

    // Fall back to old format with suffix
    std::filesystem::path suffix_file = base_dir / (pdb_id + "_pair_validation.json");

    if (std::filesystem::exists(suffix_file)) {
        return suffix_file;
    }

    // Return the old format path even if it doesn't exist (for error message)
    return suffix_file;
}

json load_json_array(const std::filesystem::path& file_path) {
    std::ifstream file_stream(file_path);
    if (!file_stream.is_open()) {
        throw std::runtime_error("Could not open file: " + file_path.string());
    }

    json data = json::parse(file_stream);

    // Handle both array format and wrapped format
    if (data.is_array()) {
        return data;
    } else if (data.contains("calculations")) {
        auto calc = data["calculations"];
        if (calc.is_array()) {
            json records;
            for (const auto& item : calc) {
                if (item.contains("type") && item["type"] == "pair_validation") {
                    records.push_back(item);
                }
            }
            return records;
        } else if (calc.is_object() && calc.contains("pair_validation")) {
            return calc["pair_validation"];
        }
    }

    // If we get here, return empty array
    return json::array();
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_id> <legacy_residue1_idx> <legacy_residue2_idx>\n";
        std::cerr << "  Note: Indices must be legacy indices (1-based) from legacy JSON files\n";
        std::cerr << "Example: " << argv[0] << " 3G8T 946 947\n";
        std::cerr << "Example: " << argv[0] << " 6CAQ 75 78\n";
        return 1;
    }

    std::string pdb_id = argv[1];
    int idx1 = std::stoi(argv[2]);
    int idx2 = std::stoi(argv[3]);

    std::filesystem::path modern_file = find_json_file(pdb_id, false);
    std::filesystem::path legacy_file = find_json_file(pdb_id, true);

    if (!std::filesystem::exists(modern_file)) {
        std::cerr << "Error: Modern JSON not found: " << modern_file << "\n";
        std::cerr << "  Tried: data/json/pair_validation/" << pdb_id << ".json\n";
        std::cerr << "         data/json/" << pdb_id << "_pair_validation.json\n";
        return 1;
    }

    if (!std::filesystem::exists(legacy_file)) {
        std::cerr << "Error: Legacy JSON not found: " << legacy_file << "\n";
        std::cerr << "  Tried: data/json_legacy/pair_validation/" << pdb_id << ".json\n";
        std::cerr << "         data/json_legacy/" << pdb_id << "_pair_validation.json\n";
        return 1;
    }

    try {
        // Load modern validation records
        json modern_data = load_json_array(modern_file);

        // Load legacy validation records
        json legacy_data = load_json_array(legacy_file);

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

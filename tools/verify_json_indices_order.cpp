/**
 * @file verify_json_indices_order.cpp
 * @brief Verify JSON files have correct legacy indices and are in correct order
 * 
 * This tool checks:
 * 1. All records have legacy_residue_idx (or base_i/base_j for pairs)
 * 2. Records are in legacy index order (1, 2, 3, ...)
 * 3. Indices match between legacy and modern JSON
 * 
 * Usage: verify_json_indices_order <pdb_id> [record_type]
 *   record_type: frame_calc, base_frame_calc, pair_validation, find_bestpair_selection (default: all)
 */

#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <vector>
#include <map>
#include <set>
#include <iomanip>

using json = nlohmann::json;

struct IndexInfo {
    int legacy_idx;
    std::string residue_name;
    char chain_id;
    int residue_seq;
    char insertion;
    size_t position; // Position in JSON array
};

bool verify_frame_records(const json& records, const std::string& source) {
    std::cout << "\n=== Verifying " << source << " frame records ===\n";
    
    std::vector<IndexInfo> indices;
    bool has_errors = false;
    
    for (size_t i = 0; i < records.size(); i++) {
        const auto& record = records[i];
        
        if (!record.contains("legacy_residue_idx")) {
            std::cerr << "ERROR: Record at position " << i << " missing legacy_residue_idx\n";
            has_errors = true;
            continue;
        }
        
        int legacy_idx = record["legacy_residue_idx"];
        IndexInfo info;
        info.legacy_idx = legacy_idx;
        info.residue_name = record.value("residue_name", "");
        std::string chain_str = record.value("chain_id", " ");
        info.chain_id = chain_str.empty() ? ' ' : chain_str[0];
        info.residue_seq = record.value("residue_seq", 0);
        std::string ins_str = record.value("insertion", " ");
        info.insertion = ins_str.empty() ? ' ' : ins_str[0];
        info.position = i;
        
        indices.push_back(info);
    }
    
    // Check order
    std::cout << "Checking order...\n";
    int prev_idx = 0;
    for (const auto& info : indices) {
        if (info.legacy_idx <= prev_idx) {
            std::cerr << "ERROR: Out of order at position " << info.position 
                      << ": legacy_idx=" << info.legacy_idx 
                      << " (previous was " << prev_idx << ")\n";
            has_errors = true;
        }
        prev_idx = info.legacy_idx;
    }
    
    if (!has_errors) {
        std::cout << "✓ All " << indices.size() << " records have legacy_residue_idx and are in order\n";
        std::cout << "  Range: " << (indices.empty() ? 0 : indices[0].legacy_idx) 
                  << " to " << (indices.empty() ? 0 : indices.back().legacy_idx) << "\n";
    }
    
    return !has_errors;
}

bool verify_pair_records(const json& records, const std::string& source) {
    std::cout << "\n=== Verifying " << source << " pair records ===\n";
    
    std::vector<std::pair<int, int>> pairs;
    bool has_errors = false;
    
    for (size_t i = 0; i < records.size(); i++) {
        const auto& record = records[i];
        
        int idx1 = -1, idx2 = -1;
        
        // Prefer base_i/base_j (always legacy indices)
        if (record.contains("base_i") && record.contains("base_j")) {
            idx1 = record["base_i"];
            idx2 = record["base_j"];
        } else if (record.contains("residue1_idx") && record.contains("residue2_idx")) {
            idx1 = record["residue1_idx"];
            idx2 = record["residue2_idx"];
        } else {
            std::cerr << "ERROR: Record at position " << i 
                      << " missing base_i/base_j or residue1_idx/residue2_idx\n";
            has_errors = true;
            continue;
        }
        
        if (idx1 <= 0 || idx2 <= 0) {
            std::cerr << "ERROR: Record at position " << i 
                      << " has invalid indices: (" << idx1 << ", " << idx2 << ")\n";
            has_errors = true;
            continue;
        }
        
        pairs.push_back({idx1, idx2});
    }
    
    if (!has_errors) {
        std::cout << "✓ All " << pairs.size() << " pair records have valid legacy indices\n";
    }
    
    return !has_errors;
}

bool compare_indices(const json& legacy_records, const json& modern_records, const std::string& record_type) {
    std::cout << "\n=== Comparing " << record_type << " indices ===\n";
    
    if (record_type == "frame_calc" || record_type == "base_frame_calc") {
        // Build maps by legacy_residue_idx
        std::map<int, json> legacy_map, modern_map;
        
        for (const auto& record : legacy_records) {
            if (record.contains("legacy_residue_idx")) {
                int idx = record["legacy_residue_idx"];
                legacy_map[idx] = record;
            }
        }
        
        for (const auto& record : modern_records) {
            if (record.contains("legacy_residue_idx")) {
                int idx = record["legacy_residue_idx"];
                modern_map[idx] = record;
            }
        }
        
        // Compare
        std::set<int> all_indices;
        for (const auto& [idx, _] : legacy_map) all_indices.insert(idx);
        for (const auto& [idx, _] : modern_map) all_indices.insert(idx);
        
        int matches = 0;
        // int mismatches = 0;  // Unused - commenting out to avoid warning
        int only_legacy = 0;
        int only_modern = 0;
        
        for (int idx : all_indices) {
            bool in_legacy = legacy_map.count(idx) > 0;
            bool in_modern = modern_map.count(idx) > 0;
            
            if (in_legacy && in_modern) {
                matches++;
            } else if (in_legacy && !in_modern) {
                only_legacy++;
                std::cerr << "WARNING: Index " << idx << " only in legacy\n";
            } else if (!in_legacy && in_modern) {
                only_modern++;
                std::cerr << "WARNING: Index " << idx << " only in modern\n";
            }
        }
        
        std::cout << "Matches: " << matches << "\n";
        if (only_legacy > 0) std::cout << "Only in legacy: " << only_legacy << "\n";
        if (only_modern > 0) std::cout << "Only in modern: " << only_modern << "\n";
        
        return (only_legacy == 0 && only_modern == 0);
    } else if (record_type == "pair_validation" || record_type == "find_bestpair_selection") {
        // Build sets of pairs
        std::set<std::pair<int, int>> legacy_pairs, modern_pairs;
        
        for (const auto& record : legacy_records) {
            int idx1 = -1, idx2 = -1;
            if (record.contains("base_i") && record.contains("base_j")) {
                idx1 = record["base_i"];
                idx2 = record["base_j"];
            } else if (record.contains("residue1_idx") && record.contains("residue2_idx")) {
                idx1 = record["residue1_idx"];
                idx2 = record["residue2_idx"];
            }
            if (idx1 > 0 && idx2 > 0) {
                legacy_pairs.insert({std::min(idx1, idx2), std::max(idx1, idx2)});
            }
        }
        
        for (const auto& record : modern_records) {
            int idx1 = -1, idx2 = -1;
            if (record.contains("base_i") && record.contains("base_j")) {
                idx1 = record["base_i"];
                idx2 = record["base_j"];
            } else if (record.contains("residue1_idx") && record.contains("residue2_idx")) {
                idx1 = record["residue1_idx"];
                idx2 = record["residue2_idx"];
            }
            if (idx1 > 0 && idx2 > 0) {
                modern_pairs.insert({std::min(idx1, idx2), std::max(idx1, idx2)});
            }
        }
        
        int matches = 0;
        for (const auto& pair : legacy_pairs) {
            if (modern_pairs.count(pair)) matches++;
        }
        
        std::cout << "Legacy pairs: " << legacy_pairs.size() << "\n";
        std::cout << "Modern pairs: " << modern_pairs.size() << "\n";
        std::cout << "Matching pairs: " << matches << "\n";
        
        return (legacy_pairs == modern_pairs);
    }
    
    return false;
}

json load_json_file(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + path.string());
    }
    json data;
    file >> data;
    return data;
}

std::vector<json> extract_records(const json& data, const std::string& record_type) {
    std::vector<json> records;
    
    if (data.is_array()) {
        for (const auto& item : data) {
            if (item.is_object() && item.value("type", "") == record_type) {
                records.push_back(item);
            }
        }
    } else if (data.contains("calculations") && data["calculations"].is_array()) {
        for (const auto& item : data["calculations"]) {
            if (item.is_object() && item.value("type", "") == record_type) {
                records.push_back(item);
            }
        }
    }
    
    return records;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <pdb_id> [record_type]\n";
        std::cerr << "  record_type: frame_calc, base_frame_calc, pair_validation, find_bestpair_selection\n";
        std::cerr << "  (default: all)\n";
        return 1;
    }
    
    std::string pdb_id = argv[1];
    std::string record_type = (argc >= 3) ? argv[2] : "all";
    
    std::filesystem::path legacy_dir = "data/json_legacy";
    std::filesystem::path modern_dir = "data/json";
    
    std::vector<std::string> types_to_check;
    if (record_type == "all") {
        types_to_check = {"base_frame_calc", "pair_validation", "find_bestpair_selection"};
    } else {
        types_to_check = {record_type};
    }
    
    bool all_ok = true;
    
    for (const auto& type : types_to_check) {
        std::filesystem::path legacy_file = legacy_dir / type / (pdb_id + ".json");
        std::filesystem::path modern_file = modern_dir / type / (pdb_id + ".json");
        
        if (!std::filesystem::exists(legacy_file)) {
            std::cerr << "WARNING: Legacy file not found: " << legacy_file << "\n";
            continue;
        }
        
        if (!std::filesystem::exists(modern_file)) {
            std::cerr << "WARNING: Modern file not found: " << modern_file << "\n";
            continue;
        }
        
        try {
            json legacy_data = load_json_file(legacy_file);
            json modern_data = load_json_file(modern_file);
            
            std::vector<json> legacy_records = extract_records(legacy_data, type);
            std::vector<json> modern_records = extract_records(modern_data, type);
            
            std::cout << "\n" << std::string(60, '=') << "\n";
            std::cout << "Checking " << type << " for " << pdb_id << "\n";
            std::cout << std::string(60, '=') << "\n";
            std::cout << "Legacy records: " << legacy_records.size() << "\n";
            std::cout << "Modern records: " << modern_records.size() << "\n";
            
            // Verify legacy
            bool legacy_ok = false;
            if (type == "base_frame_calc" || type == "frame_calc") {
                legacy_ok = verify_frame_records(legacy_records, "legacy");
            } else {
                legacy_ok = verify_pair_records(legacy_records, "legacy");
            }
            
            // Verify modern
            bool modern_ok = false;
            if (type == "base_frame_calc" || type == "frame_calc") {
                modern_ok = verify_frame_records(modern_records, "modern");
            } else {
                modern_ok = verify_pair_records(modern_records, "modern");
            }
            
            // Compare
            bool compare_ok = compare_indices(legacy_records, modern_records, type);
            
            if (!legacy_ok || !modern_ok || !compare_ok) {
                all_ok = false;
            }
            
        } catch (const std::exception& e) {
            std::cerr << "ERROR processing " << type << ": " << e.what() << "\n";
            all_ok = false;
        }
    }
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    if (all_ok) {
        std::cout << "✓ All checks passed!\n";
        return 0;
    } else {
        std::cout << "✗ Some checks failed!\n";
        return 1;
    }
}


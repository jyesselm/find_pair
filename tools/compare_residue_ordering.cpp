/**
 * @file compare_residue_ordering.cpp
 * @brief Compare residue ordering JSON files between modern and legacy
 * 
 * This tool compares two residue ordering JSON files and reports differences.
 * 
 * Usage: compare_residue_ordering <modern_json> <legacy_json>
 */

#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <algorithm>

using json = nlohmann::json;

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
    
    std::string key() const {
        return residue_name + "|" + chain_id + "|" + 
               std::to_string(residue_seq) + "|" + insertion_code;
    }
};

json load_json(const std::filesystem::path& file) {
    std::ifstream in_file(file);
    if (!in_file.is_open()) {
        throw std::runtime_error("Cannot open file: " + file.string());
    }
    json j;
    in_file >> j;
    return j;
}

std::vector<ResidueInfo> parse_residues(const json& j) {
    std::vector<ResidueInfo> residues;
    
    if (!j.contains("residues") || !j["residues"].is_array()) {
        return residues;
    }
    
    for (const auto& res_json : j["residues"]) {
        ResidueInfo info;
        info.legacy_index = res_json.value("legacy_index", 0);
        info.residue_name = res_json.value("residue_name", "");
        std::string chain_str = res_json.value("chain_id", " ");
        info.chain_id = chain_str.empty() ? " " : std::string(1, chain_str[0]);
        info.residue_seq = res_json.value("residue_seq", 0);
        std::string ins_str = res_json.value("insertion_code", " ");
        info.insertion_code = ins_str.empty() ? " " : std::string(1, ins_str[0]);
        residues.push_back(info);
    }
    
    return residues;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <modern_json> <legacy_json>\n";
        std::cerr << "Example: " << argv[0] 
                  << " data/residue_ordering/3G8T.json data/residue_ordering_legacy/3G8T.json\n";
        return 1;
    }
    
    std::filesystem::path modern_json = argv[1];
    std::filesystem::path legacy_json = argv[2];
    
    try {
        // Load JSON files
        json modern_j = load_json(modern_json);
        json legacy_j = load_json(legacy_json);
        
        // Parse residues
        std::vector<ResidueInfo> modern_residues = parse_residues(modern_j);
        std::vector<ResidueInfo> legacy_residues = parse_residues(legacy_j);
        
        std::cout << "\n" << std::string(60, '=') << "\n";
        std::cout << "Residue Ordering Comparison\n";
        std::cout << std::string(60, '=') << "\n";
        std::cout << "Modern: " << modern_json << "\n";
        std::cout << "Legacy: " << legacy_json << "\n\n";
        
        // Compare counts
        std::cout << "Total Residues:\n";
        std::cout << "  Modern: " << modern_residues.size() << "\n";
        std::cout << "  Legacy: " << legacy_residues.size() << "\n";
        
        if (modern_residues.size() != legacy_residues.size()) {
            std::cout << "  ✗ Count mismatch!\n";
        } else {
            std::cout << "  ✓ Counts match\n";
        }
        
        // Compare ordering
        size_t min_size = std::min(modern_residues.size(), legacy_residues.size());
        size_t matches = 0;
        size_t mismatches = 0;
        std::vector<size_t> mismatch_indices;
        
        for (size_t i = 0; i < min_size; i++) {
            if (modern_residues[i] == legacy_residues[i]) {
                matches++;
            } else {
                mismatches++;
                if (mismatch_indices.size() < 10) { // Show first 10 mismatches
                    mismatch_indices.push_back(i);
                }
            }
        }
        
        std::cout << "\nOrdering Comparison:\n";
        std::cout << "  Matches: " << matches << "\n";
        std::cout << "  Mismatches: " << mismatches << "\n";
        
        if (mismatches == 0 && modern_residues.size() == legacy_residues.size()) {
            std::cout << "  ✓ Perfect match!\n";
        } else {
            std::cout << "  ✗ Ordering differences found\n";
        }
        
        // Show mismatch details
        if (!mismatch_indices.empty()) {
            std::cout << "\nFirst " << mismatch_indices.size() << " mismatches:\n";
            for (size_t idx : mismatch_indices) {
                const auto& modern = modern_residues[idx];
                const auto& legacy = legacy_residues[idx];
                std::cout << "  Index " << (idx + 1) << ":\n";
                std::cout << "    Modern: " << modern.residue_name 
                          << " (chain " << modern.chain_id 
                          << ", seq " << modern.residue_seq << ")\n";
                std::cout << "    Legacy: " << legacy.residue_name 
                          << " (chain " << legacy.chain_id 
                          << ", seq " << legacy.residue_seq << ")\n";
            }
        }
        
        // Build key maps to find where residues moved
        std::map<std::string, size_t> modern_key_to_idx;
        std::map<std::string, size_t> legacy_key_to_idx;
        
        for (size_t i = 0; i < modern_residues.size(); i++) {
            modern_key_to_idx[modern_residues[i].key()] = i;
        }
        for (size_t i = 0; i < legacy_residues.size(); i++) {
            legacy_key_to_idx[legacy_residues[i].key()] = i;
        }
        
        // Find residues that are in different positions
        size_t moved_count = 0;
        for (const auto& [key, modern_idx] : modern_key_to_idx) {
            auto legacy_it = legacy_key_to_idx.find(key);
            if (legacy_it != legacy_key_to_idx.end()) {
                size_t legacy_idx = legacy_it->second;
                if (modern_idx != legacy_idx) {
                    moved_count++;
                }
            }
        }
        
        if (moved_count > 0) {
            std::cout << "\nResidues in different positions: " << moved_count << "\n";
        }
        
        std::cout << "\n" << std::string(60, '=') << "\n";
        
        // Return exit code based on match
        if (mismatches == 0 && modern_residues.size() == legacy_residues.size()) {
            return 0; // Success
        } else {
            return 1; // Mismatch
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}


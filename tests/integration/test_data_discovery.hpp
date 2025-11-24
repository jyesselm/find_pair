/**
 * @file test_data_discovery.hpp
 * @brief Utility to discover PDB/JSON pairs for integration testing
 */

#pragma once

#include <filesystem>
#include <vector>
#include <string>
#include <set>
#include <fstream>
#include <nlohmann/json.hpp>

namespace x3dna::test {

/**
 * @struct pdb_json_pair
 * @brief Represents a PDB file with its corresponding JSON file
 */
struct pdb_json_pair {
    std::filesystem::path pdb_file;
    std::filesystem::path json_file;
    std::filesystem::path globals_file;
    std::string pdb_name;
};

/**
 * @class test_data_discovery
 * @brief Discovers PDB files that have corresponding JSON files
 */
class test_data_discovery {
public:
    /**
     * @brief Discover all PDB files that have corresponding JSON files
     * @param pdb_dir Directory containing PDB files
     * @param json_dir Directory containing JSON files
     * @return Vector of PDB/JSON pairs
     */
    static std::vector<pdb_json_pair>
    discover_pairs(const std::filesystem::path& pdb_dir = "data/pdb",
                   const std::filesystem::path& json_dir = "data/json_legacy") {
        std::vector<pdb_json_pair> pairs;

        if (!std::filesystem::exists(pdb_dir) || !std::filesystem::exists(json_dir)) {
            return pairs;
        }

        // Iterate through all PDB files
        for (const auto& pdb_file : std::filesystem::directory_iterator(pdb_dir)) {
            if (pdb_file.path().extension() != ".pdb") {
                continue;
            }

            std::string pdb_name = pdb_file.path().stem().string();
            std::filesystem::path json_file = json_dir / (pdb_name + ".json");
            std::filesystem::path globals_file = json_dir / (pdb_name + "_globals.json");

            // Check if JSON file exists and doesn't contain "globals" in the name
            std::string json_filename = json_file.filename().string();
            if (json_filename.find("globals") != std::string::npos) {
                continue; // Skip files with "globals" in the name
            }

            if (std::filesystem::exists(json_file)) {
                pairs.push_back({pdb_file.path(), json_file, globals_file, pdb_name});
            }
        }

        return pairs;
    }

    /**
     * @brief Check if a PDB file has a corresponding JSON file
     * @param pdb_file Path to PDB file
     * @param json_dir Directory containing JSON files
     * @return True if JSON file exists
     */
    static bool has_json(const std::filesystem::path& pdb_file,
                         const std::filesystem::path& json_dir = "data/json_legacy") {
        std::string pdb_name = pdb_file.stem().string();
        std::filesystem::path json_file = json_dir / (pdb_name + ".json");
        return std::filesystem::exists(json_file);
    }

    /**
     * @brief Load a test set from a JSON file
     * @param test_set_file Path to test set JSON file
     * @return Set of PDB IDs in the test set, or empty set if file doesn't exist
     */
    static std::set<std::string> load_test_set(const std::filesystem::path& test_set_file) {
        std::set<std::string> pdb_ids;
        
        if (!std::filesystem::exists(test_set_file)) {
            return pdb_ids;
        }
        
        try {
            std::ifstream file(test_set_file);
            nlohmann::json json;
            file >> json;
            
            if (json.contains("pdb_ids") && json["pdb_ids"].is_array()) {
                for (const auto& id : json["pdb_ids"]) {
                    if (id.is_string()) {
                        pdb_ids.insert(id.get<std::string>());
                    }
                }
            }
        } catch (const std::exception&) {
            // If parsing fails, return empty set
        }
        
        return pdb_ids;
    }

    /**
     * @brief Discover PDB/JSON pairs filtered by a test set
     * @param test_set_size Size of test set to use (10, 50, 100, 500, or 1000)
     * @param pdb_dir Directory containing PDB files
     * @param json_dir Directory containing JSON files
     * @param test_sets_dir Directory containing test set JSON files
     * @return Vector of PDB/JSON pairs from the test set
     */
    static std::vector<pdb_json_pair>
    discover_pairs_from_test_set(int test_set_size,
                                 const std::filesystem::path& pdb_dir = "data/pdb",
                                 const std::filesystem::path& json_dir = "data/json_legacy",
                                 const std::filesystem::path& test_sets_dir = "data/test_sets") {
        std::vector<pdb_json_pair> pairs;
        
        // Load test set
        std::filesystem::path test_set_file = test_sets_dir / ("test_set_" + std::to_string(test_set_size) + ".json");
        std::set<std::string> test_set_pdb_ids = load_test_set(test_set_file);
        
        if (test_set_pdb_ids.empty()) {
            return pairs;
        }
        
        if (!std::filesystem::exists(pdb_dir) || !std::filesystem::exists(json_dir)) {
            return pairs;
        }
        
        // Iterate through test set PDB IDs
        for (const std::string& pdb_name : test_set_pdb_ids) {
            std::filesystem::path pdb_file = pdb_dir / (pdb_name + ".pdb");
            std::filesystem::path json_file = json_dir / (pdb_name + ".json");
            std::filesystem::path globals_file = json_dir / (pdb_name + "_globals.json");
            
            // Check if both PDB and JSON files exist
            if (std::filesystem::exists(pdb_file) && std::filesystem::exists(json_file)) {
                pairs.push_back({pdb_file, json_file, globals_file, pdb_name});
            }
        }
        
        return pairs;
    }

    /**
     * @brief Check if a JSON file has a pdb_atoms record (with or without atoms array)
     * @param json_file Path to JSON file
     * @return True if JSON has pdb_atoms record
     */
    static bool has_pdb_atoms_record(const std::filesystem::path& json_file) {
        if (!std::filesystem::exists(json_file)) {
            return false;
        }
        
        try {
            std::ifstream file(json_file);
            nlohmann::json json;
            file >> json;
            
            if (!json.contains("calculations") || !json["calculations"].is_array()) {
                return false;
            }
            
            // Look for pdb_atoms record (atoms array is optional - some legacy JSONs have split files)
            for (const auto& calc : json["calculations"]) {
                if (calc.contains("type") && calc["type"] == "pdb_atoms") {
                    return true;
                }
            }
        } catch (const std::exception&) {
            return false;
        }
        
        return false;
    }
};

} // namespace x3dna::test

/**
 * @file test_data_discovery.hpp
 * @brief Utility to discover PDB/JSON pairs for integration testing
 */

#pragma once

#include <filesystem>
#include <vector>
#include <string>

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

            // Check if JSON file exists
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
};

} // namespace x3dna::test

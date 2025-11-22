/**
 * @file test_json_generation_filtered.cpp
 * @brief Generate JSON files for specific PDBs only (problematic ones)
 * 
 * This is a simplified version that can be configured to process only
 * PDBs from a list file.
 */

#include <gtest/gtest.h>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include "integration_test_base.hpp"
#include "test_data_discovery.hpp"
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <sstream>
#include <thread>
#include <mutex>
#include <atomic>

namespace x3dna::test {

using namespace x3dna::core;
using namespace x3dna::io;

class JsonGenerationFilteredTest : public integration_test_base {
protected:
    void SetUp() override {
        integration_test_base::SetUp();
        output_dir_ = std::filesystem::path("data/json");
        std::filesystem::create_directories(output_dir_);
        
        // Load problematic PDBs from file
        load_problematic_pdbs();
    }

    std::filesystem::path output_dir_;
    std::set<std::string> problematic_pdbs_;

    void load_problematic_pdbs() {
        std::filesystem::path problem_file = "docs/problematic_pdbs.txt";
        if (!std::filesystem::exists(problem_file)) {
            return;
        }

        std::ifstream file(problem_file);
        std::string line;
        while (std::getline(file, line)) {
            // Skip comments and empty lines
            if (line.empty() || line[0] == '#') {
                continue;
            }
            
            // Extract PDB ID (first word)
            std::istringstream iss(line);
            std::string pdb_id;
            if (iss >> pdb_id) {
                problematic_pdbs_.insert(pdb_id);
            }
        }
    }

    void generate_json_for_pair(const pdb_json_pair& pair) {
        try {
            // Create output JSON structure matching legacy format exactly
            nlohmann::json output_json;
            output_json["pdb_file"] = pair.pdb_file.string();
            output_json["pdb_name"] = pair.pdb_name;
            output_json["calculations"] = nlohmann::json::array();

            // Add metadata
            nlohmann::json metadata;
            metadata["version"] = "X3DNA Modernized C++ Library";
            output_json["metadata"] = metadata;

            // Parse PDB file
            PdbParser parser;
            parser.set_include_hetatm(true);
            parser.set_include_waters(true);

            Structure structure = parser.parse_file(pair.pdb_file);

            // Generate pdb_atoms record using Structure::to_json_legacy()
            // This uses the actual modernized code to generate the JSON from parsed PDB
            nlohmann::json structure_json = structure.to_json_legacy();

            // Wrap in calculations array format (legacy JSON format)
            nlohmann::json pdb_atoms_record;
            pdb_atoms_record["type"] = "pdb_atoms";
            pdb_atoms_record["num_atoms"] = structure_json["num_atoms"];
            pdb_atoms_record["atoms"] = structure_json["atoms"];
            
            output_json["calculations"].push_back(pdb_atoms_record);

            // Write output JSON file
            std::filesystem::path output_file = output_dir_ / (pair.pdb_name + ".json");
            std::ofstream out_file(output_file);
            if (!out_file.is_open()) {
                throw std::runtime_error("Cannot open output file: " + output_file.string());
            }
            out_file << output_json.dump(2);
            out_file.close();

        } catch (const std::exception& e) {
            FAIL() << "Error generating JSON for " << pair.pdb_name << ": " << e.what();
        }
    }
};

TEST_F(JsonGenerationFilteredTest, GenerateProblematicPdbs) {
    // Get all PDB/JSON pairs
    auto all_pairs = test_data_discovery::discover_pairs();
    
    // Filter to only problematic PDBs
    std::vector<pdb_json_pair> filtered_pairs;
    for (const auto& pair : all_pairs) {
        if (problematic_pdbs_.empty() || problematic_pdbs_.count(pair.pdb_name) > 0) {
            filtered_pairs.push_back(pair);
        }
    }
    
    std::cout << "Processing " << filtered_pairs.size() << " problematic PDBs out of " 
              << all_pairs.size() << " total" << std::endl;
    
    size_t success_count = 0;
    size_t failure_count = 0;
    
    // Process in parallel
    std::vector<std::thread> threads;
    size_t num_threads = std::min(filtered_pairs.size(), 
                                  static_cast<size_t>(std::thread::hardware_concurrency()));
    
    std::mutex mtx;
    std::atomic<size_t> processed{0};
    std::atomic<size_t> success{0};
    std::atomic<size_t> fail{0};
    std::vector<std::string> failures;
    
    auto process_batch = [&](size_t start_idx, size_t end_idx) {
        for (size_t i = start_idx; i < end_idx && i < filtered_pairs.size(); ++i) {
            try {
                generate_json_for_pair(filtered_pairs[i]);
                success++;
                size_t current = ++processed;
                if (current % 10 == 0) {
                    std::lock_guard<std::mutex> lock(mtx);
                    std::cout << "Progress: " << current << "/" << filtered_pairs.size() << std::endl;
                }
            } catch (const std::exception& e) {
                fail++;
                processed++;
                std::lock_guard<std::mutex> lock(mtx);
                std::string error_msg = filtered_pairs[i].pdb_name + ": " + e.what();
                failures.push_back(error_msg);
                std::cerr << "Failed: " << error_msg << std::endl;
            }
        }
    };
    
    // Divide work among threads
    size_t batch_size = (filtered_pairs.size() + num_threads - 1) / num_threads;
    for (size_t t = 0; t < num_threads; ++t) {
        size_t start = t * batch_size;
        size_t end = start + batch_size;
        threads.emplace_back(process_batch, start, end);
    }
    
    // Wait for all threads
    for (auto& thread : threads) {
        thread.join();
    }
    
    success_count = success.load();
    failure_count = fail.load();
    
    std::cout << "\nJSON Generation Summary:" << std::endl;
    std::cout << "  Success: " << success_count << std::endl;
    std::cout << "  Failures: " << failure_count << std::endl;
    
    if (!failures.empty()) {
        std::cout << "\nFailures:" << std::endl;
        for (const auto& failure : failures) {
            std::cout << "  - " << failure << std::endl;
        }
    }
    
    EXPECT_GT(success_count, 0) << "No JSON files were generated successfully";
}

} // namespace x3dna::test


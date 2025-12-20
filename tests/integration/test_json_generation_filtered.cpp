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
#include <x3dna/algorithms/base_frame_calculator.hpp>
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
using namespace x3dna::algorithms;

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

            // Calculate frames and record frame calculations
            BaseFrameCalculator calculator("data/templates");
            calculator.calculate_all_frames(structure);

            // Record frame calculations for each residue
            // Legacy residue_idx is 1-based and counts ALL residues (including amino acids, etc.)
            size_t residue_idx = 1;
            for (auto& chain : structure.chains()) {
                for (auto& residue : chain.residues()) {
                    // Only process nucleotide residues
                    if (residue.residue_type() != ResidueType::UNKNOWN &&
                        residue.residue_type() != ResidueType::AMINO_ACID) {

                        // Calculate frame directly (don't check has_reference_frame)
                        // Use calculate_frame_const to avoid modifying residue during iteration
                        FrameCalculationResult frame_result = calculator.calculate_frame_const(residue);

                        if (frame_result.is_valid) {
                            // Record base_frame_calc
                            nlohmann::json base_frame_record;
                            base_frame_record["type"] = "base_frame_calc";
                            base_frame_record["residue_idx"] = residue_idx;
                            base_frame_record["base_type"] = std::string(1, residue.one_letter_code());
                            base_frame_record["residue_name"] = residue.name();
                            base_frame_record["chain_id"] = residue.chain_id();
                            base_frame_record["residue_seq"] = residue.seq_num();
                            if (!residue.insertion().empty()) {
                                base_frame_record["insertion"] = residue.insertion();
                            }
                            base_frame_record["standard_template"] = frame_result.template_file.string();
                            base_frame_record["rms_fit"] = frame_result.rms_fit;
                            base_frame_record["num_matched_atoms"] = frame_result.num_matched;
                            base_frame_record["matched_atoms"] = frame_result.matched_atoms;
                            output_json["calculations"].push_back(base_frame_record);

                            // Record ls_fitting
                            nlohmann::json ls_fitting_record;
                            ls_fitting_record["type"] = "ls_fitting";
                            ls_fitting_record["residue_idx"] = residue_idx;
                            ls_fitting_record["residue_name"] = residue.name();
                            ls_fitting_record["chain_id"] = residue.chain_id();
                            ls_fitting_record["residue_seq"] = residue.seq_num();
                            if (!residue.insertion().empty()) {
                                ls_fitting_record["insertion"] = residue.insertion();
                            }
                            ls_fitting_record["num_points"] = frame_result.num_matched;
                            ls_fitting_record["rms_fit"] = frame_result.rms_fit;

                            // Rotation matrix
                            nlohmann::json rot_array = nlohmann::json::array();
                            for (int i = 0; i < 3; ++i) {
                                nlohmann::json row = nlohmann::json::array();
                                for (int j = 0; j < 3; ++j) {
                                    row.push_back(frame_result.rotation_matrix.at(i, j));
                                }
                                rot_array.push_back(row);
                            }
                            ls_fitting_record["rotation_matrix"] = rot_array;

                            // Translation
                            ls_fitting_record["translation"] = nlohmann::json::array({frame_result.translation.x(),
                                                                                      frame_result.translation.y(),
                                                                                      frame_result.translation.z()});
                            output_json["calculations"].push_back(ls_fitting_record);
                        }
                    }

                    // Count all residues (to match legacy residue_idx behavior)
                    residue_idx++;
                }
            }

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

    std::cout << "Processing " << filtered_pairs.size() << " problematic PDBs out of " << all_pairs.size() << " total"
              << std::endl;

    size_t success_count = 0;
    size_t failure_count = 0;

    // Process in parallel
    std::vector<std::thread> threads;
    size_t num_threads = std::min(filtered_pairs.size(), static_cast<size_t>(std::thread::hardware_concurrency()));

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

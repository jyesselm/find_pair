/**
 * @file test_json_generation.cpp
 * @brief Integration test to generate JSON files using new/modernized code
 *
 * This test uses the modernized C++ code to generate JSON files in data/json/
 * directory. It only generates what is currently implemented:
 *
 * Currently implemented (Stage 3):
 * - pdb_atoms: Using PdbParser to parse PDB files, then Structure::to_json_legacy()
 *
 * Future (as algorithms are implemented):
 * - ref_frame: Using BaseFrameCalculator (Stage 4)
 * - base_pair: Using BasePairFinder (Stage 5)
 * - bpstep_params: Using ParameterCalculator (Stage 6)
 * - helical_params: Using ParameterCalculator (Stage 6)
 *
 * The generated files can then be compared with legacy files using
 * scripts/compare_json_files.py to verify correctness.
 */

#include <gtest/gtest.h>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp> // Stage 3 - PDB parser
// TODO: Add these includes as algorithms are implemented:
// #include <x3dna/algorithms/base_frame_calculator.hpp>  // Stage 4
// #include <x3dna/algorithms/base_pair_finder.hpp>  // Stage 5
// #include <x3dna/algorithms/parameter_calculator.hpp>  // Stage 6
#include "integration_test_base.hpp"
#include "test_data_discovery.hpp"
#include <filesystem>
#include <fstream>
#include <vector>
#include <mutex>
#include <thread>
#include <atomic>

namespace x3dna::test {

using namespace x3dna::core;
using namespace x3dna::io;

/**
 * @brief Test class for generating JSON files from legacy data
 */
class JsonGenerationTest : public integration_test_base {
protected:
    void SetUp() override {
        integration_test_base::SetUp();

        // Create output directory if it doesn't exist
        output_dir_ = std::filesystem::path("data/json");
        std::filesystem::create_directories(output_dir_);
    }

    std::filesystem::path output_dir_;

    /**
     * @brief Generate JSON file for a single PDB/JSON pair using new code
     *
     * This function uses the modernized code to generate JSON files.
     * It only generates what is currently implemented:
     * - pdb_atoms: Using PdbParser to parse PDB files, then Structure::to_json_legacy() (Stage 3 -
     * I/O)
     * - Other records: Will be added as algorithms are implemented
     */
    void generate_json_for_pair(const pdb_json_pair& pair) {
        try {
            // Create output JSON structure matching legacy format exactly
            nlohmann::json output_json;
            output_json["pdb_file"] = pair.pdb_file.string();
            output_json["pdb_name"] = pair.pdb_name;
            output_json["calculations"] = nlohmann::json::array();

            // Add metadata to match legacy format
            nlohmann::json metadata;
            metadata["version"] = "X3DNA Modernized C++ Library";
            output_json["metadata"] = metadata;

            // Use PdbParser to parse PDB file directly
            // Configure to match legacy JSON exactly (includes HETATM and waters)
            PdbParser parser;
            parser.set_include_hetatm(true); // Legacy JSON includes HETATM
            parser.set_include_waters(true); // Legacy JSON includes waters

            Structure structure = parser.parse_file(pair.pdb_file);

            // Optional: Compare with legacy JSON for validation (but use parsed structure)
            // This helps verify PDB parsing correctness
            try {
                auto legacy_json = load_legacy_json(pair.json_file);
                auto pdb_atoms_records = find_records_by_type(legacy_json, "pdb_atoms");
                if (!pdb_atoms_records.empty()) {
                    const auto& legacy_atoms = pdb_atoms_records[0];
                    Structure legacy_structure = Structure::from_json_legacy(legacy_atoms);

                    // Compare atom counts (basic validation)
                    if (structure.num_atoms() != legacy_structure.num_atoms()) {
                        std::cerr << "Warning: Atom count mismatch for " << pair.pdb_name
                                  << ": parsed=" << structure.num_atoms()
                                  << ", legacy=" << legacy_structure.num_atoms() << std::endl;
                    }
                }
            } catch (const std::exception& e) {
                // If legacy JSON doesn't exist or can't be loaded, continue anyway
                // We're using the parsed structure, not the legacy one
            }

            // Generate pdb_atoms record using Structure::to_json_legacy()
            // This uses the actual modernized code to generate the JSON from parsed PDB
            nlohmann::json structure_json = structure.to_json_legacy();

            // Wrap in calculations array format (legacy JSON format)
            nlohmann::json pdb_atoms_record;
            pdb_atoms_record["type"] = "pdb_atoms";
            pdb_atoms_record["num_atoms"] = structure_json["num_atoms"];
            pdb_atoms_record["atoms"] = structure_json["atoms"];
            // Note: pdb_id, num_residues, num_chains are in structure_json but not in legacy
            // format

            output_json["calculations"].push_back(pdb_atoms_record);

            // TODO: Stage 4 - Generate ref_frame records using BaseFrameCalculator
            // Once BaseFrameCalculator is implemented:
            // BaseFrameCalculator calculator;
            // calculator.calculate_all_frames(structure);
            // for (const auto* residue : structure.nucleotides()) {
            //     if (residue->reference_frame()) {
            //         nlohmann::json ref_frame_record;
            //         ref_frame_record["type"] = "ref_frame";
            //         ref_frame_record["residue_idx"] = residue->sequence_number();
            //         ref_frame_record["orien"] =
            //         residue->reference_frame()->rotation().to_json_legacy();
            //         ref_frame_record["org"] = residue->reference_frame()->origin().to_json();
            //         output_json["calculations"].push_back(ref_frame_record);
            //     }
            // }

            // TODO: Stage 5 - Generate base_pair records using BasePairFinder
            // Once BasePairFinder is implemented:
            // BasePairFinder finder;
            // auto base_pairs = finder.find_pairs(structure);
            // for (const auto& bp : base_pairs) {
            //     output_json["calculations"].push_back(bp.to_json_legacy());
            // }

            // TODO: Stage 6 - Generate bpstep_params and helical_params using ParameterCalculator
            // Once ParameterCalculator is implemented:
            // ParameterCalculator param_calc;
            // for (size_t i = 0; i < base_pairs.size() - 1; ++i) {
            //     auto step_params = param_calc.calculate_step_parameters(base_pairs[i],
            //     base_pairs[i+1]);
            //     output_json["calculations"].push_back(step_params.to_json_legacy(i, i+1));
            //     auto helical_params = param_calc.calculate_helical_parameters(base_pairs[i],
            //     base_pairs[i+1]);
            //     output_json["calculations"].push_back(helical_params.to_json_legacy(i, i+1));
            // }

            // For now, we can only generate pdb_atoms with the current implementation
            // Other record types will be added as algorithms are implemented

            // Write output JSON file
            std::filesystem::path output_file = output_dir_ / (pair.pdb_name + ".json");
            std::ofstream out_file(output_file);
            if (!out_file.is_open()) {
                throw std::runtime_error("Cannot open output file: " + output_file.string());
            }
            out_file << output_json.dump(2); // Pretty print with 2-space indent
            out_file.close();

        } catch (const std::exception& e) {
            FAIL() << "Error generating JSON for " << pair.pdb_name << ": " << e.what();
        }
    }
};

/**
 * @brief Generate JSON files for all discovered PDB/JSON pairs
 */
TEST_F(JsonGenerationTest, GenerateAllJsonFiles) {
    std::mutex mtx;
    std::atomic<size_t> success_count{0};
    std::atomic<size_t> failure_count{0};
    std::vector<std::string> failures;

    // Process pairs in parallel
    std::vector<std::thread> threads;
    size_t num_threads =
        std::min(pairs_.size(), static_cast<size_t>(std::thread::hardware_concurrency()));

    auto process_batch = [&](size_t start_idx, size_t end_idx) {
        for (size_t i = start_idx; i < end_idx && i < pairs_.size(); ++i) {
            try {
                generate_json_for_pair(pairs_[i]);
                success_count++;
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> lock(mtx);
                failures.push_back(pairs_[i].pdb_name + ": " + e.what());
                failure_count++;
            }
        }
    };

    // Divide work among threads
    size_t batch_size = (pairs_.size() + num_threads - 1) / num_threads;
    for (size_t t = 0; t < num_threads; ++t) {
        size_t start = t * batch_size;
        size_t end = start + batch_size;
        threads.emplace_back(process_batch, start, end);
    }

    // Wait for all threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Report results
    std::cout << "\nJSON Generation Summary:" << std::endl;
    std::cout << "  Success: " << success_count.load() << std::endl;
    std::cout << "  Failures: " << failure_count.load() << std::endl;

    if (!failures.empty()) {
        std::cout << "\nFailures:" << std::endl;
        for (const auto& failure : failures) {
            std::cout << "  - " << failure << std::endl;
        }
    }

    // Verify at least some files were generated
    EXPECT_GT(success_count.load(), 0) << "No JSON files were generated successfully";
}

/**
 * @brief Generate JSON file for a single test case (for debugging)
 */
TEST_F(JsonGenerationTest, GenerateSingleJsonFile) {
    if (pairs_.empty()) {
        GTEST_SKIP() << "No PDB/JSON pairs available";
    }

    // Use first pair for testing
    const auto& pair = pairs_[0];
    generate_json_for_pair(pair);

    // Verify file was created
    std::filesystem::path output_file = output_dir_ / (pair.pdb_name + ".json");
    EXPECT_TRUE(std::filesystem::exists(output_file))
        << "Output file was not created: " << output_file;

    // Verify file is valid JSON
    std::ifstream file(output_file);
    nlohmann::json json;
    EXPECT_NO_THROW(file >> json) << "Generated file is not valid JSON";
}

/**
 * @brief Generate JSON files for first 10 PDB/JSON pairs (for testing)
 */
TEST_F(JsonGenerationTest, GenerateFirstTenJsonFiles) {
    if (pairs_.empty()) {
        GTEST_SKIP() << "No PDB/JSON pairs available";
    }

    size_t num_to_process = std::min(pairs_.size(), size_t(10));
    std::cout << "\nGenerating JSON for first " << num_to_process << " PDB files..." << std::endl;

    std::mutex mtx;
    std::atomic<size_t> success_count{0};
    std::atomic<size_t> failure_count{0};
    std::vector<std::string> failures;

    // Process pairs in parallel
    std::vector<std::thread> threads;
    size_t num_threads =
        std::min(num_to_process, static_cast<size_t>(std::thread::hardware_concurrency()));

    auto process_batch = [&](size_t start_idx, size_t end_idx) {
        for (size_t i = start_idx; i < end_idx && i < num_to_process; ++i) {
            try {
                generate_json_for_pair(pairs_[i]);
                success_count++;
                std::lock_guard<std::mutex> lock(mtx);
                std::cout << "  ✓ Generated: " << pairs_[i].pdb_name << std::endl;
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> lock(mtx);
                failures.push_back(pairs_[i].pdb_name + ": " + e.what());
                failure_count++;
                std::cout << "  ✗ Failed: " << pairs_[i].pdb_name << " - " << e.what() << std::endl;
            }
        }
    };

    // Divide work among threads
    size_t batch_size = (num_to_process + num_threads - 1) / num_threads;
    for (size_t t = 0; t < num_threads; ++t) {
        size_t start = t * batch_size;
        size_t end = start + batch_size;
        threads.emplace_back(process_batch, start, end);
    }

    // Wait for all threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Report results
    std::cout << "\nJSON Generation Summary (first " << num_to_process << " files):" << std::endl;
    std::cout << "  Success: " << success_count.load() << std::endl;
    std::cout << "  Failures: " << failure_count.load() << std::endl;

    if (!failures.empty()) {
        std::cout << "\nFailures:" << std::endl;
        for (const auto& failure : failures) {
            std::cout << "  - " << failure << std::endl;
        }
    }

    // Verify at least some files were generated
    EXPECT_GT(success_count.load(), 0) << "No JSON files were generated successfully";
}

} // namespace x3dna::test

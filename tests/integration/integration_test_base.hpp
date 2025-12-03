/**
 * @file integration_test_base.hpp
 * @brief Base class for integration tests
 */

#pragma once

#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include "test_data_discovery.hpp"

namespace x3dna::test {

// pdb_json_pair is defined in test_data_discovery.hpp

/**
 * @class integration_test_base
 * @brief Base class for integration tests that compare with legacy JSON
 */
class integration_test_base : public ::testing::Test {
protected:
    void SetUp() override {
        // Try to use test set first (test_set_10 by default)
        pairs_ = test_data_discovery::discover_pairs_from_test_set(10);

        // Filter to only pairs with pdb_atoms records
        if (!pairs_.empty()) {
            std::vector<pdb_json_pair> filtered_pairs;
            for (const auto& pair : pairs_) {
                if (test_data_discovery::has_pdb_atoms_record(pair.json_file)) {
                    filtered_pairs.push_back(pair);
                }
            }
            pairs_ = filtered_pairs;
        }

        // Fall back to all pairs if test set is empty
        if (pairs_.empty()) {
            pairs_ = test_data_discovery::discover_pairs();

            // Filter to only pairs with pdb_atoms records
            if (!pairs_.empty()) {
                std::vector<pdb_json_pair> filtered_pairs;
                for (const auto& pair : pairs_) {
                    if (test_data_discovery::has_pdb_atoms_record(pair.json_file)) {
                        filtered_pairs.push_back(pair);
                    }
                }
                pairs_ = filtered_pairs;
            }
        }

        if (pairs_.empty()) {
            GTEST_SKIP() << "No PDB/JSON pairs found for testing. "
                         << "Place JSON files in data/json_legacy/ to enable tests.";
        }
    }

    std::vector<pdb_json_pair> pairs_;

    /**
     * @brief Load legacy JSON file
     * @param json_file Path to JSON file
     * @return Parsed JSON object
     */
    nlohmann::json load_legacy_json(const std::filesystem::path& json_file) {
        std::ifstream file(json_file);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open JSON file: " + json_file.string());
        }
        nlohmann::json json;
        file >> json;
        return json;
    }

    /**
     * @brief Find records by type in calculations array
     * @param json JSON object containing calculations
     * @param type Record type to find
     * @return Vector of matching records
     */
    std::vector<nlohmann::json> find_records_by_type(const nlohmann::json& json,
                                                     const std::string& type) {
        std::vector<nlohmann::json> results;

        if (!json.contains("calculations")) {
            return results;
        }

        for (const auto& calc : json["calculations"]) {
            if (calc.contains("type") && calc["type"] == type) {
                results.push_back(calc);
            }
        }

        return results;
    }
};

} // namespace x3dna::test

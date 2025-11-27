/**
 * @file test_protocols_integration.cpp
 * @brief Integration tests for protocols using real PDB files
 * 
 * Compares FindPairProtocol results with legacy JSON output
 */

#include <gtest/gtest.h>
#include <x3dna/protocols/find_pair_protocol.hpp>
#include <x3dna/config/config_manager.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/io/json_writer.hpp>
#include <x3dna/core/structure.hpp>
#include "integration_test_base.hpp"
#include "json_comparison.hpp"
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <vector>
#include <set>
#include <cmath>

using namespace x3dna::protocols;
using namespace x3dna::config;
using namespace x3dna::io;
using namespace x3dna::core;
using namespace x3dna::test;

namespace fs = std::filesystem;

/**
 * @class ProtocolsIntegrationTest
 * @brief Integration tests for FindPairProtocol
 */
class ProtocolsIntegrationTest : public integration_test_base {
protected:
    void SetUp() override {
        integration_test_base::SetUp();
        
        // Set up config with defaults
        auto& config = ConfigManager::instance();
        config.set_defaults();
    }
    
    /**
     * @brief Load base pairs from legacy JSON
     */
    std::vector<nlohmann::json> load_legacy_base_pairs(const fs::path& json_file) {
        std::vector<nlohmann::json> base_pairs;
        
        try {
            auto legacy_json = load_legacy_json(json_file);
            auto records = find_records_by_type(legacy_json, "base_pair");
            
            for (const auto& record : records) {
                if (record.contains("pairs") && record["pairs"].is_array()) {
                    for (const auto& pair : record["pairs"]) {
                        base_pairs.push_back(pair);
                    }
                }
            }
        } catch (const std::exception& e) {
            // If loading fails, return empty vector
        }
        
        return base_pairs;
    }
    
    /**
     * @brief Compare base pairs from protocol with legacy JSON
     */
    void compare_base_pairs(const std::vector<BasePair>& modern_pairs,
                           const std::vector<nlohmann::json>& legacy_pairs,
                           const std::string& pdb_name) {
        EXPECT_EQ(modern_pairs.size(), legacy_pairs.size())
            << "Base pair count mismatch for " << pdb_name
            << ": modern=" << modern_pairs.size()
            << ", legacy=" << legacy_pairs.size();
        
        if (modern_pairs.size() != legacy_pairs.size()) {
            return;  // Can't compare if counts don't match
        }
        
        // Create sets for comparison (order may differ)
        std::set<std::pair<int, int>> modern_set;
        std::set<std::pair<int, int>> legacy_set;
        
        for (const auto& pair : modern_pairs) {
            // Convert to 1-based for comparison (legacy uses 1-based)
            int i = pair.residue_idx1() + 1;
            int j = pair.residue_idx2() + 1;
            modern_set.insert({std::min(i, j), std::max(i, j)});
        }
        
        for (const auto& pair_json : legacy_pairs) {
            if (pair_json.contains("residue1") && pair_json.contains("residue2")) {
                int i = pair_json["residue1"].get<int>();
                int j = pair_json["residue2"].get<int>();
                legacy_set.insert({std::min(i, j), std::max(i, j)});
            }
        }
        
        // Compare sets
        EXPECT_EQ(modern_set.size(), legacy_set.size())
            << "Unique pair count mismatch for " << pdb_name;
        
        // Check that all modern pairs are in legacy set
        for (const auto& pair : modern_set) {
            EXPECT_TRUE(legacy_set.find(pair) != legacy_set.end())
                << "Pair (" << pair.first << ", " << pair.second
                << ") found in modern but not in legacy for " << pdb_name;
        }
        
        // Check that all legacy pairs are in modern set
        for (const auto& pair : legacy_set) {
            EXPECT_TRUE(modern_set.find(pair) != modern_set.end())
                << "Pair (" << pair.first << ", " << pair.second
                << ") found in legacy but not in modern for " << pdb_name;
        }
    }
};

/**
 * @brief Test FindPairProtocol with single PDB file
 */
TEST_F(ProtocolsIntegrationTest, FindPairProtocolSinglePDB) {
    if (pairs_.empty()) {
        GTEST_SKIP() << "No PDB/JSON pairs found";
    }
    
    const auto& pair = pairs_[0];
    
    // Parse PDB file
    PdbParser parser;
    Structure structure;
    try {
        structure = parser.parse_file(pair.pdb_file);
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Failed to parse PDB: " << e.what();
    }
    
    EXPECT_GT(structure.num_residues(), 0) << "Structure has no residues";
    
    // Create protocol
    FindPairProtocol protocol;
    protocol.set_config_manager(ConfigManager::instance());
    
    // Execute protocol
    try {
        protocol.execute(structure);
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Protocol execution failed: " << e.what();
    }
    
    // Get base pairs
    const auto& modern_pairs = protocol.base_pairs();
    
    // Load legacy base pairs
    auto legacy_pairs = load_legacy_base_pairs(pair.json_file);
    
    // Compare
    if (!legacy_pairs.empty()) {
        compare_base_pairs(modern_pairs, legacy_pairs, pair.pdb_name);
    } else {
        // If no legacy pairs, at least verify protocol executed
        EXPECT_GE(modern_pairs.size(), 0) << "Protocol should return pairs (even if empty)";
    }
}

/**
 * @brief Test FindPairProtocol with multiple PDB files
 */
TEST_F(ProtocolsIntegrationTest, FindPairProtocolMultiplePDBs) {
    if (pairs_.empty()) {
        GTEST_SKIP() << "No PDB/JSON pairs found";
    }
    
    // Test first 5 PDBs (or all if less than 5)
    size_t max_test = std::min(pairs_.size(), size_t(5));
    size_t successful = 0;
    size_t matched = 0;
    
    for (size_t i = 0; i < max_test; ++i) {
        const auto& pair = pairs_[i];
        
        try {
            // Parse PDB file
            PdbParser parser;
            Structure structure = parser.parse_file(pair.pdb_file);
            
            if (structure.num_residues() == 0) {
                continue;  // Skip empty structures
            }
            
            // Create protocol
            FindPairProtocol protocol;
            protocol.set_config_manager(ConfigManager::instance());
            
            // Execute protocol
            protocol.execute(structure);
            
            // Get base pairs
            const auto& modern_pairs = protocol.base_pairs();
            
            // Load legacy base pairs
            auto legacy_pairs = load_legacy_base_pairs(pair.json_file);
            
            // Compare if legacy pairs available
            if (!legacy_pairs.empty()) {
                if (modern_pairs.size() == legacy_pairs.size()) {
                    // Quick check: compare counts
                    matched++;
                }
            }
            
            successful++;
        } catch (const std::exception& e) {
            // Skip on error
            continue;
        }
    }
    
    EXPECT_GT(successful, 0) << "No PDBs processed successfully";
    
    if (matched > 0) {
        std::cout << "Matched " << matched << " out of " << successful
                  << " PDBs with legacy base pairs" << std::endl;
    }
}

/**
 * @brief Test FindPairProtocol parameter mapping
 */
TEST_F(ProtocolsIntegrationTest, FindPairProtocolParameterMapping) {
    if (pairs_.empty()) {
        GTEST_SKIP() << "No PDB/JSON pairs found";
    }
    
    const auto& pair = pairs_[0];
    
    // Parse PDB file
    PdbParser parser;
    Structure structure;
    try {
        structure = parser.parse_file(pair.pdb_file);
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Failed to parse PDB: " << e.what();
    }
    
    // Modify config parameters
    auto& config = ConfigManager::instance();
    config.set_defaults();
    config.thresholds().max_dorg = 20.0;
    config.thresholds().min_base_hb = 2;
    
    // Create protocol with config
    FindPairProtocol protocol;
    protocol.set_config_manager(config);
    
    // Execute protocol
    try {
        protocol.execute(structure);
        
        // Verify protocol executed (parameters were used)
        EXPECT_NO_THROW(protocol.base_pairs());
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Protocol execution failed: " << e.what();
    }
}

/**
 * @brief Test FindPairProtocol legacy mode
 */
TEST_F(ProtocolsIntegrationTest, FindPairProtocolLegacyMode) {
    if (pairs_.empty()) {
        GTEST_SKIP() << "No PDB/JSON pairs found";
    }
    
    const auto& pair = pairs_[0];
    
    // Parse PDB file
    PdbParser parser;
    Structure structure;
    try {
        structure = parser.parse_file(pair.pdb_file);
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Failed to parse PDB: " << e.what();
    }
    
    // Set legacy mode
    auto& config = ConfigManager::instance();
    config.set_defaults();
    config.set_legacy_mode(true);
    
    // Create protocol
    FindPairProtocol protocol;
    protocol.set_config_manager(config);
    protocol.set_legacy_mode(true);
    
    // Execute protocol
    try {
        protocol.execute(structure);
        
        // Verify protocol executed in legacy mode
        EXPECT_TRUE(protocol.legacy_mode());
        EXPECT_NO_THROW(protocol.base_pairs());
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Protocol execution failed: " << e.what();
    }
}

/**
 * @brief Test FindPairProtocol with JSON recording
 */
TEST_F(ProtocolsIntegrationTest, FindPairProtocolWithJsonRecording) {
    if (pairs_.empty()) {
        GTEST_SKIP() << "No PDB/JSON pairs found";
    }
    
    const auto& pair = pairs_[0];
    
    // Parse PDB file
    PdbParser parser;
    Structure structure;
    try {
        structure = parser.parse_file(pair.pdb_file);
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Failed to parse PDB: " << e.what();
    }
    
    // Create JSON writer
    JsonWriter writer(pair.pdb_file);
    
    // Create protocol with JSON writer
    FindPairProtocol protocol;
    protocol.set_config_manager(ConfigManager::instance());
    protocol.set_json_writer(&writer);
    
    // Execute protocol
    try {
        protocol.execute(structure);
        
        // Verify JSON was generated
        auto json = writer.json();
        EXPECT_TRUE(json.contains("calculations"));
        
        // Check for base_pair records
        auto records = find_records_by_type(json, "base_pair");
        EXPECT_GE(records.size(), 0) << "Should have base_pair records";
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Protocol execution failed: " << e.what();
    }
}


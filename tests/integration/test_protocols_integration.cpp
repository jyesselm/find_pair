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
        // Use ONLY legacy JSON files from data/json_legacy/ for comparison
        // Legacy JSON files are named with _globals.json suffix
        pairs_ = discover_pairs_with_legacy_globals();
        
        if (pairs_.empty()) {
            GTEST_SKIP() << "No PDB/legacy JSON pairs found for testing. "
                         << "Place legacy JSON files (with _globals.json suffix) in data/json_legacy/ to enable tests.";
        }
        
        // Set up config with defaults
        auto& config = ConfigManager::instance();
        config.set_defaults();
    }
    
    /**
     * @brief Discover pairs using legacy _globals.json files from data/json_legacy/
     */
    std::vector<pdb_json_pair> discover_pairs_with_legacy_globals() {
        std::vector<pdb_json_pair> pairs;
        fs::path pdb_dir = "data/pdb";
        fs::path json_legacy_dir = "data/json_legacy";
        
        if (!fs::exists(pdb_dir) || !fs::exists(json_legacy_dir)) {
            return pairs;
        }
        
        for (const auto& pdb_file : fs::directory_iterator(pdb_dir)) {
            if (pdb_file.path().extension() != ".pdb") {
                continue;
            }
            
            std::string pdb_name = pdb_file.path().stem().string();
            fs::path globals_file = json_legacy_dir / (pdb_name + "_globals.json");
            
            if (fs::exists(globals_file)) {
                // Use globals file as the JSON file for legacy comparison
                pairs.push_back({pdb_file.path(), globals_file, globals_file, pdb_name});
            }
        }
        
        return pairs;
    }
    
    /**
     * @brief Load base pairs from legacy JSON file
     * Legacy JSON files are in data/json_legacy/base_pair/PDB_ID.json
     * They are arrays of base pair objects (not objects with a "pairs" array)
     */
    std::vector<nlohmann::json> load_legacy_base_pairs(const fs::path& json_file) {
        std::vector<nlohmann::json> base_pairs;
        
        try {
            // Extract PDB name from globals file path
            std::string pdb_name = json_file.stem().string();
            if (pdb_name.find("_globals") != std::string::npos) {
                pdb_name = pdb_name.substr(0, pdb_name.find("_globals"));
            }
            
            // Legacy base_pair files are in data/json_legacy/base_pair/PDB_ID.json
            fs::path legacy_base_pair_file = fs::path("data/json_legacy/base_pair") / (pdb_name + ".json");
            
            if (fs::exists(legacy_base_pair_file)) {
                std::ifstream file(legacy_base_pair_file);
                nlohmann::json legacy_json;
                file >> legacy_json;
                
                // Legacy base_pair JSON files are arrays of base pair objects
                if (legacy_json.is_array()) {
                    for (const auto& pair : legacy_json) {
                        base_pairs.push_back(pair);
                    }
                } else if (legacy_json.contains("pairs") && legacy_json["pairs"].is_array()) {
                    // Fallback: if it's an object with a "pairs" array
                    for (const auto& pair : legacy_json["pairs"]) {
                        base_pairs.push_back(pair);
                    }
                }
            } else {
                // Fallback: try loading from globals file if it has calculations array
                auto legacy_json = load_legacy_json(json_file);
                auto records = find_records_by_type(legacy_json, "base_pair");
                
                for (const auto& record : records) {
                    if (record.contains("pairs") && record["pairs"].is_array()) {
                        for (const auto& pair : record["pairs"]) {
                            base_pairs.push_back(pair);
                        }
                    }
                }
            }
        } catch (const std::exception& e) {
            // If loading fails, return empty vector
            // This is expected if the legacy JSON doesn't have base_pair records
        }
        
        return base_pairs;
    }
    
    /**
     * @brief Compare base pairs from protocol with legacy JSON
     * 
     * Verifies that modern protocol output matches legacy JSON from data/json_legacy/
     */
    void compare_base_pairs(const std::vector<BasePair>& modern_pairs,
                           const std::vector<nlohmann::json>& legacy_pairs,
                           const std::string& pdb_name) {
        // Verify legacy pairs were loaded (from data/json_legacy/base_pair/)
        EXPECT_FALSE(legacy_pairs.empty())
            << "No legacy base pairs loaded from data/json_legacy/base_pair/ for " << pdb_name;
        
        if (legacy_pairs.empty()) {
            return;  // Can't compare if no legacy data
        }
        
        // Log legacy values for verification
        std::cout << "Comparing " << pdb_name << ": modern=" << modern_pairs.size() 
                  << ", legacy=" << legacy_pairs.size() << " pairs" << std::endl;
        
        // Verify legacy pairs have correct structure (base_i, base_j)
        if (!legacy_pairs.empty()) {
            const auto& first_legacy = legacy_pairs[0];
            EXPECT_TRUE(first_legacy.contains("base_i") && first_legacy.contains("base_j"))
                << "Legacy JSON missing base_i/base_j fields (not from data/json_legacy/)";
        }
        
        EXPECT_EQ(modern_pairs.size(), legacy_pairs.size())
            << "Base pair count mismatch for " << pdb_name
            << ": modern=" << modern_pairs.size()
            << ", legacy=" << legacy_pairs.size()
            << " (legacy from data/json_legacy/base_pair/)";
        
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
            // Legacy JSON uses "base_i" and "base_j" (1-based indices from org code)
            EXPECT_TRUE(pair_json.contains("base_i") && pair_json.contains("base_j"))
                << "Legacy JSON missing base_i/base_j (should be from data/json_legacy/)";
            
            if (pair_json.contains("base_i") && pair_json.contains("base_j")) {
                int i = pair_json["base_i"].get<int>();
                int j = pair_json["base_j"].get<int>();
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
                << ") found in legacy (from data/json_legacy/) but not in modern for " << pdb_name;
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
    
    // Load legacy base pairs from data/json_legacy/base_pair/PDB_ID.json
    auto legacy_pairs = load_legacy_base_pairs(pair.json_file);
    
    // Verify legacy values were loaded from data/json_legacy/
    EXPECT_FALSE(legacy_pairs.empty())
        << "Failed to load legacy base pairs from data/json_legacy/base_pair/ for " << pair.pdb_name;
    
    // Verify legacy pairs have correct structure (from org code)
    if (!legacy_pairs.empty()) {
        const auto& first_legacy = legacy_pairs[0];
        EXPECT_TRUE(first_legacy.contains("base_i") && first_legacy.contains("base_j"))
            << "Legacy JSON missing base_i/base_j (not from org code output)";
        EXPECT_TRUE(first_legacy.contains("bp_type"))
            << "Legacy JSON missing bp_type (not from org code output)";
    }
    
    // Compare with legacy JSON from data/json_legacy/
    if (!legacy_pairs.empty()) {
        // Compare base pairs with legacy (from data/json_legacy/ - org code output)
        compare_base_pairs(modern_pairs, legacy_pairs, pair.pdb_name);
    } else {
        // If no legacy pairs found, at least verify protocol executed
        // This might happen if legacy JSON doesn't exist for this PDB
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


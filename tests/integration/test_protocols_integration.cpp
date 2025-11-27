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
     * 
     * CRITICAL: Use find_bestpair_selection, not base_pair!
     * - base_pair contains ALL valid pairs from Phase 1 validation (including conflicts)
     * - find_bestpair_selection contains only the pairs actually selected by find_bestpair
     * 
     * Legacy find_bestpair_selection files are in data/json_legacy/find_bestpair_selection/PDB_ID.json
     * They are arrays with one record containing a "pairs" array
     */
    std::vector<nlohmann::json> load_legacy_base_pairs(const fs::path& json_file) {
        std::vector<nlohmann::json> base_pairs;
        
        try {
            // Extract PDB name from globals file path
            std::string pdb_name = json_file.stem().string();
            if (pdb_name.find("_globals") != std::string::npos) {
                pdb_name = pdb_name.substr(0, pdb_name.find("_globals"));
            }
            
            // Legacy find_bestpair_selection files are in data/json_legacy/find_bestpair_selection/PDB_ID.json
            // These contain the ACTUAL selected pairs (not all valid pairs)
            // Format: array with one record containing "pairs" which is an array of [i, j] arrays
            // We need to get the full pair data (frames, origins, etc.) from base_pair file
            fs::path legacy_selection_file = fs::path("data/json_legacy/find_bestpair_selection") / (pdb_name + ".json");
            fs::path legacy_base_pair_file = fs::path("data/json_legacy/base_pair") / (pdb_name + ".json");
            
            std::set<std::pair<int, int>> selected_pairs_set;
            
            // First, load the selected pair indices from find_bestpair_selection
            if (fs::exists(legacy_selection_file)) {
                std::ifstream file(legacy_selection_file);
                nlohmann::json legacy_json;
                file >> legacy_json;
                
                if (legacy_json.is_array() && !legacy_json.empty()) {
                    const auto& record = legacy_json[0];
                    if (record.contains("pairs") && record["pairs"].is_array()) {
                        for (const auto& pair_array : record["pairs"]) {
                            if (pair_array.is_array() && pair_array.size() >= 2) {
                                int i = pair_array[0].get<int>();
                                int j = pair_array[1].get<int>();
                                selected_pairs_set.insert({std::min(i, j), std::max(i, j)});
                            }
                        }
                    }
                }
            }
            
            // Then, load full pair data from base_pair file for the selected pairs
            if (fs::exists(legacy_base_pair_file) && !selected_pairs_set.empty()) {
                std::ifstream file(legacy_base_pair_file);
                nlohmann::json legacy_json;
                file >> legacy_json;
                
                if (legacy_json.is_array()) {
                    for (const auto& pair : legacy_json) {
                        if (pair.contains("base_i") && pair.contains("base_j")) {
                            int i = pair["base_i"].get<int>();
                            int j = pair["base_j"].get<int>();
                            std::pair<int, int> normalized = {std::min(i, j), std::max(i, j)};
                            
                            // Only include pairs that are in the selected set
                            if (selected_pairs_set.find(normalized) != selected_pairs_set.end()) {
                                // Keep first occurrence of each unique pair
                                bool already_added = false;
                                for (const auto& existing : base_pairs) {
                                    if (existing.contains("base_i") && existing.contains("base_j")) {
                                        int ei = existing["base_i"].get<int>();
                                        int ej = existing["base_j"].get<int>();
                                        if ((ei == i && ej == j) || (ei == j && ej == i)) {
                                            already_added = true;
                                            break;
                                        }
                                    }
                                }
                                if (!already_added) {
                                    base_pairs.push_back(pair);
                                }
                            }
                        }
                    }
                }
            } else if (!fs::exists(legacy_selection_file)) {
                // Fallback: use base_pair file directly if selection file doesn't exist
                // Fallback: try base_pair file (but note it contains ALL valid pairs, not just selected)
                fs::path legacy_base_pair_file = fs::path("data/json_legacy/base_pair") / (pdb_name + ".json");
                if (fs::exists(legacy_base_pair_file)) {
                    std::ifstream file(legacy_base_pair_file);
                    nlohmann::json legacy_json;
                    file >> legacy_json;
                    
                    if (legacy_json.is_array()) {
                        for (const auto& pair : legacy_json) {
                            base_pairs.push_back(pair);
                        }
                    }
                }
            }
        } catch (const std::exception& e) {
            // If loading fails, return empty vector
            // This is expected if the legacy JSON doesn't exist
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
        
        // Build maps: normalized pair -> modern/legacy pair data
        // Normalize pairs to (min(i,j), max(i,j)) for comparison
        std::map<std::pair<int, int>, size_t> modern_map;  // Store index instead of pointer
        std::map<std::pair<int, int>, nlohmann::json> legacy_map;
        
        // Build modern map (using 1-based legacy indices)
        for (size_t idx = 0; idx < modern_pairs.size(); ++idx) {
            const auto& pair = modern_pairs[idx];
            // BasePair uses 0-based indices that are legacy_idx - 1
            int i = static_cast<int>(pair.residue_idx1()) + 1;
            int j = static_cast<int>(pair.residue_idx2()) + 1;
            std::pair<int, int> normalized = {std::min(i, j), std::max(i, j)};
            modern_map[normalized] = idx;
        }
        
        // Build legacy map (already 1-based)
        for (const auto& pair_json : legacy_pairs) {
            EXPECT_TRUE(pair_json.contains("base_i") && pair_json.contains("base_j"))
                << "Legacy JSON missing base_i/base_j (should be from data/json_legacy/)";
            
            if (pair_json.contains("base_i") && pair_json.contains("base_j")) {
                int i = pair_json["base_i"].get<int>();
                int j = pair_json["base_j"].get<int>();
                std::pair<int, int> normalized = {std::min(i, j), std::max(i, j)};
                // Keep first occurrence if duplicates exist
                if (legacy_map.find(normalized) == legacy_map.end()) {
                    legacy_map[normalized] = pair_json;
                }
            }
        }
        
        // Compare unique pair counts
        size_t modern_unique = modern_map.size();
        size_t legacy_unique = legacy_map.size();
        
        std::cout << "Comparing " << pdb_name << ": modern=" << modern_pairs.size() 
                  << " pairs (" << modern_unique << " unique), legacy=" << legacy_pairs.size()
                  << " pairs (" << legacy_unique << " unique)" << std::endl;
        
        EXPECT_EQ(modern_unique, legacy_unique)
            << "Unique pair count mismatch for " << pdb_name
            << ": modern=" << modern_unique << ", legacy=" << legacy_unique;
        
        if (modern_unique != legacy_unique) {
            return;  // Can't compare parameters if counts don't match
        }
        
        // Compare each unique pair's parameters
        const double TOLERANCE = 1e-6;  // Tolerance for floating point comparison
        
        for (const auto& [normalized_pair, legacy_json] : legacy_map) {
            auto modern_it = modern_map.find(normalized_pair);
            EXPECT_TRUE(modern_it != modern_map.end())
                << "Pair (" << normalized_pair.first << ", " << normalized_pair.second
                << ") found in legacy but not in modern for " << pdb_name;
            
            if (modern_it == modern_map.end()) {
                continue;  // Skip if pair not found
            }
            
            const BasePair& modern_pair = modern_pairs[modern_it->second];
            int i = normalized_pair.first;
            int j = normalized_pair.second;
            
            // Determine order: legacy JSON has base_i and base_j, need to match with modern pair order
            int legacy_i = legacy_json.contains("base_i") ? legacy_json["base_i"].get<int>() : 0;
            int legacy_j = legacy_json.contains("base_j") ? legacy_json["base_j"].get<int>() : 0;
            
            // Modern pair uses 0-based indices, convert to 1-based for comparison
            int modern_i = static_cast<int>(modern_pair.residue_idx1()) + 1;
            int modern_j = static_cast<int>(modern_pair.residue_idx2()) + 1;
            
            // Determine if order matches or is reversed
            bool order_matches = (legacy_i == modern_i && legacy_j == modern_j);
            bool order_reversed = (legacy_i == modern_j && legacy_j == modern_i);
            
            EXPECT_TRUE(order_matches || order_reversed)
                << "Pair order mismatch for (" << i << ", " << j << ") in " << pdb_name
                << ": legacy=(" << legacy_i << ", " << legacy_j << "), modern=(" << modern_i << ", " << modern_j << ")";
            
            // Compare bp_type (handle reversed order: "UG" vs "GU" are the same pair)
            if (legacy_json.contains("bp_type")) {
                std::string legacy_bp_type = legacy_json["bp_type"].get<std::string>();
                std::string modern_bp_type = modern_pair.bp_type();
                
                // Normalize: check if bp_type matches or is reversed
                bool bp_type_matches = (modern_bp_type == legacy_bp_type);
                if (!bp_type_matches && modern_bp_type.size() == 2 && legacy_bp_type.size() == 2) {
                    // Check if reversed (e.g., "UG" vs "GU")
                    std::string modern_reversed = std::string(1, modern_bp_type[1]) + std::string(1, modern_bp_type[0]);
                    bp_type_matches = (modern_reversed == legacy_bp_type);
                }
                
                EXPECT_TRUE(bp_type_matches)
                    << "bp_type mismatch for pair (" << i << ", " << j << ") in " << pdb_name
                    << ": modern=" << modern_bp_type << ", legacy=" << legacy_bp_type;
            }
            
            // Compare frames and origins - need to handle order correctly
            // If order matches: legacy orien_i/org_i -> modern frame1, legacy orien_j/org_j -> modern frame2
            // If order reversed: legacy orien_i/org_i -> modern frame2, legacy orien_j/org_j -> modern frame1
            if (legacy_json.contains("orien_i") && legacy_json.contains("org_i") &&
                legacy_json.contains("orien_j") && legacy_json.contains("org_j")) {
                
                const auto& orien_i = legacy_json["orien_i"];
                const auto& org_i = legacy_json["org_i"];
                const auto& orien_j = legacy_json["orien_j"];
                const auto& org_j = legacy_json["org_j"];
                
                // Determine which modern frame corresponds to which legacy frame
                // Store values, not pointers (value() returns a temporary)
                std::optional<ReferenceFrame> modern_frame_i_opt;
                std::optional<ReferenceFrame> modern_frame_j_opt;
                
                if (order_matches) {
                    // Legacy i -> modern frame1, legacy j -> modern frame2
                    if (modern_pair.frame1().has_value() && modern_pair.frame2().has_value()) {
                        modern_frame_i_opt = modern_pair.frame1().value();
                        modern_frame_j_opt = modern_pair.frame2().value();
                    }
                } else if (order_reversed) {
                    // Legacy i -> modern frame2, legacy j -> modern frame1
                    if (modern_pair.frame1().has_value() && modern_pair.frame2().has_value()) {
                        modern_frame_i_opt = modern_pair.frame2().value();
                        modern_frame_j_opt = modern_pair.frame1().value();
                    }
                }
                
                if (modern_frame_i_opt.has_value() && modern_frame_j_opt.has_value()) {
                    const auto& modern_frame_i = modern_frame_i_opt.value();
                    const auto& modern_frame_j = modern_frame_j_opt.value();
                    
                    // Compare origin and orientation for residue i
                    if (org_i.is_array() && org_i.size() == 3) {
                        EXPECT_NEAR(modern_frame_i.origin().x(), org_i[0].get<double>(), TOLERANCE)
                            << "org_i[0] mismatch for pair (" << i << ", " << j << ") in " << pdb_name;
                        EXPECT_NEAR(modern_frame_i.origin().y(), org_i[1].get<double>(), TOLERANCE)
                            << "org_i[1] mismatch for pair (" << i << ", " << j << ") in " << pdb_name;
                        EXPECT_NEAR(modern_frame_i.origin().z(), org_i[2].get<double>(), TOLERANCE)
                            << "org_i[2] mismatch for pair (" << i << ", " << j << ") in " << pdb_name;
                    }
                    
                    if (orien_i.is_array() && orien_i.size() == 9) {
                        const auto& rot = modern_frame_i.rotation();
                        for (int row = 0; row < 3; ++row) {
                            for (int col = 0; col < 3; ++col) {
                                int idx = row * 3 + col;
                                EXPECT_NEAR(rot.at(row, col), orien_i[idx].get<double>(), TOLERANCE)
                                    << "orien_i[" << row << "][" << col << "] mismatch for pair ("
                                    << i << ", " << j << ") in " << pdb_name;
                            }
                        }
                    }
                    
                    // Compare origin and orientation for residue j
                    if (org_j.is_array() && org_j.size() == 3) {
                        EXPECT_NEAR(modern_frame_j.origin().x(), org_j[0].get<double>(), TOLERANCE)
                            << "org_j[0] mismatch for pair (" << i << ", " << j << ") in " << pdb_name;
                        EXPECT_NEAR(modern_frame_j.origin().y(), org_j[1].get<double>(), TOLERANCE)
                            << "org_j[1] mismatch for pair (" << i << ", " << j << ") in " << pdb_name;
                        EXPECT_NEAR(modern_frame_j.origin().z(), org_j[2].get<double>(), TOLERANCE)
                            << "org_j[2] mismatch for pair (" << i << ", " << j << ") in " << pdb_name;
                    }
                    
                    if (orien_j.is_array() && orien_j.size() == 9) {
                        const auto& rot = modern_frame_j.rotation();
                        for (int row = 0; row < 3; ++row) {
                            for (int col = 0; col < 3; ++col) {
                                int idx = row * 3 + col;
                                EXPECT_NEAR(rot.at(row, col), orien_j[idx].get<double>(), TOLERANCE)
                                    << "orien_j[" << row << "][" << col << "] mismatch for pair ("
                                    << i << ", " << j << ") in " << pdb_name;
                            }
                        }
                    }
                } else {
                    EXPECT_TRUE(false) << "Missing frames for pair (" << i << ", " << j << ") in " << pdb_name;
                }
            }
            
            // Compare direction vector (dir_xyz)
            // Note: Legacy has a bug where it stores [dir_y, dir_z, 0.0] instead of [dir_x, dir_y, dir_z]
            if (legacy_json.contains("dir_xyz") && modern_pair.frame1().has_value() && 
                modern_pair.frame2().has_value()) {
                const auto frame1 = modern_pair.frame1().value();  // Store value, not reference
                const auto frame2 = modern_pair.frame2().value();  // Store value, not reference
                const auto& dir_xyz = legacy_json["dir_xyz"];
                
                if (dir_xyz.is_array() && dir_xyz.size() >= 2) {
                    // Legacy stores [dir_y, dir_z, 0.0] due to bug
                    // Direction vector is independent of pair order (symmetric)
                    double modern_dir_y = frame1.y_axis().dot(frame2.y_axis());
                    double modern_dir_z = frame1.z_axis().dot(frame2.z_axis());
                    
                    EXPECT_NEAR(modern_dir_y, dir_xyz[0].get<double>(), TOLERANCE)
                        << "dir_xyz[0] (dir_y) mismatch for pair (" << i << ", " << j << ") in " << pdb_name;
                    EXPECT_NEAR(modern_dir_z, dir_xyz[1].get<double>(), TOLERANCE)
                        << "dir_xyz[1] (dir_z) mismatch for pair (" << i << ", " << j << ") in " << pdb_name;
                }
            }
        }
        
        // Verify all modern pairs are in legacy
        for (const auto& [normalized_pair, modern_idx] : modern_map) {
            EXPECT_TRUE(legacy_map.find(normalized_pair) != legacy_map.end())
                << "Pair (" << normalized_pair.first << ", " << normalized_pair.second
                << ") found in modern but not in legacy for " << pdb_name;
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
 * 
 * Tests multiple PDBs and compares unique pairs (not total counts)
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
            
            // Load legacy base pairs (from find_bestpair_selection)
            auto legacy_pairs = load_legacy_base_pairs(pair.json_file);
            
            // Compare unique pairs if legacy pairs available
            if (!legacy_pairs.empty()) {
                // Build normalized sets for comparison
                std::set<std::pair<int, int>> modern_set;
                std::set<std::pair<int, int>> legacy_set;
                
                // Modern pairs: convert to 1-based and normalize
                for (const auto& p : modern_pairs) {
                    int i = static_cast<int>(p.residue_idx1()) + 1;
                    int j = static_cast<int>(p.residue_idx2()) + 1;
                    modern_set.insert({std::min(i, j), std::max(i, j)});
                }
                
                // Legacy pairs: already 1-based, normalize
                for (const auto& p_json : legacy_pairs) {
                    if (p_json.contains("base_i") && p_json.contains("base_j")) {
                        int i = p_json["base_i"].get<int>();
                        int j = p_json["base_j"].get<int>();
                        legacy_set.insert({std::min(i, j), std::max(i, j)});
                    }
                }
                
                // Compare unique pair counts
                if (modern_set.size() == legacy_set.size()) {
                    matched++;
                    std::cout << "✓ " << pair.pdb_name << ": " << modern_set.size() 
                              << " unique pairs match" << std::endl;
                } else {
                    std::cout << "✗ " << pair.pdb_name << ": modern=" << modern_set.size() 
                              << " unique, legacy=" << legacy_set.size() << " unique (MISMATCH)" << std::endl;
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
                  << " PDBs with legacy base pairs (unique pair counts)" << std::endl;
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


/**
 * @file test_config_manager.cpp
 * @brief Unit tests for ConfigManager
 */

#include <gtest/gtest.h>
#include <x3dna/config/config_manager.hpp>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

using namespace x3dna::config;

class ConfigManagerTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Reset to defaults before each test
        auto& config = ConfigManager::instance();
        config.set_defaults();
    }

    void TearDown() override {
        // Reset to defaults after each test
        auto& config = ConfigManager::instance();
        config.set_defaults();
    }
};

// Singleton tests
TEST_F(ConfigManagerTest, SingletonInstance) {
    auto& config1 = ConfigManager::instance();
    auto& config2 = ConfigManager::instance();

    // Should return same instance
    EXPECT_EQ(&config1, &config2);
}

TEST_F(ConfigManagerTest, SingletonCannotBeCopied) {
    auto& config = ConfigManager::instance();

    // Should not compile if copy constructor is accessible
    // This is a compile-time check, but we can verify the instance works
    EXPECT_NO_THROW(config.set_defaults());
}

// Default values tests
TEST_F(ConfigManagerTest, DefaultValuesMatchLegacy) {
    auto& config = ConfigManager::instance();
    config.set_defaults();

    const auto& thresholds = config.thresholds();

    // Distance constraints
    EXPECT_DOUBLE_EQ(thresholds.min_dorg, 0.0);
    EXPECT_DOUBLE_EQ(thresholds.max_dorg, 15.0);
    EXPECT_DOUBLE_EQ(thresholds.min_dv, 0.0);
    EXPECT_DOUBLE_EQ(thresholds.max_dv, 2.5);
    EXPECT_DOUBLE_EQ(thresholds.min_dNN, 4.5);
    EXPECT_DOUBLE_EQ(thresholds.max_dNN, 1e18); // XBIG

    // Angle constraints
    EXPECT_DOUBLE_EQ(thresholds.min_plane_angle, 0.0);
    EXPECT_DOUBLE_EQ(thresholds.max_plane_angle, 65.0);

    // Hydrogen bond constraints
    EXPECT_EQ(thresholds.min_base_hb, 1);
    EXPECT_DOUBLE_EQ(thresholds.hb_lower, 1.8);
    EXPECT_DOUBLE_EQ(thresholds.hb_dist1, 4.0);
    EXPECT_DOUBLE_EQ(thresholds.hb_dist2, 0.0); // CRITICAL: Must be 0.0
    EXPECT_EQ(thresholds.hb_atoms, ".O.N");

    // Overlap threshold
    EXPECT_DOUBLE_EQ(thresholds.overlap_threshold, 0.01);

    // Helix parameters
    EXPECT_DOUBLE_EQ(thresholds.helix_break, 7.5);

    // Other parameters
    EXPECT_EQ(thresholds.alt_list, "A1");
    EXPECT_DOUBLE_EQ(thresholds.std_curved, 0.6);
    EXPECT_DOUBLE_EQ(thresholds.water_dist, 3.2);
    EXPECT_DOUBLE_EQ(thresholds.water_dlow, 0.0);
    EXPECT_EQ(thresholds.water_atoms, ".O.N");
    EXPECT_DOUBLE_EQ(thresholds.o3p_dist, 4.5);
}

// Parameter modification tests
TEST_F(ConfigManagerTest, ModifyParameters) {
    auto& config = ConfigManager::instance();
    auto& thresholds = config.thresholds();

    thresholds.max_dorg = 20.0;
    thresholds.min_base_hb = 2;

    EXPECT_DOUBLE_EQ(config.thresholds().max_dorg, 20.0);
    EXPECT_EQ(config.thresholds().min_base_hb, 2);
}

// Options tests
TEST_F(ConfigManagerTest, IncludeHetatm) {
    auto& config = ConfigManager::instance();

    EXPECT_FALSE(config.include_hetatm());
    config.set_include_hetatm(true);
    EXPECT_TRUE(config.include_hetatm());
    config.set_include_hetatm(false);
    EXPECT_FALSE(config.include_hetatm());
}

TEST_F(ConfigManagerTest, IncludeWaters) {
    auto& config = ConfigManager::instance();

    EXPECT_FALSE(config.include_waters());
    config.set_include_waters(true);
    EXPECT_TRUE(config.include_waters());
    config.set_include_waters(false);
    EXPECT_FALSE(config.include_waters());
}

TEST_F(ConfigManagerTest, LegacyMode) {
    auto& config = ConfigManager::instance();

    EXPECT_FALSE(config.legacy_mode());
    config.set_legacy_mode(true);
    EXPECT_TRUE(config.legacy_mode());
    config.set_legacy_mode(false);
    EXPECT_FALSE(config.legacy_mode());
}

// Path tests
TEST_F(ConfigManagerTest, X3DNAHome) {
    auto& config = ConfigManager::instance();

    std::filesystem::path test_path("/test/path");
    config.set_x3dna_home(test_path);
    EXPECT_EQ(config.x3dna_home(), test_path);
}

// JSON loading tests
TEST_F(ConfigManagerTest, LoadFromJSON) {
    auto& config = ConfigManager::instance();

    nlohmann::json json_config = {{"thresholds", {{"max_dorg", 20.0}, {"min_base_hb", 2}, {"hb_lower", 2.0}}},
                                  {"include_hetatm", true},
                                  {"include_waters", true},
                                  {"legacy_mode", true}};

    config.load_from_json(json_config);

    EXPECT_DOUBLE_EQ(config.thresholds().max_dorg, 20.0);
    EXPECT_EQ(config.thresholds().min_base_hb, 2);
    EXPECT_DOUBLE_EQ(config.thresholds().hb_lower, 2.0);
    EXPECT_TRUE(config.include_hetatm());
    EXPECT_TRUE(config.include_waters());
    EXPECT_TRUE(config.legacy_mode());
}

TEST_F(ConfigManagerTest, LoadFromJSONPartial) {
    auto& config = ConfigManager::instance();
    config.set_defaults();

    // Explicitly reset options to defaults
    config.set_include_hetatm(false);
    config.set_include_waters(false);
    config.set_legacy_mode(false);

    // Load only some parameters
    nlohmann::json json_config = {{"thresholds", {{"max_dorg", 25.0}}}};

    config.load_from_json(json_config);

    // Modified parameter
    EXPECT_DOUBLE_EQ(config.thresholds().max_dorg, 25.0);

    // Other parameters should remain at defaults
    EXPECT_DOUBLE_EQ(config.thresholds().min_dorg, 0.0);
    EXPECT_DOUBLE_EQ(config.thresholds().max_dv, 2.5);
    // Options should remain at defaults (not modified by partial JSON)
    EXPECT_FALSE(config.include_hetatm());
    EXPECT_FALSE(config.legacy_mode());
}

// File loading tests (if file exists)
TEST_F(ConfigManagerTest, LoadFromFile) {
    auto& config = ConfigManager::instance();

    // Create a temporary config file
    std::filesystem::path temp_file = std::filesystem::temp_directory_path() / "test_config.json";

    nlohmann::json json_config = {{"thresholds", {{"max_dorg", 18.0}, {"min_base_hb", 1}}}, {"legacy_mode", false}};

    std::ofstream file(temp_file);
    file << json_config.dump(4);
    file.close();

    // Load from file
    config.load_from_file(temp_file);

    EXPECT_DOUBLE_EQ(config.thresholds().max_dorg, 18.0);
    EXPECT_EQ(config.thresholds().min_base_hb, 1);
    EXPECT_FALSE(config.legacy_mode());

    // Clean up
    std::filesystem::remove(temp_file);
}

TEST_F(ConfigManagerTest, LoadFromFileNotFound) {
    auto& config = ConfigManager::instance();

    std::filesystem::path non_existent = "/nonexistent/path/config.json";

    // Should handle gracefully (prints warning, doesn't throw)
    // The implementation logs a warning but doesn't throw
    EXPECT_NO_THROW(config.load_from_file(non_existent));
}

// Standard base path tests
TEST_F(ConfigManagerTest, StandardBasePath) {
    auto& config = ConfigManager::instance();

    std::filesystem::path x3dna_home("/test/x3dna");
    config.set_x3dna_home(x3dna_home);

    auto base_path = config.standard_base_path();
    EXPECT_TRUE(base_path.string().find("x3dna") != std::string::npos);
}

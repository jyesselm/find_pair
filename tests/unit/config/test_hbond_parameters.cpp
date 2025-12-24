/**
 * @file test_hbond_parameters.cpp
 * @brief Unit tests for HBondParameters and HBondParametersLoader
 */

#include <gtest/gtest.h>
#include <x3dna/config/hbond_parameters.hpp>
#include <x3dna/config/hbond_parameters_loader.hpp>
#include <x3dna/config/resource_locator.hpp>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

using namespace x3dna::config;

class HBondParametersTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Clear cached parameters before each test
        HBondParametersLoader::reload();
    }
};

// === Default values tests ===

TEST_F(HBondParametersTest, DefaultsHaveCorrectValues) {
    auto params = HBondParameters::defaults();

    // Detection distance thresholds
    EXPECT_DOUBLE_EQ(params.detection.distance.min, 1.8);
    EXPECT_DOUBLE_EQ(params.detection.distance.base_base_max, 4.0);
    EXPECT_DOUBLE_EQ(params.detection.distance.base_backbone_max, 3.5);
    EXPECT_DOUBLE_EQ(params.detection.distance.conflict_filter, 4.5);

    // Detection thresholds
    EXPECT_DOUBLE_EQ(params.detection.thresholds.good_bond.min, 2.5);
    EXPECT_DOUBLE_EQ(params.detection.thresholds.good_bond.max, 3.5);
    EXPECT_DOUBLE_EQ(params.detection.thresholds.post_validation_max, 3.6);

    // Detection options
    EXPECT_FALSE(params.detection.options.enable_angle_filtering);
    EXPECT_FALSE(params.detection.options.include_backbone_backbone);

    // Geometry
    EXPECT_DOUBLE_EQ(params.geometry.donor_angle.min, 90.0);
    EXPECT_DOUBLE_EQ(params.geometry.donor_angle.ideal, 165.0);
    EXPECT_DOUBLE_EQ(params.geometry.acceptor_angle.min, 70.0);
    EXPECT_DOUBLE_EQ(params.geometry.acceptor_angle.ideal_sp2, 130.0);
    EXPECT_DOUBLE_EQ(params.geometry.acceptor_angle.ideal_sp3, 110.0);

    // Scoring
    EXPECT_DOUBLE_EQ(params.scoring.distance.ideal, 2.9);
    EXPECT_DOUBLE_EQ(params.scoring.distance.sigma, 0.3);
    EXPECT_DOUBLE_EQ(params.scoring.weights.distance, 0.45);
    EXPECT_DOUBLE_EQ(params.scoring.weights.donor_angle, 0.30);
    EXPECT_DOUBLE_EQ(params.scoring.weights.acceptor_angle, 0.25);

    // Quality tiers
    EXPECT_DOUBLE_EQ(params.quality_tiers.excellent_min, 90.0);
    EXPECT_DOUBLE_EQ(params.quality_tiers.standard_min, 70.0);
    EXPECT_DOUBLE_EQ(params.quality_tiers.acceptable_min, 50.0);
    EXPECT_DOUBLE_EQ(params.quality_tiers.questionable_min, 30.0);
}

// === JSON loading tests ===

TEST_F(HBondParametersTest, LoadFromJsonBasic) {
    nlohmann::json json = R"({
        "detection": {
            "distance": {
                "min": 2.0,
                "base_base_max": 3.8
            }
        },
        "scoring": {
            "distance": {
                "ideal": 2.85
            }
        }
    })"_json;

    auto params = HBondParametersLoader::load_from_json(json);

    EXPECT_DOUBLE_EQ(params.detection.distance.min, 2.0);
    EXPECT_DOUBLE_EQ(params.detection.distance.base_base_max, 3.8);
    EXPECT_DOUBLE_EQ(params.scoring.distance.ideal, 2.85);

    // Unspecified values should be defaults
    EXPECT_DOUBLE_EQ(params.detection.distance.base_backbone_max, 3.5);
    EXPECT_DOUBLE_EQ(params.scoring.distance.sigma, 0.3);
}

TEST_F(HBondParametersTest, LoadFromJsonNestedRanges) {
    nlohmann::json json = R"({
        "detection": {
            "thresholds": {
                "good_bond": { "min": 2.6, "max": 3.4 },
                "nonstandard": { "min": 2.7, "max": 3.1 }
            }
        }
    })"_json;

    auto params = HBondParametersLoader::load_from_json(json);

    EXPECT_DOUBLE_EQ(params.detection.thresholds.good_bond.min, 2.6);
    EXPECT_DOUBLE_EQ(params.detection.thresholds.good_bond.max, 3.4);
    EXPECT_DOUBLE_EQ(params.detection.thresholds.nonstandard.min, 2.7);
    EXPECT_DOUBLE_EQ(params.detection.thresholds.nonstandard.max, 3.1);
}

TEST_F(HBondParametersTest, LoadFromJsonOptions) {
    nlohmann::json json = R"({
        "detection": {
            "options": {
                "enable_angle_filtering": true,
                "include_backbone_backbone": true
            }
        }
    })"_json;

    auto params = HBondParametersLoader::load_from_json(json);

    EXPECT_TRUE(params.detection.options.enable_angle_filtering);
    EXPECT_TRUE(params.detection.options.include_backbone_backbone);
    // Other options should be defaults
    EXPECT_FALSE(params.detection.options.enable_quality_scoring);
}

TEST_F(HBondParametersTest, LoadFromJsonGeometry) {
    nlohmann::json json = R"({
        "geometry": {
            "donor_angle": { "min": 100.0, "ideal": 170.0 },
            "acceptor_angle": { "min": 80.0, "ideal_sp2": 140.0, "ideal_sp3": 115.0 }
        }
    })"_json;

    auto params = HBondParametersLoader::load_from_json(json);

    EXPECT_DOUBLE_EQ(params.geometry.donor_angle.min, 100.0);
    EXPECT_DOUBLE_EQ(params.geometry.donor_angle.ideal, 170.0);
    EXPECT_DOUBLE_EQ(params.geometry.acceptor_angle.min, 80.0);
    EXPECT_DOUBLE_EQ(params.geometry.acceptor_angle.ideal_sp2, 140.0);
    EXPECT_DOUBLE_EQ(params.geometry.acceptor_angle.ideal_sp3, 115.0);
}

TEST_F(HBondParametersTest, LoadFromJsonQualityTiers) {
    nlohmann::json json = R"({
        "quality_tiers": {
            "excellent": { "min_score": 95 },
            "standard": { "min_score": 75 },
            "acceptable": { "min_score": 55 },
            "questionable": { "min_score": 35 }
        }
    })"_json;

    auto params = HBondParametersLoader::load_from_json(json);

    EXPECT_DOUBLE_EQ(params.quality_tiers.excellent_min, 95.0);
    EXPECT_DOUBLE_EQ(params.quality_tiers.standard_min, 75.0);
    EXPECT_DOUBLE_EQ(params.quality_tiers.acceptable_min, 55.0);
    EXPECT_DOUBLE_EQ(params.quality_tiers.questionable_min, 35.0);
}

// === File loading tests ===

TEST_F(HBondParametersTest, LoadFromFile) {
    // Create a temporary config file
    std::filesystem::path temp_file = std::filesystem::temp_directory_path() / "test_hbond_params.json";

    nlohmann::json json = R"({
        "detection": {
            "distance": {
                "min": 1.9,
                "base_base_max": 3.9
            }
        },
        "presets": {}
    })"_json;

    std::ofstream file(temp_file);
    file << json.dump(4);
    file.close();

    auto params = HBondParametersLoader::load_from_file(temp_file);

    EXPECT_DOUBLE_EQ(params.detection.distance.min, 1.9);
    EXPECT_DOUBLE_EQ(params.detection.distance.base_base_max, 3.9);

    // Clean up
    std::filesystem::remove(temp_file);
}

TEST_F(HBondParametersTest, LoadFromFileNotFound) {
    std::filesystem::path non_existent = "/nonexistent/hbond_params.json";

    EXPECT_THROW(HBondParametersLoader::load_from_file(non_existent), std::runtime_error);
}

// === Preset tests ===

TEST_F(HBondParametersTest, LoadPresetFromJson) {
    // Create a temporary config file with presets
    std::filesystem::path temp_file = std::filesystem::temp_directory_path() / "test_hbond_presets.json";

    nlohmann::json json = R"({
        "detection": {
            "distance": {
                "min": 1.8,
                "base_base_max": 4.0
            }
        },
        "presets": {
            "test_preset": {
                "detection": {
                    "distance": {
                        "min": 2.0,
                        "base_base_max": 3.5
                    }
                }
            }
        }
    })"_json;

    std::ofstream file(temp_file);
    file << json.dump(4);
    file.close();

    // Load from file to populate cache
    HBondParametersLoader::load_from_file(temp_file);

    // Check that preset exists
    EXPECT_TRUE(HBondParametersLoader::has_preset("test_preset"));
    EXPECT_FALSE(HBondParametersLoader::has_preset("nonexistent_preset"));

    // Load preset
    auto params = HBondParametersLoader::load_preset("test_preset");
    EXPECT_DOUBLE_EQ(params.detection.distance.min, 2.0);
    EXPECT_DOUBLE_EQ(params.detection.distance.base_base_max, 3.5);

    // Clean up
    std::filesystem::remove(temp_file);
}

TEST_F(HBondParametersTest, LoadPresetNotFound) {
    // Create a minimal config file
    std::filesystem::path temp_file = std::filesystem::temp_directory_path() / "test_hbond_empty.json";

    nlohmann::json json = R"({
        "presets": {}
    })"_json;

    std::ofstream file(temp_file);
    file << json.dump(4);
    file.close();

    // Load from file
    HBondParametersLoader::load_from_file(temp_file);

    EXPECT_THROW(HBondParametersLoader::load_preset("nonexistent"), std::runtime_error);

    // Clean up
    std::filesystem::remove(temp_file);
}

// === Available presets test ===

TEST_F(HBondParametersTest, AvailablePresets) {
    // Create a temporary config file with presets
    std::filesystem::path temp_file = std::filesystem::temp_directory_path() / "test_hbond_avail.json";

    nlohmann::json json = R"({
        "presets": {
            "_description": "test presets",
            "preset1": {},
            "preset2": {}
        }
    })"_json;

    std::ofstream file(temp_file);
    file << json.dump(4);
    file.close();

    // Load from file
    HBondParametersLoader::load_from_file(temp_file);

    auto presets = HBondParametersLoader::available_presets();

    // Should not include _description
    EXPECT_EQ(presets.size(), 2u);
    EXPECT_TRUE(std::find(presets.begin(), presets.end(), "preset1") != presets.end());
    EXPECT_TRUE(std::find(presets.begin(), presets.end(), "preset2") != presets.end());

    // Clean up
    std::filesystem::remove(temp_file);
}

// === Range struct test ===

TEST_F(HBondParametersTest, RangeDefaults) {
    Range range;
    EXPECT_DOUBLE_EQ(range.min, 0.0);
    EXPECT_DOUBLE_EQ(range.max, 0.0);
}

// === Integration test with actual config file ===

TEST_F(HBondParametersTest, LoadActualConfigFile) {
    // Try to initialize resource locator
    if (!ResourceLocator::is_initialized()) {
        // Try common paths
        std::vector<std::filesystem::path> paths = {
            "resources",
            "../resources",
            "../../resources",
        };

        for (const auto& path : paths) {
            if (std::filesystem::exists(path / "config" / "hbond_parameters.json")) {
                ResourceLocator::initialize(path);
                break;
            }
        }
    }

    if (ResourceLocator::is_initialized() &&
        ResourceLocator::config_exists("hbond_parameters.json")) {
        // Load from actual config file
        auto params = HBondParametersLoader::load();

        // Verify some key values match expected defaults
        EXPECT_GE(params.detection.distance.min, 1.0);
        EXPECT_LE(params.detection.distance.min, 3.0);
        EXPECT_GE(params.detection.distance.base_base_max, 3.0);
        EXPECT_LE(params.detection.distance.base_base_max, 5.0);
    } else {
        // Skip test if resources not found
        GTEST_SKIP() << "Resource locator not initialized, skipping actual config test";
    }
}

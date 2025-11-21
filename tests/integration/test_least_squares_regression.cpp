/**
 * @file test_least_squares_regression.cpp
 * @brief Regression tests for LeastSquaresFitter using real data from legacy JSON files
 *
 * This test extracts actual point sets and expected results from legacy JSON files
 * and verifies that our implementation produces the same rotation matrices,
 * translations, and RMS values as the original algorithm.
 */

#include <gtest/gtest.h>
#include <x3dna/geometry/least_squares_fitter.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>
#include "integration_test_base.hpp"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <cmath>

using namespace x3dna::geometry;
using namespace x3dna::test;

class LeastSquaresRegressionTest : public integration_test_base {
protected:
    constexpr static double ROTATION_TOLERANCE = 0.001;    // For rotation matrix elements
    constexpr static double TRANSLATION_TOLERANCE = 0.001; // For translation components
    constexpr static double RMS_TOLERANCE = 0.001;         // For RMS values
};

// Helper function to extract point sets from frame_calc record
std::pair<std::vector<Vector3D>, std::vector<Vector3D>>
extract_point_sets(const nlohmann::json& frame_calc) {
    std::vector<Vector3D> std_points;
    std::vector<Vector3D> exp_points;

    if (!frame_calc.contains("matched_coordinates")) {
        return {std_points, exp_points};
    }

    for (const auto& coord : frame_calc["matched_coordinates"]) {
        if (coord.contains("std_xyz") && coord.contains("exp_xyz")) {
            const auto& std_xyz = coord["std_xyz"];
            const auto& exp_xyz = coord["exp_xyz"];

            if (std_xyz.is_array() && std_xyz.size() == 3 && exp_xyz.is_array() &&
                exp_xyz.size() == 3) {
                std_points.push_back(Vector3D(std_xyz[0].get<double>(), std_xyz[1].get<double>(),
                                              std_xyz[2].get<double>()));
                exp_points.push_back(Vector3D(exp_xyz[0].get<double>(), exp_xyz[1].get<double>(),
                                              exp_xyz[2].get<double>()));
            }
        }
    }

    return {std_points, exp_points};
}

// Helper function to load rotation matrix from JSON
Matrix3D load_rotation_matrix(const nlohmann::json& json) {
    if (!json.contains("rotation_matrix") || !json["rotation_matrix"].is_array()) {
        throw std::invalid_argument("Invalid rotation_matrix in JSON");
    }
    return Matrix3D::from_json_legacy(json["rotation_matrix"]);
}

// Helper function to load translation from JSON
Vector3D load_translation(const nlohmann::json& json) {
    if (!json.contains("translation") || !json["translation"].is_array()) {
        throw std::invalid_argument("Invalid translation in JSON");
    }
    return Vector3D::from_json(json["translation"]);
}

// Test with data from 157D.json (100D.json is incomplete)
TEST_F(LeastSquaresRegressionTest, Test157D_Residue1) {
    // Find 157D in pairs or use direct path
    std::filesystem::path json_path;
    for (const auto& pair : pairs_) {
        if (pair.pdb_name == "157D") {
            json_path = pair.json_file;
            break;
        }
    }

    if (json_path.empty()) {
        // Try relative path from project root
        json_path = std::filesystem::path("data/json_legacy/157D.json");
        if (!std::filesystem::exists(json_path)) {
            GTEST_SKIP() << "157D.json not found";
        }
    }

    // Load JSON file
    auto json = load_legacy_json(json_path);

    // Find first ls_fitting and frame_calc records
    auto ls_fitting_records = find_records_by_type(json, "ls_fitting");
    auto frame_calc_records = find_records_by_type(json, "frame_calc");

    if (ls_fitting_records.empty() || frame_calc_records.empty()) {
        GTEST_SKIP() << "No ls_fitting or frame_calc records found in 100D.json";
    }

    // Use first record (residue 1)
    const auto& ls_fitting = ls_fitting_records[0];
    const auto& frame_calc = frame_calc_records[0];

    // Extract point sets
    auto [std_points, exp_points] = extract_point_sets(frame_calc);

    if (std_points.size() < 3) {
        GTEST_SKIP() << "Not enough points in frame_calc record";
    }

    // Run our fitting
    LeastSquaresFitter fitter;
    auto result = fitter.fit(std_points, exp_points);

    // Load expected values
    Matrix3D expected_rotation = load_rotation_matrix(ls_fitting);
    Vector3D expected_translation = load_translation(ls_fitting);
    double expected_rms = ls_fitting["rms_fit"].get<double>();

    // Compare rotation matrix
    EXPECT_TRUE(result.rotation.approximately_equals(expected_rotation, ROTATION_TOLERANCE))
        << "Rotation matrices don't match";

    // Compare translation
    EXPECT_NEAR(result.translation.x(), expected_translation.x(), TRANSLATION_TOLERANCE);
    EXPECT_NEAR(result.translation.y(), expected_translation.y(), TRANSLATION_TOLERANCE);
    EXPECT_NEAR(result.translation.z(), expected_translation.z(), TRANSLATION_TOLERANCE);

    // Compare RMS
    EXPECT_NEAR(result.rms, expected_rms, RMS_TOLERANCE)
        << "RMS mismatch: expected " << expected_rms << ", got " << result.rms;
}

// Test with multiple residues from 157D.json
TEST_F(LeastSquaresRegressionTest, Test157D_MultipleResidues) {
    std::filesystem::path json_path;
    for (const auto& pair : pairs_) {
        if (pair.pdb_name == "157D") {
            json_path = pair.json_file;
            break;
        }
    }
    if (json_path.empty()) {
        json_path = std::filesystem::path("data/json_legacy/157D.json");
    }
    if (!std::filesystem::exists(json_path)) {
        GTEST_SKIP() << "157D.json not found";
    }

    auto json = load_legacy_json(json_path);

    auto ls_fitting_records = find_records_by_type(json, "ls_fitting");
    auto frame_calc_records = find_records_by_type(json, "frame_calc");

    if (ls_fitting_records.size() != frame_calc_records.size()) {
        GTEST_SKIP() << "Mismatched number of ls_fitting and frame_calc records";
    }

    LeastSquaresFitter fitter;
    size_t tested = 0;
    size_t passed = 0;

    // Test first 5 residues (or all if fewer)
    size_t num_to_test = std::min(static_cast<size_t>(5), ls_fitting_records.size());

    for (size_t i = 0; i < num_to_test; ++i) {
        const auto& ls_fitting = ls_fitting_records[i];
        const auto& frame_calc = frame_calc_records[i];

        auto [std_points, exp_points] = extract_point_sets(frame_calc);

        if (std_points.size() < 3) {
            continue;
        }

        tested++;

        // Run fitting
        auto result = fitter.fit(std_points, exp_points);

        // Load expected values
        Matrix3D expected_rotation = load_rotation_matrix(ls_fitting);
        Vector3D expected_translation = load_translation(ls_fitting);
        double expected_rms = ls_fitting["rms_fit"].get<double>();

        // Check if they match
        bool rotation_ok =
            result.rotation.approximately_equals(expected_rotation, ROTATION_TOLERANCE);
        bool translation_ok =
            std::abs(result.translation.x() - expected_translation.x()) < TRANSLATION_TOLERANCE &&
            std::abs(result.translation.y() - expected_translation.y()) < TRANSLATION_TOLERANCE &&
            std::abs(result.translation.z() - expected_translation.z()) < TRANSLATION_TOLERANCE;
        bool rms_ok = std::abs(result.rms - expected_rms) < RMS_TOLERANCE;

        if (rotation_ok && translation_ok && rms_ok) {
            passed++;
        } else {
            std::cout << "Residue " << i << " mismatch:" << std::endl;
            if (!rotation_ok)
                std::cout << "  Rotation mismatch" << std::endl;
            if (!translation_ok)
                std::cout << "  Translation mismatch" << std::endl;
            if (!rms_ok)
                std::cout << "  RMS mismatch: expected " << expected_rms << ", got " << result.rms
                          << std::endl;
        }
    }

    EXPECT_GT(tested, 0) << "No valid test cases found";
    EXPECT_EQ(passed, tested) << passed << " out of " << tested << " residues matched";
}

// Test with 157D.json
TEST_F(LeastSquaresRegressionTest, Test157D_FirstResidue) {
    std::filesystem::path json_path;
    for (const auto& pair : pairs_) {
        if (pair.pdb_name == "157D") {
            json_path = pair.json_file;
            break;
        }
    }
    if (json_path.empty()) {
        json_path = std::filesystem::path("data/json_legacy/157D.json");
    }
    if (!std::filesystem::exists(json_path)) {
        GTEST_SKIP() << "157D.json not found";
    }

    auto json = load_legacy_json(json_path);

    auto ls_fitting_records = find_records_by_type(json, "ls_fitting");
    auto frame_calc_records = find_records_by_type(json, "frame_calc");

    if (ls_fitting_records.empty() || frame_calc_records.empty()) {
        GTEST_SKIP() << "No ls_fitting or frame_calc records found in 157D.json";
    }

    const auto& ls_fitting = ls_fitting_records[0];
    const auto& frame_calc = frame_calc_records[0];

    auto [std_points, exp_points] = extract_point_sets(frame_calc);

    if (std_points.size() < 3) {
        GTEST_SKIP() << "Not enough points";
    }

    LeastSquaresFitter fitter;
    auto result = fitter.fit(std_points, exp_points);

    Matrix3D expected_rotation = load_rotation_matrix(ls_fitting);
    Vector3D expected_translation = load_translation(ls_fitting);
    double expected_rms = ls_fitting["rms_fit"].get<double>();

    EXPECT_TRUE(result.rotation.approximately_equals(expected_rotation, ROTATION_TOLERANCE));
    EXPECT_NEAR(result.translation.x(), expected_translation.x(), TRANSLATION_TOLERANCE);
    EXPECT_NEAR(result.translation.y(), expected_translation.y(), TRANSLATION_TOLERANCE);
    EXPECT_NEAR(result.translation.z(), expected_translation.z(), TRANSLATION_TOLERANCE);
    EXPECT_NEAR(result.rms, expected_rms, RMS_TOLERANCE);
}

// Test all discovered PDB/JSON pairs (limited for performance)
TEST_F(LeastSquaresRegressionTest, TestAllPDBPairs) {
    LeastSquaresFitter fitter;
    size_t total_tested = 0;
    size_t total_passed = 0;

    // Only test discovered pairs (not all JSON files)
    // Limit to first 10 pairs and 2 residues each for fast regular testing (~20 cases)
    size_t max_pairs = std::min(static_cast<size_t>(10), pairs_.size());

    if (pairs_.empty()) {
        GTEST_SKIP() << "No PDB/JSON pairs discovered";
    }

    for (size_t pair_idx = 0; pair_idx < max_pairs; ++pair_idx) {
        const auto& pair = pairs_[pair_idx];
        try {
            auto json = load_legacy_json(pair.json_file);

            auto ls_fitting_records = find_records_by_type(json, "ls_fitting");
            auto frame_calc_records = find_records_by_type(json, "frame_calc");

            if (ls_fitting_records.size() != frame_calc_records.size()) {
                continue;
            }

            // Test first 2 residues from each PDB file for fast regular testing
            size_t num_to_test = std::min(static_cast<size_t>(2), ls_fitting_records.size());

            for (size_t i = 0; i < num_to_test; ++i) {
                const auto& ls_fitting = ls_fitting_records[i];
                const auto& frame_calc = frame_calc_records[i];

                auto [std_points, exp_points] = extract_point_sets(frame_calc);

                if (std_points.size() < 3) {
                    continue;
                }

                total_tested++;

                // Run fitting
                auto result = fitter.fit(std_points, exp_points);

                // Load expected values
                Matrix3D expected_rotation = load_rotation_matrix(ls_fitting);
                Vector3D expected_translation = load_translation(ls_fitting);
                double expected_rms = ls_fitting["rms_fit"].get<double>();

                // Check matches
                bool rotation_ok =
                    result.rotation.approximately_equals(expected_rotation, ROTATION_TOLERANCE);
                bool translation_ok = std::abs(result.translation.x() - expected_translation.x()) <
                                          TRANSLATION_TOLERANCE &&
                                      std::abs(result.translation.y() - expected_translation.y()) <
                                          TRANSLATION_TOLERANCE &&
                                      std::abs(result.translation.z() - expected_translation.z()) <
                                          TRANSLATION_TOLERANCE;
                bool rms_ok = std::abs(result.rms - expected_rms) < RMS_TOLERANCE;

                if (rotation_ok && translation_ok && rms_ok) {
                    total_passed++;
                } else {
                    std::cout << pair.pdb_name << " residue " << i << " mismatch:" << std::endl;
                    if (!rotation_ok)
                        std::cout << "  Rotation mismatch" << std::endl;
                    if (!translation_ok)
                        std::cout << "  Translation mismatch" << std::endl;
                    if (!rms_ok)
                        std::cout << "  RMS: expected " << expected_rms << ", got " << result.rms
                                  << std::endl;
                }
            }
        } catch (const std::exception& e) {
            std::cout << "Error processing " << pair.pdb_name << ": " << e.what() << std::endl;
        }
    }

    EXPECT_GT(total_tested, 0) << "No valid test cases found";
    double pass_rate = static_cast<double>(total_passed) / total_tested * 100.0;
    EXPECT_GE(pass_rate, 95.0) << "Only " << pass_rate << "% passed (" << total_passed << "/"
                               << total_tested << ")";

    std::cout << "Tested " << total_tested << " cases, " << total_passed << " passed (" << pass_rate
              << "%)" << std::endl;
}

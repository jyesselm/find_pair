/**
 * @file test_parameters.cpp
 * @brief Unit tests for parameter structures (BasePairStepParameters, HelicalParameters)
 */

#include <gtest/gtest.h>
#include <x3dna/core/parameters.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>
#include <nlohmann/json.hpp>

using namespace x3dna::core;
using namespace x3dna::geometry;

class ParametersTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test midstep frame
        Matrix3D rotation = Matrix3D::identity();
        Vector3D origin(1.0, 2.0, 3.0);
        test_frame_ = ReferenceFrame(rotation, origin);
    }

    ReferenceFrame test_frame_;
};

// ============================================================================
// BasePairStepParameters Tests
// ============================================================================

TEST_F(ParametersTest, BasePairStepParametersDefaultConstructor) {
    BasePairStepParameters params;
    EXPECT_EQ(params.shift, 0.0);
    EXPECT_EQ(params.slide, 0.0);
    EXPECT_EQ(params.rise, 0.0);
    EXPECT_EQ(params.tilt, 0.0);
    EXPECT_EQ(params.roll, 0.0);
    EXPECT_EQ(params.twist, 0.0);
    EXPECT_FALSE(params.midstep_frame.has_value());
}

TEST_F(ParametersTest, BasePairStepParametersConstructor) {
    BasePairStepParameters params(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    EXPECT_DOUBLE_EQ(params.shift, 1.0);
    EXPECT_DOUBLE_EQ(params.slide, 2.0);
    EXPECT_DOUBLE_EQ(params.rise, 3.0);
    EXPECT_DOUBLE_EQ(params.tilt, 4.0);
    EXPECT_DOUBLE_EQ(params.roll, 5.0);
    EXPECT_DOUBLE_EQ(params.twist, 6.0);
}

TEST_F(ParametersTest, BasePairStepParametersAsArray) {
    BasePairStepParameters params(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    auto arr = params.as_array();
    EXPECT_DOUBLE_EQ(arr[0], 1.0);
    EXPECT_DOUBLE_EQ(arr[1], 2.0);
    EXPECT_DOUBLE_EQ(arr[2], 3.0);
    EXPECT_DOUBLE_EQ(arr[3], 4.0);
    EXPECT_DOUBLE_EQ(arr[4], 5.0);
    EXPECT_DOUBLE_EQ(arr[5], 6.0);
}

TEST_F(ParametersTest, BasePairStepParametersFromArray) {
    std::array<double, 6> arr = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    auto params = BasePairStepParameters::from_array(arr);
    EXPECT_DOUBLE_EQ(params.shift, 1.0);
    EXPECT_DOUBLE_EQ(params.slide, 2.0);
    EXPECT_DOUBLE_EQ(params.rise, 3.0);
    EXPECT_DOUBLE_EQ(params.tilt, 4.0);
    EXPECT_DOUBLE_EQ(params.roll, 5.0);
    EXPECT_DOUBLE_EQ(params.twist, 6.0);
}

TEST_F(ParametersTest, BasePairStepParametersEquality) {
    BasePairStepParameters params1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    BasePairStepParameters params2(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    BasePairStepParameters params3(1.1, 2.0, 3.0, 4.0, 5.0, 6.0);

    EXPECT_TRUE(params1 == params2);
    EXPECT_FALSE(params1 == params3);
}

TEST_F(ParametersTest, BasePairStepParametersApproximatelyEqual) {
    BasePairStepParameters params1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    BasePairStepParameters params2(1.0000001, 2.0, 3.0, 4.0, 5.0, 6.0);
    BasePairStepParameters params3(1.01, 2.0, 3.0, 4.0, 5.0, 6.0);

    EXPECT_TRUE(params1.approximately_equal(params2, 1e-5));
    EXPECT_FALSE(params1.approximately_equal(params3, 1e-5));
}

TEST_F(ParametersTest, BasePairStepParametersToJson) {
    BasePairStepParameters params(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    params.midstep_frame = test_frame_;
    auto j = params.to_json();

    EXPECT_DOUBLE_EQ(j["shift"], 1.0);
    EXPECT_DOUBLE_EQ(j["slide"], 2.0);
    EXPECT_DOUBLE_EQ(j["rise"], 3.0);
    EXPECT_DOUBLE_EQ(j["tilt"], 4.0);
    EXPECT_DOUBLE_EQ(j["roll"], 5.0);
    EXPECT_DOUBLE_EQ(j["twist"], 6.0);
    EXPECT_TRUE(j.contains("midstep_frame"));
}

TEST_F(ParametersTest, BasePairStepParametersFromJson) {
    nlohmann::json j;
    j["shift"] = 1.0;
    j["slide"] = 2.0;
    j["rise"] = 3.0;
    j["tilt"] = 4.0;
    j["roll"] = 5.0;
    j["twist"] = 6.0;

    auto params = BasePairStepParameters::from_json(j);
    EXPECT_DOUBLE_EQ(params.shift, 1.0);
    EXPECT_DOUBLE_EQ(params.slide, 2.0);
    EXPECT_DOUBLE_EQ(params.rise, 3.0);
    EXPECT_DOUBLE_EQ(params.tilt, 4.0);
    EXPECT_DOUBLE_EQ(params.roll, 5.0);
    EXPECT_DOUBLE_EQ(params.twist, 6.0);
}

TEST_F(ParametersTest, BasePairStepParametersToJsonLegacy) {
    BasePairStepParameters params(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    params.midstep_frame = test_frame_;
    auto j = params.to_json_legacy(5, 6);

    EXPECT_EQ(j["type"], "bpstep_params");
    EXPECT_EQ(j["bp_idx1"], 5);
    EXPECT_EQ(j["bp_idx2"], 6);
    EXPECT_DOUBLE_EQ(j["params"]["Shift"], 1.0);
    EXPECT_DOUBLE_EQ(j["params"]["Slide"], 2.0);
    EXPECT_DOUBLE_EQ(j["params"]["Rise"], 3.0);
    EXPECT_DOUBLE_EQ(j["params"]["Tilt"], 4.0);
    EXPECT_DOUBLE_EQ(j["params"]["Roll"], 5.0);
    EXPECT_DOUBLE_EQ(j["params"]["Twist"], 6.0);
    EXPECT_TRUE(j.contains("mst_org"));
    EXPECT_TRUE(j.contains("mst_orien"));
}

TEST_F(ParametersTest, BasePairStepParametersFromJsonLegacy) {
    nlohmann::json j;
    j["type"] = "bpstep_params";
    j["bp_idx1"] = 5;
    j["bp_idx2"] = 6;
    j["params"]["Shift"] = 1.0;
    j["params"]["Slide"] = 2.0;
    j["params"]["Rise"] = 3.0;
    j["params"]["Tilt"] = 4.0;
    j["params"]["Roll"] = 5.0;
    j["params"]["Twist"] = 6.0;
    j["mst_org"] = {1.0, 2.0, 3.0};
    j["mst_orien"] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    auto params = BasePairStepParameters::from_json_legacy(j);
    EXPECT_DOUBLE_EQ(params.shift, 1.0);
    EXPECT_DOUBLE_EQ(params.slide, 2.0);
    EXPECT_DOUBLE_EQ(params.rise, 3.0);
    EXPECT_DOUBLE_EQ(params.tilt, 4.0);
    EXPECT_DOUBLE_EQ(params.roll, 5.0);
    EXPECT_DOUBLE_EQ(params.twist, 6.0);
    EXPECT_TRUE(params.midstep_frame.has_value());
}

TEST_F(ParametersTest, BasePairStepParametersJsonRoundTrip) {
    BasePairStepParameters original(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    original.midstep_frame = test_frame_;

    auto j = original.to_json_legacy(5, 6);
    auto restored = BasePairStepParameters::from_json_legacy(j);

    EXPECT_TRUE(original.approximately_equal(restored, 1e-9));
    EXPECT_TRUE(restored.midstep_frame.has_value());
}

// ============================================================================
// HelicalParameters Tests
// ============================================================================

TEST_F(ParametersTest, HelicalParametersDefaultConstructor) {
    HelicalParameters params;
    EXPECT_EQ(params.x_displacement, 0.0);
    EXPECT_EQ(params.y_displacement, 0.0);
    EXPECT_EQ(params.rise, 0.0);
    EXPECT_EQ(params.inclination, 0.0);
    EXPECT_EQ(params.tip, 0.0);
    EXPECT_EQ(params.twist, 0.0);
    EXPECT_FALSE(params.midstep_frame.has_value());
}

TEST_F(ParametersTest, HelicalParametersConstructor) {
    HelicalParameters params(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    EXPECT_DOUBLE_EQ(params.x_displacement, 1.0);
    EXPECT_DOUBLE_EQ(params.y_displacement, 2.0);
    EXPECT_DOUBLE_EQ(params.rise, 3.0);
    EXPECT_DOUBLE_EQ(params.inclination, 4.0);
    EXPECT_DOUBLE_EQ(params.tip, 5.0);
    EXPECT_DOUBLE_EQ(params.twist, 6.0);
}

TEST_F(ParametersTest, HelicalParametersAsArray) {
    HelicalParameters params(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    auto arr = params.as_array();
    EXPECT_DOUBLE_EQ(arr[0], 1.0);
    EXPECT_DOUBLE_EQ(arr[1], 2.0);
    EXPECT_DOUBLE_EQ(arr[2], 3.0);
    EXPECT_DOUBLE_EQ(arr[3], 4.0);
    EXPECT_DOUBLE_EQ(arr[4], 5.0);
    EXPECT_DOUBLE_EQ(arr[5], 6.0);
}

TEST_F(ParametersTest, HelicalParametersFromArray) {
    std::array<double, 6> arr = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    auto params = HelicalParameters::from_array(arr);
    EXPECT_DOUBLE_EQ(params.x_displacement, 1.0);
    EXPECT_DOUBLE_EQ(params.y_displacement, 2.0);
    EXPECT_DOUBLE_EQ(params.rise, 3.0);
    EXPECT_DOUBLE_EQ(params.inclination, 4.0);
    EXPECT_DOUBLE_EQ(params.tip, 5.0);
    EXPECT_DOUBLE_EQ(params.twist, 6.0);
}

TEST_F(ParametersTest, HelicalParametersEquality) {
    HelicalParameters params1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    HelicalParameters params2(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    HelicalParameters params3(1.1, 2.0, 3.0, 4.0, 5.0, 6.0);

    EXPECT_TRUE(params1 == params2);
    EXPECT_FALSE(params1 == params3);
}

TEST_F(ParametersTest, HelicalParametersApproximatelyEqual) {
    HelicalParameters params1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    HelicalParameters params2(1.0000001, 2.0, 3.0, 4.0, 5.0, 6.0);
    HelicalParameters params3(1.01, 2.0, 3.0, 4.0, 5.0, 6.0);

    EXPECT_TRUE(params1.approximately_equal(params2, 1e-5));
    EXPECT_FALSE(params1.approximately_equal(params3, 1e-5));
}

TEST_F(ParametersTest, HelicalParametersToJson) {
    HelicalParameters params(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    params.midstep_frame = test_frame_;
    auto j = params.to_json();

    EXPECT_DOUBLE_EQ(j["x_displacement"], 1.0);
    EXPECT_DOUBLE_EQ(j["y_displacement"], 2.0);
    EXPECT_DOUBLE_EQ(j["rise"], 3.0);
    EXPECT_DOUBLE_EQ(j["inclination"], 4.0);
    EXPECT_DOUBLE_EQ(j["tip"], 5.0);
    EXPECT_DOUBLE_EQ(j["twist"], 6.0);
    EXPECT_TRUE(j.contains("midstep_frame"));
}

TEST_F(ParametersTest, HelicalParametersFromJson) {
    nlohmann::json j;
    j["x_displacement"] = 1.0;
    j["y_displacement"] = 2.0;
    j["rise"] = 3.0;
    j["inclination"] = 4.0;
    j["tip"] = 5.0;
    j["twist"] = 6.0;

    auto params = HelicalParameters::from_json(j);
    EXPECT_DOUBLE_EQ(params.x_displacement, 1.0);
    EXPECT_DOUBLE_EQ(params.y_displacement, 2.0);
    EXPECT_DOUBLE_EQ(params.rise, 3.0);
    EXPECT_DOUBLE_EQ(params.inclination, 4.0);
    EXPECT_DOUBLE_EQ(params.tip, 5.0);
    EXPECT_DOUBLE_EQ(params.twist, 6.0);
}

TEST_F(ParametersTest, HelicalParametersToJsonLegacy) {
    HelicalParameters params(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    params.midstep_frame = test_frame_;
    auto j = params.to_json_legacy(5, 6);

    EXPECT_EQ(j["type"], "helical_params");
    EXPECT_EQ(j["bp_idx1"], 5);
    EXPECT_EQ(j["bp_idx2"], 6);
    EXPECT_TRUE(j["params"].is_array());
    EXPECT_EQ(j["params"].size(), 6);
    EXPECT_DOUBLE_EQ(j["params"][0], 1.0);
    EXPECT_DOUBLE_EQ(j["params"][1], 2.0);
    EXPECT_DOUBLE_EQ(j["params"][2], 3.0);
    EXPECT_DOUBLE_EQ(j["params"][3], 4.0);
    EXPECT_DOUBLE_EQ(j["params"][4], 5.0);
    EXPECT_DOUBLE_EQ(j["params"][5], 6.0);
    EXPECT_TRUE(j.contains("mst_orgH"));
    EXPECT_TRUE(j.contains("mst_orienH"));
}

TEST_F(ParametersTest, HelicalParametersFromJsonLegacy) {
    nlohmann::json j;
    j["type"] = "helical_params";
    j["bp_idx1"] = 5;
    j["bp_idx2"] = 6;
    j["params"] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    j["mst_orgH"] = {1.0, 2.0, 3.0};
    j["mst_orienH"] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    auto params = HelicalParameters::from_json_legacy(j);
    EXPECT_DOUBLE_EQ(params.x_displacement, 1.0);
    EXPECT_DOUBLE_EQ(params.y_displacement, 2.0);
    EXPECT_DOUBLE_EQ(params.rise, 3.0);
    EXPECT_DOUBLE_EQ(params.inclination, 4.0);
    EXPECT_DOUBLE_EQ(params.tip, 5.0);
    EXPECT_DOUBLE_EQ(params.twist, 6.0);
    EXPECT_TRUE(params.midstep_frame.has_value());
}

TEST_F(ParametersTest, HelicalParametersJsonRoundTrip) {
    HelicalParameters original(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    original.midstep_frame = test_frame_;

    auto j = original.to_json_legacy(5, 6);
    auto restored = HelicalParameters::from_json_legacy(j);

    EXPECT_TRUE(original.approximately_equal(restored, 1e-9));
    EXPECT_TRUE(restored.midstep_frame.has_value());
}

// ============================================================================
// Real-world data tests (from legacy JSON)
// ============================================================================

TEST_F(ParametersTest, BasePairStepParametersRealWorldData) {
    // From 2Y8Y.json example
    nlohmann::json j;
    j["type"] = "bpstep_params";
    j["bp_idx1"] = 5;
    j["bp_idx2"] = 6;
    j["params"]["Shift"] = 0.398253;
    j["params"]["Slide"] = -1.454919;
    j["params"]["Rise"] = 3.129627;
    j["params"]["Tilt"] = -6.340293;
    j["params"]["Roll"] = 2.748534;
    j["params"]["Twist"] = 28.086432;
    j["mst_org"] = {14.814408, 0.063380, -9.362866};
    j["mst_orien"] = {{-0.594480, 0.480071, 0.645078},
                      {0.580022, 0.811631, -0.069494},
                      {-0.556927, 0.332847, -0.760950}};

    auto params = BasePairStepParameters::from_json_legacy(j);
    EXPECT_NEAR(params.shift, 0.398253, 1e-6);
    EXPECT_NEAR(params.slide, -1.454919, 1e-6);
    EXPECT_NEAR(params.rise, 3.129627, 1e-6);
    EXPECT_NEAR(params.tilt, -6.340293, 1e-6);
    EXPECT_NEAR(params.roll, 2.748534, 1e-6);
    EXPECT_NEAR(params.twist, 28.086432, 1e-6);
    EXPECT_TRUE(params.midstep_frame.has_value());
}

TEST_F(ParametersTest, HelicalParametersRealWorldData) {
    // From 2Y8Y.json example
    nlohmann::json j;
    j["type"] = "helical_params";
    j["bp_idx1"] = 5;
    j["bp_idx2"] = 6;
    j["params"] = {-3.477697, -2.096910, 2.821410, 5.557507, 12.820006, 28.907457};
    j["mst_orgH"] = {14.178829, 3.614285, -10.928471};
    j["mst_orienH"] = {{-0.437044, 0.411822, 0.799622},
                       {0.559742, 0.820422, -0.116600},
                       {-0.704046, 0.396623, -0.589074}};

    auto params = HelicalParameters::from_json_legacy(j);
    EXPECT_NEAR(params.x_displacement, -3.477697, 1e-6);
    EXPECT_NEAR(params.y_displacement, -2.096910, 1e-6);
    EXPECT_NEAR(params.rise, 2.821410, 1e-6);
    EXPECT_NEAR(params.inclination, 5.557507, 1e-6);
    EXPECT_NEAR(params.tip, 12.820006, 1e-6);
    EXPECT_NEAR(params.twist, 28.907457, 1e-6);
    EXPECT_TRUE(params.midstep_frame.has_value());
}


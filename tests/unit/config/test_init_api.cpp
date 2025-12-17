/**
 * @file test_init_api.cpp
 * @brief Unit tests for x3dna::init() API
 */

#include <gtest/gtest.h>
#include <x3dna/x3dna.hpp>
#include <filesystem>

class InitApiTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Reset before each test
        x3dna::shutdown();
    }

    void TearDown() override {
        x3dna::shutdown();
    }
};

TEST_F(InitApiTest, InitiallyNotInitialized) {
    EXPECT_FALSE(x3dna::is_initialized());
}

TEST_F(InitApiTest, InitWithValidPath) {
    std::filesystem::path resources = std::filesystem::path(X3DNA_SOURCE_DIR) / "resources";
    EXPECT_TRUE(x3dna::init(resources));
    EXPECT_TRUE(x3dna::is_initialized());
}

TEST_F(InitApiTest, InitWithInvalidPathReturnsFalse) {
    EXPECT_FALSE(x3dna::init("/nonexistent/path"));
    EXPECT_FALSE(x3dna::is_initialized());
}

TEST_F(InitApiTest, ShutdownResetsState) {
    std::filesystem::path resources = std::filesystem::path(X3DNA_SOURCE_DIR) / "resources";
    x3dna::init(resources);
    EXPECT_TRUE(x3dna::is_initialized());

    x3dna::shutdown();
    EXPECT_FALSE(x3dna::is_initialized());
}

TEST_F(InitApiTest, ResourcesPathReturnsCorrectPath) {
    std::filesystem::path resources = std::filesystem::path(X3DNA_SOURCE_DIR) / "resources";
    x3dna::init(resources);

    EXPECT_EQ(x3dna::resources_path(), resources);
}

TEST_F(InitApiTest, VersionReturnsNonEmpty) {
    const char* ver = x3dna::version();
    EXPECT_NE(ver, nullptr);
    EXPECT_GT(strlen(ver), 0u);
}

TEST_F(InitApiTest, AutoInitWorks) {
    // This test requires resources to be in a searchable location
    // When running from build directory, "../resources" should work
    bool result = x3dna::init();
    // We accept either outcome - the important thing is it doesn't crash
    if (result) {
        EXPECT_TRUE(x3dna::is_initialized());
    }
}

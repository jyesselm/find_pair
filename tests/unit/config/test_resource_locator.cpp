/**
 * @file test_resource_locator.cpp
 * @brief Unit tests for ResourceLocator
 */

#include <gtest/gtest.h>
#include <x3dna/config/resource_locator.hpp>
#include <filesystem>

using namespace x3dna::config;

class ResourceLocatorTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Reset ResourceLocator before each test
        ResourceLocator::reset();
    }

    void TearDown() override {
        // Clean up after each test
        ResourceLocator::reset();
    }
};

TEST_F(ResourceLocatorTest, InitiallyNotInitialized) {
    EXPECT_FALSE(ResourceLocator::is_initialized());
}

TEST_F(ResourceLocatorTest, InitializeWithValidPath) {
    // Find the resources directory relative to test execution
    std::filesystem::path resources_path = "resources";
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../resources";
    }
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../../resources";
    }

    if (std::filesystem::exists(resources_path)) {
        ASSERT_NO_THROW(ResourceLocator::initialize(resources_path));
        EXPECT_TRUE(ResourceLocator::is_initialized());
    } else {
        GTEST_SKIP() << "Resources directory not found - skipping test";
    }
}

TEST_F(ResourceLocatorTest, InitializeWithInvalidPathThrows) {
    EXPECT_THROW(ResourceLocator::initialize("/nonexistent/path/to/resources"), std::runtime_error);
    EXPECT_FALSE(ResourceLocator::is_initialized());
}

TEST_F(ResourceLocatorTest, TemplatesDirReturnsCorrectPath) {
    std::filesystem::path resources_path = "resources";
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../resources";
    }
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../../resources";
    }

    if (std::filesystem::exists(resources_path)) {
        ResourceLocator::initialize(resources_path);
        auto templates_dir = ResourceLocator::templates_dir();
        EXPECT_TRUE(templates_dir.string().find("templates") != std::string::npos);
        EXPECT_TRUE(std::filesystem::exists(templates_dir));
    } else {
        GTEST_SKIP() << "Resources directory not found - skipping test";
    }
}

TEST_F(ResourceLocatorTest, ConfigDirReturnsCorrectPath) {
    std::filesystem::path resources_path = "resources";
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../resources";
    }
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../../resources";
    }

    if (std::filesystem::exists(resources_path)) {
        ResourceLocator::initialize(resources_path);
        auto config_dir = ResourceLocator::config_dir();
        EXPECT_TRUE(config_dir.string().find("config") != std::string::npos);
        EXPECT_TRUE(std::filesystem::exists(config_dir));
    } else {
        GTEST_SKIP() << "Resources directory not found - skipping test";
    }
}

TEST_F(ResourceLocatorTest, TemplateFileReturnsCorrectPath) {
    std::filesystem::path resources_path = "resources";
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../resources";
    }
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../../resources";
    }

    if (std::filesystem::exists(resources_path)) {
        ResourceLocator::initialize(resources_path);
        auto template_file = ResourceLocator::template_file("Atomic_A.pdb");
        EXPECT_TRUE(template_file.string().find("Atomic_A.pdb") != std::string::npos);
    } else {
        GTEST_SKIP() << "Resources directory not found - skipping test";
    }
}

TEST_F(ResourceLocatorTest, ConfigFileReturnsCorrectPath) {
    std::filesystem::path resources_path = "resources";
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../resources";
    }
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../../resources";
    }

    if (std::filesystem::exists(resources_path)) {
        ResourceLocator::initialize(resources_path);
        auto config_file = ResourceLocator::config_file("atomlist.dat");
        EXPECT_TRUE(config_file.string().find("atomlist.dat") != std::string::npos);
    } else {
        GTEST_SKIP() << "Resources directory not found - skipping test";
    }
}

TEST_F(ResourceLocatorTest, TemplateExistsReturnsTrueForExistingFile) {
    std::filesystem::path resources_path = "resources";
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../resources";
    }
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../../resources";
    }

    if (std::filesystem::exists(resources_path)) {
        ResourceLocator::initialize(resources_path);
        EXPECT_TRUE(ResourceLocator::template_exists("Atomic_A.pdb"));
        EXPECT_FALSE(ResourceLocator::template_exists("nonexistent.pdb"));
    } else {
        GTEST_SKIP() << "Resources directory not found - skipping test";
    }
}

TEST_F(ResourceLocatorTest, ConfigExistsReturnsTrueForExistingFile) {
    std::filesystem::path resources_path = "resources";
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../resources";
    }
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../../resources";
    }

    if (std::filesystem::exists(resources_path)) {
        ResourceLocator::initialize(resources_path);
        EXPECT_TRUE(ResourceLocator::config_exists("atomlist.dat"));
        EXPECT_FALSE(ResourceLocator::config_exists("nonexistent.dat"));
    } else {
        GTEST_SKIP() << "Resources directory not found - skipping test";
    }
}

TEST_F(ResourceLocatorTest, ResetClearsInitialization) {
    std::filesystem::path resources_path = "resources";
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../resources";
    }
    if (!std::filesystem::exists(resources_path)) {
        resources_path = "../../resources";
    }

    if (std::filesystem::exists(resources_path)) {
        ResourceLocator::initialize(resources_path);
        EXPECT_TRUE(ResourceLocator::is_initialized());

        ResourceLocator::reset();
        EXPECT_FALSE(ResourceLocator::is_initialized());
    } else {
        GTEST_SKIP() << "Resources directory not found - skipping test";
    }
}

TEST_F(ResourceLocatorTest, AutoInitializeFromEnvironmentWorks) {
    // This test verifies that initialize_from_environment works
    // It may or may not succeed depending on the test environment
    bool result = ResourceLocator::initialize_from_environment();

    if (result) {
        EXPECT_TRUE(ResourceLocator::is_initialized());
    }
    // If it fails, that's also acceptable - just means resources weren't found
}

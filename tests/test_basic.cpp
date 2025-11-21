/**
 * @file test_basic.cpp
 * @brief Basic tests to verify test framework works
 */

#include <gtest/gtest.h>
#include <x3dna/x3dna.hpp>
#include "test_helpers.hpp"
#include "integration/test_data_discovery.hpp"

using namespace x3dna;

// Test that the library can be included
TEST(BasicTest, LibraryIncludes) {
    EXPECT_NE(version(), nullptr);
    EXPECT_STREQ(version(), "1.0.0");
}

// Test helper functions
TEST(BasicTest, TestHelpers) {
    using namespace x3dna::test;
    
    EXPECT_TRUE(approximately_equal(1.0, 1.0001, 0.001));
    EXPECT_FALSE(approximately_equal(1.0, 1.01, 0.001));
}

// Test that integration test discovery works
TEST(BasicTest, TestDataDiscovery) {
    using namespace x3dna::test;
    
    auto pairs = test_data_discovery::discover_pairs();
    
    // Should find at least 100D and 157D if JSON files exist
    if (!pairs.empty()) {
        std::cout << "Found " << pairs.size() << " PDB/JSON pairs for testing" << std::endl;
        for (const auto& pair : pairs) {
            std::cout << "  - " << pair.pdb_name << std::endl;
        }
    } else {
        std::cout << "No PDB/JSON pairs found (this is OK if JSON files haven't been generated yet)" << std::endl;
    }
}


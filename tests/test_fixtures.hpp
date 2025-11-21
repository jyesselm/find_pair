/**
 * @file test_fixtures.hpp
 * @brief Test fixtures and sample data
 */

#pragma once

#include <gtest/gtest.h>

namespace x3dna::test {

/**
 * @brief Base test fixture
 */
class test_fixture : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

} // namespace x3dna::test

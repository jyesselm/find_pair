/**
 * @file test_helpers.hpp
 * @brief Helper utilities for testing
 */

#pragma once

#include <gtest/gtest.h>
#include <cmath>
#include <string>

namespace x3dna::test {

/**
 * @brief Check if two floating point values are approximately equal
 * @param a First value
 * @param b Second value
 * @param tolerance Tolerance for comparison
 * @return True if values are within tolerance
 */
inline bool approximately_equal(double a, double b, double tolerance = 0.001) {
    return std::abs(a - b) < tolerance;
}

/**
 * @brief Assert that two values are approximately equal
 */
#define EXPECT_APPROX_EQ(a, b, tolerance)                                                          \
    EXPECT_NEAR((a), (b), (tolerance)) << "Values not approximately equal: " << (a) << " vs " << (b)

} // namespace x3dna::test

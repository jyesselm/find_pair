# Integration Tests

# Least squares regression test (uses real data from JSON files)
add_executable(test_least_squares_regression
    test_least_squares_regression.cpp
)

target_include_directories(test_least_squares_regression PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Copy test data to build directory so tests can find it
file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_BINARY_DIR})

target_link_libraries(test_least_squares_regression
    x3dna
    gtest
    gtest_main
)

gtest_discover_tests(test_least_squares_regression)

# Note: More integration tests will be added as we implement each stage


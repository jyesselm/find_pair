# Testing Configuration

# Google Test
include(GoogleTest)

# Note: Test helper utilities are header-only and included directly
# No library targets needed for header-only files

# Basic test to verify framework works
# Only add if Google Test is available
if(TARGET gtest AND TARGET gtest_main)
    add_executable(test_basic
        test_main.cpp
        test_basic.cpp
    )
    
    target_link_libraries(test_basic
        x3dna
        gtest_main
    )
    
    gtest_discover_tests(test_basic)
else()
    message(WARNING "Google Test not available, skipping test_basic")
endif()

# Unit tests
add_subdirectory(unit)

# Integration tests
add_subdirectory(integration)

# Regression tests
add_subdirectory(regression)

# Integration Tests

# Least squares regression test (uses real data from JSON files)
add_executable(test_least_squares_regression
    test_least_squares_regression.cpp
)

target_include_directories(test_least_squares_regression PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Create symlink to test data (much faster than copying 383MB every time)
# Tests run from build directory, so we need data accessible there
# Only create symlink if it doesn't already exist
if(NOT EXISTS ${CMAKE_BINARY_DIR}/data)
    if(UNIX OR APPLE)
        # Unix/macOS: use symlink (instant, no copy needed)
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${CMAKE_SOURCE_DIR}/data
            ${CMAKE_BINARY_DIR}/data
            RESULT_VARIABLE symlink_result
            ERROR_QUIET
        )
        # If symlink failed (e.g., on Windows with Unix tools), fall back to copy
        if(symlink_result AND NOT symlink_result EQUAL 0)
            file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_BINARY_DIR})
        endif()
    else()
        # Windows: copy (only once)
        file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_BINARY_DIR})
    endif()
endif()

target_link_libraries(test_least_squares_regression
    PRIVATE
    x3dna
    gtest
    gtest_main
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_least_squares_regression PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_least_squares_regression PRIVATE c++fs)
endif()

gtest_discover_tests(test_least_squares_regression)

add_executable(test_pdb_atom_validation
    test_pdb_atom_validation.cpp
)

target_include_directories(test_pdb_atom_validation PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_pdb_atom_validation
    PRIVATE
    x3dna
    gtest
    gtest_main
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_pdb_atom_validation PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_pdb_atom_validation PRIVATE c++fs)
endif()

gtest_discover_tests(test_pdb_atom_validation)

# JSON generation test
add_executable(test_json_generation
    test_json_generation.cpp
)

target_include_directories(test_json_generation PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_json_generation
    PRIVATE
    x3dna
    gtest
    gtest_main
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_json_generation PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_json_generation PRIVATE c++fs)
endif()

gtest_discover_tests(test_json_generation)

# Single PDB test tool (not a gtest, just a standalone tool)
add_executable(test_single_pdb
    test_single_pdb.cpp
)

target_include_directories(test_single_pdb PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_single_pdb
    PRIVATE
    x3dna
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_single_pdb PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_single_pdb PRIVATE c++fs)
endif()

# Note: More integration tests will be added as we implement each stage


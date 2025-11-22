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

# Filtered JSON generation test (only problematic PDBs)
add_executable(test_json_generation_filtered
    test_json_generation_filtered.cpp
)

target_include_directories(test_json_generation_filtered PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_json_generation_filtered
    PRIVATE
    x3dna
    gtest
    gtest_main
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_json_generation_filtered PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_json_generation_filtered PRIVATE c++fs)
endif()

gtest_discover_tests(test_json_generation_filtered)

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

# Frame calculation integration test
add_executable(test_frame_calculation
    test_frame_calculation.cpp
)

target_include_directories(test_frame_calculation PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_frame_calculation
    PRIVATE
    x3dna
    gtest
    gtest_main
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_frame_calculation PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_frame_calculation PRIVATE c++fs)
endif()

gtest_discover_tests(test_frame_calculation)

# Frame calculation legacy comparison test
add_executable(test_frame_calculation_legacy
    test_frame_calculation_legacy.cpp
)

target_include_directories(test_frame_calculation_legacy PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_frame_calculation_legacy
    PRIVATE
    x3dna
    gtest
    gtest_main
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_frame_calculation_legacy PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_frame_calculation_legacy PRIVATE c++fs)
endif()

gtest_discover_tests(test_frame_calculation_legacy)

# Mismatch analysis tool
add_executable(analyze_frame_mismatches
    analyze_frame_mismatches.cpp
)

target_include_directories(analyze_frame_mismatches PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(analyze_frame_mismatches
    PRIVATE
    x3dna
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(analyze_frame_mismatches PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(analyze_frame_mismatches PRIVATE c++fs)
endif()

# Check failing residues
add_executable(check_failing_residues
    check_failing_residues.cpp
)

target_include_directories(check_failing_residues PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(check_failing_residues
    PRIVATE
    x3dna
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(check_failing_residues PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(check_failing_residues PRIVATE c++fs)
endif()

# Debug all failures in detail
add_executable(debug_all_failures
    debug_all_failures.cpp
)

target_include_directories(debug_all_failures PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(debug_all_failures
    PRIVATE
    x3dna
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(debug_all_failures PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(debug_all_failures PRIVATE c++fs)
endif()

# Investigate test failures
add_executable(investigate_failures
    investigate_failures.cpp
)

target_include_directories(investigate_failures PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(investigate_failures
    PRIVATE
    x3dna
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(investigate_failures PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(investigate_failures PRIVATE c++fs)
endif()

# Compare LS fitting directly
add_executable(compare_ls_fitting
    compare_ls_fitting.cpp
)

target_include_directories(compare_ls_fitting PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(compare_ls_fitting
    PRIVATE
    x3dna
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(compare_ls_fitting PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(compare_ls_fitting PRIVATE c++fs)
endif()

# Debug atom matching
add_executable(debug_atom_matching
    debug_atom_matching.cpp
)

target_include_directories(debug_atom_matching PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(debug_atom_matching
    PRIVATE
    x3dna
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(debug_atom_matching PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(debug_atom_matching PRIVATE c++fs)
endif()

# Debug tool for frame calculation
add_executable(debug_frame_calculation
    debug_frame_calculation.cpp
)

target_include_directories(debug_frame_calculation PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(debug_frame_calculation
    PRIVATE
    x3dna
)

# Link filesystem library if needed (CMake 3.15+)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(debug_frame_calculation PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(debug_frame_calculation PRIVATE c++fs)
endif()

# Note: More integration tests will be added as we implement each stage


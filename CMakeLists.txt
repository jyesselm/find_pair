cmake_minimum_required(VERSION 3.15)
project(x3dna VERSION 1.0.0 LANGUAGES CXX)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Generate compile_commands.json for language servers (clangd, etc.)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build options
option(BUILD_TESTS "Build tests" ON)
option(BUILD_EXAMPLES "Build examples" ON)
option(BUILD_DOCS "Build documentation" OFF)

# Compiler flags
if(MSVC)
    add_compile_options(/W4 /WX-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)
endif()

# Include directories
include_directories(include)

# Dependencies (must be included before tests)
include(cmake/Dependencies.cmake)

# Library (header-only for now, will add source files as we progress)
# Changed to regular library to support source files (PdbParser, etc.)
add_library(x3dna
    src/x3dna/io/pdb_parser.cpp
    src/x3dna/io/json_writer.cpp
    src/x3dna/io/json_reader.cpp
    src/x3dna/io/pdb_writer.cpp
    src/x3dna/io/input_file_parser.cpp
    src/x3dna/algorithms/standard_base_templates.cpp
    src/x3dna/algorithms/ring_atom_matcher.cpp
    src/x3dna/algorithms/base_frame_calculator.cpp
    src/x3dna/algorithms/base_pair_validator.cpp
    src/x3dna/algorithms/base_pair_finder.cpp
    src/x3dna/algorithms/hydrogen_bond_finder.cpp
    src/x3dna/algorithms/hydrogen_bond/hydrogen_bond_counter.cpp
    src/x3dna/algorithms/hydrogen_bond/hydrogen_bond_utils.cpp
)
target_include_directories(x3dna PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Link dependencies to library
target_link_libraries(x3dna PUBLIC
    nlohmann_json::nlohmann_json
)

# Testing (after dependencies are loaded)
if(BUILD_TESTS)
    enable_testing()
    # Google Test must be available before adding tests
    add_subdirectory(tests)
endif()

# Examples
if(BUILD_EXAMPLES)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/examples/CMakeLists.txt)
        add_subdirectory(examples)
    endif()
endif()

# Tools
add_executable(generate_modern_json tools/generate_modern_json.cpp)
target_link_libraries(generate_modern_json PRIVATE x3dna)

add_executable(test_overlap_calculation tools/test_overlap_calculation.cpp)
target_link_libraries(test_overlap_calculation PRIVATE x3dna)

add_executable(compare_pdb_parsing tools/compare_pdb_parsing.cpp)
target_link_libraries(compare_pdb_parsing PRIVATE x3dna)

add_executable(compare_residue_identification tools/compare_residue_identification.cpp)
target_link_libraries(compare_residue_identification PRIVATE x3dna)

add_executable(debug_frame_calculation tools/debug_frame_calculation.cpp)
target_link_libraries(debug_frame_calculation PRIVATE x3dna)

# Documentation
if(BUILD_DOCS)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        add_subdirectory(docs)
    endif()
endif()

# Installation (optional, for packaging)
# install(TARGETS x3dna EXPORT x3dnaTargets)
# install(DIRECTORY include/ DESTINATION include)
# install(EXPORT x3dnaTargets
#     FILE x3dnaTargets.cmake
#     NAMESPACE x3dna::
#     DESTINATION lib/cmake/x3dna
# )


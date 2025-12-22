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
option(DEBUG_BP_TYPE_ID "Enable debug output for bp_type_id calculation" OFF)

# Compiler flags
if(MSVC)
    add_compile_options(/W4 /WX-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)
endif()

# Include directories
include_directories(include)

# Generate parameters header from JSON (compile-time constants)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(PARAMS_JSON ${CMAKE_SOURCE_DIR}/resources/config/parameters.json)
set(PARAMS_HPP ${CMAKE_SOURCE_DIR}/include/x3dna/config/parameters_generated.hpp)
set(PARAMS_SCRIPT ${CMAKE_SOURCE_DIR}/cmake/generate_parameters.py)

add_custom_command(
    OUTPUT ${PARAMS_HPP}
    COMMAND ${Python3_EXECUTABLE} ${PARAMS_SCRIPT} ${PARAMS_JSON} ${PARAMS_HPP}
    DEPENDS ${PARAMS_JSON} ${PARAMS_SCRIPT}
    COMMENT "Generating parameters_generated.hpp from parameters.json"
)
add_custom_target(generate_parameters DEPENDS ${PARAMS_HPP})

# Dependencies (must be included before tests)
include(cmake/Dependencies.cmake)

# Library (header-only for now, will add source files as we progress)
# Changed to regular library to support source files (PdbParser, etc.)
add_library(x3dna
    src/x3dna/io/pdb_parser.cpp
    src/x3dna/io/cif_parser.cpp
    src/x3dna/io/json_writer.cpp
    src/x3dna/io/json_reader.cpp
    src/x3dna/io/pdb_writer.cpp
    src/x3dna/io/input_file_parser.cpp
    src/x3dna/io/input_file_writer.cpp
    src/x3dna/io/frame_json_recorder.cpp
    src/x3dna/residue_tracker.cpp
    src/x3dna/core/structure_legacy_order.cpp
    src/x3dna/core/structure_legacy_order_impl.cpp
    src/x3dna/core/modified_nucleotide_registry.cpp
    src/x3dna/core/atom_symbol_registry.cpp
    src/x3dna/core/typing/atom_classification.cpp
    src/x3dna/core/typing/type_registry.cpp
    src/x3dna/core/typing/residue_classification.cpp
    src/x3dna/core/nucleotide_utils.cpp
    src/x3dna/algorithms/standard_base_templates.cpp
    src/x3dna/algorithms/ring_atom_matcher.cpp
    src/x3dna/algorithms/residue_type_detector.cpp
    src/x3dna/algorithms/base_frame_calculator.cpp
    src/x3dna/algorithms/base_pair_validator.cpp
    src/x3dna/algorithms/base_pair_finder.cpp
    src/x3dna/algorithms/parameter_calculator.cpp
    src/x3dna/algorithms/helix_organizer.cpp
    src/x3dna/algorithms/hydrogen_bond_finder.cpp
    src/x3dna/algorithms/residue_index_map.cpp
    src/x3dna/algorithms/quality_score_calculator.cpp
    src/x3dna/algorithms/pair_candidate_cache.cpp
    src/x3dna/algorithms/json_writer_observer.cpp
    src/x3dna/algorithms/pair_selection_strategy.cpp
    src/x3dna/algorithms/hydrogen_bond/hydrogen_bond_utils.cpp
    src/x3dna/algorithms/validation/overlap_calculator.cpp
    src/x3dna/config/config_manager.cpp
    src/x3dna/config/resource_locator.cpp
    src/x3dna/protocols/find_pair_protocol.cpp
    src/x3dna/protocols/analyze_protocol.cpp
    src/x3dna/apps/command_line_parser.cpp
    src/x3dna/debug/pair_validation_debugger.cpp
)
target_include_directories(x3dna PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Ensure parameters header is generated before compiling
add_dependencies(x3dna generate_parameters)

# Link dependencies to library
target_link_libraries(x3dna PUBLIC
    nlohmann_json::nlohmann_json
)
# GEMMI compiled library for PDB/CIF parsing
target_link_libraries(x3dna PRIVATE
    gemmi_cpp
)

# Note: ResourceLocator handles runtime path discovery - no compile-time paths needed

# Add compile definitions for debug flags (optional, enable with -DDEBUG_BP_TYPE_ID=ON etc)
if(DEBUG_BP_TYPE_ID)
    target_compile_definitions(x3dna PRIVATE DEBUG_BP_TYPE_ID)
endif()
if(DEBUG_FRAME_CALC)
    target_compile_definitions(x3dna PRIVATE DEBUG_FRAME_CALC)
endif()

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

# Tools - Only build tools that exist
# Commented out - tools archived or removed from repository
add_executable(generate_modern_json tools/generate_modern_json.cpp)
target_link_libraries(generate_modern_json PRIVATE x3dna)

add_executable(compare_single_pair tools/compare_single_pair.cpp)
target_link_libraries(compare_single_pair PRIVATE x3dna)
# 
# add_executable(check_residue_indices tools/check_residue_indices.cpp)
# target_link_libraries(check_residue_indices PRIVATE x3dna)
# add_executable(test_overlap_calculation tools/test_overlap_calculation.cpp)
# target_link_libraries(test_overlap_calculation PRIVATE x3dna)
# 
# add_executable(compare_pdb_parsing tools/compare_pdb_parsing.cpp)
# target_link_libraries(compare_pdb_parsing PRIVATE x3dna)
# 
# add_executable(compare_residue_identification tools/compare_residue_identification.cpp)
# target_link_libraries(compare_residue_identification PRIVATE x3dna)
# 
# add_executable(debug_frame_calculation tools/debug_frame_calculation.cpp)
# target_link_libraries(debug_frame_calculation PRIVATE x3dna)
# 
# add_executable(compare_quality_scores tools/compare_quality_scores.cpp)
# target_link_libraries(compare_quality_scores PRIVATE x3dna)
# 
# add_executable(compare_initial_hbonds tools/compare_initial_hbonds.cpp)
# target_link_libraries(compare_initial_hbonds PRIVATE x3dna)
# 
# add_executable(compare_hbond_stages tools/compare_hbond_stages.cpp)
# target_link_libraries(compare_hbond_stages PRIVATE x3dna)
# 
# add_executable(compare_atom_selection tools/compare_atom_selection.cpp)
# target_link_libraries(compare_atom_selection PRIVATE x3dna)
# 
# add_executable(find_residue_mapping tools/find_residue_mapping.cpp)
# target_link_libraries(find_residue_mapping PRIVATE x3dna)
# 
# add_executable(test_legacy_order tools/test_legacy_order.cpp)
# target_link_libraries(test_legacy_order PRIVATE x3dna)
# 
# add_executable(generate_residue_ordering_json tools/generate_residue_ordering_json.cpp)
# target_link_libraries(generate_residue_ordering_json PRIVATE x3dna)
# 
# add_executable(compare_residue_ordering tools/compare_residue_ordering.cpp)
# target_link_libraries(compare_residue_ordering PRIVATE x3dna)
# 
# add_executable(compare_hbond_detection tools/compare_hbond_detection.cpp)
# target_link_libraries(compare_hbond_detection PRIVATE x3dna)
# 
# add_executable(list_all_hbonds tools/list_all_hbonds.cpp)
# target_link_libraries(list_all_hbonds PRIVATE x3dna)
# 
# add_executable(detect_hbonds_standalone tools/detect_hbonds_standalone.cpp)
# target_link_libraries(detect_hbonds_standalone PRIVATE x3dna)
# 
# add_executable(debug_donor_acceptor tools/debug_donor_acceptor.cpp)
# target_link_libraries(debug_donor_acceptor PRIVATE x3dna)
# 
# add_executable(debug_hbond_conflict_resolution tools/debug_hbond_conflict_resolution.cpp)
# target_link_libraries(debug_hbond_conflict_resolution PRIVATE x3dna)
# 
# add_executable(debug_dorg_discrepancy tools/debug_dorg_discrepancy.cpp)
# target_link_libraries(debug_dorg_discrepancy PRIVATE x3dna)
# 
# add_executable(debug_bp_type_id_step_params tools/debug_bp_type_id_step_params.cpp)
# target_link_libraries(debug_bp_type_id_step_params PRIVATE x3dna)
# 
# add_executable(debug_frame_json tools/debug_frame_json.cpp)
# target_link_libraries(debug_frame_json PRIVATE x3dna)
# 
# add_executable(investigate_missing_pairs tools/investigate_missing_pairs.cpp)
# target_link_libraries(investigate_missing_pairs PRIVATE x3dna nlohmann_json)
# 
# add_executable(test_residue_matching_by_pdb_props tools/test_residue_matching_by_pdb_props.cpp)
# target_link_libraries(test_residue_matching_by_pdb_props PRIVATE x3dna)
# 
# add_executable(fix_residue_indices_from_json tools/fix_residue_indices_from_json.cpp)
# target_link_libraries(fix_residue_indices_from_json PRIVATE x3dna)
# 
# add_executable(debug_direction_vectors tools/debug_direction_vectors.cpp)
# target_link_libraries(debug_direction_vectors PRIVATE x3dna)
# 
# add_executable(debug_protocol_6v9q tools/debug_protocol_6v9q.cpp)
# target_link_libraries(debug_protocol_6v9q PRIVATE x3dna)
# 
# add_executable(compare_validation_discrepancy tools/compare_validation_discrepancy.cpp)
# target_link_libraries(compare_validation_discrepancy PRIVATE x3dna)
# 
# add_executable(compare_quality_score_components tools/compare_quality_score_components.cpp)
# target_link_libraries(compare_quality_score_components PRIVATE x3dna)
# 
# add_executable(trace_pair_selection tools/trace_pair_selection.cpp)
# target_link_libraries(trace_pair_selection PRIVATE x3dna)

# Applications
add_executable(find_pair_app apps/find_pair_app.cpp)
target_link_libraries(find_pair_app PRIVATE x3dna)

add_executable(analyze_app apps/analyze_app.cpp)
target_link_libraries(analyze_app PRIVATE x3dna)

# add_executable(compare_bp_type_id_calculation tools/compare_bp_type_id_calculation.cpp)
# target_link_libraries(compare_bp_type_id_calculation PRIVATE x3dna)

# Commented out - tools don't exist yet
# add_executable(compare_frames_and_step_params tools/compare_frames_and_step_params.cpp)
# target_link_libraries(compare_frames_and_step_params PRIVATE x3dna)

# add_executable(verify_json_indices_order tools/verify_json_indices_order.cpp)
# target_link_libraries(verify_json_indices_order PRIVATE x3dna nlohmann_json)

# add_executable(test_bpstep_par_equivalence tools/test_bpstep_par_equivalence.cpp)
# target_link_libraries(test_bpstep_par_equivalence PRIVATE x3dna)

# add_executable(analyze_validation_difference tools/analyze_validation_difference.cpp)
# target_link_libraries(analyze_validation_difference PRIVATE x3dna)

# Documentation
if(BUILD_DOCS)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        add_subdirectory(docs)
    endif()
endif()

# Installation and CMake Package Export
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Install the library
install(TARGETS x3dna
    EXPORT x3dnaTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install headers
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
)

# Install resources (templates and config files)
install(DIRECTORY resources/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/x3dna/resources
)

# Export targets
install(EXPORT x3dnaTargets
    FILE x3dnaTargets.cmake
    NAMESPACE x3dna::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/x3dna
)

# Create config file for find_package()
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/x3dnaConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/x3dnaConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/x3dna
)

# Create version file
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/x3dnaConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Install config files
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/x3dnaConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/x3dnaConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/x3dna
)


cmake_minimum_required(VERSION 3.15)
project(x3dna_original)

set(CMAKE_C_STANDARD 99)

# Set Release mode by default if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set optimization flags based on build type
set(CMAKE_C_FLAGS_DEBUG "-g -O0 -Wall")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wall")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g -O2 -Wall")
set(CMAKE_C_FLAGS_MINSIZEREL "-Os -DNDEBUG -Wall")

# Base flags (common to all build types)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

# Include directories
include_directories(include)

# Source files
set(SOURCES
    src/find_pair.c
    src/cmn_fncs.c
    src/ana_fncs.c
    src/app_fncs.c
    src/nrutil.c
    src/fncs_slre.c
    src/json_writer.c
)

set(ANALYZE_SOURCES
    src/analyze.c
    src/cmn_fncs.c
    src/ana_fncs.c
    src/app_fncs.c
    src/nrutil.c
    src/fncs_slre.c
    src/json_writer.c
)

# Find pair executable
add_executable(find_pair_original ${SOURCES})

# Analyze executable
add_executable(analyze_original ${ANALYZE_SOURCES})

# Combined find_pair + analyze executable
set(COMBINED_SOURCES
    src/find_pair_analyze.c
    src/find_pair.c
    src/analyze.c
    src/cmn_fncs.c
    src/ana_fncs.c
    src/app_fncs.c
    src/nrutil.c
    src/fncs_slre.c
    src/json_writer.c
)
add_executable(find_pair_analyze ${COMBINED_SOURCES})
target_compile_definitions(find_pair_analyze PRIVATE BUILD_COMBINED=1)

# Test H-bond detection executable
set(TEST_HBOND_SOURCES
    src/test_hbond_detection.c
    src/cmn_fncs.c
    src/ana_fncs.c
    src/app_fncs.c
    src/nrutil.c
    src/fncs_slre.c
    src/json_writer.c
)
add_executable(test_hbond_detection ${TEST_HBOND_SOURCES})

# Test donor/acceptor executable (simpler, no PDB needed)
set(TEST_DONOR_ACCEPTOR_SOURCES
    src/test_donor_acceptor.c
    src/cmn_fncs.c
    src/ana_fncs.c
    src/app_fncs.c
    src/nrutil.c
    src/fncs_slre.c
    src/json_writer.c
)
add_executable(test_donor_acceptor ${TEST_DONOR_ACCEPTOR_SOURCES})

# Test N3->N2 debug executable (detailed debugging for specific H-bond issue)
set(TEST_N3_N2_DEBUG_SOURCES
    src/test_n3_n2_debug.c
    src/cmn_fncs.c
    src/ana_fncs.c
    src/app_fncs.c
    src/nrutil.c
    src/fncs_slre.c
    src/json_writer.c
)
add_executable(test_n3_n2_debug ${TEST_N3_N2_DEBUG_SOURCES})

# Test initial H-bond detection (isolated - no conflict resolution or validation)
set(TEST_HBOND_INITIAL_SOURCES
    src/test_hbond_initial.c
    src/cmn_fncs.c
    src/ana_fncs.c
    src/app_fncs.c
    src/nrutil.c
    src/fncs_slre.c
    src/json_writer.c
)
add_executable(test_hbond_initial ${TEST_HBOND_INITIAL_SOURCES})

# Test conflict resolution (isolated - takes initial H-bonds, applies hb_atompair)
set(TEST_HBOND_CONFLICT_SOURCES
    src/test_hbond_conflict.c
    src/cmn_fncs.c
    src/ana_fncs.c
    src/app_fncs.c
    src/nrutil.c
    src/fncs_slre.c
    src/json_writer.c
)
add_executable(test_hbond_conflict ${TEST_HBOND_CONFLICT_SOURCES})

# Test conflict resolution with detailed debug output
set(TEST_HBOND_CONFLICT_DEBUG_SOURCES
    src/test_hbond_conflict_debug.c
    src/cmn_fncs.c
    src/ana_fncs.c
    src/app_fncs.c
    src/nrutil.c
    src/fncs_slre.c
    src/json_writer.c
)
add_executable(test_hbond_conflict_debug ${TEST_HBOND_CONFLICT_DEBUG_SOURCES})

# Test validation (isolated - takes H-bonds after conflict, applies validate_hbonds)
set(TEST_HBOND_VALIDATION_SOURCES
    src/test_hbond_validation.c
    src/cmn_fncs.c
    src/ana_fncs.c
    src/app_fncs.c
    src/nrutil.c
    src/fncs_slre.c
    src/json_writer.c
)
add_executable(test_hbond_validation ${TEST_HBOND_VALIDATION_SOURCES})

# Generate residue ordering JSON (legacy)
set(GENERATE_RESIDUE_ORDERING_SOURCES
    src/generate_residue_ordering_json.c
    src/cmn_fncs.c
    src/ana_fncs.c
    src/app_fncs.c
    src/nrutil.c
    src/fncs_slre.c
    src/json_writer.c
)
add_executable(generate_residue_ordering_json ${GENERATE_RESIDUE_ORDERING_SOURCES})

# Output directory
set_target_properties(find_pair_original analyze_original find_pair_analyze test_hbond_detection test_donor_acceptor test_n3_n2_debug test_hbond_initial test_hbond_conflict test_hbond_validation generate_residue_ordering_json
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)


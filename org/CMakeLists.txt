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

# Output directory
set_target_properties(find_pair_original analyze_original find_pair_analyze test_hbond_detection test_donor_acceptor
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)


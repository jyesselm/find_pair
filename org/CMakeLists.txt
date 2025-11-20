cmake_minimum_required(VERSION 3.15)
project(x3dna_original_debug)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -Wall")

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

# Output directory
set_target_properties(find_pair_original analyze_original find_pair_analyze
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)


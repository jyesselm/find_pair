# Dependencies Configuration

include(FetchContent)
include(GNUInstallDirs)

# nlohmann/json - Header-only JSON library
# Use the single header file approach to avoid CMake version issues
set(JSON_SINGLE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/_deps/json-single/include/nlohmann)
file(MAKE_DIRECTORY ${JSON_SINGLE_INCLUDE_DIR})

# Create interface library
add_library(nlohmann_json INTERFACE)
target_include_directories(nlohmann_json INTERFACE
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_deps/json-single/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Download the single header file
file(DOWNLOAD
    https://github.com/nlohmann/json/releases/download/v3.11.2/json.hpp
    ${JSON_SINGLE_INCLUDE_DIR}/json.hpp
    SHOW_PROGRESS
    STATUS JSON_DOWNLOAD_STATUS
)

# Check if download was successful
list(GET JSON_DOWNLOAD_STATUS 0 JSON_DOWNLOAD_CODE)
if(NOT JSON_DOWNLOAD_CODE EQUAL 0)
    message(WARNING "Failed to download nlohmann/json. You may need to download it manually.")
    message(STATUS "Please download json.hpp from https://github.com/nlohmann/json")
    message(STATUS "and place it in: ${JSON_SINGLE_INCLUDE_DIR}/json.hpp")
else()
    message(STATUS "Successfully downloaded nlohmann/json")
endif()

# Create alias for consistency
if(NOT TARGET nlohmann_json::nlohmann_json)
    add_library(nlohmann_json::nlohmann_json ALIAS nlohmann_json)
endif()

# Install nlohmann_json header for consumers
install(FILES ${JSON_SINGLE_INCLUDE_DIR}/json.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nlohmann
)

# Export nlohmann_json target
install(TARGETS nlohmann_json
    EXPORT x3dnaTargets
)

# GEMMI - Structure file parsing library (for PDB/CIF/mmCIF support)
# Build as compiled library for full functionality
FetchContent_Declare(
    gemmi
    GIT_REPOSITORY https://github.com/project-gemmi/gemmi.git
    GIT_TAG v0.6.5
)

# Configure GEMMI build options
set(BUILD_SHARED_LIBS_SAVED ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "Build static library")

FetchContent_GetProperties(gemmi)
if(NOT gemmi_POPULATED)
    FetchContent_Populate(gemmi)

    # Build gemmi_cpp library (compiled GEMMI code)
    add_subdirectory(${gemmi_SOURCE_DIR} ${gemmi_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# Restore BUILD_SHARED_LIBS
set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_SAVED})

# Find ZLIB for compressed file support (optional but recommended)
find_package(ZLIB QUIET)
if(ZLIB_FOUND)
    message(STATUS "GEMMI: ZLIB support enabled for compressed files")
else()
    message(STATUS "GEMMI: ZLIB not found, compressed files not supported")
endif()

message(STATUS "GEMMI library configured for PDB/CIF/mmCIF parsing")

# Google Test (for testing)
if(BUILD_TESTS)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
    )
    
    FetchContent_MakeAvailable(googletest)
    
    # Disable tests for googletest itself
    set(BUILD_GMOCK ON CACHE INTERNAL "")
    set(INSTALL_GTEST OFF CACHE INTERNAL "")
endif()

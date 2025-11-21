# Dependencies Configuration

include(FetchContent)

# nlohmann/json - Header-only JSON library
# Use the single header file approach to avoid CMake version issues
set(JSON_SINGLE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/_deps/json-single/include/nlohmann)
file(MAKE_DIRECTORY ${JSON_SINGLE_INCLUDE_DIR})

# Create interface library
add_library(nlohmann_json INTERFACE)
target_include_directories(nlohmann_json INTERFACE ${CMAKE_BINARY_DIR}/_deps/json-single/include)

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

# Embedding x3dna in Your Project

This guide explains how to use the x3dna library in your C++ project.

## Option 1: Add as a Subdirectory (Recommended for Development)

Copy the entire x3dna directory into your project's `external/` folder:

```
your_project/
  external/
    x3dna/
      include/
      src/
      resources/
      cmake/
      CMakeLists.txt
  src/
  CMakeLists.txt
```

In your project's `CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.15)
project(your_project)

set(CMAKE_CXX_STANDARD 17)

# Disable x3dna tests when building as subdirectory
set(BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

# Add x3dna as subdirectory
add_subdirectory(external/x3dna)

# Your executable
add_executable(your_app src/main.cpp)
target_link_libraries(your_app PRIVATE x3dna)
```

## Option 2: Install and find_package

First, build and install x3dna:

```bash
cd x3dna
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build . --config Release
cmake --install .
```

Then in your project's `CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.15)
project(your_project)

set(CMAKE_CXX_STANDARD 17)

# Find the installed x3dna package
find_package(x3dna REQUIRED)

add_executable(your_app src/main.cpp)
target_link_libraries(your_app PRIVATE x3dna::x3dna)
```

Configure with the install prefix:

```bash
cmake -DCMAKE_PREFIX_PATH=/path/to/install ..
```

## Initialization

Before using x3dna functionality, initialize the ResourceLocator in your main():

```cpp
#include <x3dna/config/resource_locator.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/protocols/find_pair_protocol.hpp>

int main() {
    // Option 1: Explicit path (most portable)
    x3dna::config::ResourceLocator::initialize("/path/to/x3dna/resources");

    // Option 2: Auto-detect from environment or relative paths
    // Searches: "resources", "../resources", etc.
    // Falls back to X3DNA_HOMEDIR environment variable
    x3dna::config::ResourceLocator::initialize_from_environment();

    // Option 3: Use CMake-provided path (when using find_package)
    // In CMake: target_compile_definitions(your_app PRIVATE
    //           X3DNA_RESOURCES="${X3DNA_RESOURCES_DIR}")
    // Then in code:
    // x3dna::config::ResourceLocator::initialize(X3DNA_RESOURCES);

    // Now you can use x3dna
    x3dna::io::PdbParser parser;
    auto structure = parser.parse_file("input.pdb");

    x3dna::protocols::FindPairProtocol protocol;
    auto pairs = protocol.run(structure);

    return 0;
}
```

## Minimal Example

```cpp
#include <x3dna/config/resource_locator.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/protocols/find_pair_protocol.hpp>
#include <x3dna/io/json_writer.hpp>
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file>\n";
        return 1;
    }

    // Initialize resources (auto-detect)
    if (!x3dna::config::ResourceLocator::initialize_from_environment()) {
        std::cerr << "Could not find x3dna resources\n";
        return 1;
    }

    // Parse PDB
    x3dna::io::PdbParser parser;
    auto structure = parser.parse_file(argv[1]);

    // Find base pairs
    x3dna::protocols::FindPairProtocol protocol;
    auto result = protocol.run(structure);

    // Write JSON output
    x3dna::io::JsonWriter writer;
    writer.write_pairs("output.json", result);

    std::cout << "Found " << result.pairs.size() << " base pairs\n";
    return 0;
}
```

## Required Dependencies

x3dna requires:
- C++17 compiler
- nlohmann/json (included when using as subdirectory)

## Resources Directory Structure

The `resources/` directory must contain:
```
resources/
  templates/      # Base frame templates (Atomic_*.pdb)
  config/         # Configuration files
    modified_nucleotides.json
    atomlist.dat
    baselist.dat
    misc_3dna.par
```

When embedding, copy the entire `resources/` directory alongside your binary or set an explicit path during initialization.

# Stage 0: Project Setup & Infrastructure

## Objectives

Set up the project structure, build system, testing framework, and development environment for the modernization effort.

## Duration

**1 week**

## Dependencies

- None (starting point)

## Tasks

### Task 0.1: Create Directory Structure
- [ ] Create `include/x3dna/` directory hierarchy
- [ ] Create `src/` directory hierarchy (mirrors include)
- [ ] Create `apps/` directory
- [ ] Create `tests/` directory with subdirectories:
  - [ ] `tests/unit/`
  - [ ] `tests/integration/`
  - [ ] `tests/regression/`
- [ ] Create `examples/` directory
- [ ] Create `docs/api/` for generated documentation

**Deliverable**: Complete directory structure matching the plan

### Task 0.2: Set Up CMake Build System
- [ ] Create root `CMakeLists.txt`
- [ ] Configure C++17/20 standard
- [ ] Set up library target structure
- [ ] Configure compiler flags (warnings, optimization)
- [ ] Set up install targets
- [ ] Configure for both Debug and Release builds

**Files to Create**:
- `CMakeLists.txt` (root)
- `cmake/FindX3DNA.cmake` (if needed)
- `.clang-format` (code formatting)
- `.clang-tidy` (static analysis)

**Deliverable**: `CMakeLists.txt` that builds empty library structure

### Task 0.3: Set Up Testing Framework
- [ ] Add Google Test via FetchContent or find_package
- [ ] Create test executable structure
- [ ] Set up test discovery
- [ ] Configure test output format
- [ ] Add code coverage tools (gcov/lcov)
- [ ] Create test helper utilities

**Files to Create**:
- `tests/CMakeLists.txt`
- `tests/test_main.cpp`
- `tests/TestHelpers.hpp`

**Deliverable**: Test framework that can run empty tests

### Task 0.4: Set Up Dependencies
- [ ] Add nlohmann/json (header-only, via FetchContent)
- [ ] Configure dependency management
- [ ] Document dependency versions
- [ ] Set up optional dependencies (if any)

**Files to Create**:
- `cmake/Dependencies.cmake`
- `docs/DEPENDENCIES.md`

**Deliverable**: All dependencies configured and documented

### Task 0.5: Set Up Development Tools
- [ ] Configure IDE project files (if needed)
- [ ] Set up pre-commit hooks (optional)
- [ ] Configure debugger settings
- [ ] Set up continuous integration (CI) configuration
- [ ] Create `.gitignore` file

**Files to Create**:
- `.gitignore`
- `.github/workflows/ci.yml` (or equivalent)
- `.editorconfig`

**Deliverable**: Development environment ready

### Task 0.6: Create Base Header Structure
- [ ] Create namespace structure in headers
- [ ] Add forward declarations header
- [ ] Create common types header
- [ ] Add version information header
- [ ] Create main library header (`x3dna.hpp`)

**Files to Create**:
- `include/x3dna/x3dna.hpp` (main include)
- `include/x3dna/version.hpp`
- `include/x3dna/forward_declarations.hpp`
- `include/x3dna/common_types.hpp`

**Deliverable**: Base header structure that compiles

### Task 0.7: Set Up Documentation Generation
- [ ] Configure Doxygen
- [ ] Create Doxygen configuration file
- [ ] Set up documentation build target
- [ ] Create initial README files
- [ ] Document build process

**Files to Create**:
- `Doxyfile`
- `README.md` (project root)
- `docs/BUILD.md`

**Deliverable**: Documentation can be generated

### Task 0.8: Create Initial Test Infrastructure
- [ ] Create test fixtures for common data
- [ ] Set up JSON comparison utilities
- [ ] Create test data directory structure
- [ ] Add helper functions for loading legacy JSON
- [ ] Create mock objects (if needed)
- [ ] Create integration test base classes
- [ ] Create PDB/JSON pair discovery utility

**Files to Create**:
- `tests/TestFixtures.hpp`
- `tests/JsonTestHelpers.hpp`
- `tests/integration/IntegrationTestBase.hpp`
- `tests/integration/TestDataDiscovery.hpp`
- `tests/integration/JsonComparison.hpp`
- `tests/data/` directory

**Deliverable**: Test infrastructure ready for use, including integration test framework

## Testing Plan

### Unit Tests
- [ ] Test that CMake configuration works
- [ ] Test that empty library compiles
- [ ] Test that test framework runs
- [ ] Test that dependencies are found

### Integration Tests
- [ ] Test full build process (clean → configure → build → test)
- [ ] Test install process
- [ ] Test documentation generation

### Validation
- [ ] Verify directory structure matches plan
- [ ] Verify all files compile without errors
- [ ] Verify tests can be run
- [ ] Verify documentation can be generated

## Deliverables

1. ✅ Complete project directory structure
2. ✅ Working CMake build system
3. ✅ Configured testing framework (Google Test)
4. ✅ All dependencies configured
5. ✅ Development tools set up
6. ✅ Base header structure
7. ✅ Documentation generation configured
8. ✅ Test infrastructure ready

## Success Criteria

- [ ] Project can be cloned and built from scratch
- [ ] All directories exist and are properly structured
- [ ] CMake generates build files successfully
- [ ] Empty library compiles without errors
- [ ] Test framework runs (even with empty tests)
- [ ] Dependencies are properly configured
- [ ] Documentation can be generated
- [ ] CI/CD pipeline runs (if configured)

## Code Structure

```
find_pair_2/
├── CMakeLists.txt
├── .gitignore
├── .clang-format
├── .clang-tidy
├── Doxyfile
├── README.md
├── include/
│   └── x3dna/
│       ├── x3dna.hpp
│       ├── version.hpp
│       ├── forward_declarations.hpp
│       └── common_types.hpp
├── src/
│   └── (empty for now)
├── apps/
│   └── (empty for now)
├── tests/
│   ├── CMakeLists.txt
│   ├── test_main.cpp
│   ├── TestHelpers.hpp
│   ├── TestFixtures.hpp
│   ├── JsonTestHelpers.hpp
│   ├── unit/
│   ├── integration/
│   ├── regression/
│   └── data/
├── examples/
│   └── (empty for now)
├── docs/
│   ├── api/ (generated)
│   └── BUILD.md
└── cmake/
    ├── Dependencies.cmake
    └── FindX3DNA.cmake
```

## Example CMakeLists.txt Structure

```cmake
cmake_minimum_required(VERSION 3.15)
project(x3dna VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build options
option(BUILD_TESTS "Build tests" ON)
option(BUILD_EXAMPLES "Build examples" ON)
option(BUILD_DOCS "Build documentation" ON)

# Include directories
include_directories(include)

# Dependencies
include(cmake/Dependencies.cmake)

# Library
add_library(x3dna INTERFACE)
target_include_directories(x3dna INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Testing
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Examples
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Documentation
if(BUILD_DOCS)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        add_subdirectory(docs)
    endif()
endif()
```

## Risks & Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| CMake configuration issues | High | Test on multiple platforms early |
| Dependency version conflicts | Medium | Pin dependency versions, document clearly |
| Test framework setup problems | Medium | Use well-established framework (Google Test) |
| CI/CD configuration issues | Low | Start simple, add complexity later |

## Next Stage

After completing Stage 0, proceed to **Stage 1: Geometry Classes** (`STAGE_01_GEOMETRY.md`)

---

*Estimated Completion: Week 1*


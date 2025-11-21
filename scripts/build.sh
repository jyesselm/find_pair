#!/bin/bash
# Build script for x3dna project with maximum CPU utilization using Ninja
# Usage: ./scripts/build.sh [Debug|Release|RelWithDebInfo]

set -e

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$PROJECT_ROOT/build"

# Check if Ninja is available
if ! command -v ninja &> /dev/null; then
    echo "Error: Ninja build system is not installed."
    echo ""
    echo "Install Ninja:"
    echo "  macOS:   brew install ninja"
    echo "  Linux:   sudo apt-get install ninja-build  # Debian/Ubuntu"
    echo "           sudo yum install ninja-build      # RHEL/CentOS"
    echo ""
    exit 1
fi

# Detect number of CPU cores
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    NUM_CPUS=$(sysctl -n hw.ncpu)
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    # Linux
    NUM_CPUS=$(nproc)
else
    # Fallback
    NUM_CPUS=4
fi

# Handle help flag
if [[ "$1" == "--help" ]] || [[ "$1" == "-h" ]]; then
    echo "Usage: $0 [Debug|Release|RelWithDebInfo|MinSizeRel]"
    echo ""
    echo "Build the x3dna project with maximum CPU utilization using Ninja."
    echo ""
    echo "Options:"
    echo "  Debug          Build with debug symbols (-g, -O0)"
    echo "  Release        Build optimized (-O3, -DNDEBUG) [default]"
    echo "  RelWithDebInfo Build optimized with debug symbols (-O2, -g)"
    echo "  MinSizeRel     Build optimized for size (-Os, -DNDEBUG)"
    echo ""
    echo "Examples:"
    echo "  $0              # Build in Release mode"
    echo "  $0 Debug        # Build in Debug mode"
    echo "  $0 RelWithDebInfo  # Build with optimizations and debug symbols"
    echo ""
    echo "Note: Requires Ninja build system to be installed."
    exit 0
fi

# Build type (default: Release)
BUILD_TYPE="${1:-Release}"

# Validate build type
if [[ ! "$BUILD_TYPE" =~ ^(Debug|Release|RelWithDebInfo|MinSizeRel)$ ]]; then
    echo "Error: Invalid build type '$BUILD_TYPE'"
    echo "Usage: $0 [Debug|Release|RelWithDebInfo|MinSizeRel]"
    echo "Use '$0 --help' for more information"
    exit 1
fi

echo "=========================================="
echo "Building x3dna project (Ninja)"
echo "=========================================="
echo "Build type: $BUILD_TYPE"
echo "CPU cores:  $NUM_CPUS"
echo "Build dir:  $BUILD_DIR"
echo "Generator:  Ninja"
echo "=========================================="
echo ""

# Create build directory if it doesn't exist
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure CMake if needed
NEED_RECONFIGURE=false

if [ ! -f "CMakeCache.txt" ]; then
    NEED_RECONFIGURE=true
elif [ "$PROJECT_ROOT/CMakeLists.txt" -nt "CMakeCache.txt" ]; then
    NEED_RECONFIGURE=true
else
    # Check if build type changed
    CACHED_BUILD_TYPE=$(grep -i "CMAKE_BUILD_TYPE:STRING" CMakeCache.txt 2>/dev/null | cut -d'=' -f2 || echo "")
    if [ "$CACHED_BUILD_TYPE" != "$BUILD_TYPE" ]; then
        NEED_RECONFIGURE=true
    fi
    # Check if generator changed (should be Ninja)
    CACHED_GENERATOR=$(grep -i "CMAKE_GENERATOR:INTERNAL" CMakeCache.txt 2>/dev/null | cut -d'=' -f2 || echo "")
    if [[ ! "$CACHED_GENERATOR" =~ [Nn]inja ]]; then
        NEED_RECONFIGURE=true
    fi
fi

if [ "$NEED_RECONFIGURE" = true ]; then
    echo "Configuring CMake with Ninja generator..."
    cmake "$PROJECT_ROOT" \
        -G Ninja \
        -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    echo ""
fi

# Build with maximum parallel jobs using Ninja
echo "Building with $NUM_CPUS parallel jobs (Ninja)..."
ninja -j "$NUM_CPUS"

echo ""
echo "=========================================="
echo "Build complete!"
echo "=========================================="


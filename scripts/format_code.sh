#!/bin/bash
# Format all C++ source files using clang-format

set -e

# Find all C++ files, excluding build directories and git
find . -type f \( -name "*.cpp" -o -name "*.hpp" -o -name "*.h" -o -name "*.cc" -o -name "*.cxx" \) \
    ! -path "./build/*" \
    ! -path "./.git/*" \
    ! -path "./org/build/*" \
    ! -path "./Testing/*" \
    -exec clang-format -i {} +

echo "Formatting complete!"


/**
 * @file x3dna.hpp
 * @brief Main header file for the X3DNA modernized library
 * @version 1.0.0
 */

#pragma once

#include <x3dna/version.hpp>
#include <x3dna/forward_declarations.hpp>
#include <x3dna/common_types.hpp>

/**
 * @namespace x3dna
 * @brief Main namespace for the X3DNA library
 */
namespace x3dna {

/**
 * @brief Get library version string
 * @return Version string (e.g., "1.0.0")
 */
inline const char* version() {
    return X3DNA_VERSION_STRING;
}

} // namespace x3dna

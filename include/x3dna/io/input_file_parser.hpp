/**
 * @file input_file_parser.hpp
 * @brief Parser for X3DNA .inp input files
 */

#pragma once

#include <string>
#include <filesystem>
#include <vector>
#include <x3dna/core/base_pair.hpp>

namespace x3dna {
namespace io {

/**
 * @struct InputData
 * @brief Data extracted from .inp file
 */
struct InputData {
    std::filesystem::path pdb_file;         // Path to PDB file
    std::string output_file;                // Output file name
    int duplex_number = 2;                  // Duplex number (usually 2)
    size_t num_base_pairs = 0;              // Number of base pairs
    int flags = 0;                          // Flags (explicit bp numbering, etc.)
    std::vector<core::BasePair> base_pairs; // Base pairs
    std::string criteria_line;              // Base-pair criteria line (if present)
    std::vector<std::string> helix_info;    // Helix information lines
};

/**
 * @class InputFileParser
 * @brief Parses X3DNA .inp input files
 */
class InputFileParser {
public:
    /**
     * @brief Parse .inp file
     * @param input_file Path to .inp file
     * @return InputData structure with parsed data
     * @throws std::runtime_error if file cannot be read or parsed
     */
    static InputData parse(const std::filesystem::path& input_file);

    /**
     * @brief Parse .inp file from stream
     * @param stream Input stream
     * @return InputData structure with parsed data
     */
    static InputData parse_stream(std::istream& stream);

private:
    /**
     * @brief Parse base pair line
     * @param line Input line
     * @param line_number Line number (for error reporting)
     * @return Base pair indices (residue_idx1, residue_idx2) - 0-based
     */
    static std::pair<size_t, size_t> parse_base_pair_line(const std::string& line,
                                                          size_t line_number);
};

} // namespace io
} // namespace x3dna

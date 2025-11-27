/**
 * @file input_file_writer.hpp
 * @brief Writer for X3DNA .inp input files
 */

#pragma once

#include <string>
#include <filesystem>
#include <vector>
#include <x3dna/core/base_pair.hpp>

namespace x3dna {
namespace io {

/**
 * @class InputFileWriter
 * @brief Writes X3DNA .inp input files
 */
class InputFileWriter {
public:
    /**
     * @brief Write .inp file
     * @param output_path Path to output .inp file
     * @param pdb_file Path to PDB file (will be written in .inp)
     * @param base_pairs Vector of base pairs to write
     * @param duplex_number Duplex number (default: 2)
     * @param flags Flags value (default: 1 for explicit bp numbering)
     */
    static void write(const std::filesystem::path& output_path,
                      const std::filesystem::path& pdb_file,
                      const std::vector<core::BasePair>& base_pairs,
                      int duplex_number = 2,
                      int flags = 1);

    /**
     * @brief Write .inp file with additional options
     * @param output_path Path to output .inp file
     * @param pdb_file Path to PDB file
     * @param output_file_name Output file name for analyze step
     * @param base_pairs Vector of base pairs to write
     * @param duplex_number Duplex number
     * @param flags Flags value
     */
    static void write(const std::filesystem::path& output_path,
                      const std::filesystem::path& pdb_file,
                      const std::string& output_file_name,
                      const std::vector<core::BasePair>& base_pairs,
                      int duplex_number = 2,
                      int flags = 1);

private:
    /**
     * @brief Generate default output file name from PDB file
     */
    static std::string default_output_filename(const std::filesystem::path& pdb_file);
};

} // namespace io
} // namespace x3dna


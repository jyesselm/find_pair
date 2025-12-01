/**
 * @file input_file_writer.hpp
 * @brief Writer for X3DNA .inp input files and ref_frames files
 */

#pragma once

#include <string>
#include <filesystem>
#include <vector>
#include <map>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/algorithms/parameter_calculator.hpp>

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

    /**
     * @brief Write ref_frames file (ref_frames_modern.dat format)
     * @param output_path Path to output file
     * @param base_pairs Vector of base pairs with reference frames
     * @param structure Structure containing residue information
     */
    static void write_ref_frames(const std::filesystem::path& output_path,
                                  const std::vector<core::BasePair>& base_pairs,
                                  const core::Structure& structure);

    /**
     * @brief Write ref_frames file with legacy pair ordering
     * @param output_path Path to output file
     * @param base_pairs Vector of base pairs with reference frames
     * @param structure Structure containing residue information
     * @param legacy_pair_ordering Map from (min_idx, max_idx) -> first residue index in legacy
     *        If first_idx == max_idx, then legacy had larger-first ordering
     */
    static void write_ref_frames(const std::filesystem::path& output_path,
                                  const std::vector<core::BasePair>& base_pairs,
                                  const core::Structure& structure,
                                  const std::map<std::pair<int,int>, int>& legacy_pair_ordering);

    /**
     * @brief Parse legacy .inp file to get pair ordering
     * @param inp_file Path to legacy .inp file
     * @return Map from (min_idx, max_idx) -> first residue index in legacy
     */
    static std::map<std::pair<int,int>, int> parse_legacy_inp_ordering(
        const std::filesystem::path& inp_file);

private:
    /**
     * @brief Generate default output file name from PDB file
     */
    static std::string default_output_filename(const std::filesystem::path& pdb_file);

    /**
     * @brief Format residue description for ref_frames output
     */
    static std::string format_residue_description(const core::Residue& residue);
};

} // namespace io
} // namespace x3dna


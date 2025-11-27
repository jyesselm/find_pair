/**
 * @file input_file_writer.cpp
 * @brief InputFileWriter implementation
 */

#include <x3dna/io/input_file_writer.hpp>
#include <fstream>
#include <iomanip>

namespace x3dna {
namespace io {

void InputFileWriter::write(const std::filesystem::path& output_path,
                             const std::filesystem::path& pdb_file,
                             const std::vector<core::BasePair>& base_pairs,
                             int duplex_number,
                             int flags) {
    std::string output_file_name = default_output_filename(pdb_file);
    write(output_path, pdb_file, output_file_name, base_pairs, duplex_number, flags);
}

void InputFileWriter::write(const std::filesystem::path& output_path,
                             const std::filesystem::path& pdb_file,
                             const std::string& output_file_name,
                             const std::vector<core::BasePair>& base_pairs,
                             int duplex_number,
                             int flags) {
    std::ofstream out(output_path);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + output_path.string());
    }

    // Line 1: PDB file path
    // Write the path as a string - use the original path string as provided
    // This preserves relative paths and matches legacy behavior
    out << pdb_file.string() << "\n";

    // Line 2: Output file name
    out << output_file_name << "\n";

    // Line 3: Duplex number
    out << "    " << duplex_number << "         # duplex\n";

    // Line 4: Number of base pairs
    out << std::setw(5) << base_pairs.size() << "         # number of base-pairs\n";

    // Line 5: Flags
    out << "    " << flags << " " << std::setw(5) << 0
        << "    # explicit bp numbering/hetero atoms\n";

    // Base pair lines
    // Format: bp_num res1 res2 flag # comment
    // Note: Convert from 0-based to 1-based for residue indices
    for (size_t i = 0; i < base_pairs.size(); ++i) {
        const auto& pair = base_pairs[i];
        size_t bp_num = i + 1;  // 1-based base pair number
        size_t res1 = pair.residue_idx1() + 1;  // Convert to 1-based
        size_t res2 = pair.residue_idx2() + 1;  // Convert to 1-based
        int flag = 0;  // Default flag value

        // Get base pair type string if available
        std::string bp_type = pair.bp_type();
        std::string comment = "";
        if (!bp_type.empty()) {
            comment = " # " + bp_type;
        }

        out << std::setw(5) << bp_num << " " << std::setw(5) << res1 << " "
            << std::setw(5) << res2 << " " << std::setw(5) << flag << comment << "\n";
    }

    out.close();
}

std::string InputFileWriter::default_output_filename(const std::filesystem::path& pdb_file) {
    // Default output file name: remove extension and add .outp
    std::string stem = pdb_file.stem().string();
    return stem + ".outp";
}

} // namespace io
} // namespace x3dna


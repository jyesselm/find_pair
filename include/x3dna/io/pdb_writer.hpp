/**
 * @file pdb_writer.hpp
 * @brief PDB file writer
 */

#pragma once

#include <string>
#include <filesystem>
#include <ostream>
#include <x3dna/core/structure.hpp>

namespace x3dna {
namespace io {

/**
 * @class PdbWriter
 * @brief Writes Structure objects to PDB file format
 */
class PdbWriter {
public:
    /**
     * @brief Write Structure to PDB file
     * @param structure Structure to write
     * @param path Output file path
     */
    void write_file(const core::Structure& structure, const std::filesystem::path& path);

    /**
     * @brief Write Structure to output stream
     * @param structure Structure to write
     * @param stream Output stream
     */
    void write_stream(const core::Structure& structure, std::ostream& stream);

    /**
     * @brief Convert Structure to PDB format string
     * @param structure Structure to convert
     * @return PDB format string
     */
    std::string to_string(const core::Structure& structure);

private:
    /**
     * @brief Format atom record line (ATOM or HETATM)
     * @param atom Atom to format
     * @param atom_serial Serial number (1-based)
     * @return Formatted PDB line
     */
    std::string format_atom_line(const core::Atom& atom, int atom_serial) const;

    /**
     * @brief Format coordinate value for PDB (8.3 format)
     * @param coord Coordinate value
     * @return Formatted string (8 characters)
     */
    static std::string format_coordinate(double coord);
};

} // namespace io
} // namespace x3dna


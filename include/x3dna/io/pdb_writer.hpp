/**
 * @file pdb_writer.hpp
 * @brief PDB file writer
 */

#pragma once

#include <string>
#include <filesystem>
#include <ostream>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/structure/residue.hpp>  // Polymorphic types

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

    // === Polymorphic Structure writing ===

    /**
     * @brief Write polymorphic Structure to PDB file
     * @param structure Structure to write
     * @param path Output file path
     */
    void write_file(const core::structure::Structure& structure, const std::filesystem::path& path);

    /**
     * @brief Write polymorphic Structure to output stream
     * @param structure Structure to write
     * @param stream Output stream
     */
    void write_stream(const core::structure::Structure& structure, std::ostream& stream);

    /**
     * @brief Convert polymorphic Structure to PDB format string
     * @param structure Structure to convert
     * @return PDB format string
     */
    std::string to_string(const core::structure::Structure& structure);

private:
    /**
     * @brief Format atom record line using Residue context
     * @param atom Atom to format
     * @param residue Residue containing the atom (provides residue-level fields)
     * @param structure Structure containing record type map
     * @param atom_serial Serial number (1-based)
     * @return Formatted PDB line
     */
    std::string format_atom_line(const core::Atom& atom, const core::Residue& residue, const core::Structure& structure,
                                 int atom_serial) const;

    /**
     * @brief Format coordinate value for PDB (8.3 format)
     * @param coord Coordinate value
     * @return Formatted string (8 characters)
     */
    static std::string format_coordinate(double coord);

    /**
     * @brief Format atom record line using polymorphic Residue context
     * @param atom Atom to format
     * @param residue Polymorphic residue containing the atom
     * @param structure Polymorphic structure containing record type map
     * @param atom_serial Serial number (1-based)
     * @return Formatted PDB line
     */
    std::string format_atom_line(const core::Atom& atom, const core::structure::IResidue& residue,
                                 const core::structure::Structure& structure, int atom_serial) const;
};

} // namespace io
} // namespace x3dna

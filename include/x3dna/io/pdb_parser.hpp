/**
 * @file pdb_parser.hpp
 * @brief PDB file parser for reading PDB format files
 *
 * This parser reads PDB files and creates Structure objects containing
 * atoms, residues, and chains. It handles ATOM and HETATM records.
 */

#pragma once

#include <filesystem>
#include <istream>
#include <string>
#include <memory>
#include <stdexcept>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>

namespace x3dna {
namespace io {

/**
 * @class PdbParser
 * @brief Parser for PDB (Protein Data Bank) format files
 *
 * Parses PDB files and creates Structure objects. Supports:
 * - ATOM records (standard atoms)
 * - HETATM records (heteroatoms, optional)
 * - Chain identification
 * - Residue numbering
 * - Alternate conformations
 */
class PdbParser {
public:
    /**
     * @brief Constructor
     */
    PdbParser() = default;

    /**
     * @brief Destructor
     */
    ~PdbParser() = default;

    // Copy and move semantics
    PdbParser(const PdbParser&) = default;
    PdbParser& operator=(const PdbParser&) = default;
    PdbParser(PdbParser&&) noexcept = default;
    PdbParser& operator=(PdbParser&&) noexcept = default;

    /**
     * @brief Parse PDB file from filesystem path
     * @param path Path to PDB file
     * @return Structure object containing parsed data
     * @throws ParseError if file cannot be read or parsed
     */
    core::Structure parse_file(const std::filesystem::path& path);

    /**
     * @brief Parse PDB file from input stream
     * @param stream Input stream containing PDB data
     * @return Structure object containing parsed data
     * @throws ParseError if stream cannot be read or parsed
     */
    core::Structure parse_stream(std::istream& stream);

    /**
     * @brief Parse PDB file from string content
     * @param content String containing PDB file content
     * @return Structure object containing parsed data
     * @throws ParseError if content cannot be parsed
     */
    core::Structure parse_string(const std::string& content);

    /**
     * @brief Set whether to include HETATM records
     * @param value True to include HETATM, false to skip
     */
    void set_include_hetatm(bool value) { include_hetatm_ = value; }

    /**
     * @brief Get whether HETATM records are included
     * @return True if HETATM records are included
     */
    bool include_hetatm() const { return include_hetatm_; }

    /**
     * @brief Set whether to include water molecules (HOH)
     * @param value True to include waters, false to skip
     */
    void set_include_waters(bool value) { include_waters_ = value; }

    /**
     * @brief Get whether water molecules are included
     * @return True if waters are included
     */
    bool include_waters() const { return include_waters_; }

    /**
     * @brief Exception class for parsing errors
     */
    class ParseError : public std::runtime_error {
    public:
        /**
         * @brief Constructor
         * @param message Error message
         * @param line_number Line number where error occurred (0 if unknown)
         */
        ParseError(const std::string& message, size_t line_number = 0);

        /**
         * @brief Get line number where error occurred
         * @return Line number (0 if unknown)
         */
        size_t line_number() const { return line_number_; }

    private:
        size_t line_number_;
    };

private:
    bool include_hetatm_ = false;  // Include HETATM records
    bool include_waters_ = false;  // Include water molecules (HOH)
    bool filter_by_occupancy_ = false;  // Filter atoms by occupancy (default: false, match legacy)

    /**
     * @brief Parse a single ATOM or HETATM line
     * @param line Line from PDB file
     * @param line_number Line number for error reporting
     * @return Atom object parsed from line
     * @throws ParseError if line cannot be parsed
     */
    core::Atom parse_atom_line(const std::string& line, size_t line_number);

    /**
     * @brief Parse a single HETATM line
     * @param line Line from PDB file
     * @param line_number Line number for error reporting
     * @return Atom object parsed from line
     * @throws ParseError if line cannot be parsed
     */
    core::Atom parse_hetatm_line(const std::string& line, size_t line_number);

    /**
     * @brief Parse atom name from PDB line (columns 13-16)
     * @param line PDB line
     * @return Atom name (normalized)
     */
    std::string parse_atom_name(const std::string& line);

    /**
     * @brief Parse residue name from PDB line (columns 18-20)
     * @param line PDB line
     * @return Residue name (normalized)
     */
    std::string parse_residue_name(const std::string& line);

    /**
     * @brief Parse chain ID from PDB line (column 22)
     * @param line PDB line
     * @return Chain ID character (' ' if not specified)
     */
    char parse_chain_id(const std::string& line);

    /**
     * @brief Parse alternate location indicator from PDB line (column 17)
     * @param line PDB line
     * @return Alternate location character (' ' if not specified)
     */
    char parse_alt_loc(const std::string& line);

    /**
     * @brief Parse insertion code from PDB line (column 27)
     * @param line PDB line
     * @return Insertion code character (' ' if not specified)
     */
    char parse_insertion(const std::string& line);

    /**
     * @brief Parse occupancy from PDB line (columns 55-60)
     * @param line PDB line
     * @return Occupancy value (1.0 if not specified or parsing fails)
     */
    double parse_occupancy(const std::string& line);

    /**
     * @brief Parse residue sequence number from PDB line (columns 23-26)
     * @param line PDB line
     * @return Residue sequence number
     * @throws ParseError if cannot parse number
     */
    int parse_residue_seq(const std::string& line);

    /**
     * @brief Parse coordinates from PDB line (columns 31-54)
     * @param line PDB line
     * @return Vector3D with x, y, z coordinates
     * @throws ParseError if coordinates cannot be parsed
     */
    geometry::Vector3D parse_coordinates(const std::string& line);

    /**
     * @brief Check if residue is a water molecule
     * @param residue_name Residue name
     * @return True if residue is water (HOH, WAT, etc.)
     */
    bool is_water(const std::string& residue_name) const;

    /**
     * @brief Normalize atom name (add leading space if needed for 4-char format)
     * @param name Atom name from PDB
     * @return Normalized atom name
     */
    std::string normalize_atom_name(const std::string& name) const;

    /**
     * @brief Normalize residue name (ensure 3 characters with leading space if needed)
     * @param name Residue name from PDB
     * @return Normalized residue name
     */
    std::string normalize_residue_name(const std::string& name) const;

    /**
     * @brief Internal parsing method
     * @param stream Input stream
     * @param pdb_id PDB identifier (extracted from file or provided)
     * @return Structure object
     */
    core::Structure parse_impl(std::istream& stream, const std::string& pdb_id = "");
};

} // namespace io
} // namespace x3dna


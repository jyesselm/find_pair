/**
 * @file pdb_parser.hpp
 * @brief PDB file parser using GEMMI library
 *
 * This parser reads PDB files using GEMMI and creates Structure objects
 * containing atoms, residues, and chains.
 */

#pragma once

#include <filesystem>
#include <istream>
#include <string>
#include <memory>
#include <stdexcept>
#include <map>
#include <vector>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/io/residue_key.hpp>

// Forward declare GEMMI types to avoid header inclusion
namespace gemmi {
struct Structure;
}

namespace x3dna {
namespace io {

/**
 * @class PdbParser
 * @brief Parser for PDB (Protein Data Bank) format files using GEMMI
 *
 * Parses PDB files using GEMMI library and creates Structure objects. Supports:
 * - ATOM records (standard atoms)
 * - HETATM records (heteroatoms, optional)
 * - Chain identification
 * - Residue numbering
 * - Alternate conformations
 * - Compressed files (.pdb.gz)
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
     * @param path Path to PDB file (.pdb or .pdb.gz)
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
    void set_include_hetatm(bool value) {
        include_hetatm_ = value;
    }

    /**
     * @brief Get whether HETATM records are included
     * @return True if HETATM records are included
     */
    bool include_hetatm() const {
        return include_hetatm_;
    }

    /**
     * @brief Set whether to include water molecules (HOH)
     * @param value True to include waters, false to skip
     */
    void set_include_waters(bool value) {
        include_waters_ = value;
    }

    /**
     * @brief Get whether water molecules are included
     * @return True if waters are included
     */
    bool include_waters() const {
        return include_waters_;
    }

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
        size_t line_number() const {
            return line_number_;
        }

    private:
        size_t line_number_;
    };

private:
    bool include_hetatm_ = false; // Include HETATM records
    bool include_waters_ = false; // Include water molecules (HOH)

    /**
     * @brief Convert GEMMI Structure to our Structure
     * @param gemmi_struct GEMMI structure
     * @param pdb_id PDB identifier
     * @return Structure object
     */
    core::Structure convert_gemmi_structure(const gemmi::Structure& gemmi_struct, const std::string& pdb_id);

    /**
     * @brief Normalize atom name from GEMMI format to PDB 4-character format
     * @param name Atom name from GEMMI
     * @return Normalized atom name
     */
    std::string normalize_atom_name_from_gemmi(const std::string& name) const;

    /**
     * @brief Normalize residue name from GEMMI
     * @param name Residue name from GEMMI
     * @return Normalized residue name
     */
    std::string normalize_residue_name_from_gemmi(const std::string& name) const;

    /**
     * @brief Check if residue is a water molecule
     * @param residue_name Residue name
     * @return True if residue is water (HOH, WAT, etc.)
     */
    bool is_water(const std::string& residue_name) const;

    /**
     * @brief Check if residue name is a modified nucleotide
     * @param residue_name Residue name
     * @return True if residue is a known modified nucleotide
     */
    bool is_modified_nucleotide_name(const std::string& residue_name) const;

    /**
     * @brief Check alternate location filter
     * @param alt_loc Alternate location character
     * @return True if this alt_loc should be kept
     */
    bool check_alt_loc_filter(char alt_loc) const;

    /**
     * @brief Apply atom name exact matches (phosphate atoms, etc.)
     * @param name Atom name
     * @return Transformed atom name
     */
    std::string apply_atom_name_exact_matches(const std::string& name) const;

    /**
     * @brief Ensure atom name is exactly 4 characters
     * @param name Atom name
     * @return Padded/truncated atom name
     */
    std::string ensure_atom_name_length(const std::string& name) const;

    /**
     * @brief Normalize atom name (legacy compatibility)
     * @param name Atom name
     * @return Normalized atom name
     */
    std::string normalize_atom_name(const std::string& name) const;

    /**
     * @brief Normalize residue name (legacy compatibility)
     * @param name Residue name
     * @return Normalized residue name
     */
    std::string normalize_residue_name(const std::string& name) const;

    /**
     * @brief Build Structure from grouped residue atoms
     * @param pdb_id Structure identifier
     * @param residue_atoms Map of residue key to atoms
     * @param legacy_idx_map Map of residue key to legacy index (preserves PDB file order)
     * @param chain_order Vector of chain IDs in file encounter order
     * @return Structure object
     */
    core::Structure build_structure_from_residues(const std::string& pdb_id,
                                                  const std::map<ResidueKey, std::vector<core::Atom>>& residue_atoms,
                                                  const std::map<ResidueKey, int>& legacy_idx_map,
                                                  const std::vector<std::string>& chain_order) const;
};

} // namespace io
} // namespace x3dna

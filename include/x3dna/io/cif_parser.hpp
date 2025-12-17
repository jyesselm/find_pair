/**
 * @file cif_parser.hpp
 * @brief mmCIF/CIF file parser using GEMMI library
 *
 * This parser reads mmCIF files and creates Structure objects containing
 * atoms, residues, and chains. It handles coordinate data from the
 * _atom_site loop.
 */

#pragma once

#include <filesystem>
#include <string>
#include <stdexcept>
#include <map>
#include <tuple>
#include <vector>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>

// Forward declare GEMMI types to avoid header inclusion
namespace gemmi {
    struct Structure;
}

namespace x3dna {
namespace io {

/**
 * @class CifParser
 * @brief Parser for mmCIF (macromolecular Crystallographic Information File) format
 *
 * Parses mmCIF files using GEMMI and creates Structure objects. Supports:
 * - ATOM records (group_PDB = "ATOM")
 * - HETATM records (group_PDB = "HETATM", optional)
 * - Chain identification (auth_asym_id preferred, label_asym_id fallback)
 * - Residue numbering (auth_seq_id preferred)
 * - Alternate conformations (label_alt_id)
 * - Compressed files (.cif.gz)
 */
class CifParser {
public:
    /**
     * @brief Constructor
     */
    CifParser() = default;

    /**
     * @brief Destructor
     */
    ~CifParser() = default;

    // Copy and move semantics
    CifParser(const CifParser&) = default;
    CifParser& operator=(const CifParser&) = default;
    CifParser(CifParser&&) noexcept = default;
    CifParser& operator=(CifParser&&) noexcept = default;

    /**
     * @brief Parse CIF file from filesystem path
     * @param path Path to CIF file (.cif or .cif.gz)
     * @return Structure object containing parsed data
     * @throws ParseError if file cannot be read or parsed
     */
    core::Structure parse_file(const std::filesystem::path& path);

    /**
     * @brief Parse CIF file from string content
     * @param content String containing CIF file content
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
     * @brief Set whether to use auth_* fields (default) or label_* fields
     * @param use_auth True to use auth_asym_id/auth_seq_id (PDB-style), false for label_*
     *
     * auth_* fields correspond to original PDB identifiers and are preferred
     * for compatibility with legacy X3DNA output.
     */
    void set_use_auth_fields(bool use_auth) { use_auth_fields_ = use_auth; }

    /**
     * @brief Get whether auth_* fields are used
     * @return True if using auth_* fields
     */
    bool use_auth_fields() const { return use_auth_fields_; }

    /**
     * @brief Exception class for parsing errors
     */
    class ParseError : public std::runtime_error {
    public:
        /**
         * @brief Constructor
         * @param message Error message
         */
        explicit ParseError(const std::string& message);
    };

private:
    bool include_hetatm_ = false;     // Include HETATM records
    bool include_waters_ = false;     // Include water molecules
    bool use_auth_fields_ = true;     // Use auth_* fields for PDB compatibility

    /**
     * @brief Convert GEMMI Structure to our Structure
     * @param gemmi_struct GEMMI structure
     * @param pdb_id PDB identifier
     * @return Structure object
     */
    core::Structure convert_gemmi_structure(const gemmi::Structure& gemmi_struct,
                                            const std::string& pdb_id);

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
     * @brief Normalize atom name to PDB 4-character format
     * @param name Atom name from CIF
     * @return Normalized atom name
     */
    std::string normalize_atom_name(const std::string& name) const;

    /**
     * @brief Normalize residue name (trim whitespace)
     * @param name Residue name from CIF
     * @return Normalized residue name
     */
    std::string normalize_residue_name(const std::string& name) const;

    /**
     * @brief Check if atom should be kept based on filtering rules
     * @param is_hetatm True if HETATM record
     * @param alt_loc Alternate location indicator
     * @param residue_name Residue name
     * @return True if atom should be kept
     */
    bool should_keep_atom(bool is_hetatm, char alt_loc, const std::string& residue_name) const;

    /**
     * @brief Check alternate location filter
     * @param alt_loc Alternate location character
     * @return True if this alt_loc should be kept
     */
    bool check_alt_loc_filter(char alt_loc) const;

    /**
     * @brief Apply atom name formatting rules (PDB-style spacing)
     * @param name Atom name
     * @return Formatted atom name
     */
    std::string apply_atom_name_formatting_rules(const std::string& name) const;

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
     * @brief Build Structure from grouped residue atoms
     * @param pdb_id Structure identifier
     * @param residue_atoms Map of residue key to atoms
     * @return Structure object
     */
    core::Structure build_structure_from_residues(
        const std::string& pdb_id,
        const std::map<std::tuple<std::string, char, int, char>, std::vector<core::Atom>>& residue_atoms) const;
};

} // namespace io
} // namespace x3dna

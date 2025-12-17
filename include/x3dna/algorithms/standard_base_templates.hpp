/**
 * @file standard_base_templates.hpp
 * @brief Standard base template loader for frame calculation
 */

#pragma once

#include <string>
#include <filesystem>
#include <map>
#include <memory>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/io/pdb_parser.hpp>

namespace x3dna {
namespace algorithms {

/**
 * @class StandardBaseTemplates
 * @brief Loads and caches standard base PDB template files
 *
 * Standard templates are ideal base geometries used for frame calculation.
 * Template files are named: Atomic_A.pdb, Atomic_C.pdb, Atomic_G.pdb, etc.
 */
class StandardBaseTemplates {
public:
    /**
     * @brief Default constructor
     *
     * Uses X3DNA_HOMEDIR environment variable if set, otherwise looks in
     * common installation paths or current directory.
     */
    StandardBaseTemplates();

    /**
     * @brief Constructor with explicit template path
     * @param template_path Base path to template directory (e.g., "/path/to/x3dna/config")
     */
    explicit StandardBaseTemplates(const std::filesystem::path& template_path);

    /**
     * @brief Load standard base template for a residue type
     * @param type Residue type (ADENINE, CYTOSINE, GUANINE, THYMINE, URACIL)
     * @param is_modified If true, use lowercase template (Atomic.x.pdb) for modified nucleotides
     * @return Structure containing the standard base atoms
     * @throws std::runtime_error if template file cannot be found or loaded
     */
    core::Structure load_template(core::ResidueType type, bool is_modified);

    /**
     * @brief Load standard base template for a residue type (standard nucleotide)
     * @param type Residue type
     * @return Structure containing the standard base atoms
     */
    core::Structure load_template(core::ResidueType type);

    /**
     * @brief Get template file path for a residue type
     * @param type Residue type
     * @param is_modified If true, use lowercase template (Atomic.x.pdb) for modified nucleotides
     * @return Path to template file
     */
    std::filesystem::path get_template_path(core::ResidueType type, bool is_modified) const;

    /**
     * @brief Get template file path for a residue type (standard nucleotide)
     * @param type Residue type
     * @return Path to template file
     */
    std::filesystem::path get_template_path(core::ResidueType type) const;

    /**
     * @brief Set base template directory path
     * @param template_path Base path to template directory
     */
    void set_template_path(const std::filesystem::path& template_path);

    /**
     * @brief Get current template path
     * @return Current template base path
     */
    std::filesystem::path template_path() const {
        return template_path_;
    }

    /**
     * @brief Clear cached templates (force reload on next access)
     */
    void clear_cache();

    /**
     * @brief Check if template exists for a residue type
     * @param type Residue type
     * @return true if template file exists
     */
    bool template_exists(core::ResidueType type) const;

private:
    std::filesystem::path template_path_;
    std::map<core::ResidueType, std::shared_ptr<core::Structure>> cache_;
    io::PdbParser parser_;

    /**
     * @brief Convert ResidueType to template filename
     * @param type Residue type
     * @param is_modified If true, use lowercase template (Atomic.x.pdb) for modified nucleotides
     * @return Template filename (e.g., "Atomic_A.pdb" or "Atomic.a.pdb")
     */
    static std::string type_to_filename(core::ResidueType type, bool is_modified);

    /**
     * @brief Convert ResidueType to template filename (standard nucleotide)
     * @param type Residue type
     * @return Template filename (e.g., "Atomic_A.pdb")
     */
    static std::string type_to_filename(core::ResidueType type);

    /**
     * @brief Check if a file exists
     * @param path File path
     * @return true if file exists
     */
    static bool file_exists(const std::filesystem::path& path);
};

} // namespace algorithms
} // namespace x3dna

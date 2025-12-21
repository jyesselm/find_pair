#pragma once

#include <string>
#include <map>
#include <optional>

namespace x3dna {
namespace core {

/**
 * @brief Registry for PDB atom name to element symbol mapping
 *
 * Maps 4-character PDB atom name patterns to element symbols (H, C, N, O, S, P).
 * Pattern format uses '.' as wildcard for non-alphabetic characters.
 * Examples: ".H.A" -> "H", ".C.." -> "C", ".N.." -> "N"
 *
 * This replaces the duplicate atom_list loading code in hydrogen_bond_utils.cpp
 * and base_pair_validator.cpp with a centralized, data-driven approach.
 *
 * Data is loaded from resources/config/atomlist.json on first access.
 */
class AtomSymbolRegistry {
public:
    /**
     * @brief Get element symbol for an atom name
     * @param atom_name Atom name (trimmed or padded, will be normalized)
     * @return Element symbol (e.g., "H", "C", "N", "O", "S", "P") or "XX" if unknown
     */
    [[nodiscard]] static std::string get_symbol(const std::string& atom_name);

    /**
     * @brief Get element index for an atom name
     * @param atom_name Atom name (trimmed or padded, will be normalized)
     * @return Index: 1=C, 2=O, 3=H, 4=N, 5=S, 6=P, 0=unknown
     *
     * This matches the legacy asym_idx values used in hydrogen bond detection.
     */
    [[nodiscard]] static int get_atom_idx(const std::string& atom_name);

    /**
     * @brief Check if a pattern exists in the registry
     * @param pattern 4-character pattern with '.' wildcards
     * @return true if pattern is registered
     */
    [[nodiscard]] static bool contains_pattern(const std::string& pattern);

private:
    // Symbol to index mapping (matches legacy asym_idx)
    static const std::map<std::string, int> SYMBOL_TO_IDX;

    // Get the lazy-loaded pattern registry
    static const std::map<std::string, std::string>& get_patterns();

    // Convert atom name to pattern (replace non-alpha with '.')
    static std::string atom_name_to_pattern(const std::string& atom_name);

    // Pad atom name to 4 characters
    static std::string pad_atom_name(const std::string& atom_name);
};

} // namespace core
} // namespace x3dna

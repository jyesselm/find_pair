/**
 * @file hydrogen_bond_utils.hpp
 * @brief Utilities for hydrogen bond calculations
 */

#pragma once

#include <string>
#include <x3dna/core/atom_symbol_registry.hpp>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

/**
 * @brief Atom list utilities - delegates to AtomSymbolRegistry
 *
 * This class provides backward compatibility with existing code.
 * All functionality is now provided by core::AtomSymbolRegistry.
 */
class AtomListUtils {
public:
    /**
     * @brief Load atom list (no-op, registry is lazy-loaded)
     * @param x3dna_home Ignored - kept for API compatibility
     */
    static void load_atom_list(const std::string& x3dna_home = "") {
        (void)x3dna_home;  // Registry is lazy-loaded from atomlist.json
    }

    /**
     * @brief Get atom index from atom name
     * @param atom_name Atom name (e.g., "N1", " N1 ")
     * @return Atom type index (0=UNK, 1=C, 2=O, 3=H, 4=N, 5=S, 6=P)
     */
    [[nodiscard]] static int get_atom_idx(const std::string& atom_name) {
        return core::AtomSymbolRegistry::get_atom_idx(atom_name);
    }

    /**
     * @brief Check if atom list is loaded (always true - registry is lazy-loaded)
     */
    [[nodiscard]] static bool is_loaded() {
        return true;
    }
};

/**
 * @brief Check if atom is a base atom (matches legacy is_baseatom)
 * @param atom_name Atom name
 * @return True if base atom
 */
[[nodiscard]] bool is_base_atom(const std::string& atom_name);

/**
 * @brief Check if two atoms can form a hydrogen bond (matches legacy good_hbatoms)
 * @param atom1 First atom name
 * @param atom2 Second atom name
 * @param hb_atoms H-bond atom list (default ".O.N")
 * @param include_backbone_backbone If true, allow backbone-backbone bonds (default: false)
 * @return True if atoms can form H-bond
 */
[[nodiscard]] bool good_hb_atoms(const std::string& atom1, const std::string& atom2, const std::string& hb_atoms,
                                  bool include_backbone_backbone = false);

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna

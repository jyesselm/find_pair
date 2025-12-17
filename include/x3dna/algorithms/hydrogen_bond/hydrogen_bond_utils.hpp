/**
 * @file hydrogen_bond_utils.hpp
 * @brief Utilities for hydrogen bond calculations (atom list management, etc.)
 */

#pragma once

#include <string>
#include <map>
#include <vector>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

/**
 * @brief Atom list management utilities
 *
 * Provides shared atom list loading and lookup for hydrogen bond calculations
 */
class AtomListUtils {
public:
    /**
     * @brief Load atom list from file
     * @param x3dna_home X3DNA home directory (empty to use environment variable)
     */
    static void load_atom_list(const std::string& x3dna_home = "");

    /**
     * @brief Get atom index from atom name (matches legacy aname2asym + asym_idx)
     * @param atom_name Atom name (e.g., " N1 ")
     * @return Atom type index (0=UNK, 1=C, 2=O, 3=H, 4=N, 5=S, 6=P, ...)
     */
    [[nodiscard]] static int get_atom_idx(const std::string& atom_name);

    /**
     * @brief Check if atom list is loaded
     */
    [[nodiscard]] static bool is_loaded() {
        return atom_list_loaded_;
    }

private:
    static std::map<std::string, std::string> atom_list_;
    static bool atom_list_loaded_;
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
 * @return True if atoms can form H-bond
 */
[[nodiscard]] bool good_hb_atoms(const std::string& atom1, const std::string& atom2, const std::string& hb_atoms);

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna

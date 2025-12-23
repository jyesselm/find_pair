/**
 * @file hydrogen_bond_utils.cpp
 * @brief Implementation of hydrogen bond utilities
 *
 * Note: AtomListUtils is now fully inline in the header and delegates to
 * core::AtomSymbolRegistry. This file only contains the free functions.
 */

#include <x3dna/algorithms/hydrogen_bond/hydrogen_bond_utils.hpp>
#include <algorithm>
#include <cctype>
#include <vector>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

bool is_base_atom(const std::string& atom_name) {
    // Matches legacy is_baseatom (line 4652 in cmn_fncs.c)
    // Base atoms: C5M or atoms matching pattern " XD " where X is not H or P and D is digit
    // Input is already trimmed from Atom class

    if (atom_name == "C5M") {
        return true;
    }

    // For trimmed names like "N1", "C2", "N9", etc.
    // Pattern: exactly 2 chars - letter (not H or P) followed by digit
    // This excludes sugar atoms like "C5'" and backbone atoms like "O1P"
    if (atom_name.length() == 2 && atom_name[0] != 'H' && atom_name[0] != 'P' &&
        std::isdigit(static_cast<unsigned char>(atom_name[1]))) {
        return true;
    }

    return false;
}

bool good_hb_atoms(const std::string& atom1, const std::string& atom2, const std::string& hb_atoms,
                   bool include_backbone_backbone) {
    // Match legacy good_hbatoms() logic EXACTLY (lines 3864-3877 in cmn_fncs.c)
    // Input atom names are already trimmed from Atom class

    // Step 1: PO list check (matches legacy lines 3866-3870)
    // Skip this check if include_backbone_backbone is true (for DSSR-like detection)
    if (!include_backbone_backbone) {
        // Include both old (O1P/O2P) and new (OP1/OP2) atom naming conventions
        static const std::vector<std::string> PO = {"O1P", "O2P", "OP1", "OP2", "O3'", "O4'", "O5'", "N7"};
        bool atom1_in_po = std::find(PO.begin(), PO.end(), atom1) != PO.end();
        bool atom2_in_po = std::find(PO.begin(), PO.end(), atom2) != PO.end();
        if (atom1_in_po && atom2_in_po) {
            return false;
        }
    }

    // Step 2: idx-based check (matches legacy lines 3871-3874)
    // Build hb_idx array from hb_atoms string (format: ".O.N" means O and N)
    std::vector<int> hb_idx;
    std::string hb_atoms_upper = hb_atoms;
    std::transform(hb_atoms_upper.begin(), hb_atoms_upper.end(), hb_atoms_upper.begin(), ::toupper);

    // Element symbol to index mapping (matches legacy asym_idx)
    static const std::map<std::string, int> atoms_list_idx = {{"C", 1}, {"O", 2}, {"H", 3},
                                                              {"N", 4}, {"S", 5}, {"P", 6}};

    // Parse hb_atoms format: ".O.N" - extract element chars after each '.'
    for (size_t i = 0; i < hb_atoms_upper.length(); ++i) {
        if (hb_atoms_upper[i] == '.' && i + 1 < hb_atoms_upper.length()) {
            std::string sym(1, hb_atoms_upper[i + 1]);
            auto it = atoms_list_idx.find(sym);
            if (it != atoms_list_idx.end()) {
                hb_idx.push_back(it->second);
            }
        }
    }

    // Get atom indices via the registry (delegates to AtomSymbolRegistry)
    int idx1 = AtomListUtils::get_atom_idx(atom1);
    int idx2 = AtomListUtils::get_atom_idx(atom2);

    // Legacy check: (idx1 == 2 || idx1 == 4 || idx2 == 2 || idx2 == 4) AND both in hb_idx
    bool at_least_one_on = (idx1 == 2 || idx1 == 4 || idx2 == 2 || idx2 == 4);

    if (at_least_one_on) {
        bool idx1_in_hb = std::find(hb_idx.begin(), hb_idx.end(), idx1) != hb_idx.end();
        bool idx2_in_hb = std::find(hb_idx.begin(), hb_idx.end(), idx2) != hb_idx.end();
        if (idx1_in_hb && idx2_in_hb) {
            return true;
        }
    }

    return false;
}

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna

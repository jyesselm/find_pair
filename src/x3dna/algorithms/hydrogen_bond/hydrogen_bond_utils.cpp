/**
 * @file hydrogen_bond_utils.cpp
 * @brief Implementation of hydrogen bond utilities
 */

#include <x3dna/algorithms/hydrogen_bond/hydrogen_bond_utils.hpp>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstdlib>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

// Static members
std::map<std::string, std::string> AtomListUtils::atom_list_;
bool AtomListUtils::atom_list_loaded_ = false;

void AtomListUtils::load_atom_list(const std::string& x3dna_home) {
    if (atom_list_loaded_) {
        return; // Already loaded
    }

    // Priority 1: Check local resources/config directory (self-contained)
    std::vector<std::string> search_paths = {
        "resources/config/atomlist.dat",
        "../resources/config/atomlist.dat",
        "../../resources/config/atomlist.dat",
#ifdef X3DNA_SOURCE_DIR
        std::string(X3DNA_SOURCE_DIR) + "/resources/config/atomlist.dat",
#endif
    };

    std::string atomlist_file;
    for (const auto& path : search_paths) {
        std::ifstream test(path);
        if (test.good()) {
            atomlist_file = path;
            break;
        }
    }

    // Priority 2: Fall back to X3DNA environment variable
    if (atomlist_file.empty()) {
        std::string home_dir = x3dna_home;
        if (home_dir.empty()) {
            const char* env_home = std::getenv("X3DNA");
            if (env_home) {
                home_dir = env_home;
            }
        }
        if (!home_dir.empty()) {
            atomlist_file = home_dir + "/config/atomlist.dat";
        }
    }
    std::ifstream file(atomlist_file);
    if (!file.is_open()) {
        // If file not found, use fallback logic only
        atom_list_loaded_ = true;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Parse line: atom_name_pattern atomic_symbol
        std::istringstream iss(line);
        std::string aname4, asym;
        if (!(iss >> aname4 >> asym)) {
            continue;
        }

        // Skip if starts with #
        if (aname4[0] == '#' || asym[0] == '#') {
            continue;
        }

        // Validate format (matches legacy get_atomlist)
        if (aname4.length() != 4) {
            continue;
        }
        if (asym.length() != 1 && asym.length() != 2) {
            continue;
        }

        // Normalize: convert to uppercase and pad atomic symbol to 2 chars
        std::transform(aname4.begin(), aname4.end(), aname4.begin(), ::toupper);
        std::transform(asym.begin(), asym.end(), asym.begin(), ::toupper);

        // Pad atomic symbol to 2 chars with space if needed
        if (asym.length() == 1) {
            asym = " " + asym;
        }

        // Store in map: pattern -> atomic symbol
        atom_list_[aname4] = asym;
    }

    file.close();
    atom_list_loaded_ = true;
}

int AtomListUtils::get_atom_idx(const std::string& atom_name) {
    // Load atom list if not already loaded
    if (!atom_list_loaded_) {
        load_atom_list();
    }

    if (atom_name.length() < 4) {
        return 0; // UNK - need at least 4 chars for PDB format
    }

    // Match legacy aname2asym logic EXACTLY (lines 4160-4193 in cmn_fncs.c)
    // Step 1: Create pattern by replacing non-alphabetic with '.'
    std::string aname_pattern = atom_name;
    for (size_t i = 0; i < 4 && i < aname_pattern.length(); i++) {
        if (!std::isalpha(static_cast<unsigned char>(aname_pattern[i]))) {
            aname_pattern[i] = '.';
        }
    }

    // Step 2: Try to match against atomlist (matches legacy lines 4168-4170)
    std::string my_asym;
    bool found_match = false;
    for (const auto& [pattern, sym] : atom_list_) {
        // str_pmatch(atomlist[i], aname) checks if atomlist[i] starts with aname
        if (pattern.length() >= aname_pattern.length() && pattern.substr(0, aname_pattern.length()) == aname_pattern) {
            my_asym = sym;
            found_match = true;
            break;
        }
    }

    // Step 3: If no match, use fallback logic (matches legacy lines 4171-4184)
    if (!found_match) {
        bool unknown = (aname_pattern == ".UNK");

        // Fallback case 1: aname[0]!='.' && aname[1]!='.' && aname[2]=='.' && aname[3]=='.'
        if (aname_pattern.length() >= 4 && aname_pattern[0] != '.' && aname_pattern[1] != '.' &&
            aname_pattern[2] == '.' && aname_pattern[3] == '.') {
            my_asym = std::string(" ") + aname_pattern[1]; // Space + second char
        }
        // Fallback case 2: aname[0]=='.' && aname[1]!='.' && !unknown
        else if (aname_pattern.length() >= 2 && aname_pattern[0] == '.' && aname_pattern[1] != '.' && !unknown) {
            my_asym = std::string(" ") + aname_pattern[1];
        }
        // Fallback case 3: aname[0]=='H' -> " H"
        else if (aname_pattern.length() >= 1 && aname_pattern[0] == 'H') {
            my_asym = " H";
        }
        // No match - use UNK
        else {
            my_asym = "XX"; // UNKATM
        }
    }

    // Step 4: Map atomic symbol to idx (matches legacy asym_idx)
    static const std::map<std::string, int> atoms_list_idx = {{" C", 1}, {" O", 2}, {" H", 3},
                                                              {" N", 4}, {" S", 5}, {" P", 6}};
    auto it = atoms_list_idx.find(my_asym);
    if (it != atoms_list_idx.end()) {
        return it->second;
    }
    return 0;
}

bool is_base_atom(const std::string& atom_name) {
    // Matches legacy is_baseatom (line 4652 in cmn_fncs.c)
    // Base atoms: C5M or atoms matching pattern " HP" where H is not H or P
    if (atom_name == " C5M") {
        return true;
    }

    // Pattern: space, character (not H or P), digit, space
    if (atom_name.length() >= 4 && atom_name[0] == ' ' && atom_name[1] != 'H' && atom_name[1] != 'P' &&
        std::isdigit(static_cast<unsigned char>(atom_name[2])) && atom_name[3] == ' ') {
        return true;
    }

    return false;
}

bool good_hb_atoms(const std::string& atom1, const std::string& atom2, const std::string& hb_atoms) {
    // Match legacy good_hbatoms() logic EXACTLY (lines 3864-3877 in cmn_fncs.c)

    // Step 1: PO list check (matches legacy lines 3866-3870)
    static const std::vector<std::string> PO = {" O1P", " O2P", " O3'", " O4'", " O5'", " N7 "};
    bool atom1_in_po = std::find(PO.begin(), PO.end(), atom1) != PO.end();
    bool atom2_in_po = std::find(PO.begin(), PO.end(), atom2) != PO.end();
    if (atom1_in_po && atom2_in_po) {
        return false;
    }

    // Step 2: idx-based check (matches legacy lines 3871-3874)
    // Build hb_idx array from hb_atoms string
    std::vector<int> hb_idx;
    std::string hb_atoms_upper = hb_atoms;
    std::transform(hb_atoms_upper.begin(), hb_atoms_upper.end(), hb_atoms_upper.begin(), ::toupper);
    std::replace(hb_atoms_upper.begin(), hb_atoms_upper.end(), '.', ' ');

    static const std::map<std::string, int> atoms_list_idx = {{" C", 1}, {" O", 2}, {" H", 3},
                                                              {" N", 4}, {" S", 5}, {" P", 6}};

    for (size_t i = 0; i + 1 < hb_atoms_upper.length(); i += 2) {
        std::string sym = hb_atoms_upper.substr(i, 2);
        auto it = atoms_list_idx.find(sym);
        if (it != atoms_list_idx.end()) {
            hb_idx.push_back(it->second);
        }
    }

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

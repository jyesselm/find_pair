/**
 * @file atom_capacity.cpp
 * @brief Implementation of atom capacity lookup
 */

#include <x3dna/algorithms/hydrogen_bond/slot/atom_capacity.hpp>
#include <algorithm>
#include <cctype>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {
namespace slot {

// Key type for capacity maps: (base_type, atom_name) - internal to this file
struct CapacityKey {
    char base_type;
    std::string atom_name;

    bool operator==(const CapacityKey& other) const {
        return base_type == other.base_type && atom_name == other.atom_name;
    }
};

struct CapacityKeyHash {
    size_t operator()(const CapacityKey& key) const {
        return std::hash<char>()(key.base_type) ^
               (std::hash<std::string>()(key.atom_name) << 1);
    }
};

using CapacityMap = std::unordered_map<CapacityKey, int, CapacityKeyHash>;

namespace {

// Helper to create capacity key
CapacityKey make_key(char base, const std::string& atom) {
    return {base, atom};
}

// Standard bases
const std::string STANDARD_BASES = "AGCUTP I";

// Donor capacity map (static singleton)
const CapacityMap& get_donor_capacity_map() {
    static const CapacityMap map = {
        // NH2 amino groups - 2 hydrogens
        {make_key('A', "N6"), 2},
        {make_key('C', "N4"), 2},
        {make_key('G', "N2"), 2},

        // Imino NH - 1 hydrogen
        {make_key('G', "N1"), 1},
        {make_key('U', "N3"), 1},
        {make_key('T', "N3"), 1},

        // Ribose O2' hydroxyl - can donate 1
        {make_key('A', "O2'"), 1},
        {make_key('G', "O2'"), 1},
        {make_key('C', "O2'"), 1},
        {make_key('U', "O2'"), 1},
        {make_key('T', "O2'"), 1},
        {make_key('P', "O2'"), 1},
        {make_key('I', "O2'"), 1},

        // Ribose O3' hydroxyl - can donate 1
        {make_key('A', "O3'"), 1},
        {make_key('G', "O3'"), 1},
        {make_key('C', "O3'"), 1},
        {make_key('U', "O3'"), 1},
        {make_key('T', "O3'"), 1},
        {make_key('P', "O3'"), 1},
        {make_key('I', "O3'"), 1},

        // Ribose O5' hydroxyl - can donate 1
        {make_key('A', "O5'"), 1},
        {make_key('G', "O5'"), 1},
        {make_key('C', "O5'"), 1},
        {make_key('U', "O5'"), 1},
        {make_key('T', "O5'"), 1},
        {make_key('P', "O5'"), 1},
        {make_key('I', "O5'"), 1},

        // Pseudouridine (P)
        {make_key('P', "N1"), 1},
        {make_key('P', "N3"), 1},

        // Inosine (I)
        {make_key('I', "N1"), 1},
    };
    return map;
}

// Acceptor capacity map (static singleton)
const CapacityMap& get_acceptor_capacity_map() {
    static const CapacityMap map = {
        // sp2 carbonyl oxygens - 2 lone pairs
        {make_key('G', "O6"), 2},
        {make_key('U', "O2"), 2},
        {make_key('U', "O4"), 2},
        {make_key('C', "O2"), 2},
        {make_key('T', "O2"), 2},
        {make_key('T', "O4"), 2},

        // sp2 ring nitrogens - 1 lone pair
        {make_key('A', "N1"), 1},
        {make_key('A', "N3"), 1},
        {make_key('A', "N7"), 1},
        {make_key('G', "N3"), 1},
        {make_key('G', "N7"), 1},
        {make_key('C', "N3"), 1},

        // Ribose O2' - 2 lone pairs
        {make_key('A', "O2'"), 2},
        {make_key('G', "O2'"), 2},
        {make_key('C', "O2'"), 2},
        {make_key('U', "O2'"), 2},
        {make_key('T', "O2'"), 2},
        {make_key('P', "O2'"), 2},
        {make_key('I', "O2'"), 2},

        // Ribose O4' ring - 1 lone pair accessible
        {make_key('A', "O4'"), 1},
        {make_key('G', "O4'"), 1},
        {make_key('C', "O4'"), 1},
        {make_key('U', "O4'"), 1},
        {make_key('T', "O4'"), 1},
        {make_key('P', "O4'"), 1},
        {make_key('I', "O4'"), 1},

        // Ribose O3' - 2 lone pairs
        {make_key('A', "O3'"), 2},
        {make_key('G', "O3'"), 2},
        {make_key('C', "O3'"), 2},
        {make_key('U', "O3'"), 2},
        {make_key('T', "O3'"), 2},
        {make_key('P', "O3'"), 2},
        {make_key('I', "O3'"), 2},

        // Ribose O5' - 2 lone pairs
        {make_key('A', "O5'"), 2},
        {make_key('G', "O5'"), 2},
        {make_key('C', "O5'"), 2},
        {make_key('U', "O5'"), 2},
        {make_key('T', "O5'"), 2},
        {make_key('P', "O5'"), 2},
        {make_key('I', "O5'"), 2},

        // Phosphate oxygens - 3 lone pairs (both naming conventions)
        {make_key('A', "OP1"), 3}, {make_key('A', "O1P"), 3},
        {make_key('G', "OP1"), 3}, {make_key('G', "O1P"), 3},
        {make_key('C', "OP1"), 3}, {make_key('C', "O1P"), 3},
        {make_key('U', "OP1"), 3}, {make_key('U', "O1P"), 3},
        {make_key('T', "OP1"), 3}, {make_key('T', "O1P"), 3},
        {make_key('A', "OP2"), 3}, {make_key('A', "O2P"), 3},
        {make_key('G', "OP2"), 3}, {make_key('G', "O2P"), 3},
        {make_key('C', "OP2"), 3}, {make_key('C', "O2P"), 3},
        {make_key('U', "OP2"), 3}, {make_key('U', "O2P"), 3},
        {make_key('T', "OP2"), 3}, {make_key('T', "O2P"), 3},

        // Pseudouridine (P)
        {make_key('P', "O2"), 2},
        {make_key('P', "O4"), 2},
        {make_key('P', "OP1"), 3}, {make_key('P', "O1P"), 3},
        {make_key('P', "OP2"), 3}, {make_key('P', "O2P"), 3},

        // Inosine (I)
        {make_key('I', "O6"), 2},
        {make_key('I', "N3"), 1},
        {make_key('I', "N7"), 1},
        {make_key('I', "OP1"), 3}, {make_key('I', "O1P"), 3},
        {make_key('I', "OP2"), 3}, {make_key('I', "O2P"), 3},
    };
    return map;
}

} // anonymous namespace

std::string AtomCapacity::normalize_atom_name(const std::string& atom_name) {
    std::string result = atom_name;
    // Trim whitespace
    result.erase(0, result.find_first_not_of(" \t"));
    result.erase(result.find_last_not_of(" \t") + 1);
    return result;
}

bool AtomCapacity::is_backbone_atom(const std::string& atom_name) {
    std::string name = normalize_atom_name(atom_name);
    return name == "P" || name == "OP1" || name == "OP2" ||
           name == "O1P" || name == "O2P" || name == "O3'" || name == "O5'";
}

std::optional<char> AtomCapacity::get_parent_base_type(const std::string& residue_code) {
    if (residue_code.empty()) {
        return std::nullopt;
    }

    // Standard single-letter codes
    if (residue_code.length() == 1) {
        char c = std::toupper(residue_code[0]);
        if (STANDARD_BASES.find(c) != std::string::npos) {
            return c;
        }
    }

    // DNA variants
    if (residue_code == "DA") return 'A';
    if (residue_code == "DG") return 'G';
    if (residue_code == "DC") return 'C';
    if (residue_code == "DT") return 'T';

    // Common modified bases
    std::string upper = residue_code;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

    // Guanine derivatives
    if (upper.find("GUA") != std::string::npos) return 'G';
    if (upper[0] == 'G' && upper.length() <= 3) return 'G';

    // Adenine derivatives
    if (upper.find("ADE") != std::string::npos) return 'A';
    if (upper[0] == 'A' && upper.length() <= 3) return 'A';
    if (upper.find("1MA") != std::string::npos) return 'A';

    // Cytosine derivatives
    if (upper.find("CYT") != std::string::npos) return 'C';
    if (upper[0] == 'C' && upper.length() <= 3) return 'C';
    if (upper == "5MC" || upper == "OMC") return 'C';

    // Uracil derivatives
    if (upper.find("URA") != std::string::npos) return 'U';
    if (upper[0] == 'U' && upper.length() <= 3) return 'U';
    if (upper.find("PSU") != std::string::npos) return 'P';
    if (upper.find("H2U") != std::string::npos) return 'U';

    // Thymine derivatives
    if (upper.find("THY") != std::string::npos) return 'T';
    if (upper[0] == 'T' && upper.length() <= 3) return 'T';

    // Pattern-based matching for modified bases
    if (upper.find("GN") != std::string::npos || upper.find("RG") != std::string::npos) return 'G';
    if (upper.find("RU") != std::string::npos || upper.find("UR") != std::string::npos) return 'U';

    return std::nullopt;
}

int AtomCapacity::get_donor_capacity(const std::string& residue_code,
                                      const std::string& atom_name) {
    std::string norm_atom = normalize_atom_name(atom_name);
    const auto& map = get_donor_capacity_map();

    // Try direct lookup with single-letter code
    if (residue_code.length() == 1) {
        char base = std::toupper(residue_code[0]);
        CapacityKey key{base, norm_atom};
        auto it = map.find(key);
        if (it != map.end()) {
            return it->second;
        }
    }

    // Try parent base type fallback
    auto parent = get_parent_base_type(residue_code);
    if (parent) {
        CapacityKey key{*parent, norm_atom};
        auto it = map.find(key);
        if (it != map.end()) {
            return it->second;
        }
    }

    return 0;
}

int AtomCapacity::get_acceptor_capacity(const std::string& residue_code,
                                         const std::string& atom_name) {
    std::string norm_atom = normalize_atom_name(atom_name);
    const auto& map = get_acceptor_capacity_map();

    // Try direct lookup with single-letter code
    if (residue_code.length() == 1) {
        char base = std::toupper(residue_code[0]);
        CapacityKey key{base, norm_atom};
        auto it = map.find(key);
        if (it != map.end()) {
            return it->second;
        }
    }

    // Try parent base type fallback
    auto parent = get_parent_base_type(residue_code);
    if (parent) {
        CapacityKey key{*parent, norm_atom};
        auto it = map.find(key);
        if (it != map.end()) {
            return it->second;
        }
    }

    return 0;
}

} // namespace slot
} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna

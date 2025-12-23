/**
 * @file edge_classifier.cpp
 * @brief Implementation of Leontis-Westhof edge classification
 */

#include <x3dna/algorithms/hydrogen_bond/edge_classifier.hpp>
#include <unordered_map>
#include <unordered_set>
#include <cctype>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

namespace {

// Edge assignments for each base type
// Key: atom name, Value: edge
// Based on Leontis-Westhof 2001 definitions

const std::unordered_map<std::string, core::BaseEdge> ADENINE_EDGES = {
    // Watson edge - canonical WC face
    {"N1", core::BaseEdge::WATSON},
    {"C2", core::BaseEdge::WATSON},
    {"N6", core::BaseEdge::WATSON}, // Primary W, but also participates in H
    // Hoogsteen edge - major groove face
    {"N7", core::BaseEdge::HOOGSTEEN},
    {"C8", core::BaseEdge::HOOGSTEEN},
    // Sugar edge - minor groove / sugar face
    {"N3", core::BaseEdge::SUGAR},
    {"C4", core::BaseEdge::SUGAR},
    {"O2'", core::BaseEdge::SUGAR},
};

const std::unordered_map<std::string, core::BaseEdge> GUANINE_EDGES = {
    // Watson edge
    {"N1", core::BaseEdge::WATSON},
    {"C2", core::BaseEdge::WATSON},
    {"O6", core::BaseEdge::WATSON}, // Primary W, participates in tWH too
    // Hoogsteen edge
    {"N7", core::BaseEdge::HOOGSTEEN},
    {"C8", core::BaseEdge::HOOGSTEEN},
    // Sugar edge
    {"N2", core::BaseEdge::SUGAR}, // Exocyclic amino - sugar edge
    {"N3", core::BaseEdge::SUGAR},
    {"C4", core::BaseEdge::SUGAR},
    {"O2'", core::BaseEdge::SUGAR},
};

const std::unordered_map<std::string, core::BaseEdge> CYTOSINE_EDGES = {
    // Watson edge
    {"N3", core::BaseEdge::WATSON},
    {"C4", core::BaseEdge::WATSON},
    {"N4", core::BaseEdge::WATSON}, // Exocyclic amino
    // Hoogsteen edge (less common for pyrimidines)
    {"C5", core::BaseEdge::HOOGSTEEN},
    {"C6", core::BaseEdge::HOOGSTEEN},
    // Sugar edge
    {"O2", core::BaseEdge::SUGAR},
    {"N1", core::BaseEdge::SUGAR}, // Glycosidic N
    {"O2'", core::BaseEdge::SUGAR},
};

const std::unordered_map<std::string, core::BaseEdge> URACIL_EDGES = {
    // Watson edge
    {"N3", core::BaseEdge::WATSON},
    {"C4", core::BaseEdge::WATSON},
    {"O4", core::BaseEdge::WATSON},
    // Hoogsteen edge
    {"C5", core::BaseEdge::HOOGSTEEN},
    {"C6", core::BaseEdge::HOOGSTEEN},
    // Sugar edge
    {"O2", core::BaseEdge::SUGAR},
    {"N1", core::BaseEdge::SUGAR},
    {"O2'", core::BaseEdge::SUGAR},
};

const std::unordered_map<std::string, core::BaseEdge> THYMINE_EDGES = {
    // Same as uracil
    {"N3", core::BaseEdge::WATSON},
    {"C4", core::BaseEdge::WATSON},
    {"O4", core::BaseEdge::WATSON},
    {"C5", core::BaseEdge::HOOGSTEEN},
    {"C6", core::BaseEdge::HOOGSTEEN},
    {"O2", core::BaseEdge::SUGAR},
    {"N1", core::BaseEdge::SUGAR},
    // DNA has no O2'
};

// Base atoms (excludes backbone and most sugar atoms)
const std::unordered_set<std::string> BASE_ATOMS = {
    "N1", "C2", "N2", "N3", "C4", "N4", "C5", "C6", "N6", "O6",
    "N7", "C8", "N9", "O2", "O4", "C5M", // C5M is thymine methyl
    "O2'" // Include O2' as it participates in sugar edge H-bonds
};

// Modified base parent types
const std::unordered_map<std::string, char> MODIFIED_BASE_PARENTS = {
    // Modified adenines
    {"1MA", 'A'}, {"MIA", 'A'}, {"I", 'A'}, {"M2A", 'A'},
    // Modified guanines
    {"2MG", 'G'}, {"7MG", 'G'}, {"M2G", 'G'}, {"OMG", 'G'}, {"YG", 'G'},
    // Modified cytosines
    {"5MC", 'C'}, {"OMC", 'C'}, {"4AC", 'C'},
    // Modified uracils
    {"PSU", 'U'}, {"H2U", 'U'}, {"5MU", 'U'}, {"4SU", 'U'}, {"DHU", 'U'},
    // Modified thymines
    {"5HT", 'T'},
    // DNA bases
    {"DA", 'A'}, {"DC", 'C'}, {"DG", 'G'}, {"DT", 'T'},
};

} // anonymous namespace

core::BaseEdge EdgeClassifier::classify(const std::string& atom_name, char base_type) {
    const std::unordered_map<std::string, core::BaseEdge>* edge_map = nullptr;

    switch (std::toupper(base_type)) {
        case 'A':
            edge_map = &ADENINE_EDGES;
            break;
        case 'G':
            edge_map = &GUANINE_EDGES;
            break;
        case 'C':
            edge_map = &CYTOSINE_EDGES;
            break;
        case 'U':
            edge_map = &URACIL_EDGES;
            break;
        case 'T':
            edge_map = &THYMINE_EDGES;
            break;
        default:
            return core::BaseEdge::UNKNOWN;
    }

    auto it = edge_map->find(atom_name);
    if (it != edge_map->end()) {
        return it->second;
    }

    return core::BaseEdge::UNKNOWN;
}

core::BaseEdge EdgeClassifier::classify_from_residue(const std::string& atom_name,
                                                      const std::string& residue_name) {
    char base_type = get_base_type(residue_name);
    if (base_type == '?') {
        return core::BaseEdge::UNKNOWN;
    }
    return classify(atom_name, base_type);
}

std::vector<std::string> EdgeClassifier::atoms_on_edge(char base_type, core::BaseEdge edge) {
    std::vector<std::string> result;
    const std::unordered_map<std::string, core::BaseEdge>* edge_map = nullptr;

    switch (std::toupper(base_type)) {
        case 'A':
            edge_map = &ADENINE_EDGES;
            break;
        case 'G':
            edge_map = &GUANINE_EDGES;
            break;
        case 'C':
            edge_map = &CYTOSINE_EDGES;
            break;
        case 'U':
            edge_map = &URACIL_EDGES;
            break;
        case 'T':
            edge_map = &THYMINE_EDGES;
            break;
        default:
            return result;
    }

    for (const auto& [atom, atom_edge] : *edge_map) {
        if (atom_edge == edge) {
            result.push_back(atom);
        }
    }

    return result;
}

bool EdgeClassifier::is_base_atom(const std::string& atom_name) {
    return BASE_ATOMS.count(atom_name) > 0;
}

char EdgeClassifier::get_base_type(const std::string& residue_name) {
    // Single letter names
    if (residue_name.size() == 1) {
        char c = std::toupper(residue_name[0]);
        if (c == 'A' || c == 'C' || c == 'G' || c == 'U' || c == 'T') {
            return c;
        }
    }

    // Check modified base lookup
    auto it = MODIFIED_BASE_PARENTS.find(residue_name);
    if (it != MODIFIED_BASE_PARENTS.end()) {
        return it->second;
    }

    // Try uppercase version
    std::string upper = residue_name;
    for (char& c : upper) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    it = MODIFIED_BASE_PARENTS.find(upper);
    if (it != MODIFIED_BASE_PARENTS.end()) {
        return it->second;
    }

    return '?';
}

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna

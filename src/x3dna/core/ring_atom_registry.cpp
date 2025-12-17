/**
 * @file ring_atom_registry.cpp
 * @brief Implementation of RingAtomRegistry
 */

#include <x3dna/core/ring_atom_registry.hpp>
#include <algorithm>
#include <cctype>

namespace x3dna {
namespace core {

// Purine ring atoms: fused 6+5 ring system (A, G, I)
const std::vector<std::string> RingAtomRegistry::PURINE_RING_ATOMS = {
    "N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"
};

// Pyrimidine ring atoms: single 6-membered ring (C, U, T, P)
const std::vector<std::string> RingAtomRegistry::PYRIMIDINE_RING_ATOMS = {
    "N1", "C2", "N3", "C4", "C5", "C6"
};

namespace {

/**
 * @brief Trim leading and trailing whitespace from a string_view
 */
[[nodiscard]] std::string_view trim(std::string_view sv) {
    while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.front()))) {
        sv.remove_prefix(1);
    }
    while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.back()))) {
        sv.remove_suffix(1);
    }
    return sv;
}

} // anonymous namespace

bool RingAtomRegistry::is_ring_atom(std::string_view atom_name) {
    const auto trimmed = trim(atom_name);

    // Check against purine atoms (superset includes all pyrimidine atoms)
    for (const auto& ring_atom : PURINE_RING_ATOMS) {
        if (trimmed == ring_atom) {
            return true;
        }
    }
    return false;
}

const std::vector<std::string>& RingAtomRegistry::purine_atoms() {
    return PURINE_RING_ATOMS;
}

const std::vector<std::string>& RingAtomRegistry::pyrimidine_atoms() {
    return PYRIMIDINE_RING_ATOMS;
}

const std::vector<std::string>& RingAtomRegistry::atoms_for_type(ResidueType type) {
    if (is_purine(type)) {
        return PURINE_RING_ATOMS;
    }
    return PYRIMIDINE_RING_ATOMS;
}

bool RingAtomRegistry::is_purine(ResidueType type) {
    return type == ResidueType::ADENINE ||
           type == ResidueType::GUANINE ||
           type == ResidueType::INOSINE;
}

} // namespace core
} // namespace x3dna

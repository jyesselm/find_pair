/**
 * @file hydrogen_bond.hpp
 * @brief HydrogenBond struct representing a hydrogen bond in a base pair
 */

#pragma once

#include <string>
#include <optional>

namespace x3dna {
namespace core {

/**
 * @brief Represents a hydrogen bond between two atoms in a base pair
 *
 * Hydrogen bonds are detected between donor and acceptor atoms.
 * The type field indicates whether the bond is standard ('-') or non-standard (' ').
 */
struct HydrogenBond {
    std::string donor_atom;              // Donor atom name (e.g., " N6 ")
    std::string acceptor_atom;           // Acceptor atom name (e.g., " O4 ")
    double distance = 0.0;               // Bond distance in Angstroms
    char type = ' ';                     // '-' for standard, ' ' for non-standard
    std::optional<size_t> hbond_idx;     // Optional index for tracking

    /**
     * @brief Check if this is a standard hydrogen bond
     * @return true if type is '-'
     */
    [[nodiscard]] bool is_standard() const {
        return type == '-';
    }

    /**
     * @brief Default constructor
     */
    HydrogenBond() = default;

    /**
     * @brief Constructor with all fields
     */
    HydrogenBond(std::string donor, std::string acceptor, double dist, char t)
        : donor_atom(std::move(donor))
        , acceptor_atom(std::move(acceptor))
        , distance(dist)
        , type(t) {}
};

} // namespace core
} // namespace x3dna

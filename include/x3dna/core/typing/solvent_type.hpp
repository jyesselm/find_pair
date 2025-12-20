/**
 * @file solvent_type.hpp
 * @brief Solvent and ion type classifications
 */

#pragma once

namespace x3dna {
namespace core {
namespace typing {

/**
 * @enum SolventType
 * @brief Classification of solvent molecules
 */
enum class SolventType {
    UNKNOWN,
    WATER,     ///< Water molecules (HOH, WAT, H2O, OH2, SOL)
    ORGANIC    ///< Organic solvents (DMSO, ethanol, etc.)
};

/**
 * @enum IonType
 * @brief Classification of ion species
 */
enum class IonType {
    UNKNOWN,
    // Alkali metals
    LITHIUM,    ///< Li+
    SODIUM,     ///< Na+
    POTASSIUM,  ///< K+
    RUBIDIUM,   ///< Rb+
    CESIUM,     ///< Cs+
    // Alkaline earth metals
    MAGNESIUM,  ///< Mg2+
    CALCIUM,    ///< Ca2+
    STRONTIUM,  ///< Sr2+
    BARIUM,     ///< Ba2+
    // Transition metals
    MANGANESE,  ///< Mn2+
    IRON,       ///< Fe2+/Fe3+
    COBALT,     ///< Co2+
    NICKEL,     ///< Ni2+
    COPPER,     ///< Cu+/Cu2+
    ZINC,       ///< Zn2+
    CADMIUM,    ///< Cd2+
    // Halogens
    FLUORIDE,   ///< F-
    CHLORIDE,   ///< Cl-
    BROMIDE,    ///< Br-
    IODIDE      ///< I-
};

/**
 * @brief Convert SolventType to string
 */
[[nodiscard]] inline const char* to_string(SolventType type) {
    switch (type) {
        case SolventType::UNKNOWN: return "UNKNOWN";
        case SolventType::WATER: return "WATER";
        case SolventType::ORGANIC: return "ORGANIC";
    }
    return "UNKNOWN";
}

/**
 * @brief Convert IonType to string
 */
[[nodiscard]] inline const char* to_string(IonType type) {
    switch (type) {
        case IonType::UNKNOWN: return "UNKNOWN";
        case IonType::LITHIUM: return "LITHIUM";
        case IonType::SODIUM: return "SODIUM";
        case IonType::POTASSIUM: return "POTASSIUM";
        case IonType::RUBIDIUM: return "RUBIDIUM";
        case IonType::CESIUM: return "CESIUM";
        case IonType::MAGNESIUM: return "MAGNESIUM";
        case IonType::CALCIUM: return "CALCIUM";
        case IonType::STRONTIUM: return "STRONTIUM";
        case IonType::BARIUM: return "BARIUM";
        case IonType::MANGANESE: return "MANGANESE";
        case IonType::IRON: return "IRON";
        case IonType::COBALT: return "COBALT";
        case IonType::NICKEL: return "NICKEL";
        case IonType::COPPER: return "COPPER";
        case IonType::ZINC: return "ZINC";
        case IonType::CADMIUM: return "CADMIUM";
        case IonType::FLUORIDE: return "FLUORIDE";
        case IonType::CHLORIDE: return "CHLORIDE";
        case IonType::BROMIDE: return "BROMIDE";
        case IonType::IODIDE: return "IODIDE";
    }
    return "UNKNOWN";
}

/**
 * @brief Check if an ion type is a cation (positive)
 */
[[nodiscard]] inline bool is_cation(IonType type) {
    switch (type) {
        case IonType::LITHIUM:
        case IonType::SODIUM:
        case IonType::POTASSIUM:
        case IonType::RUBIDIUM:
        case IonType::CESIUM:
        case IonType::MAGNESIUM:
        case IonType::CALCIUM:
        case IonType::STRONTIUM:
        case IonType::BARIUM:
        case IonType::MANGANESE:
        case IonType::IRON:
        case IonType::COBALT:
        case IonType::NICKEL:
        case IonType::COPPER:
        case IonType::ZINC:
        case IonType::CADMIUM:
            return true;
        default:
            return false;
    }
}

/**
 * @brief Check if an ion type is an anion (negative)
 */
[[nodiscard]] inline bool is_anion(IonType type) {
    switch (type) {
        case IonType::FLUORIDE:
        case IonType::CHLORIDE:
        case IonType::BROMIDE:
        case IonType::IODIDE:
            return true;
        default:
            return false;
    }
}

} // namespace typing

// Expose in core namespace for backwards compatibility
using typing::SolventType;
using typing::IonType;

} // namespace core
} // namespace x3dna

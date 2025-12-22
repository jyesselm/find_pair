/**
 * @file molecule_type.hpp
 * @brief Hierarchical molecule and residue type enums
 */

#pragma once

namespace x3dna {
namespace core {

/**
 * @enum MoleculeType
 * @brief Top-level classification of molecular entities
 */
enum class MoleculeType {
    UNKNOWN,
    NUCLEIC_ACID,  ///< RNA or DNA nucleotides
    PROTEIN,       ///< Amino acids
    LIPID,         ///< Lipid molecules
    WATER,         ///< Water molecules (HOH, WAT)
    ION,           ///< Metal ions and small charged species
    LIGAND         ///< Other small molecules, drugs, cofactors
};

/**
 * @enum NucleicAcidType
 * @brief Classification of nucleic acid type
 */
enum class NucleicAcidType {
    UNKNOWN,
    RNA,   ///< Ribonucleic acid (has 2'-OH)
    DNA    ///< Deoxyribonucleic acid (no 2'-OH)
};

/**
 * @enum BaseType
 * @brief Nucleobase identity (canonical or modified maps to this)
 */
enum class BaseType {
    UNKNOWN,
    ADENINE,
    GUANINE,
    CYTOSINE,
    THYMINE,
    URACIL,
    INOSINE,       ///< Hypoxanthine base
    PSEUDOURIDINE  ///< Isomer of uridine
};

/**
 * @enum BaseCategory
 * @brief Purine vs pyrimidine classification
 */
enum class BaseCategory {
    UNKNOWN,
    PURINE,      ///< Two-ring bases: A, G, I
    PYRIMIDINE   ///< Single-ring bases: C, T, U, PSU
};

/**
 * @brief Convert MoleculeType to string
 */
[[nodiscard]] inline const char* to_string(MoleculeType type) {
    switch (type) {
        case MoleculeType::UNKNOWN: return "UNKNOWN";
        case MoleculeType::NUCLEIC_ACID: return "NUCLEIC_ACID";
        case MoleculeType::PROTEIN: return "PROTEIN";
        case MoleculeType::LIPID: return "LIPID";
        case MoleculeType::WATER: return "WATER";
        case MoleculeType::ION: return "ION";
        case MoleculeType::LIGAND: return "LIGAND";
    }
    return "UNKNOWN";
}

/**
 * @brief Convert NucleicAcidType to string
 */
[[nodiscard]] inline const char* to_string(NucleicAcidType type) {
    switch (type) {
        case NucleicAcidType::UNKNOWN: return "UNKNOWN";
        case NucleicAcidType::RNA: return "RNA";
        case NucleicAcidType::DNA: return "DNA";
    }
    return "UNKNOWN";
}

/**
 * @brief Convert BaseType to string
 */
[[nodiscard]] inline const char* to_string(BaseType type) {
    switch (type) {
        case BaseType::UNKNOWN: return "UNKNOWN";
        case BaseType::ADENINE: return "ADENINE";
        case BaseType::GUANINE: return "GUANINE";
        case BaseType::CYTOSINE: return "CYTOSINE";
        case BaseType::THYMINE: return "THYMINE";
        case BaseType::URACIL: return "URACIL";
        case BaseType::INOSINE: return "INOSINE";
        case BaseType::PSEUDOURIDINE: return "PSEUDOURIDINE";
    }
    return "UNKNOWN";
}

/**
 * @brief Convert BaseCategory to string
 */
[[nodiscard]] inline const char* to_string(BaseCategory cat) {
    switch (cat) {
        case BaseCategory::UNKNOWN: return "UNKNOWN";
        case BaseCategory::PURINE: return "PURINE";
        case BaseCategory::PYRIMIDINE: return "PYRIMIDINE";
    }
    return "UNKNOWN";
}

/**
 * @brief Get BaseCategory for a BaseType
 */
[[nodiscard]] inline BaseCategory get_base_category(BaseType type) {
    switch (type) {
        case BaseType::ADENINE:
        case BaseType::GUANINE:
        case BaseType::INOSINE:
            return BaseCategory::PURINE;
        case BaseType::CYTOSINE:
        case BaseType::THYMINE:
        case BaseType::URACIL:
        case BaseType::PSEUDOURIDINE:
            return BaseCategory::PYRIMIDINE;
        default:
            return BaseCategory::UNKNOWN;
    }
}

} // namespace core
} // namespace x3dna

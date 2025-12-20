/**
 * @file protein_type.hpp
 * @brief Protein/amino acid type classifications
 */

#pragma once

#include <string>

namespace x3dna {
namespace core {
namespace typing {

/**
 * @enum AminoAcidType
 * @brief Classification of amino acid identity
 */
enum class AminoAcidType {
    UNKNOWN,
    // Standard 20 amino acids
    ALA,  ///< Alanine
    ARG,  ///< Arginine
    ASN,  ///< Asparagine
    ASP,  ///< Aspartic acid
    CYS,  ///< Cysteine
    GLN,  ///< Glutamine
    GLU,  ///< Glutamic acid
    GLY,  ///< Glycine
    HIS,  ///< Histidine
    ILE,  ///< Isoleucine
    LEU,  ///< Leucine
    LYS,  ///< Lysine
    MET,  ///< Methionine
    PHE,  ///< Phenylalanine
    PRO,  ///< Proline
    SER,  ///< Serine
    THR,  ///< Threonine
    TRP,  ///< Tryptophan
    TYR,  ///< Tyrosine
    VAL,  ///< Valine
    // Non-standard amino acids
    SEC,  ///< Selenocysteine (21st amino acid)
    PYL,  ///< Pyrrolysine (22nd amino acid)
    // Ambiguous codes
    ASX,  ///< Asparagine or Aspartic acid
    GLX,  ///< Glutamine or Glutamic acid
    XLE,  ///< Leucine or Isoleucine
    UNK   ///< Unknown amino acid
};

/**
 * @enum AminoAcidCategory
 * @brief Chemical property classification of amino acids
 */
enum class AminoAcidCategory {
    UNKNOWN,
    HYDROPHOBIC,  ///< A, V, L, I, M, F, W, P, G
    POLAR,        ///< S, T, C, Y, N, Q
    POSITIVE,     ///< K, R, H (basic)
    NEGATIVE      ///< D, E (acidic)
};

/**
 * @brief Convert AminoAcidType to string
 */
[[nodiscard]] inline const char* to_string(AminoAcidType type) {
    switch (type) {
        case AminoAcidType::UNKNOWN: return "UNKNOWN";
        case AminoAcidType::ALA: return "ALA";
        case AminoAcidType::ARG: return "ARG";
        case AminoAcidType::ASN: return "ASN";
        case AminoAcidType::ASP: return "ASP";
        case AminoAcidType::CYS: return "CYS";
        case AminoAcidType::GLN: return "GLN";
        case AminoAcidType::GLU: return "GLU";
        case AminoAcidType::GLY: return "GLY";
        case AminoAcidType::HIS: return "HIS";
        case AminoAcidType::ILE: return "ILE";
        case AminoAcidType::LEU: return "LEU";
        case AminoAcidType::LYS: return "LYS";
        case AminoAcidType::MET: return "MET";
        case AminoAcidType::PHE: return "PHE";
        case AminoAcidType::PRO: return "PRO";
        case AminoAcidType::SER: return "SER";
        case AminoAcidType::THR: return "THR";
        case AminoAcidType::TRP: return "TRP";
        case AminoAcidType::TYR: return "TYR";
        case AminoAcidType::VAL: return "VAL";
        case AminoAcidType::SEC: return "SEC";
        case AminoAcidType::PYL: return "PYL";
        case AminoAcidType::ASX: return "ASX";
        case AminoAcidType::GLX: return "GLX";
        case AminoAcidType::XLE: return "XLE";
        case AminoAcidType::UNK: return "UNK";
    }
    return "UNKNOWN";
}

/**
 * @brief Convert AminoAcidCategory to string
 */
[[nodiscard]] inline const char* to_string(AminoAcidCategory cat) {
    switch (cat) {
        case AminoAcidCategory::UNKNOWN: return "UNKNOWN";
        case AminoAcidCategory::HYDROPHOBIC: return "HYDROPHOBIC";
        case AminoAcidCategory::POLAR: return "POLAR";
        case AminoAcidCategory::POSITIVE: return "POSITIVE";
        case AminoAcidCategory::NEGATIVE: return "NEGATIVE";
    }
    return "UNKNOWN";
}

/**
 * @brief Get AminoAcidCategory for an AminoAcidType
 */
[[nodiscard]] inline AminoAcidCategory get_amino_acid_category(AminoAcidType type) {
    switch (type) {
        // Hydrophobic
        case AminoAcidType::ALA:
        case AminoAcidType::VAL:
        case AminoAcidType::LEU:
        case AminoAcidType::ILE:
        case AminoAcidType::MET:
        case AminoAcidType::PHE:
        case AminoAcidType::TRP:
        case AminoAcidType::PRO:
        case AminoAcidType::GLY:
            return AminoAcidCategory::HYDROPHOBIC;
        // Polar
        case AminoAcidType::SER:
        case AminoAcidType::THR:
        case AminoAcidType::CYS:
        case AminoAcidType::TYR:
        case AminoAcidType::ASN:
        case AminoAcidType::GLN:
        case AminoAcidType::SEC:
            return AminoAcidCategory::POLAR;
        // Positive (basic)
        case AminoAcidType::LYS:
        case AminoAcidType::ARG:
        case AminoAcidType::HIS:
        case AminoAcidType::PYL:
            return AminoAcidCategory::POSITIVE;
        // Negative (acidic)
        case AminoAcidType::ASP:
        case AminoAcidType::GLU:
            return AminoAcidCategory::NEGATIVE;
        default:
            return AminoAcidCategory::UNKNOWN;
    }
}

/**
 * @brief Get one-letter code for an AminoAcidType
 */
[[nodiscard]] inline char get_one_letter_code(AminoAcidType type) {
    switch (type) {
        case AminoAcidType::ALA: return 'A';
        case AminoAcidType::ARG: return 'R';
        case AminoAcidType::ASN: return 'N';
        case AminoAcidType::ASP: return 'D';
        case AminoAcidType::CYS: return 'C';
        case AminoAcidType::GLN: return 'Q';
        case AminoAcidType::GLU: return 'E';
        case AminoAcidType::GLY: return 'G';
        case AminoAcidType::HIS: return 'H';
        case AminoAcidType::ILE: return 'I';
        case AminoAcidType::LEU: return 'L';
        case AminoAcidType::LYS: return 'K';
        case AminoAcidType::MET: return 'M';
        case AminoAcidType::PHE: return 'F';
        case AminoAcidType::PRO: return 'P';
        case AminoAcidType::SER: return 'S';
        case AminoAcidType::THR: return 'T';
        case AminoAcidType::TRP: return 'W';
        case AminoAcidType::TYR: return 'Y';
        case AminoAcidType::VAL: return 'V';
        case AminoAcidType::SEC: return 'U';
        case AminoAcidType::PYL: return 'O';
        case AminoAcidType::ASX: return 'B';
        case AminoAcidType::GLX: return 'Z';
        case AminoAcidType::XLE: return 'J';
        case AminoAcidType::UNK: return 'X';
        default: return '?';
    }
}

/**
 * @brief Check if an AminoAcidType is one of the standard 20
 */
[[nodiscard]] inline bool is_standard_amino_acid(AminoAcidType type) {
    switch (type) {
        case AminoAcidType::ALA:
        case AminoAcidType::ARG:
        case AminoAcidType::ASN:
        case AminoAcidType::ASP:
        case AminoAcidType::CYS:
        case AminoAcidType::GLN:
        case AminoAcidType::GLU:
        case AminoAcidType::GLY:
        case AminoAcidType::HIS:
        case AminoAcidType::ILE:
        case AminoAcidType::LEU:
        case AminoAcidType::LYS:
        case AminoAcidType::MET:
        case AminoAcidType::PHE:
        case AminoAcidType::PRO:
        case AminoAcidType::SER:
        case AminoAcidType::THR:
        case AminoAcidType::TRP:
        case AminoAcidType::TYR:
        case AminoAcidType::VAL:
            return true;
        default:
            return false;
    }
}

} // namespace typing

// Expose in core namespace for backwards compatibility
using typing::AminoAcidType;
using typing::AminoAcidCategory;
using typing::get_amino_acid_category;

} // namespace core
} // namespace x3dna

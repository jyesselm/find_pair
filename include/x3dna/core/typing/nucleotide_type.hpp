/**
 * @file nucleotide_type.hpp
 * @brief Nucleotide-specific type classifications
 */

#pragma once

namespace x3dna {
namespace core {
namespace typing {

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

/**
 * @brief Check if a BaseType is a purine
 */
[[nodiscard]] inline bool is_purine(BaseType type) {
    return type == BaseType::ADENINE ||
           type == BaseType::GUANINE ||
           type == BaseType::INOSINE;
}

/**
 * @brief Check if a BaseType is a pyrimidine
 */
[[nodiscard]] inline bool is_pyrimidine(BaseType type) {
    return type == BaseType::CYTOSINE ||
           type == BaseType::THYMINE ||
           type == BaseType::URACIL ||
           type == BaseType::PSEUDOURIDINE;
}

} // namespace typing

// Expose in core namespace for backwards compatibility
using typing::NucleicAcidType;
using typing::BaseType;
using typing::BaseCategory;
using typing::get_base_category;

} // namespace core
} // namespace x3dna

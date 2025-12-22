/**
 * @file residue_type.hpp
 * @brief Enum for residue type classification
 */

#pragma once

namespace x3dna {
namespace core {

/**
 * @brief Classification of residue types
 */
enum class ResidueType {
    UNKNOWN = -2,
    AMINO_ACID = -1,
    NUCLEOTIDE = 0,
    ADENINE = 1,
    CYTOSINE = 2,
    GUANINE = 3,
    THYMINE = 4,
    URACIL = 5,
    NONCANONICAL_RNA = 6, // Modified nucleotides with ring atoms
    WATER = 7,            // Water molecules (HOH, WAT)
    ION = 8,              // Ions (MG, NA, CL, etc.)
    LIGAND = 9,           // Other small molecules/ligands
    PSEUDOURIDINE = 10,   // Pseudouridine (PSU) - C1' bonded to C5 instead of N1
    INOSINE = 11          // Inosine (I)
};

} // namespace core
} // namespace x3dna

/**
 * @file typing.hpp
 * @brief Convenience header for the molecular typing system
 *
 * This header provides unified access to all type classification
 * utilities for molecular entities (nucleotides, amino acids, atoms, etc.)
 */

#pragma once

// Core type enumerations
#include <x3dna/core/typing/molecule_type.hpp>
#include <x3dna/core/typing/nucleotide_type.hpp>
#include <x3dna/core/typing/protein_type.hpp>
#include <x3dna/core/typing/solvent_type.hpp>
#include <x3dna/core/typing/atom_type.hpp>

// Classification utilities
#include <x3dna/core/typing/atom_classification.hpp>
#include <x3dna/core/typing/residue_classification.hpp>

// Unified registry
#include <x3dna/core/typing/type_registry.hpp>

namespace x3dna {
namespace core {

// Re-export typing namespace contents for convenience
using namespace typing;

} // namespace core
} // namespace x3dna

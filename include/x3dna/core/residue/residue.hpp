/**
 * @file residue.hpp
 * @brief Convenience header including all polymorphic residue types
 *
 * All types are in the x3dna::core::poly namespace to avoid conflicts
 * with the legacy types in x3dna::core.
 */

#pragma once

// Interface
#include <x3dna/core/residue/iresidue.hpp>

// Concrete types
#include <x3dna/core/residue/rna.hpp>
#include <x3dna/core/residue/dna.hpp>
#include <x3dna/core/residue/protein.hpp>
#include <x3dna/core/residue/ligand.hpp>

// Factory
#include <x3dna/core/residue/residue_factory.hpp>

// Chain
#include <x3dna/core/residue/chain.hpp>

// Structure
#include <x3dna/core/residue/structure.hpp>

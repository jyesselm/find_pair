/**
 * @file residue.hpp
 * @brief Convenience header including all polymorphic residue types
 *
 * All types are in the x3dna::core::structure namespace to avoid conflicts
 * with the legacy types in x3dna::core.
 */

#pragma once

// Interface
#include <x3dna/core/structure/iresidue.hpp>

// Concrete types
#include <x3dna/core/structure/rna.hpp>
#include <x3dna/core/structure/dna.hpp>
#include <x3dna/core/structure/protein.hpp>
#include <x3dna/core/structure/ligand.hpp>

// Factory
#include <x3dna/core/structure/residue_factory.hpp>

// Chain
#include <x3dna/core/structure/chain.hpp>

// Structure
#include <x3dna/core/structure/structure.hpp>

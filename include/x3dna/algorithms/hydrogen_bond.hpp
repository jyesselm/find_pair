/**
 * @file hydrogen_bond.hpp
 * @brief Convenience header for the H-bond detection module
 *
 * This header provides a unified include for all H-bond related functionality:
 * - Core types (HBond, HBondClassification, etc.)
 * - Detection parameters with presets
 * - Role classifier for donor/acceptor determination
 * - Geometry utilities
 * - HBondDetector for general-purpose detection
 * - InteractionFilter for filtering by interaction type
 *
 * Example usage:
 * @code
 * #include <x3dna/algorithms/hydrogen_bond.hpp>
 *
 * using namespace x3dna::algorithms::hydrogen_bond;
 *
 * // Detect base-base H-bonds (legacy compatible)
 * HBondDetector detector(HBondDetectionParams::legacy_compatible());
 * auto hbonds = detector.detect_base_hbonds(residue1, residue2);
 *
 * // Filter for specific interaction types
 * auto base_hbonds = InteractionFilter::filter(hbonds, HBondInteractionType::BASE_BASE);
 *
 * // Detect all RNA H-bonds including backbone
 * HBondDetector general_detector(HBondDetectionParams::general());
 * auto all_hbonds = general_detector.detect_all_hbonds_between(
 *     residue1, residue2, MoleculeType::NUCLEIC_ACID, MoleculeType::NUCLEIC_ACID);
 * @endcode
 */

#pragma once

// Core types
#include <x3dna/algorithms/hydrogen_bond/hbond.hpp>
#include <x3dna/algorithms/hydrogen_bond/hbond_types.hpp>

// Module components
#include <x3dna/algorithms/hydrogen_bond/detection_params.hpp>
#include <x3dna/algorithms/hydrogen_bond/role_classifier.hpp>
#include <x3dna/algorithms/hydrogen_bond/geometry.hpp>
#include <x3dna/algorithms/hydrogen_bond/detector.hpp>
#include <x3dna/algorithms/hydrogen_bond/interaction_filter.hpp>

// Legacy utilities (for backward compatibility)
#include <x3dna/algorithms/hydrogen_bond/hydrogen_bond_utils.hpp>

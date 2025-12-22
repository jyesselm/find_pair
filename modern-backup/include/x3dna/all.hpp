/**
 * @file all.hpp
 * @brief Convenience header that includes all commonly used x3dna components
 *
 * For external projects, simply include this header to get access to:
 * - Core types (Structure, Residue, Atom, BasePair, ReferenceFrame)
 * - I/O (PdbParser, CifParser, JsonWriter)
 * - Protocols (FindPairProtocol, AnalyzeProtocol)
 * - Geometry types (Vector3D, Matrix3D)
 *
 * Example:
 * @code
 *   #include <x3dna/all.hpp>
 *
 *   int main() {
 *       // Initialize library (auto-detects resources)
 *       x3dna::init();
 *
 *       // Parse a PDB file
 *       auto structure = x3dna::load_structure("1ehz.pdb");
 *
 *       // Find base pairs
 *       auto pairs = x3dna::find_base_pairs(structure);
 *
 *       // Calculate step parameters
 *       auto params = x3dna::calculate_parameters(pairs);
 *
 *       return 0;
 *   }
 * @endcode
 */

#pragma once

// Main library header (init, version, etc.)
#include <x3dna/x3dna.hpp>

// Core types
#include <x3dna/core/atom.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/parameters.hpp>

// Geometry
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>

// I/O
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/io/cif_parser.hpp>
#include <x3dna/io/pdb_writer.hpp>
#include <x3dna/io/json_writer.hpp>

// Protocols
#include <x3dna/protocols/find_pair_protocol.hpp>
#include <x3dna/protocols/analyze_protocol.hpp>

// Algorithms (commonly used)
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/parameter_calculator.hpp>

namespace x3dna {

// ============================================================================
// High-level convenience functions for common use cases
// ============================================================================

/**
 * @brief Load a structure from a PDB or CIF file
 *
 * Automatically detects file format based on extension.
 *
 * @param file_path Path to structure file (.pdb, .cif, .mmcif)
 * @return Parsed structure
 * @throws std::runtime_error if file cannot be parsed
 */
[[nodiscard]] inline core::Structure load_structure(const std::filesystem::path& file_path) {
    std::string ext = file_path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if (ext == ".pdb" || ext == ".ent") {
        io::PdbParser parser;
        return parser.parse_file(file_path);
    } else if (ext == ".cif" || ext == ".mmcif") {
        io::CifParser parser;
        return parser.parse_file(file_path);
    } else {
        throw std::runtime_error("Unknown file extension: " + ext + ". Supported: .pdb, .cif, .mmcif");
    }
}

/**
 * @brief Find base pairs in a structure
 *
 * Runs the complete find_pair protocol including:
 * - Frame calculation for all residues
 * - Base pair detection and validation
 * - Hydrogen bond analysis
 *
 * @param structure Structure to analyze (will be modified with frames)
 * @return Vector of detected base pairs
 */
[[nodiscard]] inline std::vector<core::BasePair> find_base_pairs(core::Structure& structure) {
    protocols::FindPairConfig config;
    config.legacy_mode = false;

    protocols::FindPairProtocol protocol(config::ResourceLocator::templates_dir(), config);
    protocol.execute(structure);

    return protocol.base_pairs();
}

/**
 * @brief Find base pairs with custom configuration
 *
 * @param structure Structure to analyze
 * @param config Protocol configuration
 * @return Vector of detected base pairs
 */
[[nodiscard]] inline std::vector<core::BasePair> find_base_pairs(core::Structure& structure,
                                                                 const protocols::FindPairConfig& config) {
    protocols::FindPairProtocol protocol(config::ResourceLocator::templates_dir(), config);
    protocol.execute(structure);
    return protocol.base_pairs();
}

/**
 * @brief Calculate step parameters for base pairs
 *
 * Calculates the 6 base pair step parameters:
 * - Shift, Slide, Rise (translational)
 * - Tilt, Roll, Twist (rotational)
 *
 * @param pairs Vector of base pairs (must have valid frames)
 * @return Vector of step parameter sets
 */
[[nodiscard]] inline std::vector<core::BasePairStepParameters> calculate_step_parameters(
    const std::vector<core::BasePair>& pairs) {
    algorithms::ParameterCalculator calc;
    return calc.calculate_all_step_parameters(pairs);
}

/**
 * @brief Calculate helical parameters for base pairs
 *
 * Calculates the 6 helical parameters:
 * - x-displacement, y-displacement, rise (translational)
 * - inclination, tip, twist (rotational)
 *
 * @param pairs Vector of base pairs (must have valid frames)
 * @return Vector of helical parameter sets
 */
[[nodiscard]] inline std::vector<core::HelicalParameters> calculate_helical_parameters(
    const std::vector<core::BasePair>& pairs) {
    algorithms::ParameterCalculator calc;
    std::vector<core::HelicalParameters> results;

    for (size_t i = 0; i + 1 < pairs.size(); ++i) {
        results.push_back(calc.calculate_helical_parameters(pairs[i], pairs[i + 1]));
    }

    return results;
}

/**
 * @brief Write structure to PDB file
 *
 * @param structure Structure to write
 * @param file_path Output file path
 */
inline void save_structure(const core::Structure& structure, const std::filesystem::path& file_path) {
    io::PdbWriter writer;
    writer.write_file(structure, file_path);
}

} // namespace x3dna

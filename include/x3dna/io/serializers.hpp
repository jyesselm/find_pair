/**
 * @file serializers.hpp
 * @brief JSON serializers for core objects
 *
 * Provides dedicated serialization classes following Single Responsibility Principle.
 * Each serializer handles conversion to/from both legacy and modern JSON formats.
 *
 * Note: JSON serialization is intentionally separated from core classes to:
 * 1. Allow multiple JSON formats without cluttering core classes
 * 2. Keep core classes focused on domain logic
 * 3. Make it easier to add new serialization formats
 */

#pragma once

#include <nlohmann/json.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/nucleotide_utils.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/reference_frame.hpp>

namespace x3dna {
namespace io {

/**
 * @class AtomSerializer
 * @brief Serializes Atom objects to/from JSON
 *
 * Uses original (padded) names for JSON output to match legacy format.
 * When deserializing, the Atom constructor will trim names automatically.
 */
class AtomSerializer {
public:
    /**
     * @brief Convert Atom to legacy JSON format (pdb_atoms record)
     * @param atom The atom to serialize
     * @param residue Residue context (required for residue-level fields)
     */
    [[nodiscard]] static nlohmann::json to_legacy_json(const core::Atom& atom, const core::Residue& residue) {
        nlohmann::json j;
        // Use atom name (already in PDB 4-character format)
        j["atom_name"] = atom.name();
        j["xyz"] = {atom.position().x(), atom.position().y(), atom.position().z()};

        // Use residue context for residue-level fields
        j["residue_name"] = residue.name();
        j["chain_id"] = residue.chain_id();
        j["residue_seq"] = residue.seq_num();

        // Determine record_type from residue classification
        char olc = core::one_letter_code(residue);
        bool is_hetatm = !residue.is_nucleotide() || (olc >= 'a' && olc <= 'z');
        j["record_type"] = is_hetatm ? "H" : "A";

        if (atom.alt_loc() != ' ' && atom.alt_loc() != '\0') {
            j["alt_loc"] = std::string(1, atom.alt_loc());
        }

        // Use residue insertion
        const std::string& insertion = residue.insertion();
        if (!insertion.empty()) {
            j["insertion"] = insertion;
        }

        // Add debug information
        j["occupancy"] = atom.occupancy();
        if (atom.atom_serial() > 0) {
            j["atom_serial"] = atom.atom_serial();
        }
        if (atom.model_number() > 0) {
            j["model_number"] = atom.model_number();
        }
        if (atom.b_factor() != 0.0) {
            j["b_factor"] = atom.b_factor();
        }
        if (!atom.element().empty()) {
            j["element"] = atom.element();
        }

        return j;
    }

    /**
     * @brief Convert Atom to modern JSON format
     * @param atom The atom to serialize
     * @param residue Residue context (required for residue-level fields)
     */
    [[nodiscard]] static nlohmann::json to_json(const core::Atom& atom, const core::Residue& residue) {
        nlohmann::json j;

        // Use legacy_atom_idx for atom_idx to match legacy exactly (1-based)
        int legacy_atom_idx = atom.legacy_atom_idx();
        j["atom_idx"] = legacy_atom_idx > 0 ? legacy_atom_idx : 0;

        // Core fields - use atom name (already in PDB 4-character format)
        j["atom_name"] = atom.name();

        // Use residue context for residue-level fields
        j["residue_name"] = residue.name();
        j["chain_id"] = residue.chain_id();
        j["residue_seq"] = residue.seq_num();

        // Insertion code (only if non-empty)
        const std::string& insertion = residue.insertion();
        if (!insertion.empty()) {
            j["insertion"] = insertion;
        }

        // Unique residue identifier
        j["res_id"] = residue.res_id();

        // Coordinates
        j["xyz"] = nlohmann::json::array({atom.position().x(), atom.position().y(), atom.position().z()});

        // Determine record_type from residue classification
        char olc = core::one_letter_code(residue);
        bool is_hetatm = !residue.is_nucleotide() || (olc >= 'a' && olc <= 'z');
        j["record_type"] = is_hetatm ? "H" : "A";

        return j;
    }

    /**
     * @brief Create Atom from legacy JSON format
     * Note: Residue-level fields (residue_name, chain_id, residue_seq, insertion) are ignored.
     *       These fields should be obtained from the parent Residue object.
     */
    [[nodiscard]] static core::Atom from_legacy_json(const nlohmann::json& j) {
        std::string name = j.value("atom_name", "");
        std::vector<double> xyz = j.value("xyz", std::vector<double>{0.0, 0.0, 0.0});
        geometry::Vector3D position(xyz[0], xyz[1], xyz[2]);

        // Create atom with per-atom fields only
        // Note: record_type is now determined from residue classification, not stored on atom
        return core::Atom(name, position);
    }

    /**
     * @brief Create Atom from modern JSON format
     * Note: Residue-level fields (residue_name, chain_id, residue_seq, insertion) are ignored.
     *       These fields should be obtained from the parent Residue object.
     */
    [[nodiscard]] static core::Atom from_json(const nlohmann::json& j) {
        return from_legacy_json(j); // Same format for now
    }
};

/**
 * @class ReferenceFrameSerializer
 * @brief Serializes ReferenceFrame objects to/from JSON
 */
class ReferenceFrameSerializer {
public:
    /**
     * @brief Convert ReferenceFrame to legacy JSON format
     */
    [[nodiscard]] static nlohmann::json to_legacy_json(const core::ReferenceFrame& frame) {
        return frame.to_json_legacy();
    }

    /**
     * @brief Convert ReferenceFrame to modern JSON format
     */
    [[nodiscard]] static nlohmann::json to_json(const core::ReferenceFrame& frame) {
        return frame.to_json();
    }

    /**
     * @brief Create ReferenceFrame from legacy JSON format
     */
    [[nodiscard]] static core::ReferenceFrame from_legacy_json(const nlohmann::json& j) {
        return core::ReferenceFrame::from_json_legacy(j);
    }

    /**
     * @brief Create ReferenceFrame from modern JSON format
     */
    [[nodiscard]] static core::ReferenceFrame from_json(const nlohmann::json& j) {
        return core::ReferenceFrame::from_json(j);
    }
};

/**
 * @class ResidueSerializer
 * @brief Serializes Residue objects to/from JSON
 */
class ResidueSerializer {
public:
    /**
     * @brief Convert Residue to legacy JSON format
     */
    [[nodiscard]] static nlohmann::json to_legacy_json(const core::Residue& residue) {
        nlohmann::json j;
        j["residue_name"] = residue.name();
        j["residue_seq"] = residue.seq_num();
        j["chain_id"] = residue.chain_id();
        j["atoms"] = nlohmann::json::array();
        // Pass residue context to ensure atoms use residue's fields
        for (const auto& atom : residue.atoms()) {
            j["atoms"].push_back(AtomSerializer::to_legacy_json(atom, residue));
        }
        if (residue.reference_frame().has_value()) {
            j["reference_frame"] = ReferenceFrameSerializer::to_legacy_json(residue.reference_frame().value());
        }
        return j;
    }

    /**
     * @brief Convert Residue to modern JSON format
     */
    [[nodiscard]] static nlohmann::json to_json(const core::Residue& residue) {
        nlohmann::json j;
        j["name"] = residue.name();
        j["seq_num"] = residue.seq_num();
        j["chain_id"] = residue.chain_id();
        j["res_id"] = residue.res_id();
        j["atoms"] = nlohmann::json::array();
        // Pass residue context to ensure atoms use residue's fields
        for (const auto& atom : residue.atoms()) {
            j["atoms"].push_back(AtomSerializer::to_json(atom, residue));
        }
        if (residue.reference_frame().has_value()) {
            j["reference_frame"] = ReferenceFrameSerializer::to_json(residue.reference_frame().value());
        }
        return j;
    }

    /**
     * @brief Create Residue from legacy JSON format
     */
    [[nodiscard]] static core::Residue from_legacy_json(const nlohmann::json& j) {
        std::string name = j.value("residue_name", "");
        int seq_num = j.value("residue_seq", 0);
        std::string chain_id = j.value("chain_id", "");
        std::string insertion = j.value("insertion", "");

        // Trim name for classification lookup
        auto trim = [](const std::string& str) -> std::string {
            size_t start = str.find_first_not_of(" \t");
            if (start == std::string::npos)
                return "";
            size_t end = str.find_last_not_of(" \t");
            return str.substr(start, end - start + 1);
        };
        std::string trimmed = trim(name);

        // Get classification from trimmed name
        auto classification = core::typing::TypeRegistry::instance().classify_residue(trimmed);

        // Collect atoms
        std::vector<core::Atom> atoms;
        if (j.contains("atoms") && j["atoms"].is_array()) {
            for (const auto& atom_json : j["atoms"]) {
                atoms.push_back(AtomSerializer::from_legacy_json(atom_json));
            }
        }

        // Build residue preserving exact name from JSON
        // Note: one_letter_code is now computed from name, not stored
        core::Residue residue = core::Residue::create(name, seq_num, chain_id)
                                    .insertion(insertion)
                                    .classification(classification)
                                    .atoms(atoms)
                                    .build();

        // Add reference frame if present
        if (j.contains("reference_frame")) {
            residue.set_reference_frame(ReferenceFrameSerializer::from_legacy_json(j["reference_frame"]));
        }

        return residue;
    }

    /**
     * @brief Create Residue from modern JSON format
     */
    [[nodiscard]] static core::Residue from_json(const nlohmann::json& j) {
        std::string name = j.value("name", "");
        int seq_num = j.value("seq_num", 0);
        std::string chain_id = j.value("chain_id", "");
        std::string insertion = j.value("insertion", "");

        // Trim name for classification lookup
        auto trim = [](const std::string& str) -> std::string {
            size_t start = str.find_first_not_of(" \t");
            if (start == std::string::npos)
                return "";
            size_t end = str.find_last_not_of(" \t");
            return str.substr(start, end - start + 1);
        };
        std::string trimmed = trim(name);

        // Get classification from trimmed name
        auto classification = core::typing::TypeRegistry::instance().classify_residue(trimmed);

        // Collect atoms
        std::vector<core::Atom> atoms;
        if (j.contains("atoms") && j["atoms"].is_array()) {
            for (const auto& atom_json : j["atoms"]) {
                atoms.push_back(AtomSerializer::from_json(atom_json));
            }
        }

        // Build residue preserving exact name from JSON
        // Note: one_letter_code is now computed from name, not stored
        core::Residue residue = core::Residue::create(name, seq_num, chain_id)
                                    .insertion(insertion)
                                    .classification(classification)
                                    .atoms(atoms)
                                    .build();

        // Add reference frame if present
        if (j.contains("reference_frame")) {
            residue.set_reference_frame(ReferenceFrameSerializer::from_json(j["reference_frame"]));
        }

        return residue;
    }
};

/**
 * @class ChainSerializer
 * @brief Serializes Chain objects to/from JSON
 */
class ChainSerializer {
public:
    /**
     * @brief Convert Chain to legacy JSON format
     */
    [[nodiscard]] static nlohmann::json to_legacy_json(const core::Chain& chain) {
        nlohmann::json j;
        j["chain_id"] = chain.chain_id();
        j["num_residues"] = chain.num_residues();
        j["residues"] = nlohmann::json::array();
        for (const auto& residue : chain.residues()) {
            j["residues"].push_back(ResidueSerializer::to_legacy_json(residue));
        }
        return j;
    }

    /**
     * @brief Convert Chain to modern JSON format
     */
    [[nodiscard]] static nlohmann::json to_json(const core::Chain& chain) {
        nlohmann::json j;
        j["chain_id"] = chain.chain_id();
        j["residues"] = nlohmann::json::array();
        for (const auto& residue : chain.residues()) {
            j["residues"].push_back(ResidueSerializer::to_json(residue));
        }
        return j;
    }

    /**
     * @brief Create Chain from legacy JSON format
     */
    [[nodiscard]] static core::Chain from_legacy_json(const nlohmann::json& j) {
        std::string chain_id = j.value("chain_id", "");
        core::Chain chain(chain_id);

        if (j.contains("residues") && j["residues"].is_array()) {
            for (const auto& residue_json : j["residues"]) {
                chain.add_residue(ResidueSerializer::from_legacy_json(residue_json));
            }
        }

        return chain;
    }

    /**
     * @brief Create Chain from modern JSON format
     */
    [[nodiscard]] static core::Chain from_json(const nlohmann::json& j) {
        std::string chain_id = j.value("chain_id", "");
        core::Chain chain(chain_id);

        if (j.contains("residues") && j["residues"].is_array()) {
            for (const auto& residue_json : j["residues"]) {
                chain.add_residue(ResidueSerializer::from_json(residue_json));
            }
        }

        return chain;
    }
};

/**
 * @class StructureSerializer
 * @brief Serializes Structure objects to/from JSON
 */
class StructureSerializer {
public:
    /**
     * @brief Convert Structure to legacy JSON format
     */
    [[nodiscard]] static nlohmann::json to_legacy_json(const core::Structure& structure) {
        nlohmann::json j;
        j["pdb_id"] = structure.pdb_id();
        j["num_atoms"] = structure.num_atoms();
        j["atoms"] = nlohmann::json::array();

        // Flatten atoms from all chains/residues
        for (const auto& chain : structure.chains()) {
            for (const auto& residue : chain.residues()) {
                for (const auto& atom : residue.atoms()) {
                    j["atoms"].push_back(AtomSerializer::to_legacy_json(atom, residue));
                }
            }
        }

        return j;
    }

    /**
     * @brief Convert Structure to modern JSON format
     */
    [[nodiscard]] static nlohmann::json to_json(const core::Structure& structure) {
        nlohmann::json j;
        j["pdb_id"] = structure.pdb_id();
        j["chains"] = nlohmann::json::array();
        for (const auto& chain : structure.chains()) {
            j["chains"].push_back(ChainSerializer::to_json(chain));
        }
        return j;
    }

    /**
     * @brief Create Structure from legacy JSON format
     */
    [[nodiscard]] static core::Structure from_legacy_json(const nlohmann::json& j) {
        std::string pdb_id = j.value("pdb_id", "");
        core::Structure structure(pdb_id);

        if (j.contains("atoms") && j["atoms"].is_array()) {
            // Group atoms by chain and residue using JSON metadata
            std::map<std::string, std::map<std::tuple<std::string, int, std::string>, std::vector<core::Atom>>>
                chain_residue_atoms;

            for (const auto& atom_json : j["atoms"]) {
                core::Atom atom = AtomSerializer::from_legacy_json(atom_json);

                // Extract residue-level fields from JSON (not from atom)
                std::string chain_id = atom_json.value("chain_id", "");
                std::string residue_name = atom_json.value("residue_name", "");
                int seq_num = atom_json.value("residue_seq", 0);
                std::string insertion = atom_json.value("insertion", "");

                auto residue_key = std::make_tuple(residue_name, seq_num, insertion);
                chain_residue_atoms[chain_id][residue_key].push_back(atom);
            }

            // Build structure from grouped atoms
            for (const auto& [chain_id, residue_atoms] : chain_residue_atoms) {
                core::Chain chain(chain_id);
                for (const auto& [residue_key, atoms] : residue_atoms) {
                    if (!atoms.empty()) {
                        std::string residue_name = std::get<0>(residue_key);
                        int seq_num = std::get<1>(residue_key);
                        std::string insertion = std::get<2>(residue_key);
                        core::Residue residue = core::Residue::create_from_atoms(residue_name, seq_num, chain_id,
                                                                                 insertion, atoms);
                        chain.add_residue(residue);
                    }
                }
                structure.add_chain(chain);
            }
        }

        return structure;
    }

    /**
     * @brief Create Structure from modern JSON format
     */
    [[nodiscard]] static core::Structure from_json(const nlohmann::json& j) {
        std::string pdb_id = j.value("pdb_id", "");
        core::Structure structure(pdb_id);

        if (j.contains("chains") && j["chains"].is_array()) {
            for (const auto& chain_json : j["chains"]) {
                structure.add_chain(ChainSerializer::from_json(chain_json));
            }
        }

        return structure;
    }
};

/**
 * @class BasePairSerializer
 * @brief Serializes BasePair objects to/from JSON
 */
class BasePairSerializer {
public:
    /**
     * @brief Convert BasePair to legacy JSON format (base_pair record)
     */
    [[nodiscard]] static nlohmann::json to_legacy_json(const core::BasePair& bp) {
        return bp.to_json_legacy();
    }

    /**
     * @brief Convert BasePair to modern JSON format
     */
    [[nodiscard]] static nlohmann::json to_json(const core::BasePair& bp) {
        return bp.to_json();
    }

    /**
     * @brief Create BasePair from legacy JSON format
     */
    [[nodiscard]] static core::BasePair from_legacy_json(const nlohmann::json& j) {
        return core::BasePair::from_json_legacy(j);
    }

    /**
     * @brief Create BasePair from modern JSON format
     */
    [[nodiscard]] static core::BasePair from_json(const nlohmann::json& j) {
        return core::BasePair::from_json(j);
    }
};

} // namespace io
} // namespace x3dna

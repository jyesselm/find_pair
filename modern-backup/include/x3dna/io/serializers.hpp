/**
 * @file serializers.hpp
 * @brief JSON serializers for core objects
 *
 * Provides dedicated serialization classes following Single Responsibility Principle.
 * Each serializer handles conversion to/from both legacy and modern JSON formats.
 */

#pragma once

#include <nlohmann/json.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/reference_frame.hpp>

namespace x3dna {
namespace io {

/**
 * @class AtomSerializer
 * @brief Serializes Atom objects to/from JSON
 */
class AtomSerializer {
public:
    /**
     * @brief Convert Atom to legacy JSON format (pdb_atoms record)
     */
    [[nodiscard]] static nlohmann::json to_legacy_json(const core::Atom& atom) {
        return atom.to_json_legacy();
    }

    /**
     * @brief Convert Atom to modern JSON format
     */
    [[nodiscard]] static nlohmann::json to_json(const core::Atom& atom) {
        return atom.to_json();
    }

    /**
     * @brief Create Atom from legacy JSON format
     */
    [[nodiscard]] static core::Atom from_legacy_json(const nlohmann::json& j) {
        return core::Atom::from_json_legacy(j);
    }

    /**
     * @brief Create Atom from modern JSON format
     */
    [[nodiscard]] static core::Atom from_json(const nlohmann::json& j) {
        return core::Atom::from_json(j);
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
        return residue.to_json_legacy();
    }

    /**
     * @brief Convert Residue to modern JSON format
     */
    [[nodiscard]] static nlohmann::json to_json(const core::Residue& residue) {
        return residue.to_json();
    }

    /**
     * @brief Create Residue from legacy JSON format
     */
    [[nodiscard]] static core::Residue from_legacy_json(const nlohmann::json& j) {
        return core::Residue::from_json_legacy(j);
    }

    /**
     * @brief Create Residue from modern JSON format
     */
    [[nodiscard]] static core::Residue from_json(const nlohmann::json& j) {
        return core::Residue::from_json(j);
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
        return chain.to_json_legacy();
    }

    /**
     * @brief Convert Chain to modern JSON format
     */
    [[nodiscard]] static nlohmann::json to_json(const core::Chain& chain) {
        return chain.to_json();
    }

    /**
     * @brief Create Chain from legacy JSON format
     */
    [[nodiscard]] static core::Chain from_legacy_json(const nlohmann::json& j) {
        return core::Chain::from_json_legacy(j);
    }

    /**
     * @brief Create Chain from modern JSON format
     */
    [[nodiscard]] static core::Chain from_json(const nlohmann::json& j) {
        return core::Chain::from_json(j);
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
        return structure.to_json_legacy();
    }

    /**
     * @brief Convert Structure to modern JSON format
     */
    [[nodiscard]] static nlohmann::json to_json(const core::Structure& structure) {
        return structure.to_json();
    }

    /**
     * @brief Create Structure from legacy JSON format
     */
    [[nodiscard]] static core::Structure from_legacy_json(const nlohmann::json& j) {
        return core::Structure::from_json_legacy(j);
    }

    /**
     * @brief Create Structure from modern JSON format
     */
    [[nodiscard]] static core::Structure from_json(const nlohmann::json& j) {
        return core::Structure::from_json(j);
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

} // namespace io
} // namespace x3dna

/**
 * @file json_reader.hpp
 * @brief JSON file reader for modern and legacy formats
 */

#pragma once

#include <string>
#include <filesystem>
#include <vector>
#include <nlohmann/json.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/reference_frame.hpp>

namespace x3dna {
namespace io {

/**
 * @class JsonReader
 * @brief Reads Structure and calculation records from JSON format
 */
class JsonReader {
public:
    /**
     * @brief Read Structure from modern JSON format file
     * @param path Path to JSON file
     * @return Structure object
     */
    static core::Structure read_structure(const std::filesystem::path& path);

    /**
     * @brief Read Structure from modern JSON object
     * @param json JSON object
     * @return Structure object
     */
    static core::Structure read_structure(const nlohmann::json& json);

    /**
     * @brief Read Structure from legacy JSON format file
     * @param path Path to legacy JSON file
     * @return Structure object
     */
    static core::Structure read_structure_legacy(const std::filesystem::path& path);

    /**
     * @brief Read Structure from legacy JSON object
     * @param json JSON object
     * @return Structure object
     */
    static core::Structure read_structure_legacy(const nlohmann::json& json);

    /**
     * @brief Read base pairs from JSON
     * @param json JSON object
     * @return Vector of base pairs
     */
    static std::vector<core::BasePair> read_base_pairs(const nlohmann::json& json);

    /**
     * @brief Read reference frames from JSON
     * @param json JSON object
     * @return Vector of reference frames with residue indices
     */
    static std::vector<std::pair<size_t, core::ReferenceFrame>>
    read_ref_frames(const nlohmann::json& json);

    /**
     * @brief Find calculation records by type
     * @param json JSON object
     * @param record_type Type string (e.g., "pdb_atoms", "base_pair")
     * @return Vector of matching records
     */
    static std::vector<nlohmann::json> find_records_by_type(const nlohmann::json& json,
                                                            const std::string& record_type);

private:
    /**
     * @brief Load JSON from file
     * @param path Path to JSON file
     * @return JSON object
     */
    static nlohmann::json load_json_file(const std::filesystem::path& path);
};

} // namespace io
} // namespace x3dna

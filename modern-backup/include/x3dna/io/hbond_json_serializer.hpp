/**
 * @file hbond_json_serializer.hpp
 * @brief Serializes H-bonds to/from legacy JSON format
 */

#pragma once

#include <utility>
#include <vector>
#include <nlohmann/json.hpp>
#include <x3dna/core/hbond.hpp>

namespace x3dna {
namespace io {

/**
 * @brief Serializes H-bonds to/from legacy JSON format
 *
 * Ensures exact compatibility with legacy X3DNA JSON output.
 */
class HBondJsonSerializer {
public:
    /**
     * @brief Create legacy hbond_list JSON record
     * @param first_residue_idx 0-based index (output as 1-based)
     * @param second_residue_idx 0-based index (output as 1-based)
     * @param bonds H-bonds to serialize
     * @return JSON matching legacy format exactly
     */
    [[nodiscard]] static nlohmann::json to_hbond_list_record(size_t first_residue_idx, size_t second_residue_idx,
                                                             const std::vector<core::HBond>& bonds);

    /**
     * @brief Serialize single H-bond to legacy JSON
     * @param bond H-bond to serialize
     * @param one_based_index 1-based index for hbond_idx field
     */
    [[nodiscard]] static nlohmann::json bond_to_json(const core::HBond& bond, size_t one_based_index);

    /**
     * @brief Deserialize H-bond from legacy JSON
     */
    [[nodiscard]] static core::HBond bond_from_json(const nlohmann::json& j);

    /**
     * @brief Deserialize hbond_list record
     * @return Pair of (residue indices as pair, bonds vector)
     */
    [[nodiscard]] static std::pair<std::pair<size_t, size_t>, std::vector<core::HBond>> from_hbond_list_record(
        const nlohmann::json& j);

private:
    // Legacy JSON field names (MUST NOT CHANGE)
    static constexpr const char* FIELD_TYPE = "type";
    static constexpr const char* FIELD_BASE_I = "base_i";
    static constexpr const char* FIELD_BASE_J = "base_j";
    static constexpr const char* FIELD_NUM_HBONDS = "num_hbonds";
    static constexpr const char* FIELD_HBONDS = "hbonds";
    static constexpr const char* FIELD_HBOND_IDX = "hbond_idx";
    static constexpr const char* FIELD_DONOR_ATOM = "donor_atom";
    static constexpr const char* FIELD_ACCEPTOR_ATOM = "acceptor_atom";
    static constexpr const char* FIELD_DISTANCE = "distance";
    static constexpr const char* FIELD_HBOND_TYPE = "type";
    static constexpr const char* TYPE_VALUE = "hbond_list";
};

} // namespace io
} // namespace x3dna

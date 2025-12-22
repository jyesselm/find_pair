/**
 * @file hbond_json_serializer.cpp
 * @brief Implementation of H-bond JSON serialization
 */

#include <x3dna/io/hbond_json_serializer.hpp>
#include <cmath>
#include <stdexcept>

namespace x3dna {
namespace io {

nlohmann::json HBondJsonSerializer::to_hbond_list_record(size_t first_residue_idx, size_t second_residue_idx,
                                                         const std::vector<core::HBond>& bonds) {
    nlohmann::json record = nlohmann::json::object();
    record[FIELD_TYPE] = TYPE_VALUE;
    record[FIELD_BASE_I] = static_cast<long>(first_residue_idx + 1);  // Convert to 1-based
    record[FIELD_BASE_J] = static_cast<long>(second_residue_idx + 1); // Convert to 1-based
    record[FIELD_NUM_HBONDS] = bonds.size();

    nlohmann::json hbonds_array = nlohmann::json::array();
    for (size_t i = 0; i < bonds.size(); ++i) {
        hbonds_array.push_back(bond_to_json(bonds[i], i + 1)); // Use 1-based index
    }
    record[FIELD_HBONDS] = hbonds_array;

    return record;
}

nlohmann::json HBondJsonSerializer::bond_to_json(const core::HBond& bond, size_t one_based_index) {
    nlohmann::json j = nlohmann::json::object();

    // Use detection_index if set, otherwise use provided one_based_index
    if (bond.detection_index.has_value()) {
        j[FIELD_HBOND_IDX] = static_cast<long>(bond.detection_index.value());
    } else {
        j[FIELD_HBOND_IDX] = static_cast<long>(one_based_index);
    }

    // Atom names (exact 4-char format from HBond)
    j[FIELD_DONOR_ATOM] = bond.donor_atom_name;
    j[FIELD_ACCEPTOR_ATOM] = bond.acceptor_atom_name;

    // Distance with proper formatting
    if (std::isnan(bond.distance) || std::isinf(bond.distance)) {
        j[FIELD_DISTANCE] = nlohmann::json(nullptr);
    } else {
        j[FIELD_DISTANCE] = std::round(bond.distance * 1000000.0) / 1000000.0;
    }

    // Legacy type character
    j[FIELD_HBOND_TYPE] = std::string(1, bond.legacy_type_char());

    return j;
}

core::HBond HBondJsonSerializer::bond_from_json(const nlohmann::json& j) {
    if (!j.is_object()) {
        throw std::invalid_argument("HBond JSON must be an object");
    }

    core::HBond bond;

    // Parse atom names
    if (j.contains(FIELD_DONOR_ATOM) && j[FIELD_DONOR_ATOM].is_string()) {
        bond.donor_atom_name = j[FIELD_DONOR_ATOM].get<std::string>();
    }

    if (j.contains(FIELD_ACCEPTOR_ATOM) && j[FIELD_ACCEPTOR_ATOM].is_string()) {
        bond.acceptor_atom_name = j[FIELD_ACCEPTOR_ATOM].get<std::string>();
    }

    // Parse distance
    if (j.contains(FIELD_DISTANCE)) {
        if (j[FIELD_DISTANCE].is_null()) {
            bond.distance = 0.0;
        } else if (j[FIELD_DISTANCE].is_number()) {
            bond.distance = j[FIELD_DISTANCE].get<double>();
        }
    }

    // Parse hbond_idx
    if (j.contains(FIELD_HBOND_IDX) && j[FIELD_HBOND_IDX].is_number()) {
        bond.detection_index = j[FIELD_HBOND_IDX].get<size_t>();
    }

    // Parse type and set classification
    if (j.contains(FIELD_HBOND_TYPE) && j[FIELD_HBOND_TYPE].is_string()) {
        std::string type_str = j[FIELD_HBOND_TYPE].get<std::string>();
        if (!type_str.empty()) {
            char type_char = type_str[0];
            if (type_char == '-') {
                bond.classification = core::HBondClassification::STANDARD;
            } else if (type_char == '*') {
                bond.classification = core::HBondClassification::NON_STANDARD;
            } else if (type_char == ' ') {
                bond.classification = core::HBondClassification::INVALID;
            } else {
                bond.classification = core::HBondClassification::UNKNOWN;
            }
        }
    }

    return bond;
}

std::pair<std::pair<size_t, size_t>, std::vector<core::HBond>> HBondJsonSerializer::from_hbond_list_record(
    const nlohmann::json& j) {
    if (!j.is_object()) {
        throw std::invalid_argument("hbond_list record must be an object");
    }

    if (!j.contains(FIELD_TYPE) || j[FIELD_TYPE] != TYPE_VALUE) {
        throw std::invalid_argument("Invalid hbond_list record: missing or incorrect type field");
    }

    // Parse residue indices (convert from 1-based to 0-based)
    size_t base_i = 0;
    size_t base_j = 0;

    if (j.contains(FIELD_BASE_I) && j[FIELD_BASE_I].is_number()) {
        base_i = j[FIELD_BASE_I].get<size_t>() - 1; // Convert to 0-based
    }

    if (j.contains(FIELD_BASE_J) && j[FIELD_BASE_J].is_number()) {
        base_j = j[FIELD_BASE_J].get<size_t>() - 1; // Convert to 0-based
    }

    // Parse bonds array
    std::vector<core::HBond> bonds;
    if (j.contains(FIELD_HBONDS) && j[FIELD_HBONDS].is_array()) {
        for (const auto& bond_json : j[FIELD_HBONDS]) {
            bonds.push_back(bond_from_json(bond_json));
        }
    }

    return std::make_pair(std::make_pair(base_i, base_j), bonds);
}

} // namespace io
} // namespace x3dna

/**
 * @file json_reader.cpp
 * @brief Implementation of JSON reader
 */

#include <x3dna/io/json_reader.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/base_pair.hpp>
#include <fstream>
#include <stdexcept>

namespace x3dna {
namespace io {

nlohmann::json JsonReader::load_json_file(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open JSON file: " + path.string());
    }
    
    nlohmann::json json;
    try {
        file >> json;
    } catch (const nlohmann::json::parse_error& e) {
        throw std::runtime_error("JSON parse error in " + path.string() + ": " + e.what());
    }
    
    return json;
}

core::Structure JsonReader::read_structure(const std::filesystem::path& path) {
    nlohmann::json json = load_json_file(path);
    return read_structure(json);
}

core::Structure JsonReader::read_structure(const nlohmann::json& json) {
    // Modern format: direct Structure JSON
    if (json.contains("pdb_id") && json.contains("chains")) {
        return core::Structure::from_json(json);
    }
    
    throw std::invalid_argument("JSON does not contain valid Structure data");
}

core::Structure JsonReader::read_structure_legacy(const std::filesystem::path& path) {
    nlohmann::json json = load_json_file(path);
    return read_structure_legacy(json);
}

core::Structure JsonReader::read_structure_legacy(const nlohmann::json& json) {
    // Legacy format: look for pdb_atoms record in calculations array
    if (!json.contains("calculations")) {
        throw std::invalid_argument("Legacy JSON does not contain calculations array");
    }
    
    // Find pdb_atoms record
    auto pdb_atoms_records = find_records_by_type(json, "pdb_atoms");
    if (pdb_atoms_records.empty()) {
        throw std::invalid_argument("No pdb_atoms record found in legacy JSON");
    }
    
    // Use first pdb_atoms record to build structure
    const auto& atoms_record = pdb_atoms_records[0];
    
    // Use Structure::from_json_legacy which handles the pdb_atoms format
    // But we need to add pdb_id from the top-level JSON
    nlohmann::json structure_json = atoms_record;
    std::string pdb_name = json.value("pdb_name", "");
    if (!pdb_name.empty()) {
        structure_json["pdb_id"] = pdb_name;
    }
    
    return core::Structure::from_json_legacy(structure_json);
}

std::vector<core::BasePair> JsonReader::read_base_pairs(const nlohmann::json& json) {
    std::vector<core::BasePair> pairs;
    
    if (!json.contains("calculations")) {
        return pairs;
    }
    
    auto base_pair_records = find_records_by_type(json, "base_pair");
    for (const auto& record : base_pair_records) {
        pairs.push_back(core::BasePair::from_json_legacy(record));
    }
    
    return pairs;
}

std::vector<std::pair<size_t, core::ReferenceFrame>> JsonReader::read_ref_frames(const nlohmann::json& json) {
    std::vector<std::pair<size_t, core::ReferenceFrame>> frames;
    
    if (!json.contains("calculations")) {
        return frames;
    }
    
    // Look for frame_calc or all_ref_frames records
    auto frame_records = find_records_by_type(json, "frame_calc");
    for (const auto& record : frame_records) {
        if (record.contains("residue_idx") && record.contains("frame")) {
            size_t residue_idx = record["residue_idx"].get<size_t>();
            core::ReferenceFrame frame = core::ReferenceFrame::from_json_legacy(record["frame"]);
            frames.emplace_back(residue_idx, frame);
        }
    }
    
    // Also check all_ref_frames record
    auto all_frames_records = find_records_by_type(json, "all_ref_frames");
    for (const auto& record : all_frames_records) {
        if (record.contains("frames") && record["frames"].is_array()) {
            for (const auto& frame_json : record["frames"]) {
                if (frame_json.contains("residue_idx") && frame_json.contains("frame")) {
                    size_t residue_idx = frame_json["residue_idx"].get<size_t>();
                    core::ReferenceFrame frame = core::ReferenceFrame::from_json_legacy(frame_json["frame"]);
                    frames.emplace_back(residue_idx, frame);
                }
            }
        }
    }
    
    return frames;
}

std::vector<nlohmann::json> JsonReader::find_records_by_type(const nlohmann::json& json,
                                                             const std::string& record_type) {
    std::vector<nlohmann::json> records;
    
    if (!json.contains("calculations") || !json["calculations"].is_array()) {
        return records;
    }
    
    for (const auto& record : json["calculations"]) {
        if (record.contains("type") && record["type"] == record_type) {
            records.push_back(record);
        }
    }
    
    return records;
}

} // namespace io
} // namespace x3dna


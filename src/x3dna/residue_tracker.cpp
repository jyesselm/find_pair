#include "x3dna/residue_tracker.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

using json = nlohmann::json;

namespace x3dna {

void ResidueTracker::add_residue(const std::string& chain_id, int residue_seq, const std::string& insertion,
                                 const std::string& residue_name) {
    int read_idx = static_cast<int>(residues_.size());
    residues_.emplace_back(read_idx, chain_id, residue_seq, insertion, residue_name);
}

void ResidueTracker::mark_filtered(int read_index, const std::string& reason) {
    if (read_index >= 0 && read_index < static_cast<int>(residues_.size())) {
        residues_[read_index].filtered = true;
        residues_[read_index].filter_reason = reason;
    }
}

void ResidueTracker::assign_modern_index(int read_index, int modern_index) {
    if (read_index >= 0 && read_index < static_cast<int>(residues_.size())) {
        residues_[read_index].modern_index = modern_index;
    }
}

bool ResidueTracker::load_legacy_indices(const std::string& legacy_json_path) {
    // Try to open file
    std::ifstream file(legacy_json_path);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open legacy JSON: " << legacy_json_path << std::endl;
        return false;
    }

    // Parse JSON
    json j;
    try {
        file >> j;
    } catch (const json::exception& e) {
        std::cerr << "Error parsing legacy JSON: " << e.what() << std::endl;
        return false;
    }

    // Check if this is array format or object format
    json records;
    if (j.is_array()) {
        // Array format - could be records with "type" field or direct records
        records = json::array();
        for (const auto& item : j) {
            // Check if has type field and it's base_frame_calc
            if (item.contains("type") && item["type"] == "base_frame_calc") {
                records.push_back(item);
            }
            // Or if it has residue_idx field (legacy direct format)
            else if (item.contains("residue_idx") && item.contains("base_type")) {
                records.push_back(item);
            }
        }
    } else if (j.is_object() && j.contains("base_frame_calc")) {
        // Object format with "base_frame_calc" key
        records = j["base_frame_calc"];
    } else {
        std::cerr << "Unexpected JSON format - no base_frame_calc records found" << std::endl;
        return false;
    }

    if (records.empty()) {
        std::cerr << "No base_frame_calc records in legacy JSON" << std::endl;
        return false;
    }

    // Load indices from records
    int loaded_count = 0;
    for (const auto& record : records) {
        // Match by PDB properties
        std::string chain = record.value("chain_id", "");
        int seq = record.value("residue_seq", 0);
        std::string ins = record.value("insertion", "");
        int legacy_idx = record.value("residue_idx", -1);

        // Find matching residue
        int read_idx = find_by_pdb_props(chain, seq, ins);
        if (read_idx >= 0) {
            residues_[read_idx].legacy_index = legacy_idx;
            loaded_count++;
        } else {
            std::cerr << "Warning: Legacy residue " << chain << seq << ins << " not found in read residues"
                      << std::endl;
        }
    }

    std::cout << "Loaded " << loaded_count << " legacy indices from " << legacy_json_path << std::endl;

    return loaded_count > 0;
}

int ResidueTracker::find_by_pdb_props(const std::string& chain, int seq, const std::string& ins) const {
    for (size_t i = 0; i < residues_.size(); ++i) {
        const auto& r = residues_[i];
        if (r.chain_id == chain && r.residue_seq == seq && r.insertion == ins) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

ResidueTracker::ValidationResult ResidueTracker::validate() const {
    ValidationResult result;
    result.success = true;
    result.num_residues_read = static_cast<int>(residues_.size());
    result.num_filtered = 0;
    result.num_legacy = 0;
    result.num_modern = 0;
    result.num_matched = 0;
    result.num_unmatched = 0;

    // Count residues by category
    for (const auto& r : residues_) {
        if (r.filtered) {
            result.num_filtered++;
        }
        if (r.legacy_index >= 0) {
            result.num_legacy++;
        }
        if (r.modern_index >= 0) {
            result.num_modern++;
        }
        if (r.legacy_index >= 0 && r.modern_index >= 0) {
            result.num_matched++;
        }
    }

    // Check 1: num_modern should equal num_legacy for perfect match
    if (result.num_modern != result.num_legacy) {
        result.success = false;
        std::ostringstream oss;
        oss << "Count mismatch: modern=" << result.num_modern << " legacy=" << result.num_legacy;
        result.errors.push_back(oss.str());
    }

    // Check 2: all non-filtered residues should have both indices
    for (const auto& r : residues_) {
        if (!r.filtered) {
            if (r.modern_index < 0) {
                result.success = false;
                std::ostringstream oss;
                oss << "Non-filtered residue " << r.chain_id << r.residue_seq << r.insertion << " has no modern index";
                result.errors.push_back(oss.str());
                // Only report first 10 errors to avoid spam
                if (result.errors.size() >= 10)
                    break;
            }
            if (r.legacy_index < 0) {
                result.success = false;
                std::ostringstream oss;
                oss << "Non-filtered residue " << r.chain_id << r.residue_seq << r.insertion << " has no legacy index";
                result.errors.push_back(oss.str());
                // Only report first 10 errors to avoid spam
                if (result.errors.size() >= 10)
                    break;
            }
        }
    }

    // Check 3: filtered residues should NOT have modern index
    for (const auto& r : residues_) {
        if (r.filtered && r.modern_index >= 0) {
            result.success = false;
            std::ostringstream oss;
            oss << "Filtered residue " << r.chain_id << r.residue_seq << r.insertion << " has modern index "
                << r.modern_index << " (reason: " << r.filter_reason << ")";
            result.errors.push_back(oss.str());
            // Only report first 10 errors to avoid spam
            if (result.errors.size() >= 10)
                break;
        }
    }

    result.num_unmatched = result.num_modern - result.num_matched;

    return result;
}

std::string ResidueTracker::ValidationResult::to_string() const {
    std::ostringstream oss;
    oss << "\n=== Residue Index Validation ===\n";
    oss << "Status: " << (success ? "✅ PASS" : "❌ FAIL") << "\n";
    oss << "Residues read:    " << num_residues_read << "\n";
    oss << "Filtered out:     " << num_filtered << "\n";
    oss << "Modern indices:   " << num_modern << "\n";
    oss << "Legacy indices:   " << num_legacy << "\n";
    oss << "Matched:          " << num_matched << "\n";
    oss << "Unmatched:        " << num_unmatched << "\n";

    if (!errors.empty()) {
        oss << "\nErrors:\n";
        for (const auto& err : errors) {
            oss << "  - " << err << "\n";
        }
        if (errors.size() >= 10) {
            oss << "  ... (additional errors suppressed)\n";
        }
    }

    return oss.str();
}

void ResidueTracker::export_mapping(const std::string& output_path) const {
    json j = json::array();

    for (const auto& r : residues_) {
        json record;
        record["read_index"] = r.read_index;
        record["legacy_index"] = r.legacy_index;
        record["modern_index"] = r.modern_index;
        record["filtered"] = r.filtered;
        record["filter_reason"] = r.filter_reason;
        record["chain_id"] = r.chain_id;
        record["residue_seq"] = r.residue_seq;
        record["insertion"] = r.insertion;
        record["residue_name"] = r.residue_name;
        j.push_back(record);
    }

    std::ofstream file(output_path);
    if (file.is_open()) {
        file << std::setw(2) << j << std::endl;
        std::cout << "Exported index mapping to: " << output_path << std::endl;
    } else {
        std::cerr << "Error: Could not write to " << output_path << std::endl;
    }
}

std::optional<int> ResidueTracker::get_legacy_index(int modern_index) const {
    for (const auto& r : residues_) {
        if (r.modern_index == modern_index) {
            return r.legacy_index;
        }
    }
    return std::nullopt;
}

std::optional<int> ResidueTracker::get_modern_index(int legacy_index) const {
    for (const auto& r : residues_) {
        if (r.legacy_index == legacy_index) {
            return r.modern_index;
        }
    }
    return std::nullopt;
}

} // namespace x3dna

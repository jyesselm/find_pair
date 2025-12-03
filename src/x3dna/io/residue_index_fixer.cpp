/**
 * @file residue_index_fixer.cpp
 * @brief Implementation of residue index fixing from JSON
 */

#include <x3dna/io/residue_index_fixer.hpp>
#include <fstream>
#include <nlohmann/json.hpp>

using namespace x3dna;
using json = nlohmann::json;

// Residue key: (residue_name, chain_id, residue_seq, insertion)
using ResidueKey = std::tuple<std::string, char, int, char>;

namespace x3dna {
namespace io {

int fix_residue_indices_from_json(core::Structure& structure, const std::string& legacy_json_file) {
    // Step 1: Build map by PDB properties
    std::map<ResidueKey, core::Residue*> residues_by_pdb_props;

    for (auto& chain : structure.chains()) {
        for (auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                ResidueKey key = std::make_tuple(residue.name(), residue.chain_id(),
                                                 residue.seq_num(), residue.insertion());
                residues_by_pdb_props[key] = &residue;
            }
        }
    }

    // Step 2: Load legacy JSON
    std::ifstream json_file(legacy_json_file);
    if (!json_file.is_open()) {
        return -1; // Error opening file
    }

    json legacy_data;
    try {
        json_file >> legacy_data;
    } catch (const json::parse_error& e) {
        // JSON file may be corrupted/incomplete (e.g., missing closing bracket)
        // Try to read as string and fix common issues
        json_file.clear();
        json_file.seekg(0);
        std::string content((std::istreambuf_iterator<char>(json_file)),
                            std::istreambuf_iterator<char>());

        // Try to fix common issues: missing closing bracket
        // Remove trailing whitespace
        while (!content.empty() && (content.back() == ' ' || content.back() == '\n' ||
                                    content.back() == '\r' || content.back() == '\t')) {
            content.pop_back();
        }

        // Remove trailing comma if present
        while (!content.empty() && content.back() == ',') {
            content.pop_back();
            // Remove any whitespace after comma
            while (!content.empty() && (content.back() == ' ' || content.back() == '\n' ||
                                        content.back() == '\r' || content.back() == '\t')) {
                content.pop_back();
            }
        }

        // Add closing bracket if array is not closed
        if (!content.empty() && content[0] == '[' && content.back() != ']') {
            content += "\n]";
        }

        try {
            legacy_data = json::parse(content);
        } catch (const json::parse_error& e2) {
            // Still can't parse - return error
            return -3; // JSON parse error
        }
    }

    if (!legacy_data.is_array()) {
        return -2; // Not an array
    }

    // Step 3: Build legacy index map
    std::map<ResidueKey, int> legacy_idx_by_pdb_props;

    for (const auto& rec : legacy_data) {
        bool is_base_frame_calc = false;
        if (rec.contains("type")) {
            is_base_frame_calc = (rec["type"] == "base_frame_calc");
        } else {
            is_base_frame_calc = rec.contains("residue_idx");
        }

        if (is_base_frame_calc) {
            std::string residue_name = "";
            if (rec.contains("residue_name")) {
                residue_name = rec["residue_name"];
            } else if (rec.contains("base_type")) {
                std::string base_type = rec["base_type"];
                if (base_type == "A")
                    residue_name = "  A";
                else if (base_type == "C")
                    residue_name = "  C";
                else if (base_type == "G")
                    residue_name = "  G";
                else if (base_type == "U")
                    residue_name = "  U";
                else if (base_type == "T")
                    residue_name = "  T";
            }

            std::string chain_str = rec.value("chain_id", "");
            char chain_id = chain_str.empty() ? ' ' : chain_str[0];
            int residue_seq = rec.value("residue_seq", 0);
            std::string ins_str = rec.value("insertion", "");
            char insertion = ins_str.empty() ? ' ' : ins_str[0];
            int legacy_idx = rec.value("residue_idx", 0);

            if (legacy_idx > 0 && !residue_name.empty()) {
                ResidueKey key = std::make_tuple(residue_name, chain_id, residue_seq, insertion);
                legacy_idx_by_pdb_props[key] = legacy_idx;
            }
        }
    }

    // Step 4: Match and fix
    int matched_count = 0;

    for (const auto& [key, legacy_idx] : legacy_idx_by_pdb_props) {
        auto it = residues_by_pdb_props.find(key);
        if (it != residues_by_pdb_props.end()) {
            core::Residue* residue = it->second;

            // Set legacy_residue_idx on all atoms in this residue
            for (auto& atom : residue->atoms()) {
                atom.set_legacy_residue_idx(legacy_idx);
            }

            matched_count++;
        }
    }

    return matched_count;
}

} // namespace io
} // namespace x3dna

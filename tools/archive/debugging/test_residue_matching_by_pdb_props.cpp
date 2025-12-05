/**
 * @file test_residue_matching_by_pdb_props.cpp
 * @brief Test residue matching by PDB properties, then assign legacy indices from JSON
 *
 * This tool demonstrates a two-step approach:
 * 1. Match residues by PDB properties (chain_id, residue_seq, insertion, residue_name)
 * 2. Assign legacy indices from legacy JSON output
 *
 * This decouples residue matching from legacy index assignment, making debugging easier.
 *
 * Usage: test_residue_matching_by_pdb_props <pdb_file> <legacy_json_file>
 */

#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <nlohmann/json.hpp>

using namespace x3dna;
using json = nlohmann::json;

// Residue key: (residue_name, chain_id, residue_seq, insertion)
using ResidueKey = std::tuple<std::string, char, int, char>;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <legacy_json_file>\n";
        std::cerr << "Example: " << argv[0]
                  << " data/pdb/6CAQ.pdb data/json_legacy/base_frame_calc/6CAQ.json\n";
        return 1;
    }

    std::string pdb_file = argv[1];
    std::string legacy_json_file = argv[2];

    std::cout << "Testing Residue Matching by PDB Properties\n";
    std::cout << "=" << std::string(60, '=') << "\n";
    std::cout << "PDB file: " << pdb_file << "\n";
    std::cout << "Legacy JSON: " << legacy_json_file << "\n\n";

    // Step 1: Parse PDB and build map by PDB properties
    std::cout << "STEP 1: Parse PDB and match by PDB properties\n";
    std::cout << "-" << std::string(60, '-') << "\n";

    io::PdbParser parser;
    core::Structure structure = parser.parse_file(pdb_file);

    // Build map: (residue_name, chain_id, residue_seq, insertion) -> Residue*
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

    std::cout << "Parsed " << residues_by_pdb_props.size() << " residues from PDB\n\n";

    // Step 2: Load legacy JSON and assign legacy indices
    std::cout << "STEP 2: Load legacy JSON and assign legacy indices\n";
    std::cout << "-" << std::string(60, '-') << "\n";

    std::ifstream json_file(legacy_json_file);
    if (!json_file.is_open()) {
        std::cerr << "ERROR: Could not open JSON file: " << legacy_json_file << "\n";
        return 1;
    }

    json legacy_data;
    json_file >> legacy_data;

    if (!legacy_data.is_array()) {
        std::cerr << "ERROR: JSON file is not an array\n";
        return 1;
    }

    // Build map: (residue_name, chain_id, residue_seq, insertion) -> legacy_idx
    std::map<ResidueKey, int> legacy_idx_by_pdb_props;
    int matched_count = 0;
    int unmatched_count = 0;

    for (const auto& rec : legacy_data) {
        // Handle both records with 'type' field and records without
        bool is_base_frame_calc = false;
        if (rec.contains("type")) {
            is_base_frame_calc = (rec["type"] == "base_frame_calc");
        } else {
            // If no type field, assume it's a base_frame_calc record if it has residue_idx
            is_base_frame_calc = rec.contains("residue_idx");
        }

        if (is_base_frame_calc) {
            // Try different field names for residue_name
            std::string residue_name = "";
            if (rec.contains("residue_name")) {
                residue_name = rec["residue_name"];
            } else if (rec.contains("base_type")) {
                // base_type is single letter, need to convert to 3-letter
                std::string base_type = rec["base_type"];
                // Simple mapping (could be expanded)
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

    std::cout << "Loaded " << legacy_idx_by_pdb_props.size() << " legacy residue indices\n\n";

    // Step 3: Match and assign
    std::cout << "STEP 3: Match residues and assign legacy indices\n";
    std::cout << "-" << std::string(60, '-') << "\n";

    std::map<int, core::Residue*> residues_by_legacy_idx;

    for (const auto& [key, legacy_idx] : legacy_idx_by_pdb_props) {
        auto it = residues_by_pdb_props.find(key);
        if (it != residues_by_pdb_props.end()) {
            // Found match - assign legacy index
            core::Residue* residue = it->second;

            // Set legacy_residue_idx on all atoms in this residue
            for (auto& atom : residue->atoms()) {
                atom.set_legacy_residue_idx(legacy_idx);
            }

            residues_by_legacy_idx[legacy_idx] = residue;
            matched_count++;
        } else {
            unmatched_count++;
            auto [resname, chain, seq, ins] = key;
            std::cout << "⚠️  No match for legacy residue: " << resname << " Chain " << chain
                      << " Seq " << seq;
            if (ins != ' ')
                std::cout << " Ins '" << ins << "'";
            std::cout << " (legacy_idx=" << legacy_idx << ")\n";
        }
    }

    std::cout << "\nMatched: " << matched_count << " residues\n";
    std::cout << "Unmatched: " << unmatched_count << " residues\n\n";

    // Step 4: Test lookup by legacy index
    std::cout << "STEP 4: Test lookup by legacy index\n";
    std::cout << "-" << std::string(60, '-') << "\n";

    // Test a few specific indices
    std::vector<int> test_indices = {1102, 1127, 1, 100, 500};

    for (int idx : test_indices) {
        auto it = residues_by_legacy_idx.find(idx);
        if (it != residues_by_legacy_idx.end()) {
            core::Residue* res = it->second;
            std::cout << "Index " << idx << ": " << res->name() << " Chain " << res->chain_id()
                      << " Seq " << res->seq_num();
            if (res->insertion() != ' ') {
                std::cout << " Ins '" << res->insertion() << "'";
            }
            std::cout << "\n";
        } else {
            std::cout << "Index " << idx << ": Not found\n";
        }
    }

    std::cout << "\n✅ Test complete!\n";
    std::cout << "\nThis approach:\n";
    std::cout << "  1. Matches residues by PDB properties (reliable)\n";
    std::cout << "  2. Assigns legacy indices from JSON (decoupled)\n";
    std::cout << "  3. Makes debugging easier (clear separation)\n";

    return 0;
}

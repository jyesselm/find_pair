/**
 * @file fix_residue_indices_from_json.cpp
 * @brief Fix residue legacy indices by matching with legacy JSON
 * 
 * This tool uses the PDB properties matching approach to fix legacy_residue_idx
 * on all atoms in the structure by matching with legacy JSON output.
 * 
 * Usage: fix_residue_indices_from_json <pdb_file> <legacy_json_file> <output_pdb>
 */

#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/io/pdb_writer.hpp>
#include <nlohmann/json.hpp>

using namespace x3dna;
using json = nlohmann::json;

// Residue key: (residue_name, chain_id, residue_seq, insertion)
using ResidueKey = std::tuple<std::string, char, int, char>;

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <legacy_json_file> <output_pdb>\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/6CAQ.pdb data/json_legacy/base_frame_calc/6CAQ.json data/pdb/6CAQ_fixed.pdb\n";
        return 1;
    }

    std::string pdb_file = argv[1];
    std::string legacy_json_file = argv[2];
    std::string output_pdb = argv[3];

    std::cout << "Fixing Residue Indices from Legacy JSON\n";
    std::cout << "=" << std::string(60, '=') << "\n";
    std::cout << "PDB file: " << pdb_file << "\n";
    std::cout << "Legacy JSON: " << legacy_json_file << "\n";
    std::cout << "Output PDB: " << output_pdb << "\n\n";

    // Step 1: Parse PDB
    std::cout << "STEP 1: Parse PDB\n";
    std::cout << "-" << std::string(60, '-') << "\n";
    
    io::PdbParser parser;
    core::Structure structure = parser.parse_file(pdb_file);
    
    // Build map: (residue_name, chain_id, residue_seq, insertion) -> Residue*
    std::map<ResidueKey, core::Residue*> residues_by_pdb_props;
    
    for (auto& chain : structure.chains()) {
        for (auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                ResidueKey key = std::make_tuple(
                    residue.name(),
                    residue.chain_id(),
                    residue.seq_num(),
                    residue.insertion()
                );
                residues_by_pdb_props[key] = &residue;
            }
        }
    }
    
    std::cout << "Parsed " << residues_by_pdb_props.size() << " residues\n\n";

    // Step 2: Load legacy JSON
    std::cout << "STEP 2: Load legacy JSON\n";
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
                if (base_type == "A") residue_name = "  A";
                else if (base_type == "C") residue_name = "  C";
                else if (base_type == "G") residue_name = "  G";
                else if (base_type == "U") residue_name = "  U";
                else if (base_type == "T") residue_name = "  T";
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

    // Step 3: Match and fix indices
    std::cout << "STEP 3: Match residues and fix legacy indices\n";
    std::cout << "-" << std::string(60, '-') << "\n";
    
    int matched_count = 0;
    int unmatched_count = 0;
    
    for (const auto& [key, legacy_idx] : legacy_idx_by_pdb_props) {
        auto it = residues_by_pdb_props.find(key);
        if (it != residues_by_pdb_props.end()) {
            core::Residue* residue = it->second;
            
            // Set legacy_residue_idx on all atoms in this residue
            for (auto& atom : residue->atoms()) {
                atom.set_legacy_residue_idx(legacy_idx);
            }
            
            matched_count++;
        } else {
            unmatched_count++;
        }
    }
    
    std::cout << "Fixed indices for " << matched_count << " residues\n";
    if (unmatched_count > 0) {
        std::cout << "Warning: " << unmatched_count << " residues from JSON not found in PDB\n";
    }
    std::cout << "\n";

    // Step 4: Verify a few indices
    std::cout << "STEP 4: Verify fixed indices\n";
    std::cout << "-" << std::string(60, '-') << "\n";
    
    std::map<int, core::Residue*> residues_by_legacy_idx;
    for (auto& chain : structure.chains()) {
        for (auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx > 0) {
                    residues_by_legacy_idx[legacy_idx] = &residue;
                }
            }
        }
    }
    
    std::vector<int> test_indices = {1102, 1127, 1, 100};
    for (int idx : test_indices) {
        auto it = residues_by_legacy_idx.find(idx);
        if (it != residues_by_legacy_idx.end()) {
            core::Residue* res = it->second;
            std::cout << "Index " << idx << ": " << res->name() 
                      << " Chain " << res->chain_id() 
                      << " Seq " << res->seq_num() << "\n";
        }
    }
    std::cout << "\n";

    // Step 5: Write output (optional - for now just report)
    std::cout << "STEP 5: Structure ready with fixed indices\n";
    std::cout << "-" << std::string(60, '-') << "\n";
    std::cout << "Structure has " << structure.num_atoms() << " atoms\n";
    std::cout << "Structure has " << residues_by_legacy_idx.size() << " residues with legacy indices\n";
    std::cout << "\n✅ Indices fixed! Structure is ready for use.\n";
    std::cout << "\nNote: To save the structure, uncomment the PDB writer code.\n";
    
    // Uncomment to write output PDB:
    // io::PdbWriter writer;
    // writer.write_file(structure, output_pdb);
    // std::cout << "Written to: " << output_pdb << "\n";
    
    return 0;
}


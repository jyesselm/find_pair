/**
 * @file generate_residue_ordering_json.cpp
 * @brief Generate JSON file with residue ordering for comparison with legacy
 *
 * This tool generates a JSON file containing the residue ordering information
 * that can be directly compared with legacy output.
 *
 * Usage: generate_residue_ordering_json <pdb_file> <output_json>
 */

#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/structure.hpp>
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <iomanip>

using namespace x3dna::core;
using namespace x3dna::io;
using json = nlohmann::json;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <output_json>\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/3G8T.pdb data/residue_ordering/3G8T.json\n";
        return 1;
    }

    std::filesystem::path pdb_file = argv[1];
    std::filesystem::path output_json = argv[2];

    if (!std::filesystem::exists(pdb_file)) {
        std::cerr << "Error: PDB file not found: " << pdb_file << "\n";
        return 1;
    }

    try {
        // Parse PDB with legacy-compatible settings
        PdbParser parser;
        parser.set_include_hetatm(true);
        parser.set_include_waters(true);
        Structure structure = parser.parse_file(pdb_file);

        // Get residues in legacy order (using Structure's built-in method)
        auto residues = structure.residues_in_legacy_order();

        // Build JSON output
        json output;
        output["pdb_id"] = pdb_file.stem().string();
        output["total_residues"] = residues.size();
        output["residues"] = json::array();

        // Add each residue with its index and properties
        for (size_t i = 0; i < residues.size(); i++) {
            const Residue* res = residues[i];
            json residue_info;
            residue_info["legacy_index"] = static_cast<int>(i + 1); // 1-based
            residue_info["residue_name"] = res->name();
            residue_info["chain_id"] = std::string(1, res->chain_id());
            residue_info["residue_seq"] = res->seq_num();
            residue_info["insertion_code"] = std::string(1, res->insertion());
            residue_info["num_atoms"] = res->num_atoms();

            // Add first and last atom names for verification
            if (res->num_atoms() > 0) {
                residue_info["first_atom"] = res->atoms()[0].name();
                residue_info["last_atom"] = res->atoms()[res->num_atoms() - 1].name();
            }

            output["residues"].push_back(residue_info);
        }

        // Create output directory if needed
        std::filesystem::create_directories(output_json.parent_path());

        // Write JSON file
        std::ofstream out_file(output_json);
        if (!out_file.is_open()) {
            std::cerr << "Error: Cannot open output file: " << output_json << "\n";
            return 1;
        }

        out_file << std::setw(2) << output << "\n";

        std::cout << "Generated residue ordering JSON: " << output_json << "\n";
        std::cout << "Total residues: " << residues.size() << "\n";

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

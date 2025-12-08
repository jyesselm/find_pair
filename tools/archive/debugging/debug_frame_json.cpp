/**
 * @file debug_frame_json.cpp
 * @brief Generate debug JSON files for specific residues to compare frame calculations
 *
 * This tool generates JSON files for specific residues showing:
 * - base_frame_calc records
 * - ls_fitting records
 * - frame_calc records
 * - Frame origins and rotations
 *
 * Usage: debug_frame_json <pdb_file> <legacy_idx1> [legacy_idx2] [pdb_id]
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/io/json_writer.hpp>
#include <nlohmann/json.hpp>

using namespace x3dna;
using json = nlohmann::json;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <legacy_idx1> [legacy_idx2] [pdb_id]\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/6CAQ.pdb 1101 1127 6CAQ\n";
        std::cerr << "         " << argv[0] << " data/pdb/6CAQ.pdb 1101 6CAQ\n";
        return 1;
    }

    std::string pdb_file = argv[1];
    int target_idx1 = std::stoi(argv[2]);
    int target_idx2 = (argc > 3 && std::isdigit(argv[3][0])) ? std::stoi(argv[3]) : -1;
    std::string pdb_id = (argc > 3 && !std::isdigit(argv[argc - 1][0])) ? argv[argc - 1] : (argc > 4) ? argv[4] : "";

    if (pdb_id.empty()) {
        // Extract PDB ID from filename
        size_t last_slash = pdb_file.find_last_of("/\\");
        size_t last_dot = pdb_file.find_last_of(".");
        if (last_dot != std::string::npos && last_dot > last_slash) {
            pdb_id = pdb_file.substr(last_slash + 1, last_dot - last_slash - 1);
        }
    }

    std::cout << "Generating debug JSON for residues\n";
    std::cout << "PDB file: " << pdb_file << "\n";
    std::cout << "PDB ID: " << pdb_id << "\n";
    std::cout << "Target residue 1 (legacy_idx): " << target_idx1 << "\n";
    if (target_idx2 > 0) {
        std::cout << "Target residue 2 (legacy_idx): " << target_idx2 << "\n";
    }

    // Step 1: Parse PDB
    std::cout << "\nParsing PDB file...\n";
    io::PdbParser parser;
    core::Structure structure = parser.parse_file(pdb_file);

    // Build mapping from legacy_residue_idx to residue
    std::map<int, core::Residue*> residue_by_legacy_idx;

    for (auto& chain : structure.chains()) {
        for (auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx > 0) {
                    residue_by_legacy_idx[legacy_idx] = &residue;
                }
            }
        }
    }

    // Find target residues
    std::vector<int> target_indices = {target_idx1};
    if (target_idx2 > 0) {
        target_indices.push_back(target_idx2);
    }

    std::vector<core::Residue*> target_residues;
    for (int idx : target_indices) {
        auto it = residue_by_legacy_idx.find(idx);
        if (it == residue_by_legacy_idx.end()) {
            std::cerr << "ERROR: Residue at legacy_idx " << idx << " not found!\n";
            return 1;
        }
        target_residues.push_back(it->second);
        std::cout << "Found residue " << idx << ": " << it->second->name() << " (Chain " << it->second->chain_id()
                  << ", Seq " << it->second->seq_num() << ")\n";
    }

    // Step 2: Calculate frames
    std::cout << "\nCalculating frames...\n";
    std::filesystem::path template_path = "resources/templates";
    if (!std::filesystem::exists(template_path)) {
        template_path = "data/templates";
    }

    algorithms::BaseFrameCalculator calculator(template_path.string());
    std::vector<algorithms::FrameCalculationResult> frame_results;

    for (auto* residue : target_residues) {
        auto frame_result = calculator.calculate_frame(*residue);
        frame_results.push_back(frame_result);

        if (frame_result.is_valid) {
            std::cout << "Frame calculated for " << residue->name() << " (legacy_idx "
                      << residue->atoms()[0].legacy_residue_idx() << ")\n";
            std::cout << "  Origin: [" << frame_result.frame.origin().x() << ", " << frame_result.frame.origin().y()
                      << ", " << frame_result.frame.origin().z() << "]\n";
            std::cout << "  RMS fit: " << frame_result.rms_fit << "\n";
            std::cout << "  Matched atoms: " << frame_result.num_matched << "\n";
        } else {
            std::cerr << "ERROR: Frame calculation failed for " << residue->name() << "\n";
            return 1;
        }
    }

    // Step 3: Generate JSON
    std::cout << "\nGenerating JSON output...\n";

    json output;
    output["pdb_file"] = pdb_file;
    output["pdb_name"] = pdb_id;
    output["calculations"] = json::array();

    io::JsonWriter writer(pdb_file);

    for (size_t i = 0; i < target_residues.size(); ++i) {
        auto* residue = target_residues[i];
        const auto& frame_result = frame_results[i];

        int legacy_residue_idx = residue->atoms()[0].legacy_residue_idx();
        size_t record_idx = static_cast<size_t>(legacy_residue_idx - 1); // Convert to 0-based

        char base_type = residue->one_letter_code();

        // Record base_frame_calc
        writer.record_base_frame_calc(record_idx, base_type, frame_result.template_file, frame_result.rms_fit,
                                      frame_result.matched_atoms, residue->name(), residue->chain_id(),
                                      residue->seq_num(), residue->insertion());

        // Record ls_fitting
        writer.record_ls_fitting(record_idx, frame_result.num_matched, frame_result.rms_fit,
                                 frame_result.rotation_matrix, frame_result.translation, residue->name(),
                                 residue->chain_id(), residue->seq_num(), residue->insertion());

        // Create frame_calc record manually (need matched coordinates)
        // For now, we'll create a simplified version
        json frame_calc_record;
        frame_calc_record["type"] = "frame_calc";
        frame_calc_record["residue_idx"] = record_idx;
        frame_calc_record["legacy_residue_idx"] = legacy_residue_idx;
        frame_calc_record["base_type"] = std::string(1, base_type);
        frame_calc_record["residue_name"] = residue->name();
        frame_calc_record["chain_id"] = std::string(1, residue->chain_id());
        frame_calc_record["residue_seq"] = residue->seq_num();
        if (residue->insertion() != ' ') {
            frame_calc_record["insertion"] = std::string(1, residue->insertion());
        }
        frame_calc_record["template_file"] = frame_result.template_file.string();
        frame_calc_record["rms_fit"] = frame_result.rms_fit;
        frame_calc_record["num_matched_atoms"] = frame_result.num_matched;

        // Frame origin and rotation
        frame_calc_record["frame_origin"] = json::array(
            {frame_result.frame.origin().x(), frame_result.frame.origin().y(), frame_result.frame.origin().z()});

        auto rot = frame_result.frame.rotation();
        frame_calc_record["rotation_matrix"] = json::array({json::array({rot.at(0, 0), rot.at(0, 1), rot.at(0, 2)}),
                                                            json::array({rot.at(1, 0), rot.at(1, 1), rot.at(1, 2)}),
                                                            json::array({rot.at(2, 0), rot.at(2, 1), rot.at(2, 2)})});

        output["calculations"].push_back(frame_calc_record);
    }

    // Write split files to record-type-specific directories
    std::filesystem::path output_dir = "data/json/debug_" + pdb_id;
    std::filesystem::create_directories(output_dir);
    writer.write_split_files(output_dir);

    // Also write a separate debug JSON file with detailed frame info
    std::string debug_file = "data/json/debug_" + pdb_id + "_frames_detailed.json";
    std::ofstream out(debug_file);
    out << output.dump(2) << std::endl;
    std::cout << "\nDebug JSON written to: " << debug_file << "\n";
    std::cout << "Standard JSON written to: " << output_dir << "/base_frame_calc/" << pdb_id << ".json\n";

    return 0;
}

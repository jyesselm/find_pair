/**
 * @file generate_modern_json.cpp
 * @brief Standalone tool to generate modern JSON for a single PDB file
 */

#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/io/json_writer.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/io/frame_json_recorder.hpp>
#include <x3dna/algorithms/base_pair_finder.hpp>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;

// Helper function: Detect RNA structure
bool detect_rna_structure(const Structure& structure) {
    return BaseFrameCalculator::detect_rna(structure);
}

// Helper function: Setup frame calculator with RNA detection
BaseFrameCalculator setup_frame_calculator(const std::filesystem::path& template_path,
                                           const Structure& structure) {
    BaseFrameCalculator calculator(template_path);

    bool is_rna = detect_rna_structure(structure);
    calculator.set_is_rna(is_rna);

    if (is_rna) {
        std::cout << "Detected RNA structure (O2' atoms found)\n";
    } else {
        std::cout << "Detected DNA structure (no O2' atoms)\n";
    }

    return calculator;
}

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_pdb_file> <output_json_dir> [--stage=STAGE]\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/1ABC.pdb data/json\n";
        std::cerr << "         " << argv[0] << " data/pdb/1ABC.pdb data/json --stage=atoms\n";
        std::cerr << "\n";
        std::cerr << "Stages:\n";
        std::cerr << "  atoms       - Only generate atoms JSON (fast, for testing)\n";
        std::cerr << "  frames      - Only generate frame calculations\n";
        std::cerr << "  distances   - Only generate distance checks\n";
        std::cerr << "  hbonds      - Only generate H-bonds\n";
        std::cerr << "  validation  - Only generate pair validation\n";
        std::cerr << "  selection   - Only generate pair selection\n";
        std::cerr << "  steps       - Only generate step parameters\n";
        std::cerr << "  helical     - Only generate helical parameters\n";
        std::cerr << "  all         - Generate all stages (default)\n";
        std::cerr << "\n";
        std::cerr << "Output: Creates segmented JSON files in record-type-specific directories:\n";
        std::cerr << "  <output_json_dir>/pdb_atoms/<PDB_ID>.json\n";
        std::cerr << "  <output_json_dir>/base_frame_calc/<PDB_ID>.json\n";
        std::cerr << "  <output_json_dir>/base_pair/<PDB_ID>.json\n";
        std::cerr << "  etc.\n";
        return 1;
    }

    std::filesystem::path pdb_file = argv[1];
    std::filesystem::path json_output_dir = argv[2]; // Now expects a directory, not a file

    // Parse stage flag (default: "all")
    std::string stage = "all";

    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        if (arg.find("--stage=") == 0) {
            stage = arg.substr(8); // Extract stage name after "="
            if (stage.empty()) {
                std::cerr << "Error: --stage= requires a stage name\n";
                return 1;
            }
        } else {
            std::cerr << "Error: Unknown option: " << arg << "\n";
            std::cerr << "Use --stage=STAGE\n";
            return 1;
        }
    }

    // Validate stage
    const std::set<std::string> valid_stages = {
        "atoms",      "residue_indices", "ls_fitting", "frames",  "distances", "hbonds",
        "validation", "selection",       "steps",      "helical", "all"};
    if (valid_stages.find(stage) == valid_stages.end()) {
        std::cerr << "Error: Invalid stage: " << stage << "\n";
        std::cerr << "Valid stages: atoms, residue_indices, frames, distances, hbonds, validation, "
                     "selection, "
                     "steps, helical, all\n";
        return 1;
    }

    if (!std::filesystem::exists(pdb_file)) {
        std::cerr << "Error: PDB file not found: " << pdb_file << "\n";
        return 1;
    }

    try {
        // Create output directory if needed
        std::filesystem::create_directories(json_output_dir);

        // Extract PDB name from file
        std::string pdb_name = std::filesystem::path(pdb_file).stem().string();

        // Parse PDB file
        PdbParser parser;
        parser.set_include_hetatm(true);
        parser.set_include_waters(true);

        // Parse PDB file (atoms will have legacy indices assigned sequentially during parsing)
        Structure structure = parser.parse_file(pdb_file);

        // Set PDB ID on structure (needed for write_atoms_json)
        structure.set_pdb_id(pdb_name);

        std::cout << "Processing: " << pdb_name << " (stage: " << stage << ")\n";
        std::cout << "Input: " << pdb_file << "\n";
        std::cout << "Output: " << json_output_dir << "\n\n";

        // Stage 1: Atoms
        if (stage == "atoms" || stage == "all") {
            std::cout << "Stage 1: Writing atoms...\n";
            structure.write_atoms_json(json_output_dir);
            std::cout << "  ✅ pdb_atoms/" << pdb_name << ".json\n";
            std::cout << "     " << structure.num_atoms() << " atoms written\n\n";
        }

        // Stage 2: Residue indices (maps residues to atom ranges)
        if (stage == "residue_indices" || stage == "all") {
            // Create JSON writer
            JsonWriter writer(pdb_file);

            // Record residue indices (seidx) - maps residues to atom ranges
            std::cout << "Stage 2: Writing residue indices...\n";
            writer.record_residue_indices(structure);

            // Write the JSON files
            writer.write_split_files(json_output_dir, true);

            std::cout << "  ✅ residue_indices/" << pdb_name << ".json\n";
            std::cout << "     " << structure.num_residues() << " residues mapped\n\n";
        }

        // If only atoms or residue_indices stage, we're done
        if (stage == "atoms" || stage == "residue_indices") {
            std::cout << "✅ Success! Generated JSON for: " << pdb_name << " (stage: " << stage
                      << ")\n";
            return 0;
        }

        // Stage 3: LS Fitting (least-squares fitting data only)
        if (stage == "ls_fitting" || stage == "all") {
            std::cout << "Stage 3: Writing ls_fitting...\n";

            // Create JSON writer (no legacy JSON for PDB line caching)
            JsonWriter writer(pdb_file);

            // Record residue indices (needed for frames)
            writer.record_residue_indices(structure);

            // Setup frame calculator with RNA detection
            BaseFrameCalculator calculator = setup_frame_calculator("data/templates", structure);

            FrameJsonRecorder recorder(calculator);
            size_t records_count = recorder.record_ls_fitting(structure, writer);

            // Write split JSON files
            writer.write_split_files(json_output_dir, true);

            std::cout << "  ✅ ls_fitting/" << pdb_name << ".json\n";
            std::cout << "     " << records_count << " ls_fitting records calculated\n\n";
        }

        // Stage 4: Frames (base_frame_calc and frame_calc, but not ls_fitting)
        if (stage == "frames" || stage == "all") {
            std::cout << "Stage 4: Writing frames...\n";

            // Create JSON writer (no legacy JSON for PDB line caching)
            JsonWriter writer(pdb_file);

            // Record residue indices (needed for frames)
            writer.record_residue_indices(structure);

            // Setup frame calculator with RNA detection
            BaseFrameCalculator calculator = setup_frame_calculator("data/templates", structure);

            FrameJsonRecorder recorder(calculator);
            size_t base_frame_count = recorder.record_base_frame_calc(structure, writer);
            size_t frame_calc_count = recorder.record_frame_calc(structure, writer);

            // Write split JSON files
            writer.write_split_files(json_output_dir, true);

            std::cout << "  ✅ base_frame_calc/" << pdb_name << ".json\n";
            std::cout << "  ✅ frame_calc/" << pdb_name << ".json\n";
            std::cout << "     " << base_frame_count << " base_frame_calc records, "
                      << frame_calc_count << " frame_calc records\n\n";
        }

        // If only atoms, residue_indices, ls_fitting, or frames stage, we're done
        if (stage == "atoms" || stage == "residue_indices" || stage == "ls_fitting" ||
            stage == "frames") {
            std::cout << "✅ Success! Generated JSON for: " << pdb_name << " (stage: " << stage
                      << ")\n";
            return 0;
        }

        // Stages 4-8: Only process if not just atoms, residue_indices, ls_fitting, or frames stage
        // Note: For "all" stage, frames are already handled by stages 3 and 4 above
        if (stage != "atoms" && stage != "residue_indices" && stage != "ls_fitting" &&
            stage != "frames") {
            // Create JSON writer
            JsonWriter writer(pdb_file);

            // Record residue indices (seidx) - maps residues to atom ranges (needed for later
            // stages)
            writer.record_residue_indices(structure);

            // Calculate frames (needed for base pair finding and other stages)
            // Setup frame calculator with RNA detection
            BaseFrameCalculator calculator = setup_frame_calculator("data/templates", structure);
            calculator.calculate_all_frames(structure);

            // Find and record base pairs (with validation recording)
            // Note: Legacy residue indices are already set on atoms during PDB parsing
            BasePairFinder finder;
            auto base_pairs = finder.find_pairs_with_recording(structure, &writer);
            // Base pairs are already recorded via find_pairs_with_recording
            // But we also record the final base_pair records
            for (const auto& pair : base_pairs) {
                writer.record_base_pair(pair);
            }
            std::cout << "  Base pairs found: " << base_pairs.size() << "\n";

            // Write split JSON files to record-type-specific directories
            writer.write_split_files(json_output_dir, true);

            std::cout << "Successfully generated JSON files in: " << json_output_dir << "\n";
            std::cout << "  Atoms: " << structure.num_atoms() << "\n";
            std::cout << "  Total residues: " << structure.num_residues() << "\n";
            std::cout << "  Base pairs: " << base_pairs.size() << "\n";
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

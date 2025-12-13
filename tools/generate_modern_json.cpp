/**
 * @file generate_modern_json.cpp
 * @brief Standalone tool to generate modern JSON for PDB files
 *
 * Supports single PDB or batch processing with automatic progress saving.
 */

#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <set>
#include <nlohmann/json.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/io/json_writer.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/io/frame_json_recorder.hpp>
#include <x3dna/algorithms/base_pair_finder.hpp>
#include <x3dna/algorithms/parameter_calculator.hpp>
#include <x3dna/algorithms/helix_organizer.hpp>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;
using json = nlohmann::json;

// Progress tracking structure
struct Progress {
    std::string stage;
    std::string start_time;
    std::string last_update;
    int total_pdbs = 0;
    int processed = 0;
    int succeeded = 0;
    int failed = 0;
    std::vector<std::string> completed_pdbs;
    std::vector<std::string> failed_pdbs;
    std::vector<std::string> pending_pdbs;

    json to_json() const {
        return {{"stage", stage},
                {"start_time", start_time},
                {"last_update", last_update},
                {"total_pdbs", total_pdbs},
                {"processed", processed},
                {"succeeded", succeeded},
                {"failed", failed},
                {"completed_pdbs", completed_pdbs},
                {"failed_pdbs", failed_pdbs},
                {"pending_pdbs", pending_pdbs}};
    }

    static Progress from_json(const json& j) {
        Progress p;
        p.stage = j.value("stage", "all");
        p.start_time = j.value("start_time", "");
        p.last_update = j.value("last_update", "");
        p.total_pdbs = j.value("total_pdbs", 0);
        p.processed = j.value("processed", 0);
        p.succeeded = j.value("succeeded", 0);
        p.failed = j.value("failed", 0);
        p.completed_pdbs = j.value("completed_pdbs", std::vector<std::string>{});
        p.failed_pdbs = j.value("failed_pdbs", std::vector<std::string>{});
        p.pending_pdbs = j.value("pending_pdbs", std::vector<std::string>{});
        return p;
    }
};

// Get current timestamp
std::string get_timestamp() {
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

// Save progress to file
void save_progress(const Progress& progress, const std::filesystem::path& progress_file) {
    std::ofstream f(progress_file);
    if (f.is_open()) {
        f << std::setw(2) << progress.to_json() << std::endl;
    }
}

// Load progress from file
Progress load_progress(const std::filesystem::path& progress_file) {
    if (!std::filesystem::exists(progress_file)) {
        return Progress{};
    }
    std::ifstream f(progress_file);
    if (f.is_open()) {
        json j = json::parse(f);
        return Progress::from_json(j);
    }
    return Progress{};
}

// Helper function: Detect RNA structure
bool detect_rna_structure(const Structure& structure) {
    return BaseFrameCalculator::detect_rna(structure);
}

// Helper function: Extract backbone O3' and P coordinates for 5'→3' direction checking
BackboneData extract_backbone_data(const Structure& structure) {
    BackboneData backbone;

    // Get residues in legacy order
    auto residues = structure.residues_in_legacy_order();

    for (size_t i = 0; i < residues.size(); ++i) {
        const Residue* residue = residues[i];
        if (!residue)
            continue;

        // Get legacy residue index (1-based)
        size_t legacy_idx = i + 1;

        BackboneAtoms atoms;

        // Find O3' atom
        auto o3_prime = residue->find_atom(" O3'");
        if (o3_prime.has_value()) {
            atoms.O3_prime = o3_prime->position();
        }

        // Find P atom
        auto p_atom = residue->find_atom(" P  ");
        if (p_atom.has_value()) {
            atoms.P = p_atom->position();
        }

        // Only add if we have at least one backbone atom
        if (atoms.O3_prime.has_value() || atoms.P.has_value()) {
            backbone[legacy_idx] = atoms;
        }
    }

    return backbone;
}

// Helper function: Setup frame calculator with RNA detection
BaseFrameCalculator setup_frame_calculator(const std::filesystem::path& template_path, const Structure& structure,
                                           bool verbose = true) {
    BaseFrameCalculator calculator(template_path);

    bool is_rna = detect_rna_structure(structure);
    calculator.set_is_rna(is_rna);

    if (verbose) {
        if (is_rna) {
            std::cout << "  Detected RNA structure (O2' atoms found)\n";
        } else {
            std::cout << "  Detected DNA structure (no O2' atoms)\n";
        }
    }

    return calculator;
}

// Process a single PDB file
bool process_single_pdb(const std::filesystem::path& pdb_file, const std::filesystem::path& json_output_dir,
                        const std::string& stage, bool verbose = true) {
    try {
        // Create output directory if needed
        std::filesystem::create_directories(json_output_dir);

        // Extract PDB name from file
        std::string pdb_name = std::filesystem::path(pdb_file).stem().string();

        // Parse PDB file
        PdbParser parser;
        parser.set_include_hetatm(true);
        parser.set_include_waters(true);

        Structure structure = parser.parse_file(pdb_file);
        structure.set_pdb_id(pdb_name);

        if (verbose) {
            std::cout << "Processing: " << pdb_name << " (stage: " << stage << ")\n";
        }

        // Stage 1: Atoms
        if (stage == "atoms" || stage == "all") {
            structure.write_atoms_json(json_output_dir);
            if (verbose) {
                std::cout << "  ✅ pdb_atoms/" << pdb_name << ".json (" << structure.num_atoms() << " atoms)\n";
            }
        }

        // Stage 2: Residue indices
        if (stage == "residue_indices" || stage == "all") {
            JsonWriter writer(pdb_file);
            writer.record_residue_indices(structure);
            writer.write_split_files(json_output_dir, true);
            if (verbose) {
                std::cout << "  ✅ residue_indices/" << pdb_name << ".json (" << structure.num_residues()
                          << " residues)\n";
            }
        }

        if (stage == "atoms" || stage == "residue_indices") {
            return true;
        }

        // Stage 3: LS Fitting
        if (stage == "ls_fitting" || stage == "all") {
            JsonWriter writer(pdb_file);
            writer.record_residue_indices(structure);
            BaseFrameCalculator calculator = setup_frame_calculator("data/templates", structure, verbose);
            FrameJsonRecorder recorder(calculator);
            size_t records_count = recorder.record_ls_fitting(structure, writer);
            writer.write_split_files(json_output_dir, true);
            if (verbose) {
                std::cout << "  ✅ ls_fitting/" << pdb_name << ".json (" << records_count << " records)\n";
            }
        }

        // Stage 4: Frames
        if (stage == "frames" || stage == "all") {
            JsonWriter writer(pdb_file);
            writer.record_residue_indices(structure);
            BaseFrameCalculator calculator = setup_frame_calculator("data/templates", structure, verbose);
            FrameJsonRecorder recorder(calculator);
            size_t base_frame_count = recorder.record_base_frame_calc(structure, writer);
            size_t frame_calc_count = recorder.record_frame_calc(structure, writer);
            writer.write_split_files(json_output_dir, true);
            if (verbose) {
                std::cout << "  ✅ base_frame_calc/" << pdb_name << ".json (" << base_frame_count << " records)\n";
                std::cout << "  ✅ frame_calc/" << pdb_name << ".json (" << frame_calc_count << " records)\n";
            }
        }

        if (stage == "atoms" || stage == "residue_indices" || stage == "ls_fitting" || stage == "frames") {
            return true;
        }

        // Stages 4-10: Full pair finding
        if (stage != "atoms" && stage != "residue_indices" && stage != "ls_fitting" && stage != "frames") {
            JsonWriter writer(pdb_file);
            writer.record_residue_indices(structure);
            BaseFrameCalculator calculator = setup_frame_calculator("data/templates", structure, verbose);
            calculator.calculate_all_frames(structure);

            BasePairFinder finder;
            auto base_pairs = finder.find_pairs_with_recording(structure, &writer);
            for (const auto& pair : base_pairs) {
                writer.record_base_pair(pair);
            }

            // Stages 11-12: Step and helical parameters
            // Legacy iterates through pairs in backbone connectivity order (5' to 3'),
            // NOT sequential base pair order. HelixOrganizer determines this order.
            if (stage == "all" || stage == "steps" || stage == "helical") {
                if (base_pairs.size() >= 2) {
                    // Get helix order from HelixOrganizer
                    BackboneData backbone = extract_backbone_data(structure);
                    HelixOrganizer organizer;
                    auto helix_order = organizer.organize(base_pairs, backbone);

                    ParameterCalculator param_calc;
                    size_t valid_steps = 0;

                    // Calculate step params following backbone connectivity order
                    // This matches legacy behavior - pairs are in five2three order
                    for (size_t i = 0; i + 1 < helix_order.pair_order.size(); ++i) {
                        size_t idx1 = helix_order.pair_order[i];
                        size_t idx2 = helix_order.pair_order[i + 1];
                        const auto& pair1 = base_pairs[idx1];
                        const auto& pair2 = base_pairs[idx2];

                        // Verify both pairs have frames
                        if (!pair1.frame1().has_value() || !pair1.frame2().has_value() ||
                            !pair2.frame1().has_value() || !pair2.frame2().has_value()) {
                            continue;
                        }

                        // Legacy uses frame1 (org_i, orien_i) directly for step params
                        // Verified: legacy mst_org = (bp1.org_i + bp2.org_i) / 2
                        auto bp1_frame = pair1.frame1().value();
                        auto bp2_frame = pair2.frame1().value();

                        // Calculate step parameters between org_i frames
                        auto step_params = param_calc.calculate_step_parameters(bp1_frame, bp2_frame);
                        // Use 1-based helix position indices (matching legacy)
                        size_t bp_idx1 = i + 1;
                        size_t bp_idx2 = i + 2;
                        writer.record_bpstep_params(bp_idx1, bp_idx2, step_params);

                        // Calculate helical parameters using same frames
                        auto helical_params = param_calc.calculate_helical_parameters_impl(bp1_frame, bp2_frame);
                        writer.record_helical_params(bp_idx1, bp_idx2, helical_params);

                        valid_steps++;
                    }

                    if (verbose) {
                        size_t total_steps = base_pairs.size() - 1;
                        std::cout << "  ✅ Generated step/helical params (" << valid_steps << "/" << total_steps
                                  << " steps)\n";
                    }
                }
            }

            writer.write_split_files(json_output_dir, true);

            if (verbose) {
                std::cout << "  ✅ Generated all JSON files (" << base_pairs.size() << " base pairs)\n";
            }
        }

        return true;
    } catch (const std::exception& e) {
        std::cerr << "  ❌ Error: " << e.what() << "\n";
        return false;
    }
}

// Load PDB list from file
std::vector<std::string> load_pdb_list(const std::filesystem::path& list_file) {
    std::vector<std::string> pdbs;
    std::ifstream f(list_file);
    std::string line;
    while (std::getline(f, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\n\r"));
        line.erase(line.find_last_not_of(" \t\n\r") + 1);
        if (!line.empty() && line[0] != '#') {
            pdbs.push_back(line);
        }
    }
    return pdbs;
}

// Get PDB files from directory
std::vector<std::string> get_pdbs_from_dir(const std::filesystem::path& pdb_dir) {
    std::vector<std::string> pdbs;
    for (const auto& entry : std::filesystem::directory_iterator(pdb_dir)) {
        if (entry.path().extension() == ".pdb") {
            pdbs.push_back(entry.path().stem().string());
        }
    }
    std::sort(pdbs.begin(), pdbs.end());
    return pdbs;
}

void print_usage(const char* prog_name) {
    std::cerr << "Usage:\n";
    std::cerr << "  Single PDB:\n";
    std::cerr << "    " << prog_name << " <input.pdb> <output_dir> [--stage=STAGE]\n\n";
    std::cerr << "  Multiple PDBs:\n";
    std::cerr << "    " << prog_name << " --pdb-list=<file.txt> --pdb-dir=<pdb_dir> <output_dir> [options]\n";
    std::cerr << "    " << prog_name << " --all-pdbs --pdb-dir=<pdb_dir> <output_dir> [options]\n\n";
    std::cerr << "Options:\n";
    std::cerr << "  --stage=STAGE       Stage to generate (atoms, frames, all, etc.)\n";
    std::cerr << "  --pdb-list=FILE     File with PDB IDs (one per line)\n";
    std::cerr << "  --pdb-dir=DIR       Directory containing PDB files\n";
    std::cerr << "  --all-pdbs          Process all PDBs in pdb-dir\n";
    std::cerr << "  --progress=FILE     Progress file (default: <output_dir>/progress.json)\n";
    std::cerr << "  --resume            Resume from progress file\n";
    std::cerr << "  --max=N             Maximum PDBs to process\n";
    std::cerr << "  --quiet             Less verbose output\n\n";
    std::cerr << "Stages:\n";
    std::cerr << "  atoms, residue_indices, ls_fitting, frames, distances,\n";
    std::cerr << "  hbonds, validation, selection, steps, helical, all\n\n";
    std::cerr << "Examples:\n";
    std::cerr << "  " << prog_name << " data/pdb/1EHZ.pdb data/json --stage=atoms\n";
    std::cerr << "  " << prog_name << " --pdb-list=fast_pdbs.txt --pdb-dir=data/pdb data/json --stage=frames\n";
    std::cerr << "  " << prog_name << " --all-pdbs --pdb-dir=data/pdb data/json --resume --max=100\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    // Parse arguments
    std::string stage = "all";
    std::string pdb_list_file;
    std::string pdb_dir = "data/pdb";
    std::string output_dir;
    std::string progress_file;
    bool all_pdbs = false;
    bool resume = false;
    bool quiet = false;
    int max_pdbs = -1;
    std::string single_pdb_file;

    std::vector<std::string> positional_args;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg.find("--stage=") == 0) {
            stage = arg.substr(8);
        } else if (arg.find("--pdb-list=") == 0) {
            pdb_list_file = arg.substr(11);
        } else if (arg.find("--pdb-dir=") == 0) {
            pdb_dir = arg.substr(10);
        } else if (arg.find("--progress=") == 0) {
            progress_file = arg.substr(11);
        } else if (arg.find("--max=") == 0) {
            max_pdbs = std::stoi(arg.substr(6));
        } else if (arg == "--all-pdbs") {
            all_pdbs = true;
        } else if (arg == "--resume") {
            resume = true;
        } else if (arg == "--quiet" || arg == "-q") {
            quiet = true;
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if (arg[0] != '-') {
            positional_args.push_back(arg);
        } else {
            std::cerr << "Error: Unknown option: " << arg << "\n";
            return 1;
        }
    }

    // Validate stage
    const std::set<std::string> valid_stages = {"atoms",     "residue_indices", "ls_fitting", "frames",
                                                "distances", "hbonds",          "validation", "selection",
                                                "steps",     "helical",         "all"};
    if (valid_stages.find(stage) == valid_stages.end()) {
        std::cerr << "Error: Invalid stage: " << stage << "\n";
        return 1;
    }

    // Determine mode: single PDB or batch
    bool batch_mode = all_pdbs || !pdb_list_file.empty();

    if (!batch_mode) {
        // Single PDB mode
        if (positional_args.size() < 2) {
            std::cerr << "Error: Missing arguments for single PDB mode\n";
            print_usage(argv[0]);
            return 1;
        }
        single_pdb_file = positional_args[0];
        output_dir = positional_args[1];

        if (!std::filesystem::exists(single_pdb_file)) {
            std::cerr << "Error: PDB file not found: " << single_pdb_file << "\n";
            return 1;
        }

        std::cout << "Processing: " << single_pdb_file << " (stage: " << stage << ")\n";
        std::cout << "Input: " << single_pdb_file << "\n";
        std::cout << "Output: " << output_dir << "\n\n";

        bool success = process_single_pdb(single_pdb_file, output_dir, stage, !quiet);

        if (success) {
            std::cout << "\n✅ Success!\n";
            return 0;
        } else {
            return 1;
        }
    }

    // Batch mode
    if (positional_args.empty()) {
        std::cerr << "Error: Missing output directory\n";
        return 1;
    }
    output_dir = positional_args[0];

    // Set default progress file
    if (progress_file.empty()) {
        progress_file = output_dir + "/progress.json";
    }

    // Get list of PDBs to process
    std::vector<std::string> pdb_ids;
    if (!pdb_list_file.empty()) {
        pdb_ids = load_pdb_list(pdb_list_file);
    } else if (all_pdbs) {
        pdb_ids = get_pdbs_from_dir(pdb_dir);
    }

    if (pdb_ids.empty()) {
        std::cerr << "Error: No PDBs found to process\n";
        return 1;
    }

    // Apply max limit
    if (max_pdbs > 0 && static_cast<int>(pdb_ids.size()) > max_pdbs) {
        pdb_ids.resize(max_pdbs);
    }

    // Initialize or load progress
    Progress progress;
    std::set<std::string> completed_set;

    if (resume && std::filesystem::exists(progress_file)) {
        progress = load_progress(progress_file);
        completed_set.insert(progress.completed_pdbs.begin(), progress.completed_pdbs.end());
        std::cout << "Resuming from progress file: " << progress_file << "\n";
        std::cout << "  Previously completed: " << progress.completed_pdbs.size() << "\n";
        std::cout << "  Previously failed: " << progress.failed_pdbs.size() << "\n\n";
    } else {
        progress.stage = stage;
        progress.start_time = get_timestamp();
        progress.total_pdbs = pdb_ids.size();
        progress.pending_pdbs = pdb_ids;
    }

    // Create output directory
    std::filesystem::create_directories(output_dir);

    std::cout << "Batch processing: " << pdb_ids.size() << " PDBs (stage: " << stage << ")\n";
    std::cout << "PDB directory: " << pdb_dir << "\n";
    std::cout << "Output directory: " << output_dir << "\n";
    std::cout << "Progress file: " << progress_file << "\n\n";

    int processed = 0;
    int succeeded = 0;
    int failed = 0;
    int skipped = 0;

    for (size_t i = 0; i < pdb_ids.size(); i++) {
        const auto& pdb_id = pdb_ids[i];

        // Skip if already completed
        if (completed_set.count(pdb_id)) {
            skipped++;
            continue;
        }

        std::filesystem::path pdb_path = std::filesystem::path(pdb_dir) / (pdb_id + ".pdb");

        if (!std::filesystem::exists(pdb_path)) {
            if (!quiet) {
                std::cout << "[" << (i + 1) << "/" << pdb_ids.size() << "] " << pdb_id << ": SKIP (file not found)\n";
            }
            progress.failed_pdbs.push_back(pdb_id);
            failed++;
            continue;
        }

        if (!quiet) {
            std::cout << "[" << (i + 1) << "/" << pdb_ids.size() << "] " << pdb_id << "...\n";
        }

        bool success = process_single_pdb(pdb_path, output_dir, stage, !quiet);

        processed++;
        if (success) {
            succeeded++;
            progress.completed_pdbs.push_back(pdb_id);
            completed_set.insert(pdb_id);
            if (!quiet) {
                std::cout << "  ✅ Done\n";
            }
        } else {
            failed++;
            progress.failed_pdbs.push_back(pdb_id);
        }

        // Update and save progress
        progress.processed = progress.completed_pdbs.size() + progress.failed_pdbs.size();
        progress.succeeded = progress.completed_pdbs.size();
        progress.failed = progress.failed_pdbs.size();
        progress.last_update = get_timestamp();

        // Remove from pending
        progress.pending_pdbs.erase(std::remove(progress.pending_pdbs.begin(), progress.pending_pdbs.end(), pdb_id),
                                    progress.pending_pdbs.end());

        // Save progress after each PDB
        save_progress(progress, progress_file);
    }

    // Final summary
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "BATCH PROCESSING COMPLETE\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "Total PDBs: " << pdb_ids.size() << "\n";
    std::cout << "Processed: " << processed << "\n";
    std::cout << "Succeeded: " << succeeded << "\n";
    std::cout << "Failed: " << failed << "\n";
    std::cout << "Skipped (already done): " << skipped << "\n";
    std::cout << "Progress saved to: " << progress_file << "\n";

    return (failed > 0) ? 1 : 0;
}

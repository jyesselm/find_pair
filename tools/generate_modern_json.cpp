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

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " <input_pdb_file> <output_json_file> [--legacy]\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/1ABC.pdb data/json/1ABC.json\n";
        std::cerr << "Options:\n";
        std::cerr << "  --legacy    Exclude C4 atom from matching (matches legacy behavior)\n";
        return 1;
    }

    std::filesystem::path pdb_file = argv[1];
    std::filesystem::path json_file = argv[2];

    // Check for legacy mode flag
    bool legacy_mode = false;
    if (argc == 4) {
        std::string flag = argv[3];
        if (flag == "--legacy") {
            legacy_mode = true;
            std::cout << "Legacy compatibility mode enabled (C4 atom excluded)\n";
        } else {
            std::cerr << "Error: Unknown option: " << flag << "\n";
            std::cerr << "Use --legacy to enable legacy compatibility mode\n";
            return 1;
        }
    }

    if (!std::filesystem::exists(pdb_file)) {
        std::cerr << "Error: PDB file not found: " << pdb_file << "\n";
        return 1;
    }

    try {
        // Create output directory if needed
        std::filesystem::create_directories(json_file.parent_path());

        // Parse PDB file
        PdbParser parser;
        parser.set_include_hetatm(true);
        parser.set_include_waters(true);

        // Parse PDB file (atoms will have legacy indices assigned sequentially during parsing)
        Structure structure = parser.parse_file(pdb_file);

        // Try to find legacy JSON file for PDB line caching (optional)
        std::filesystem::path legacy_json_file;
        std::filesystem::path legacy_dir = json_file.parent_path().parent_path() / "json_legacy";
        std::filesystem::path legacy_file = legacy_dir / json_file.filename();
        if (std::filesystem::exists(legacy_file)) {
            legacy_json_file = legacy_file;
        }

        // Create JSON writer (with optional legacy JSON for PDB line caching)
        JsonWriter writer(pdb_file, legacy_json_file);

        // Record PDB atoms
        writer.record_pdb_atoms(structure);

        // Calculate frames and record frame calculations
        BaseFrameCalculator calculator("data/templates");

        // Set legacy mode if requested
        calculator.set_legacy_mode(legacy_mode);

        // Detect RNA by checking for O2' atoms (RNA indicator)
        bool is_rna = false;
        for (const auto& chain : structure.chains()) {
            for (const auto& residue : chain.residues()) {
                if (residue.find_atom(" O2'").has_value() ||
                    residue.find_atom(" O2*").has_value()) {
                    is_rna = true;
                    break;
                }
            }
            if (is_rna)
                break;
        }
        calculator.set_is_rna(is_rna);

        if (is_rna) {
            std::cout << "Detected RNA structure (O2' atoms found)\n";
        } else {
            std::cout << "Detected DNA structure (no O2' atoms)\n";
        }

        // DEBUG: Verify template path
        std::filesystem::path template_path = calculator.template_path();
        std::cout << "Template path: " << template_path << "\n";

        // Verify templates exist (need to access through templates_ member, but it's private)
        // For now, just check if path exists
        std::filesystem::path template_dir = calculator.template_path();
        std::cout << "Checking templates in: " << template_dir << "\n";
        if (!std::filesystem::exists(template_dir)) {
            std::cerr << "WARNING: Template directory does not exist: " << template_dir << "\n";
        } else {
            std::filesystem::path template_a = template_dir / "Atomic_A.pdb";
            std::filesystem::path template_c = template_dir / "Atomic_C.pdb";
            std::filesystem::path template_g = template_dir / "Atomic_G.pdb";
            std::filesystem::path template_t = template_dir / "Atomic_T.pdb";

            if (!std::filesystem::exists(template_a)) {
                std::cerr << "WARNING: Template Atomic_A.pdb not found at: " << template_a << "\n";
            }
            if (!std::filesystem::exists(template_c)) {
                std::cerr << "WARNING: Template Atomic_C.pdb not found at: " << template_c << "\n";
            }
            if (!std::filesystem::exists(template_g)) {
                std::cerr << "WARNING: Template Atomic_G.pdb not found at: " << template_g << "\n";
            }
            if (!std::filesystem::exists(template_t)) {
                std::cerr << "WARNING: Template Atomic_T.pdb not found at: " << template_t << "\n";
            }
        }

        calculator.calculate_all_frames(structure);

        // Record frame calculations for each residue
        // Calculate frames directly for all nucleotide residues (don't rely on has_reference_frame
        // check) Legacy residue_idx is 1-based and counts ALL residues (including amino acids,
        // etc.)
        size_t residue_idx = 1;
        size_t frames_recorded = 0;
        size_t nucleotides_found = 0;
        size_t template_load_failures = 0;
        size_t atom_match_failures = 0;
        size_t unknown_residue_types = 0;

        for (const auto& chain : structure.chains()) {
            for (const auto& residue : chain.residues()) {
                // Check residue type
                ResidueType res_type = residue.residue_type();
                char one_letter = residue.one_letter_code();

                if (res_type == ResidueType::UNKNOWN) {
                    unknown_residue_types++;
                    // DEBUG: Log unknown residues that look like nucleotides
                    if (one_letter == 'A' || one_letter == 'C' || one_letter == 'G' ||
                        one_letter == 'T' || one_letter == 'U') {
                        std::cerr << "WARNING: Residue " << residue.name()
                                  << " has one_letter_code=" << one_letter
                                  << " but residue_type=UNKNOWN\n";
                    }
                }

                // Only process nucleotide residues
                // Include standard nucleotides and modified nucleotides (detected by ring atoms)
                bool is_nucleotide =
                    (res_type != ResidueType::UNKNOWN && res_type != ResidueType::AMINO_ACID);

                // Check for modified nucleotides that have ring atoms but aren't in standard list
                if (!is_nucleotide && res_type == ResidueType::UNKNOWN) {
                    // Check for ring atoms (similar to legacy residue_ident)
                    static const std::vector<std::string> common_ring_atoms = {
                        " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "};
                    int ring_atom_count = 0;
                    for (const auto& atom_name : common_ring_atoms) {
                        for (const auto& atom : residue.atoms()) {
                            if (atom.name() == atom_name) {
                                ring_atom_count++;
                                break;
                            }
                        }
                    }
                    // If has >= 3 ring atoms, treat as nucleotide (like legacy)
                    if (ring_atom_count >= 3) {
                        is_nucleotide = true;
                    }
                }

                if (is_nucleotide) {
                    nucleotides_found++;

                    // Calculate frame directly (don't check has_reference_frame)
                    FrameCalculationResult frame_result = calculator.calculate_frame_const(residue);

                    if (frame_result.is_valid) {
                        char base_type = residue.one_letter_code();

                        // Record base_frame_calc
                        // JsonWriter expects 0-based index, but legacy uses 1-based
                        writer.record_base_frame_calc(
                            residue_idx - 1, // Convert 1-based to 0-based
                            base_type, frame_result.template_file, frame_result.rms_fit,
                            frame_result.matched_atoms, residue.name(), residue.chain_id(),
                            residue.seq_num(), residue.insertion());

                        // Record ls_fitting
                        writer.record_ls_fitting(
                            residue_idx - 1, // Convert 1-based to 0-based
                            frame_result.num_matched, frame_result.rms_fit,
                            frame_result.rotation_matrix, frame_result.translation, residue.name(),
                            residue.chain_id(), residue.seq_num(), residue.insertion());

                        // Record frame_calc
                        // Extract coordinates from matched atoms - we need to get them from the
                        // frame result For now, use empty vectors since we don't have direct access
                        // to matched atom coords The frame_calc record is less critical than
                        // base_frame_calc and ls_fitting
                        std::vector<Vector3D> standard_coords, experimental_coords;
                        // TODO: Extract coordinates from matched atoms if needed

                        writer.record_frame_calc(
                            residue_idx - 1, // Convert 1-based to 0-based
                            base_type, frame_result.template_file, frame_result.rms_fit,
                            standard_coords, experimental_coords, residue.name(),
                            residue.chain_id(), residue.seq_num(), residue.insertion());

                        frames_recorded++;
                    } else {
                        // DEBUG: Log why it failed
                        std::cerr << "Frame calculation failed for " << residue.name() << " "
                                  << residue.chain_id() << ":" << residue.seq_num()
                                  << " (type=" << static_cast<int>(res_type)
                                  << ", one_letter=" << one_letter << ")\n";
                        if (frame_result.template_file.empty()) {
                            template_load_failures++;
                        } else if (frame_result.num_matched < 3) {
                            atom_match_failures++;
                            std::cerr << "  Only matched " << frame_result.num_matched
                                      << " atoms\n";
                        }
                    }
                }

                // Count all residues (to match legacy residue_idx behavior)
                residue_idx++;
            }
        }

        // Write JSON to file
        writer.write_to_file(json_file);

        std::cout << "Successfully generated JSON: " << json_file << "\n";
        std::cout << "  Atoms: " << structure.num_atoms() << "\n";
        std::cout << "  Total residues: " << residue_idx - 1 << "\n";
        std::cout << "  Nucleotides found: " << nucleotides_found << "\n";
        std::cout << "  Frames calculated: " << frames_recorded << "\n";
        std::cout << "  Unknown residue types: " << unknown_residue_types << "\n";
        if (template_load_failures > 0) {
            std::cout << "  Template load failures: " << template_load_failures << "\n";
        }
        if (atom_match_failures > 0) {
            std::cout << "  Atom match failures: " << atom_match_failures << "\n";
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

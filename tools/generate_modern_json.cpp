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
#include <x3dna/io/residue_index_fixer.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/base_pair_finder.hpp>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 5) {
        std::cerr << "Usage: " << argv[0] << " <input_pdb_file> <output_json_dir> [--legacy] [--fix-indices]\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/1ABC.pdb data/json\n";
        std::cerr << "Options:\n";
        std::cerr << "  --legacy       Exclude C4 atom from matching (matches legacy behavior)\n";
        std::cerr << "  --fix-indices  Fix residue indices from legacy JSON (for comparison)\n";
        std::cerr << "\n";
        std::cerr << "Output: Creates segmented JSON files in record-type-specific directories:\n";
        std::cerr << "  <output_json_dir>/pdb_atoms/<PDB_ID>.json\n";
        std::cerr << "  <output_json_dir>/base_frame_calc/<PDB_ID>.json\n";
        std::cerr << "  <output_json_dir>/base_pair/<PDB_ID>.json\n";
        std::cerr << "  etc.\n";
        return 1;
    }

    std::filesystem::path pdb_file = argv[1];
    std::filesystem::path json_output_dir = argv[2];  // Now expects a directory, not a file

    // Check for options
    bool legacy_mode = false;
    bool fix_indices = false;
    for (int i = 3; i < argc; i++) {
        std::string flag = argv[i];
        if (flag == "--legacy") {
            legacy_mode = true;
            std::cout << "Legacy compatibility mode enabled (C4 atom excluded)\n";
        } else if (flag == "--fix-indices") {
            fix_indices = true;
            std::cout << "Fix indices mode enabled (will match legacy JSON indices)\n";
        } else {
            std::cerr << "Error: Unknown option: " << flag << "\n";
            std::cerr << "Use --legacy or --fix-indices\n";
            return 1;
        }
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

        // Fix indices from legacy JSON if requested
        if (fix_indices) {
            std::filesystem::path legacy_base_frame = 
                json_output_dir.parent_path() / "json_legacy" / "base_frame_calc" / (pdb_name + ".json");
            if (std::filesystem::exists(legacy_base_frame)) {
                int fixed = x3dna::io::fix_residue_indices_from_json(structure, legacy_base_frame.string());
                if (fixed > 0) {
                    std::cout << "[INFO] Fixed " << fixed << " residue indices from: " 
                              << legacy_base_frame << "\n";
                } else {
                    std::cerr << "[WARNING] Could not fix indices from: " << legacy_base_frame << "\n";
                }
            } else {
                std::cerr << "[WARNING] Legacy JSON not found for --fix-indices: " 
                          << legacy_base_frame << "\n";
            }
        }

        // Try to find legacy JSON file for PDB line caching (optional)
        std::filesystem::path legacy_json_file;
        std::filesystem::path legacy_dir = json_output_dir.parent_path() / "json_legacy";
        std::filesystem::path legacy_file = legacy_dir / "pdb_atoms" / (pdb_name + ".json");
        if (std::filesystem::exists(legacy_file)) {
            legacy_json_file = legacy_file;
        }

        // Create JSON writer (with optional legacy JSON for PDB line caching)
        JsonWriter writer(pdb_file, legacy_json_file);

        // Record PDB atoms
        // Note: Legacy indices are already set on atoms during PDB parsing
        // We do NOT load legacy indices from legacy JSON files - they are generated fresh
        writer.record_pdb_atoms(structure);
        
        // Record residue indices (seidx) - maps residues to atom ranges
        writer.record_residue_indices(structure);

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
        // Track different residue type counts
        size_t water_count = 0;
        size_t ion_count = 0;
        size_t noncanonical_rna_count = 0;
        size_t ligand_count = 0;
        size_t unknown_count = 0;

        for (auto& chain : structure.chains()) {
            for (auto& residue : chain.residues()) {
                // Check residue type
                ResidueType res_type = residue.residue_type();
                char one_letter = residue.one_letter_code();

                // Count by type
                switch (res_type) {
                    case ResidueType::WATER:
                        water_count++;
                        break;
                    case ResidueType::ION:
                        ion_count++;
                        break;
                    case ResidueType::NONCANONICAL_RNA:
                        noncanonical_rna_count++;
                        break;
                    case ResidueType::UNKNOWN:
                        unknown_count++;
                        // DEBUG: Log unknown residues that look like nucleotides
                        if (one_letter == 'A' || one_letter == 'C' || one_letter == 'G' ||
                            one_letter == 'T' || one_letter == 'U') {
                            std::cerr << "WARNING: Residue " << residue.name()
                                      << " has one_letter_code=" << one_letter
                                      << " but residue_type=UNKNOWN\n";
                        }
                        break;
                    default:
                        break;
                }

                // Only process nucleotide residues
                // Include standard nucleotides, modified nucleotides, and noncanonical RNA
                bool is_nucleotide =
                    (res_type != ResidueType::UNKNOWN && res_type != ResidueType::AMINO_ACID &&
                     res_type != ResidueType::WATER && res_type != ResidueType::ION &&
                     res_type != ResidueType::LIGAND);

                // Check for modified nucleotides that have ring atoms but aren't in standard list
                // (This is now handled by residue_type() returning NONCANONICAL_RNA, but keep for compatibility)
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

                    // Get legacy_residue_idx from atoms (set during PDB parsing)
                    // This matches what base_pair_finder uses
                    int legacy_residue_idx = 0;
                    if (!residue.atoms().empty()) {
                        legacy_residue_idx = residue.atoms()[0].legacy_residue_idx();
                    }

                    // Calculate frame and set it on the residue (needed for find_pairs)
                    FrameCalculationResult frame_result = calculator.calculate_frame(residue);

                    if (frame_result.is_valid) {
                        // CRITICAL: Store frame on residue (required for validation in find_pairs)
                        residue.set_reference_frame(frame_result.frame);
                        
                        char base_type = residue.one_letter_code();

                        // Record base_frame_calc
                        // Use legacy_residue_idx if available, otherwise use the counter
                        // IMPORTANT: Keep as 1-based to match legacy JSON output format
                        size_t record_idx = (legacy_residue_idx > 0) ? 
                            static_cast<size_t>(legacy_residue_idx) : 
                            residue_idx;
                        writer.record_base_frame_calc(
                            record_idx, // Keep 1-based to match legacy
                            base_type, frame_result.template_file, frame_result.rms_fit,
                            frame_result.matched_atoms, residue.name(), residue.chain_id(),
                            residue.seq_num(), residue.insertion());

                        // Record ls_fitting
                        writer.record_ls_fitting(
                            record_idx, // Keep 1-based to match legacy
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
                            record_idx, // Keep 1-based to match legacy
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
        std::cout << "  Total residues: " << residue_idx - 1 << "\n";
        std::cout << "  Nucleotides found: " << nucleotides_found << "\n";
        std::cout << "  Frames calculated: " << frames_recorded << "\n";
        // Report residue type breakdown
        size_t total_non_nucleotide = water_count + ion_count + noncanonical_rna_count + ligand_count + unknown_count;
        if (total_non_nucleotide > 0) {
            std::cout << "  Non-nucleotide residues: " << total_non_nucleotide;
            if (water_count > 0) std::cout << " (water: " << water_count;
            if (ion_count > 0) std::cout << (water_count > 0 ? ", ions: " : " (ions: ") << ion_count;
            if (noncanonical_rna_count > 0) std::cout << (water_count > 0 || ion_count > 0 ? ", noncanonical RNA: " : " (noncanonical RNA: ") << noncanonical_rna_count;
            if (ligand_count > 0) std::cout << (water_count > 0 || ion_count > 0 || noncanonical_rna_count > 0 ? ", ligands: " : " (ligands: ") << ligand_count;
            if (unknown_count > 0) std::cout << (water_count > 0 || ion_count > 0 || noncanonical_rna_count > 0 || ligand_count > 0 ? ", other: " : " (other: ") << unknown_count;
            std::cout << ")\n";
        }
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

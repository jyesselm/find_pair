/**
 * @file base_frame_calculator.cpp
 * @brief Implementation of base frame calculator
 */

#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/io/json_writer.hpp>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <cctype>

namespace x3dna {
namespace algorithms {

namespace {
// Standard nucleotide ring geometry (from legacy xyz_ring array)
// Matches RA_LIST order: " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
constexpr std::array<std::array<double, 3>, 9> STANDARD_RING_GEOMETRY = {{
    {{-1.265, 3.177, 0.000}}, // C4
    {{-2.342, 2.364, 0.001}}, // N3
    {{-1.999, 1.087, 0.000}}, // C2
    {{-0.700, 0.641, 0.000}}, // N1
    {{0.424, 1.460, 0.000}},  // C6
    {{0.071, 2.833, 0.000}},  // C5
    {{0.870, 3.969, 0.000}},  // N7 (purine)
    {{0.023, 4.962, 0.000}},  // C8 (purine)
    {{-1.289, 4.551, 0.000}}, // N9 (purine)
}};

// Legacy RA_LIST order for ring atoms
constexpr std::array<const char*, 9> RING_ATOM_NAMES = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ",
                                                        " C5 ", " N7 ", " C8 ", " N9 "};

/**
 * @brief Check nucleotide type by RMSD (matches legacy check_nt_type_by_rmsd)
 * @param residue Residue to check
 * @return RMSD value if calculable, or nullopt if not enough atoms
 */
std::optional<double> check_nt_type_by_rmsd(const core::Residue& residue) {
    // Find ring atoms in residue
    std::vector<geometry::Vector3D> experimental_coords;
    std::vector<geometry::Vector3D> standard_coords;
    int nN = 0; // Count of nitrogen atoms (N1, N3, N7, N9)
    bool has_c1_prime = false;

    for (size_t i = 0; i < RING_ATOM_NAMES.size(); ++i) {
        const char* atom_name = RING_ATOM_NAMES[i];

        // Find this atom in residue
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == atom_name) {
                const auto& pos = atom.position();
                experimental_coords.push_back(geometry::Vector3D(pos.x(), pos.y(), pos.z()));

                // Use corresponding standard geometry
                standard_coords.push_back(geometry::Vector3D(STANDARD_RING_GEOMETRY[i][0],
                                                             STANDARD_RING_GEOMETRY[i][1],
                                                             STANDARD_RING_GEOMETRY[i][2]));

                // Count nitrogen atoms (indices 1=N3, 3=N1, 6=N7, 8=N9)
                if (i == 1 || i == 3 || i == 6 || i == 8) {
                    nN++;
                }
                break;
            }
        }
    }

    // Check for C1' atom (required by legacy)
    for (const auto& atom : residue.atoms()) {
        if (atom.name() == " C1'") {
            has_c1_prime = true;
            break;
        }
    }

    // Legacy requires: (!nN && !C1_prime) -> return DUMMY
    if (nN == 0 && !has_c1_prime) {
        return std::nullopt;
    }

    // Need at least 3 atoms for RMSD calculation
    if (experimental_coords.size() < 3) {
        return std::nullopt;
    }

    // Perform least-squares fitting (matches legacy ls_fitting)
    geometry::LeastSquaresFitter fitter;
    try {
        auto fit_result = fitter.fit(standard_coords, experimental_coords);
        return fit_result.rms;
    } catch (const std::exception&) {
        return std::nullopt;
    }
}
} // namespace

BaseFrameCalculator::BaseFrameCalculator(const std::filesystem::path& template_path)
    : templates_(template_path) {}

FrameCalculationResult BaseFrameCalculator::calculate_frame(core::Residue& residue) {
    FrameCalculationResult result = calculate_frame_impl(residue);

    // Store frame in residue if valid
    if (result.is_valid) {
        residue.set_reference_frame(result.frame);
    }

    return result;
}

FrameCalculationResult
BaseFrameCalculator::calculate_frame_const(const core::Residue& residue) const {
    return calculate_frame_impl(residue);
}

FrameCalculationResult
BaseFrameCalculator::calculate_frame_impl(const core::Residue& residue) const {
    FrameCalculationResult result;
    result.is_valid = false;

    // Get residue type
    core::ResidueType residue_type = residue.residue_type();
    std::string res_name = residue.name();
    while (!res_name.empty() && res_name[0] == ' ')
        res_name.erase(0, 1);
    while (!res_name.empty() && res_name.back() == ' ')
        res_name.pop_back();

    // Check if this is a modified nucleotide not in NT_LIST (requires RMSD check)
    // Legacy: All residues not in NT_LIST go through RMSD check, even if they're recognized as
    // nucleotides Standard NT_LIST: A, C, G, T, U, PSU (pseudouridine), I (inosine) H2U is NOT in
    // NT_LIST, so it needs RMSD check regardless of residue_type
    static const std::vector<std::string> NT_LIST = {
        "A", "C", "G", "T", "U", "PSU", "P5P", "PU", "I", "DI", "ADP", "GDP", "CDP", "UDP", "TDP"};
    bool is_in_nt_list = false;
    std::string res_upper = res_name;
    for (char& c : res_upper) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    for (const auto& nt : NT_LIST) {
        std::string nt_upper = nt;
        for (char& c : nt_upper) {
            c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }
        if (res_upper == nt_upper) {
            is_in_nt_list = true;
            break;
        }
    }
    bool needs_rmsd_check = !is_in_nt_list; // All non-NT_LIST residues need RMSD check

// DEBUG: Log residue info (can be enabled with -DDEBUG_FRAME_CALC)
#ifdef DEBUG_FRAME_CALC
    std::cerr << "DEBUG: Calculating frame for residue: " << residue.name() << " "
              << residue.chain_id() << ":" << residue.seq_num()
              << " (type=" << static_cast<int>(residue_type)
              << ", one_letter=" << residue.one_letter_code() << ", needs_rmsd=" << needs_rmsd_check
              << ")\n";
#endif

    // Check if valid nucleotide
    // Legacy includes modified nucleotides that have ring atoms but aren't in NT_LIST
    // Check if residue has ring atoms as a fallback (similar to legacy residue_ident)
    bool has_ring_atoms = false;
    int ring_atom_count = 0;
    bool has_purine_atoms = false;

    if (residue_type == core::ResidueType::UNKNOWN ||
        residue_type == core::ResidueType::AMINO_ACID ||
        residue_type == core::ResidueType::NONCANONICAL_RNA ||
        (needs_rmsd_check && (residue_type == core::ResidueType::ADENINE ||
                              residue_type == core::ResidueType::CYTOSINE ||
                              residue_type == core::ResidueType::GUANINE ||
                              residue_type == core::ResidueType::THYMINE ||
                              residue_type == core::ResidueType::URACIL))) {
        // Check for ring atoms (C4, N3, C2, N1, C6, C5 are common to all)
        static const std::vector<std::string> common_ring_atoms = {" C4 ", " N3 ", " C2 ",
                                                                   " N1 ", " C6 ", " C5 "};
        static const std::vector<std::string> purine_ring_atoms = {" N7 ", " C8 ", " N9 "};
        static const std::vector<std::string> nitrogen_atoms = {" N1 ", " N3 "};

        for (const auto& atom_name : common_ring_atoms) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == atom_name) {
                    ring_atom_count++;
                    break;
                }
            }
        }

        // Check for purine-specific atoms
        for (const auto& atom_name : purine_ring_atoms) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == atom_name) {
                    has_purine_atoms = true;
#ifdef DEBUG_FRAME_CALC
                    std::cerr << "DEBUG: Found purine atom: " << atom_name << "\n";
#endif
                    break;
                }
            }
        }

        // CRITICAL: Require nitrogen atoms (N1 or N3) to prevent non-nucleotides
        // like glucose (GLC) from being treated as nucleotides
        // This matches the fix in BasePairFinder::is_nucleotide()
        bool has_nitrogen = false;
        for (const auto& atom_name : nitrogen_atoms) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == atom_name) {
                    has_nitrogen = true;
                    break;
                }
            }
            if (has_nitrogen) {
                break;
            }
        }

        // Require >= 3 ring atoms AND at least one nitrogen atom (N1 or N3)
        has_ring_atoms = (ring_atom_count >= 3 && has_nitrogen);
    } else {
        // For standard nucleotides (A, C, G, T, U) that ARE in NT_LIST, skip RMSD check
        // But if they're modified nucleotides not in NT_LIST (like H2U), they need RMSD check
        if (needs_rmsd_check) {
            // Even though it's recognized as a standard nucleotide type, it's not in NT_LIST
            // So we need to check for ring atoms and do RMSD check
            static const std::vector<std::string> common_ring_atoms = {" C4 ", " N3 ", " C2 ",
                                                                       " N1 ", " C6 ", " C5 "};
            static const std::vector<std::string> purine_ring_atoms = {" N7 ", " C8 ", " N9 "};
            static const std::vector<std::string> nitrogen_atoms = {" N1 ", " N3 "};

            for (const auto& atom_name : common_ring_atoms) {
                for (const auto& atom : residue.atoms()) {
                    if (atom.name() == atom_name) {
                        ring_atom_count++;
                        break;
                    }
                }
            }

            for (const auto& atom_name : purine_ring_atoms) {
                for (const auto& atom : residue.atoms()) {
                    if (atom.name() == atom_name) {
                        has_purine_atoms = true;
                        break;
                    }
                }
            }

            bool has_nitrogen = false;
            for (const auto& atom_name : nitrogen_atoms) {
                for (const auto& atom : residue.atoms()) {
                    if (atom.name() == atom_name) {
                        has_nitrogen = true;
                        break;
                    }
                }
                if (has_nitrogen) {
                    break;
                }
            }

            has_ring_atoms = (ring_atom_count >= 3 && has_nitrogen);
        } else {
            // Standard nucleotide in NT_LIST - accept without RMSD check
            has_ring_atoms = true;
        }
    }

    // CRITICAL: ALWAYS check RMSD for ALL residues with ring atoms (even standard nucleotides)
    // This matches legacy residue_ident() which does RMSD check for ALL residues
    // Legacy rejects residues where RMSD > NT_CUTOFF (0.2618), even standard A/C/G/T/U
    // Previous code incorrectly skipped RMSD for standard nucleotides - FIXED!
    if (has_ring_atoms) {
        auto rmsd_result = check_nt_type_by_rmsd(residue);
        // Debug output for 1TTT residue 16 (check with trimmed name)
        std::string res_name_debug = residue.name();
        while (!res_name_debug.empty() && res_name_debug[0] == ' ')
            res_name_debug.erase(0, 1);
        while (!res_name_debug.empty() && res_name_debug.back() == ' ')
            res_name_debug.pop_back();
        if (res_name_debug == "H2U" && residue.chain_id() == 'D' && residue.seq_num() == 16) {
            std::cerr << "[RMSD DEBUG] Residue 16 (H2U, Chain D): ";
            if (rmsd_result.has_value()) {
                std::cerr << "RMSD=" << *rmsd_result << ", threshold=0.2618, ";
                if (*rmsd_result > 0.2618) {
                    std::cerr << "FAILED - will reject\n";
                } else {
                    std::cerr << "PASSED - will accept\n";
                }
            } else {
                std::cerr << "Could not calculate RMSD - will reject\n";
            }
        }
        if (!rmsd_result.has_value() || *rmsd_result > 0.2618) {
            // RMSD check failed - reject this residue (matches legacy behavior)
            // This correctly rejects distorted residues like H2U residue 16 in 1TTT
#ifdef DEBUG_FRAME_CALC
            std::cerr << "DEBUG: RMSD check failed (rmsd="
                      << (rmsd_result.has_value() ? std::to_string(*rmsd_result) : "N/A")
                      << ") - rejecting residue\n";
#endif
            return result; // Cannot calculate frame - RMSD check failed
        }
#ifdef DEBUG_FRAME_CALC
        std::cerr << "DEBUG: RMSD check passed (rmsd=" << *rmsd_result << ")\n";
#endif
    }

    // If no ring atoms, cannot calculate frame
    if (!has_ring_atoms) {
#ifdef DEBUG_FRAME_CALC
        std::cerr << "DEBUG: Skipping - not a nucleotide (type=" << static_cast<int>(residue_type)
                  << ", no ring atoms)\n";
#endif
        return result; // Cannot calculate frame for non-nucleotides
    }

#ifdef DEBUG_FRAME_CALC
    if (has_ring_atoms) {
        std::cerr << "DEBUG: Residue has " << ring_atom_count
                  << " ring atoms - treating as nucleotide\n";
        std::cerr << "DEBUG: Has purine atoms: " << (has_purine_atoms ? "yes" : "no") << "\n";
        if (!has_purine_atoms && ring_atom_count >= 6) {
            // Debug: list all atom names to see why purine detection failed
            std::cerr << "DEBUG: Available atom names in residue: ";
            for (const auto& atom : residue.atoms()) {
                std::cerr << atom.name() << " ";
            }
            std::cerr << "\n";
        }
    }
#endif

    // Determine residue type from atoms (for modified nucleotides or if type is still UNKNOWN)
    // If residue_type is already a standard nucleotide (A, C, G, T, U) and in NT_LIST, use it
    // Otherwise, determine type from atoms
    if (residue_type == core::ResidueType::UNKNOWN ||
        residue_type == core::ResidueType::AMINO_ACID ||
        residue_type == core::ResidueType::NONCANONICAL_RNA || needs_rmsd_check) {
        // Determine type from atoms
        if (has_ring_atoms) {
            if (has_purine_atoms) {
                // Determine purine type (A vs G) by checking for characteristic atoms
                bool has_o6 = false, has_n6 = false, has_n2 = false;
                for (const auto& atom : residue.atoms()) {
                    if (atom.name() == " O6 ")
                        has_o6 = true;
                    if (atom.name() == " N6 ")
                        has_n6 = true;
                    if (atom.name() == " N2 ")
                        has_n2 = true;
                }
                // G has O6, or N2 without N6; A has N6 without O6
                residue_type = (has_o6 || (!has_n6 && has_n2)) ? core::ResidueType::GUANINE
                                                               : core::ResidueType::ADENINE;
            } else {
                // Determine pyrimidine type (C vs T vs U vs P) by checking for characteristic atoms
                bool has_n4 = false, has_c5m = false;
                bool is_pseudouridine = false;

                // Check for characteristic atoms
                for (const auto& atom : residue.atoms()) {
                    if (atom.name() == " N4 ")
                        has_n4 = true;
                    if (atom.name() == " C5M" || atom.name() == " C7 ")
                        has_c5m = true;
                }

                // Check for pseudouridine (PSU): C1' bonded to C5 instead of N1
                // Legacy: p1p2_dist(xyz[c1p], xyz[c5]) <= 2.0 && p1p2_dist(xyz[c1p], xyz[n1]) > 2.0
                auto c1p_opt = residue.find_atom(" C1'");
                auto n1_opt = residue.find_atom(" N1 ");
                auto c5_opt = residue.find_atom(" C5 ");

                if (c1p_opt.has_value() && n1_opt.has_value() && c5_opt.has_value()) {
                    double dist_c1p_n1 = (c1p_opt->position() - n1_opt->position()).length();
                    double dist_c1p_c5 = (c1p_opt->position() - c5_opt->position()).length();
                    // Pseudouridine: C1'-C5 bond (~1.5 Å) instead of C1'-N1 bond
                    if (dist_c1p_c5 <= 2.0 && dist_c1p_n1 > 2.0) {
                        is_pseudouridine = true;
#ifdef DEBUG_FRAME_CALC
                        std::cerr << "DEBUG: Detected pseudouridine: C1'-C5=" << dist_c1p_c5
                                  << ", C1'-N1=" << dist_c1p_n1 << "\n";
#endif
                    }
                }

                // Assign type based on detection
                if (is_pseudouridine) {
                    residue_type = core::ResidueType::PSEUDOURIDINE;
                } else if (has_n4) {
                    residue_type = core::ResidueType::CYTOSINE;
                } else if (has_c5m) {
                    residue_type = core::ResidueType::THYMINE;
                } else {
                    residue_type = core::ResidueType::URACIL; // Default for pyrimidines
                }
            }
        } else {
#ifdef DEBUG_FRAME_CALC
            std::cerr << "DEBUG: Skipping - not a nucleotide (type="
                      << static_cast<int>(residue_type) << ", no ring atoms)\n";
#endif
            return result; // Cannot calculate frame for non-nucleotides
        }
    }

    // Load standard template
    core::Structure standard_template;
    try {
        standard_template = templates_.load_template(residue_type);
        result.template_file = templates_.get_template_path(residue_type);

#ifdef DEBUG_FRAME_CALC
        std::cerr << "DEBUG: Template loaded: " << result.template_file << "\n";
        std::cerr << "DEBUG: Template has " << standard_template.num_atoms() << " atoms\n";
#endif
    } catch (const std::exception& e) {
// Template not found or failed to load
#ifdef DEBUG_FRAME_CALC
        std::cerr << "DEBUG: Template loading failed: " << e.what() << "\n";
        std::cerr << "DEBUG: Template path: " << templates_.template_path() << "\n";
        std::cerr << "DEBUG: Residue type: " << static_cast<int>(residue_type) << "\n";
#endif
        return result;
    }

    // Match ring atoms
    // Pass detected residue_type to ensure correct atom list is used
    // (residue.residue_type() might still be UNKNOWN for modified nucleotides)
    MatchedAtoms matched = RingAtomMatcher::match(residue, standard_template,
                                                  std::optional<core::ResidueType>(residue_type));

#ifdef DEBUG_FRAME_CALC
    std::cerr << "DEBUG: Matched " << matched.num_matched << " atoms\n";
    std::cerr << "DEBUG: Matched atom names: ";
    for (const auto& name : matched.atom_names) {
        std::cerr << name << " ";
    }
    std::cerr << "\n";
#endif

    if (!matched.is_valid()) {
#ifdef DEBUG_FRAME_CALC
        std::cerr << "DEBUG: Not enough matched atoms (need >= 3, got " << matched.num_matched
                  << ")\n";
        std::cerr << "DEBUG: Residue has " << residue.num_atoms() << " total atoms\n";
#endif
        // Not enough matched atoms (need at least 3)
        return result;
    }

    result.num_matched = matched.num_matched;
    result.matched_atoms = matched.atom_names;

    // Extract coordinates from matched atoms
    std::vector<geometry::Vector3D> standard_coords;
    std::vector<geometry::Vector3D> experimental_coords;

    for (size_t i = 0; i < matched.num_matched; ++i) {
        standard_coords.push_back(matched.standard[i].position());
        experimental_coords.push_back(matched.experimental[i].position());
    }

    // Perform least-squares fitting
    geometry::LeastSquaresFitter fitter;
    geometry::LeastSquaresFitter::FitResult fit_result =
        fitter.fit(standard_coords, experimental_coords);

    // Extract rotation and translation
    result.rotation_matrix = fit_result.rotation;
    result.translation = fit_result.translation;
    result.rms_fit = fit_result.rms;

    // Create reference frame from rotation and translation
    result.frame = core::ReferenceFrame(result.rotation_matrix, result.translation);
    result.is_valid = true;

    return result;
}

void BaseFrameCalculator::calculate_all_frames(core::Structure& structure) {
    // CRITICAL: Iterate in legacy index order (PDB file order) to match legacy behavior
    // This ensures frames are calculated in the same order as legacy code
    // Get residues in legacy order (PDB file order)
    std::vector<core::Residue*> residues_in_order;
    for (const auto* residue_ptr : structure.residues_in_legacy_order()) {
        // Need non-const pointer for modification
        residues_in_order.push_back(const_cast<core::Residue*>(residue_ptr));
    }

    // Iterate through residues in legacy order, calculating frames
    // Note: We don't skip UNKNOWN residues here because calculate_frame_impl
    // has logic to detect modified nucleotides by ring atoms (like XGR, XCR, etc.)
    // Only skip amino acids explicitly
    for (auto* residue : residues_in_order) {
        // Skip only amino acids - let calculate_frame handle UNKNOWN residues
        // (it will check for ring atoms to detect modified nucleotides)
        if (residue->residue_type() == core::ResidueType::AMINO_ACID) {
            continue;
        }

        // Calculate frame (modifies residue to store frame)
        // calculate_frame_impl will check for ring atoms if residue_type is UNKNOWN
        calculate_frame(*residue);
    }
}

void BaseFrameCalculator::set_template_path(const std::filesystem::path& template_path) {
    templates_.set_template_path(template_path);
}

bool BaseFrameCalculator::detect_rna(const core::Structure& structure) {
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            if (residue.find_atom(" O2'").has_value() || residue.find_atom(" O2*").has_value()) {
                return true;
            }
        }
    }
    return false;
}

size_t BaseFrameCalculator::calculate_and_record_frames(core::Structure& structure,
                                                        io::JsonWriter& writer) {
    // Auto-detect RNA vs DNA
    bool is_rna = detect_rna(structure);
    set_is_rna(is_rna);

    // Calculate and record frames in one pass
    // Get residues in legacy order (PDB file order)
    auto residues = structure.residues_in_legacy_order();
    size_t frames_recorded = 0;

    for (const auto* residue_ptr : residues) {
        auto* residue = const_cast<core::Residue*>(residue_ptr);

        // Skip amino acids (calculate_frame handles other types including UNKNOWN)
        if (residue->residue_type() == core::ResidueType::AMINO_ACID) {
            continue;
        }

        // Calculate frame (stores frame on residue and returns full result)
        FrameCalculationResult frame_result = calculate_frame(*residue);

        if (!frame_result.is_valid) {
            continue;
        }

        // Get legacy_residue_idx from atoms
        int legacy_residue_idx = 0;
        if (!residue->atoms().empty()) {
            legacy_residue_idx = residue->atoms()[0].legacy_residue_idx();
        }

        if (legacy_residue_idx <= 0) {
            continue;
        }

        char base_type = residue->one_letter_code();
        size_t record_idx = static_cast<size_t>(legacy_residue_idx);

        // Record base_frame_calc
        writer.record_base_frame_calc(record_idx, base_type, frame_result.template_file,
                                      frame_result.rms_fit, frame_result.matched_atoms,
                                      residue->name(), residue->chain_id(), residue->seq_num(),
                                      residue->insertion());

        // Record ls_fitting
        writer.record_ls_fitting(record_idx, frame_result.num_matched, frame_result.rms_fit,
                                 frame_result.rotation_matrix, frame_result.translation,
                                 residue->name(), residue->chain_id(), residue->seq_num(),
                                 residue->insertion());

        // Record frame_calc
        std::vector<geometry::Vector3D> standard_coords, experimental_coords;
        writer.record_frame_calc(record_idx, base_type, frame_result.template_file,
                                 frame_result.rms_fit, standard_coords, experimental_coords,
                                 residue->name(), residue->chain_id(), residue->seq_num(),
                                 residue->insertion());

        frames_recorded++;
    }

    return frames_recorded;
}

} // namespace algorithms
} // namespace x3dna

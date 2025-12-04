/**
 * @file base_frame_calculator.cpp
 * @brief Implementation of base frame calculator
 */

#include <x3dna/algorithms/base_frame_calculator.hpp>
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
 * @brief Result of RMSD check
 */
struct RmsdCheckResult {
    std::optional<double> rmsd;
    bool found_purine_atoms; // Whether any purine atoms (N7, C8, N9) were found
};

/**
 * @brief Check nucleotide type by RMSD (matches legacy check_nt_type_by_rmsd)
 * @param residue Residue to check
 * @return RMSD value if calculable, or nullopt if not enough atoms, plus purine atom flag
 */
RmsdCheckResult check_nt_type_by_rmsd(const core::Residue& residue) {
    // Find ring atoms in residue
    // LEGACY BEHAVIOR: Try ALL 9 ring atoms first (matches legacy residue_ident)
    // Then if RMSD fails and purine atoms were found, retry with pyrimidine-only
    std::vector<geometry::Vector3D> experimental_coords;
    std::vector<geometry::Vector3D> standard_coords;
    int nN = 0; // Count of nitrogen atoms (N1, N3, N7, N9)
    bool has_c1_prime = false;
    int purine_atom_count = 0; // Count of purine atoms found (N7, C8, N9)
    
#ifdef DEBUG_FRAME_CALC
    std::string res_name_check = residue.name();
    while (!res_name_check.empty() && res_name_check[0] == ' ')
        res_name_check.erase(0, 1);
    while (!res_name_check.empty() && res_name_check.back() == ' ')
        res_name_check.pop_back();
    bool is_cvc = (res_name_check == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7);
    if (is_cvc) {
        std::cerr << "DEBUG: check_nt_type_by_rmsd for CVC B7, residue has " << residue.num_atoms() << " atoms\n";
    }
#endif
    
    // Try to match ALL ring atoms (like legacy does)
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
                
                // Count purine atoms (indices 6=N7, 7=C8, 8=N9) - for two-try fallback
                if (i >= 6) {
                    purine_atom_count++;
                }
#ifdef DEBUG_FRAME_CALC
                if (is_cvc) {
                    std::cerr << "DEBUG: Matched ring atom " << i << ": " << atom_name << "\n";
                }
#endif
                break;
            }
#ifdef DEBUG_FRAME_CALC
            if (is_cvc && i < 3) {  // Only show first few to avoid spam
                std::cerr << "DEBUG: Checking atom '" << atom.name() << "' vs '" << atom_name << "'\n";
            }
#endif
        }
    }
    
#ifdef DEBUG_FRAME_CALC
    if (is_cvc) {
        std::cerr << "DEBUG: check_nt_type_by_rmsd result: matched=" << experimental_coords.size() 
                  << ", nN=" << nN << ", purine_count=" << purine_atom_count << "\n";
    }
#endif

    // Check for C1' or C1R atom (required by legacy)
    // Some nucleotides like NMN use C1R instead of C1'
    for (const auto& atom : residue.atoms()) {
        if (atom.name() == " C1'" || atom.name() == " C1R") {
            has_c1_prime = true;
            break;
        }
    }

    // Legacy requires: (!nN && !C1_prime) -> return DUMMY
    if (nN == 0 && !has_c1_prime) {
        return {std::nullopt, purine_atom_count > 0};
    }

    // Need at least 3 atoms for RMSD calculation
    if (experimental_coords.size() < 3) {
        return {std::nullopt, purine_atom_count > 0};
    }

    // Perform least-squares fitting (matches legacy ls_fitting)
    geometry::LeastSquaresFitter fitter;
    try {
        auto fit_result = fitter.fit(standard_coords, experimental_coords);
        return {fit_result.rms, purine_atom_count > 0};
    } catch (const std::exception&) {
        return {std::nullopt, purine_atom_count > 0};
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
    
    // DEBUG: Print for any residue with B:7 to see what we get
    if (residue.chain_id() == 'B' && residue.seq_num() == 7) {
        std::cerr << "DEBUG: B7 residue - name='" << res_name << "' (len=" << res_name.length() 
                  << "), raw_name='" << residue.name() << "' (len=" << residue.name().length() 
                  << "), type=" << static_cast<int>(residue_type) 
                  << ", num_atoms=" << residue.num_atoms() << "\n";
        if (residue.num_atoms() > 0) {
            std::cerr << "DEBUG: B7 atom names:\n";
            for (const auto& atom : residue.atoms()) {
                std::cerr << "  '" << atom.name() << "'\n";
            }
        }
    }

#ifdef DEBUG_FRAME_CALC
    // Debug: show all atom names for CVC or any UNKNOWN residue
    std::string residue_name_raw = residue.name();
    while (!residue_name_raw.empty() && residue_name_raw[0] == ' ')
        residue_name_raw.erase(0, 1);
    while (!residue_name_raw.empty() && residue_name_raw.back() == ' ')
        residue_name_raw.pop_back();
    
    if (residue_name_raw == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7) {
        std::cerr << "DEBUG: CVC B7 has " << residue.num_atoms() << " atoms:\n";
        for (const auto& atom : residue.atoms()) {
            std::cerr << "  Atom: '" << atom.name() << "' (len=" << atom.name().length() << ")\n";
        }
    }
#endif

    if (residue_type == core::ResidueType::UNKNOWN ||
        residue_type == core::ResidueType::AMINO_ACID ||
        residue_type == core::ResidueType::NONCANONICAL_RNA ||
        (needs_rmsd_check && (residue_type == core::ResidueType::ADENINE ||
                              residue_type == core::ResidueType::CYTOSINE ||
                              residue_type == core::ResidueType::GUANINE ||
                              residue_type == core::ResidueType::THYMINE ||
                              residue_type == core::ResidueType::URACIL))) {
        // ALWAYS PRINT FOR DEBUGGING
        std::cerr << "DEBUG: Entered ring atom counting block for " << res_name << " " 
                  << residue.chain_id() << ":" << residue.seq_num() 
                  << " (type=" << static_cast<int>(residue_type) << ")\n";
        
        // Check for ring atoms (C4, N3, C2, N1, C6, C5 are common to all)
        static const std::vector<std::string> common_ring_atoms = {" C4 ", " N3 ", " C2 ",
                                                                   " N1 ", " C6 ", " C5 "};
        static const std::vector<std::string> purine_ring_atoms = {" N7 ", " C8 ", " N9 "};
        static const std::vector<std::string> nitrogen_atoms = {" N1 ", " N3 "};

#ifdef DEBUG_FRAME_CALC
        // Always show for UNKNOWN residues to debug CVC
        if (residue_type == core::ResidueType::UNKNOWN) {
            std::string res_name_debug = res_name;
            std::cerr << "DEBUG: UNKNOWN residue " << res_name_debug << " " << residue.chain_id() 
                      << ":" << residue.seq_num() << " has " << residue.num_atoms() << " atoms\n";
            if (res_name_debug == "CVC") {
                std::cerr << "DEBUG: CVC atom names:\n";
                for (const auto& atom : residue.atoms()) {
                    std::cerr << "  '" << atom.name() << "'\n";
                }
            }
        }
#endif

        // UNCONDITIONAL DEBUG: Print all atom names for UNKNOWN residues to diagnose
        if (residue_type == core::ResidueType::UNKNOWN) {
            std::string res_name_debug = res_name;
            std::cerr << "DEBUG: UNKNOWN residue '" << res_name_debug << "' " << residue.chain_id() 
                      << ":" << residue.seq_num() << " has " << residue.num_atoms() << " atoms\n";
            if (residue.num_atoms() > 0) {
                std::cerr << "DEBUG: Atom names: ";
                for (const auto& atom : residue.atoms()) {
                    std::cerr << "'" << atom.name() << "' ";
                }
                std::cerr << "\n";
            }
        }
        
        for (const auto& atom_name : common_ring_atoms) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == atom_name) {
                    ring_atom_count++;
#ifdef DEBUG_FRAME_CALC
                    std::cerr << "DEBUG: Found common ring atom: " << atom_name << "\n";
#endif
                    break;
                }
            }
        }
        
        if (residue_type == core::ResidueType::UNKNOWN) {
            std::cerr << "DEBUG: After common ring atoms, count=" << ring_atom_count << "\n";
        }

        // Also count purine atoms for ring_atom_count (legacy tries ALL 9 atoms)
        // This allows residues like CVC (C4, C8, N9) to be detected
        for (const auto& atom_name : purine_ring_atoms) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == atom_name) {
                    ring_atom_count++;
#ifdef DEBUG_FRAME_CALC
                    std::cerr << "DEBUG: Found purine ring atom: " << atom_name << "\n";
#endif
                    break;
                }
#ifdef DEBUG_FRAME_CALC
                // Debug: show what atom names we're comparing
                if (res_name == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7) {
                    std::cerr << "DEBUG: Comparing '" << atom.name() << "' with '" << atom_name << "'\n";
                }
#endif
            }
        }

        // Check for purine-specific atoms
            // CRITICAL FIX: Only count as purine if we have BOTH N7 and C8
            // Some modified pyrimidines (like 70U) have C8 in side chains, not as ring atom
            // Legacy checks for this by requiring both N7 and C8 for purine detection
            bool has_n7 = false, has_c8 = false;
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == " N7 ") has_n7 = true;
                if (atom.name() == " C8 ") has_c8 = true;
            }
            
            // Require both N7 and C8 for purine detection
            // This prevents false purine detection for modified pyrimidines with C8 in modifications
            has_purine_atoms = (has_n7 && has_c8);
#ifdef DEBUG_FRAME_CALC
            if (has_purine_atoms) {
                std::cerr << "DEBUG: Found purine atoms: N7=" << has_n7 
                          << ", C8=" << has_c8 << "\n";
            } else if (has_c8 && !has_n7) {
                std::cerr << "DEBUG: Has C8 but not N7 - NOT a purine (likely modified pyrimidine)\n";
            }
#endif

        // Legacy only requires >= 3 ring atoms (by name) and C1' presence
        // The RMSD check (check_nt_type_by_rmsd) rejects non-nucleotides via:
        // 1. No C1' (e.g., glucose) → rejected
        // 2. No nitrogen atoms AND no C1' → rejected (nN=0 && !C1_prime → DUMMY)
        // So we don't need explicit nitrogen check here
        has_ring_atoms = (ring_atom_count >= 3);
        
        if (residue_type == core::ResidueType::UNKNOWN) {
            std::cerr << "DEBUG: Final ring_atom_count=" << ring_atom_count << ", has_ring_atoms=" << has_ring_atoms << "\n";
        }
#ifdef DEBUG_FRAME_CALC
        std::string residue_name_debug = residue.name();
        while (!residue_name_debug.empty() && residue_name_debug[0] == ' ')
            residue_name_debug.erase(0, 1);
        while (!residue_name_debug.empty() && residue_name_debug.back() == ' ')
            residue_name_debug.pop_back();
        std::cerr << "DEBUG: Residue " << residue_name_debug << " " << residue.chain_id() << ":" << residue.seq_num()
                  << " - Ring atom count: " << ring_atom_count << ", has_ring_atoms: " << has_ring_atoms << "\n";
        if (residue_name_debug == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7) {
            std::cerr << "DEBUG: CVC B7 final check - ring_atom_count=" << ring_atom_count << ", has_ring_atoms=" << has_ring_atoms << "\n";
        }
#endif
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

            // Also count purine atoms for ring_atom_count (legacy tries ALL 9 atoms)
            for (const auto& atom_name : purine_ring_atoms) {
                for (const auto& atom : residue.atoms()) {
                    if (atom.name() == atom_name) {
                        ring_atom_count++;
                        has_purine_atoms = true;
                        break;
                    }
                }
            }

            // Legacy only requires >= 3 ring atoms
            has_ring_atoms = (ring_atom_count >= 3);
        } else {
            // Standard nucleotide in NT_LIST - accept without RMSD check
            has_ring_atoms = true;
        }
    }

    // CRITICAL: ALWAYS check RMSD for ALL residues with ring atoms (even standard nucleotides)
    // This matches legacy residue_ident() which does RMSD check for ALL residues
    // Legacy rejects residues where RMSD > NT_CUTOFF (0.2618), even standard A/C/G/T/U
    
    if (has_ring_atoms) {
        // DEBUG: Always print for CVC
        if (res_name == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7) {
            std::cerr << "DEBUG: CVC B7 - Calling check_nt_type_by_rmsd\n";
        }
        auto rmsd_check = check_nt_type_by_rmsd(residue);
        auto rmsd_result = rmsd_check.rmsd;
        bool found_purine_atoms = rmsd_check.found_purine_atoms;
        
        // DEBUG: Always print for CVC
        if (res_name == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7) {
            std::cerr << "DEBUG: CVC B7 - RMSD result: " 
                      << (rmsd_result.has_value() ? std::to_string(*rmsd_result) : "nullopt")
                      << ", found_purine_atoms=" << found_purine_atoms << "\n";
        }
        
        // Use strict threshold (0.2618) for all bases
        double rmsd_threshold = 0.2618;
        
        // Debug output disabled
        if (!rmsd_result.has_value() || *rmsd_result > rmsd_threshold) {
            // RMSD check failed with all atoms
            // Legacy's TWO-TRY approach: If purine atoms were found, retry with ONLY pyrimidine atoms
            // This handles cases where purine ring is distorted but pyrimidine core is fine
            // Note: Legacy uses kr > 0 (any purine atoms found), not requiring both N7 and C8
            if (found_purine_atoms) {
                // Retry with only 6 pyrimidine atoms (C4, N3, C2, N1, C6, C5)
                // Calculate RMSD using ONLY pyrimidine ring atoms (first 6)
                std::vector<geometry::Vector3D> experimental_coords;
                std::vector<geometry::Vector3D> standard_coords;
                
                for (size_t i = 0; i < 6; ++i) {  // Only pyrimidine atoms (indices 0-5)
                    const char* atom_name = RING_ATOM_NAMES[i];
                    for (const auto& atom : residue.atoms()) {
                        if (atom.name() == atom_name) {
                            const auto& pos = atom.position();
                            experimental_coords.push_back(geometry::Vector3D(pos.x(), pos.y(), pos.z()));
                            standard_coords.push_back(geometry::Vector3D(STANDARD_RING_GEOMETRY[i][0],
                                                                         STANDARD_RING_GEOMETRY[i][1],
                                                                         STANDARD_RING_GEOMETRY[i][2]));
                            break;
                        }
                    }
                }
                
                std::optional<double> pyrimidine_rmsd;
                if (experimental_coords.size() >= 3) {
                    geometry::LeastSquaresFitter fitter;
                    try {
                        auto fit_result = fitter.fit(standard_coords, experimental_coords);
                        pyrimidine_rmsd = fit_result.rms;
                    } catch (const std::exception&) {
                        pyrimidine_rmsd = std::nullopt;
                    }
                }
                
#ifdef DEBUG_FRAME_CALC
                std::cerr << "DEBUG: First RMSD check failed (rmsd="
                          << (rmsd_result.has_value() ? std::to_string(*rmsd_result) : "N/A")
                          << "), retrying as pyrimidine-only\n";
#endif
                
                // If still fails, reject
                if (!pyrimidine_rmsd.has_value() || *pyrimidine_rmsd > rmsd_threshold) {
#ifdef DEBUG_FRAME_CALC
                    std::cerr << "DEBUG: Second try also failed - rejecting residue\n";
#endif
                    return result; // Cannot calculate frame - both attempts failed
                }
                
                // Second try passed - accept as pyrimidine (force type to pyrimidine)
                has_purine_atoms = false; // Treat as pyrimidine
                rmsd_result = pyrimidine_rmsd;
            } else {
                // No purine atoms and first try failed - reject
#ifdef DEBUG_FRAME_CALC
                std::cerr << "DEBUG: RMSD check failed (rmsd="
                          << (rmsd_result.has_value() ? std::to_string(*rmsd_result) : "N/A")
                          << ", threshold=" << rmsd_threshold << ") - rejecting residue\n";
#endif
                return result; // Cannot calculate frame - RMSD check failed
            }
        }
#ifdef DEBUG_FRAME_CALC
        std::cerr << "DEBUG: RMSD check passed (rmsd=" << *rmsd_result
                  << ", threshold=" << rmsd_threshold << ")\n";
#endif
        
        // DEBUG: Always print for CVC
        if (res_name == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7) {
            std::cerr << "DEBUG: CVC B7 - RMSD check passed, continuing to frame calculation\n";
        }
    }

    // If no ring atoms, cannot calculate frame
    if (!has_ring_atoms) {
#ifdef DEBUG_FRAME_CALC
        std::string residue_name_final = residue.name();
        while (!residue_name_final.empty() && residue_name_final[0] == ' ')
            residue_name_final.erase(0, 1);
        while (!residue_name_final.empty() && residue_name_final.back() == ' ')
            residue_name_final.pop_back();
        std::cerr << "DEBUG: Skipping - not a nucleotide (type=" << static_cast<int>(residue_type)
                  << ", no ring atoms, ring_count=" << ring_atom_count << ", residue=" << residue_name_final
                  << " " << residue.chain_id() << ":" << residue.seq_num() << ")\n";
        if (residue_name_final == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7) {
            std::cerr << "DEBUG: CVC B7 REJECTED - ring_atom_count=" << ring_atom_count << ", has_ring_atoms=" << has_ring_atoms << "\n";
            std::cerr << "DEBUG: CVC B7 has " << residue.num_atoms() << " atoms total\n";
            for (const auto& atom : residue.atoms()) {
                std::cerr << "  Atom: '" << atom.name() << "'\n";
            }
        }
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

    // DEBUG: Always print for CVC
    if (res_name == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7) {
        std::cerr << "DEBUG: CVC B7 - After RMSD check, determining residue type. Current type=" 
                  << static_cast<int>(residue_type) << ", has_purine_atoms=" << has_purine_atoms << "\n";
    }
    
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

    // DEBUG: Always print for CVC
    if (res_name == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7) {
        std::cerr << "DEBUG: CVC B7 - Loading template for type=" << static_cast<int>(residue_type) << "\n";
    }
    
    // Load standard template
    core::Structure standard_template;
    try {
        standard_template = templates_.load_template(residue_type);
        result.template_file = templates_.get_template_path(residue_type);
        
        // DEBUG: Always print for CVC
        if (res_name == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7) {
            std::cerr << "DEBUG: CVC B7 - Template loaded: " << result.template_file << "\n";
        }

#ifdef DEBUG_FRAME_CALC
        std::cerr << "DEBUG: Template loaded: " << result.template_file << "\n";
        std::cerr << "DEBUG: Template has " << standard_template.num_atoms() << " atoms\n";
#endif
    } catch (const std::exception& e) {
// Template not found or failed to load
        // DEBUG: Always print for CVC
        if (res_name == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7) {
            std::cerr << "DEBUG: CVC B7 - Template loading failed: " << e.what() << "\n";
            std::cerr << "DEBUG: CVC B7 - Template path: " << templates_.template_path() << "\n";
            std::cerr << "DEBUG: CVC B7 - Residue type: " << static_cast<int>(residue_type) << "\n";
        }
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

    // DEBUG: Always print for CVC
    if (res_name == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7) {
        std::cerr << "DEBUG: CVC B7 - Matched " << matched.num_matched << " atoms\n";
        std::cerr << "DEBUG: CVC B7 - Matched atom names: ";
        for (const auto& name : matched.atom_names) {
            std::cerr << name << " ";
        }
        std::cerr << "\n";
        std::cerr << "DEBUG: CVC B7 - matched.is_valid()=" << matched.is_valid() << "\n";
    }

#ifdef DEBUG_FRAME_CALC
    std::cerr << "DEBUG: Matched " << matched.num_matched << " atoms\n";
    std::cerr << "DEBUG: Matched atom names: ";
    for (const auto& name : matched.atom_names) {
        std::cerr << name << " ";
    }
    std::cerr << "\n";
#endif

    if (!matched.is_valid()) {
        // DEBUG: Always print for CVC
        if (res_name == "CVC" && residue.chain_id() == 'B' && residue.seq_num() == 7) {
            std::cerr << "DEBUG: CVC B7 - REJECTED: Not enough matched atoms (need >= 3, got " 
                      << matched.num_matched << ")\n";
            std::cerr << "DEBUG: CVC B7 - Residue has " << residue.num_atoms() << " total atoms\n";
        }
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


} // namespace algorithms
} // namespace x3dna

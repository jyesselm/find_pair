/**
 * @file debug_hbond_conflict_resolution.cpp
 * @brief Debug tool to trace conflict resolution step-by-step for modern code
 * 
 * This tool traces through all phases of conflict resolution:
 * - Phase 1: Initial conflict detection
 * - Phase 2: idx2 population and linkage type calculation
 * - Phase 3: Additional conflict marking
 * 
 * Outputs detailed information at each step for comparison with legacy.
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <filesystem>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/hydrogen_bond_finder.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;

struct Phase1Debug {
    size_t num_iter;
    std::string donor;
    std::string acceptor;
    double distance;
    double dtmp1;
    double dtmp2;
    size_t ddidx1;
    size_t ddidx2;
    bool is_conflict;
};

struct Phase2Debug {
    size_t k;
    std::string donor;
    std::string acceptor;
    double distance;
    bool is_conflict_before;
    int idx2_0;
    int idx2_1;
    std::vector<std::pair<size_t, std::string>> shared_with; // (hbond_idx, which_atom)
};

struct Phase3Debug {
    size_t k;
    std::string donor;
    std::string acceptor;
    double distance_before;
    double distance_after;
    int linkage_type;
    bool is_conflict_before;
    bool is_conflict_after;
    bool should_mark_conflict; // linkage_type != 18 && in_range
};

// Modified version of resolve_conflicts with debug output
std::vector<Phase1Debug> trace_phase1(std::vector<HydrogenBondResult>& hbonds, double /*hb_lower*/, double /*hb_dist2*/) {
    std::vector<Phase1Debug> phase1_trace;
    
    if (hbonds.empty()) {
        return phase1_trace;
    }
    
    const size_t num_hbonds = hbonds.size();
    std::vector<bool> matched_idx(num_hbonds, false);
    
    size_t num_iter = 0;
    size_t m = 0;
    
    while (true) {
        Phase1Debug debug;
        
        // Find next unmatched H-bond
        while (num_iter < num_hbonds && matched_idx[num_iter]) {
            num_iter++;
        }
        
        if (num_iter >= num_hbonds) {
            break;
        }
        
        debug.num_iter = num_iter;
        debug.donor = hbonds[num_iter].donor_atom;
        debug.acceptor = hbonds[num_iter].acceptor_atom;
        debug.distance = hbonds[num_iter].distance;
        
        // Find shortest H-bonds for donor and acceptor atoms
        double dtmp[3] = {0.0, 0.0, 0.0};
        size_t ddidx[3] = {num_hbonds, num_hbonds, num_hbonds};
        
        double current_dist = std::abs(hbonds[num_iter].distance);
        dtmp[1] = current_dist;
        ddidx[1] = num_iter;
        dtmp[2] = current_dist;
        ddidx[2] = num_iter;
        
        // Find shorter H-bonds for same donor or acceptor
        for (size_t n = 0; n < num_hbonds; ++n) {
            if (n == num_iter || matched_idx[n]) {
                continue;
            }
            
            double dist_n = std::abs(hbonds[n].distance);
            
            if (hbonds[n].donor_atom == hbonds[num_iter].donor_atom && dist_n < dtmp[1]) {
                dtmp[1] = dist_n;
                ddidx[1] = n;
            }
            
            if (hbonds[n].acceptor_atom == hbonds[num_iter].acceptor_atom && dist_n < dtmp[2]) {
                dtmp[2] = dist_n;
                ddidx[2] = n;
            }
        }
        
        debug.dtmp1 = dtmp[1];
        debug.dtmp2 = dtmp[2];
        debug.ddidx1 = ddidx[1];
        debug.ddidx2 = ddidx[2];
        
        // Check if conflict
        debug.is_conflict = (ddidx[1] == ddidx[2] && ddidx[1] < num_hbonds);
        
        if (debug.is_conflict) {
            size_t k = ddidx[1];
            hbonds[k].distance = -hbonds[k].distance;
            
            num_iter = 0;
            for (size_t n = 0; n < num_hbonds; ++n) {
                if (matched_idx[n]) {
                    continue;
                }
                if (hbonds[n].donor_atom == hbonds[k].donor_atom ||
                    hbonds[n].acceptor_atom == hbonds[k].acceptor_atom) {
                    matched_idx[n] = true;
                    m++;
                }
            }
            
            if (m >= num_hbonds) {
                break;
            }
        } else {
            num_iter++;
        }
        
        phase1_trace.push_back(debug);
    }
    
    return phase1_trace;
}

std::vector<Phase2Debug> trace_phase2(const std::vector<HydrogenBondResult>& hbonds, 
                                       std::vector<std::vector<int>>& idx2) {
    std::vector<Phase2Debug> phase2_trace;
    const size_t num_hbonds = hbonds.size();
    
    for (size_t k = 0; k < num_hbonds; ++k) {
        Phase2Debug debug;
        debug.k = k;
        debug.donor = hbonds[k].donor_atom;
        debug.acceptor = hbonds[k].acceptor_atom;
        debug.distance = hbonds[k].distance;
        debug.is_conflict_before = (hbonds[k].distance < 0.0);
        
        if (hbonds[k].distance > 0.0) {
            continue; // Not a conflict
        }
        
        // Mark as conflict (9, 9)
        idx2[k][0] = 9;
        idx2[k][1] = 9;
        debug.idx2_0 = 9;
        debug.idx2_1 = 9;
        
        // Find all non-conflicted H-bonds that share atoms
        for (size_t m_idx = 0; m_idx < num_hbonds; ++m_idx) {
            if (m_idx == k || hbonds[m_idx].distance < 0.0) {
                continue;
            }
            
            if (hbonds[m_idx].donor_atom == hbonds[k].donor_atom) {
                idx2[m_idx][0] = 1;
                debug.shared_with.push_back({m_idx, "donor"});
            }
            if (hbonds[m_idx].acceptor_atom == hbonds[k].acceptor_atom) {
                idx2[m_idx][1] = 1;
                debug.shared_with.push_back({m_idx, "acceptor"});
            }
        }
        
        phase2_trace.push_back(debug);
    }
    
    // Also trace non-conflicted H-bonds to show their idx2 values
    for (size_t k = 0; k < num_hbonds; ++k) {
        if (hbonds[k].distance < 0.0) {
            continue; // Already traced
        }
        
        Phase2Debug debug;
        debug.k = k;
        debug.donor = hbonds[k].donor_atom;
        debug.acceptor = hbonds[k].acceptor_atom;
        debug.distance = hbonds[k].distance;
        debug.is_conflict_before = false;
        debug.idx2_0 = idx2[k][0];
        debug.idx2_1 = idx2[k][1];
        phase2_trace.push_back(debug);
    }
    
    return phase2_trace;
}

std::vector<Phase3Debug> trace_phase3(std::vector<HydrogenBondResult>& hbonds,
                                       const std::vector<std::vector<int>>& idx2,
                                       double hb_lower, double hb_dist2) {
    std::vector<Phase3Debug> phase3_trace;
    const size_t num_hbonds = hbonds.size();
    
    for (size_t k = 0; k < num_hbonds; ++k) {
        Phase3Debug debug;
        debug.k = k;
        debug.donor = hbonds[k].donor_atom;
        debug.acceptor = hbonds[k].acceptor_atom;
        debug.distance_before = hbonds[k].distance;
        debug.is_conflict_before = (hbonds[k].distance < 0.0);
        
        int linkage_sum = idx2[k][0] + idx2[k][1];
        hbonds[k].linkage_type = linkage_sum;
        debug.linkage_type = linkage_sum;
        
        debug.should_mark_conflict = (linkage_sum != 18 && hbonds[k].distance > 0.0 &&
                                      hbonds[k].distance >= hb_lower && hbonds[k].distance <= hb_dist2);
        
        if (debug.should_mark_conflict) {
            hbonds[k].distance = -hbonds[k].distance;
        }
        
        debug.distance_after = hbonds[k].distance;
        debug.is_conflict_after = (hbonds[k].distance < 0.0);
        
        phase3_trace.push_back(debug);
    }
    
    return phase3_trace;
}

void print_phase1_trace(const std::vector<Phase1Debug>& trace) {
    std::cout << "\n========================================\n";
    std::cout << "PHASE 1: Initial Conflict Detection\n";
    std::cout << "========================================\n";
    
    for (const auto& debug : trace) {
        std::cout << "\nIteration " << debug.num_iter << ":\n";
        std::cout << "  H-bond: " << debug.donor << " -> " << debug.acceptor 
                  << " (dist=" << std::fixed << std::setprecision(6) << debug.distance << ")\n";
        std::cout << "  Shortest for donor: idx=" << debug.ddidx1 << " dist=" << debug.dtmp1 << "\n";
        std::cout << "  Shortest for acceptor: idx=" << debug.ddidx2 << " dist=" << debug.dtmp2 << "\n";
        std::cout << "  Conflict detected: " << (debug.is_conflict ? "YES" : "NO") << "\n";
        if (debug.is_conflict) {
            std::cout << "    -> Marking H-bond " << debug.ddidx1 << " as conflict\n";
        }
    }
}

void print_phase2_trace(const std::vector<Phase2Debug>& trace) {
    std::cout << "\n========================================\n";
    std::cout << "PHASE 2: idx2 Population\n";
    std::cout << "========================================\n";
    
    for (const auto& debug : trace) {
        std::cout << "\nH-bond " << debug.k << ": " << debug.donor << " -> " << debug.acceptor << "\n";
        std::cout << "  Distance: " << std::fixed << std::setprecision(6) << debug.distance 
                  << (debug.is_conflict_before ? " (CONFLICT)" : " (positive)") << "\n";
        std::cout << "  idx2[0] = " << debug.idx2_0 << ", idx2[1] = " << debug.idx2_1 << "\n";
        if (!debug.shared_with.empty()) {
            std::cout << "  Shared atoms with:\n";
            for (const auto& shared : debug.shared_with) {
                std::cout << "    H-bond " << shared.first << " (" << shared.second << ")\n";
            }
        }
    }
}

void print_phase3_trace(const std::vector<Phase3Debug>& trace) {
    std::cout << "\n========================================\n";
    std::cout << "PHASE 3: Linkage Type & Additional Conflicts\n";
    std::cout << "========================================\n";
    
    for (const auto& debug : trace) {
        std::cout << "\nH-bond " << debug.k << ": " << debug.donor << " -> " << debug.acceptor << "\n";
        std::cout << "  Distance before: " << std::fixed << std::setprecision(6) << debug.distance_before 
                  << (debug.is_conflict_before ? " (CONFLICT)" : " (positive)") << "\n";
        std::cout << "  Linkage type: " << debug.linkage_type << "\n";
        std::cout << "  Should mark conflict: " << (debug.should_mark_conflict ? "YES" : "NO") << "\n";
        std::cout << "  Distance after: " << debug.distance_after 
                  << (debug.is_conflict_after ? " (CONFLICT)" : " (positive)") << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <residue_i> <residue_j> <output_file>\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/1VBY.pdb 45 62 debug_output.txt\n";
        return 1;
    }
    
    std::string pdb_file = argv[1];
    int residue_i = std::stoi(argv[2]);
    int residue_j = std::stoi(argv[3]);
    std::string output_file = argv[4];
    
    // Open output file
    std::ofstream out(output_file);
    if (!out) {
        std::cerr << "Error: Cannot open output file: " << output_file << "\n";
        return 1;
    }
    
    // Redirect stdout to file
    std::streambuf* cout_buf = std::cout.rdbuf();
    std::cout.rdbuf(out.rdbuf());
    
    std::cout << "========================================\n";
    std::cout << "H-bond Conflict Resolution Debug\n";
    std::cout << "========================================\n";
    std::cout << "PDB: " << pdb_file << "\n";
    std::cout << "Pair: (" << residue_i << ", " << residue_j << ")\n";
    
    // Load structure
    PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(true);
    Structure structure = parser.parse_file(pdb_file);
    
    // Find residues by legacy index
    const Residue* res1 = nullptr;
    const Residue* res2 = nullptr;
    
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx == residue_i) {
                    res1 = &residue;
                }
                if (legacy_idx == residue_j) {
                    res2 = &residue;
                }
            }
        }
    }
    
    if (!res1 || !res2) {
        std::cerr << "Error: Could not find residues " << residue_i << " and " << residue_j << "\n";
        return 1;
    }
    
    std::cout << "\nResidue " << residue_i << ": " << res1->name() << " (chain " << res1->chain_id() << ")\n";
    std::cout << "Residue " << residue_j << ": " << res2->name() << " (chain " << res2->chain_id() << ")\n";
    
    // Get base types (use one_letter_code)
    char base1 = res1->one_letter_code();
    char base2 = res2->one_letter_code();
    if (base1 == '?') {
        // For modified nucleotides, check residue_type
        auto type1 = res1->residue_type();
        if (type1 == ResidueType::ADENINE) base1 = 'A';
        else if (type1 == ResidueType::CYTOSINE) base1 = 'C';
        else if (type1 == ResidueType::GUANINE) base1 = 'G';
        else if (type1 == ResidueType::THYMINE) base1 = 'T';
        else if (type1 == ResidueType::URACIL) base1 = 'U';
    }
    if (base2 == '?') {
        auto type2 = res2->residue_type();
        if (type2 == ResidueType::ADENINE) base2 = 'A';
        else if (type2 == ResidueType::CYTOSINE) base2 = 'C';
        else if (type2 == ResidueType::GUANINE) base2 = 'G';
        else if (type2 == ResidueType::THYMINE) base2 = 'T';
        else if (type2 == ResidueType::URACIL) base2 = 'U';
    }
    std::cout << "Base types: " << base1 << " - " << base2 << "\n";
    
    // Find H-bonds
    ValidationParameters params = ValidationParameters::defaults();
    double hb_lower = params.hb_lower;
    double hb_dist1 = params.hb_dist1;
    double hb_dist2 = 4.5;
    
    DetailedHBondResult detailed = HydrogenBondFinder::find_hydrogen_bonds_detailed(
        *res1, *res2, hb_lower, hb_dist1, hb_dist2);
    
    std::cout << "\n========================================\n";
    std::cout << "INITIAL H-BONDS (before conflict resolution)\n";
    std::cout << "========================================\n";
    for (size_t i = 0; i < detailed.initial_hbonds.size(); ++i) {
        const auto& hb = detailed.initial_hbonds[i];
        std::cout << "  " << (i+1) << ". " << hb.donor_atom << " -> " << hb.acceptor_atom 
                  << " dist=" << std::fixed << std::setprecision(6) << hb.distance << "\n";
    }
    
    // Now manually trace through conflict resolution
    std::vector<HydrogenBondResult> hbonds = detailed.initial_hbonds;
    std::vector<std::vector<int>> idx2(hbonds.size(), std::vector<int>(2, 0));
    
    // Phase 1: Initial conflict detection
    std::vector<Phase1Debug> phase1 = trace_phase1(hbonds, hb_lower, hb_dist2);
    print_phase1_trace(phase1);
    
    // Phase 2: idx2 population
    std::vector<Phase2Debug> phase2 = trace_phase2(hbonds, idx2);
    print_phase2_trace(phase2);
    
    // Phase 3: Linkage type and additional conflicts
    std::vector<Phase3Debug> phase3 = trace_phase3(hbonds, idx2, hb_lower, hb_dist2);
    print_phase3_trace(phase3);
    
    // Final state
    std::cout << "\n========================================\n";
    std::cout << "FINAL STATE (after conflict resolution)\n";
    std::cout << "========================================\n";
    for (size_t i = 0; i < hbonds.size(); ++i) {
        const auto& hb = hbonds[i];
        std::cout << "  " << (i+1) << ". " << hb.donor_atom << " -> " << hb.acceptor_atom 
                  << " dist=" << std::fixed << std::setprecision(6) << hb.distance 
                  << (hb.distance < 0.0 ? " (CONFLICT)" : " (positive)")
                  << " linkage=" << hb.linkage_type << "\n";
    }
    
    // Restore stdout
    std::cout.rdbuf(cout_buf);
    
    std::cout << "Debug output written to: " << output_file << "\n";
    return 0;
}


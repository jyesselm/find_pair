/**
 * @file compare_hbond_detectors.cpp
 * @brief Compare new HBondDetector against existing HydrogenBondFinder
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/hydrogen_bond_finder.hpp>
#include <x3dna/algorithms/hydrogen_bond.hpp>
#include <x3dna/core/typing/atom_classification.hpp>

using namespace x3dna;
using namespace x3dna::core;
using namespace x3dna::core::typing;
using namespace x3dna::algorithms;
using namespace x3dna::algorithms::hydrogen_bond;

// Filter to only base-base H-bonds
std::vector<HydrogenBondResult> filter_base_only(const std::vector<HydrogenBondResult>& hbonds) {
    std::vector<HydrogenBondResult> result;
    for (const auto& hb : hbonds) {
        bool atom1_base = AtomClassifier::is_nucleobase_atom(hb.donor_atom);
        bool atom2_base = AtomClassifier::is_nucleobase_atom(hb.acceptor_atom);
        if (atom1_base && atom2_base) {
            result.push_back(hb);
        }
    }
    return result;
}

void compare_pair(const Residue& res1, const Residue& res2, int pair_num, int& matches, int& mismatches) {
    // Use existing HydrogenBondFinder (baseline) with production params
    // hb_dist1=4.0, hb_dist2=0.0 (Phase 3 disabled in production)
    auto legacy_result = HydrogenBondFinder::find_hydrogen_bonds_detailed(res1, res2, 2.0, 4.0, 0.0);

    // Filter legacy to base-only
    auto legacy_base_only = filter_base_only(legacy_result.final_hbonds);

    // Use new HBondDetector with legacy_compatible params
    HBondDetector detector(HBondDetectionParams::legacy_compatible());
    auto new_result = detector.detect_base_hbonds_detailed(res1, res2);

    // Compare counts
    const size_t legacy_count = legacy_base_only.size();
    const size_t new_count = new_result.final_bonds.size();

    if (legacy_count != new_count) {
        std::cout << "Pair " << pair_num << " (" << res1.name() << "-" << res2.name() << "): "
                  << "COUNT MISMATCH - legacy=" << legacy_count << ", new=" << new_count << "\n";

        std::cout << "  Legacy base H-bonds:\n";
        for (const auto& hb : legacy_base_only) {
            std::cout << "    " << hb.donor_atom << " - " << hb.acceptor_atom
                      << " dist=" << std::fixed << std::setprecision(3) << hb.distance
                      << " type='" << hb.type << "'\n";
        }

        std::cout << "  New H-bonds:\n";
        for (const auto& hb : new_result.final_bonds) {
            std::cout << "    " << hb.donor_atom_name << " - " << hb.acceptor_atom_name
                      << " dist=" << std::fixed << std::setprecision(3) << hb.distance
                      << " type='" << to_legacy_char(hb.classification) << "'\n";
        }
        mismatches++;
        return;
    }

    // Compare individual H-bonds
    bool all_match = true;
    for (size_t i = 0; i < legacy_count; ++i) {
        const auto& leg = legacy_base_only[i];

        // Find matching H-bond in new results
        bool found = false;
        for (const auto& newb : new_result.final_bonds) {
            if (leg.donor_atom == newb.donor_atom_name &&
                leg.acceptor_atom == newb.acceptor_atom_name) {
                found = true;

                // Check type matches
                char new_type = to_legacy_char(newb.classification);
                if (leg.type != new_type) {
                    std::cout << "Pair " << pair_num << ": TYPE MISMATCH for "
                              << leg.donor_atom << "-" << leg.acceptor_atom
                              << " legacy='" << leg.type << "' new='" << new_type << "'\n";
                    all_match = false;
                }

                // Check distance matches
                if (std::abs(leg.distance - newb.distance) > 0.001) {
                    std::cout << "Pair " << pair_num << ": DISTANCE MISMATCH for "
                              << leg.donor_atom << "-" << leg.acceptor_atom
                              << " legacy=" << leg.distance << " new=" << newb.distance << "\n";
                    all_match = false;
                }
                break;
            }
        }

        if (!found) {
            std::cout << "Pair " << pair_num << ": MISSING in new - "
                      << leg.donor_atom << "-" << leg.acceptor_atom << "\n";
            all_match = false;
        }
    }

    if (all_match) {
        if (legacy_count > 0) {
            matches++;
        }
    } else {
        mismatches++;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file>\n";
        return 1;
    }

    std::string pdb_path = argv[1];

    // Parse PDB
    io::PdbParser parser;
    auto structure = parser.parse_file(pdb_path);

    std::cout << "Comparing H-bond detectors for " << pdb_path << "\n";
    std::cout << "Structure has " << structure.chains().size() << " chains\n\n";

    // Collect all residues
    std::vector<const Residue*> residues;
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            if (residue.is_nucleotide()) {
                residues.push_back(&residue);
            }
        }
    }

    std::cout << "Found " << residues.size() << " nucleotide residues\n\n";

    // Compare all pairs within distance cutoff
    int pair_num = 0;
    int matches = 0;
    int mismatches = 0;

    for (size_t i = 0; i < residues.size(); ++i) {
        for (size_t j = i + 1; j < residues.size(); ++j) {
            const auto* res1 = residues[i];
            const auto* res2 = residues[j];

            // Quick distance check using C1' atoms
            auto c1_1 = res1->find_atom("C1'");
            auto c1_2 = res2->find_atom("C1'");

            if (!c1_1 || !c1_2) continue;

            double dist = c1_1->position().distance_to(c1_2->position());
            if (dist > 15.0) continue;  // Skip pairs too far apart

            ++pair_num;
            compare_pair(*res1, *res2, pair_num, matches, mismatches);
        }
    }

    std::cout << "\n============================================================\n";
    std::cout << "BASE-ONLY H-BOND COMPARISON SUMMARY\n";
    std::cout << "============================================================\n";
    std::cout << "Total pairs checked: " << pair_num << "\n";
    std::cout << "Pairs with H-bonds matching: " << matches << "\n";
    std::cout << "Pairs with mismatches: " << mismatches << "\n";
    if (matches + mismatches > 0) {
        std::cout << "Match rate: " << std::fixed << std::setprecision(1)
                  << (100.0 * matches / (matches + mismatches)) << "%\n";
    }
    std::cout << "============================================================\n";

    return 0;
}

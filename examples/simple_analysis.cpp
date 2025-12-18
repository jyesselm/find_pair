/**
 * @file simple_analysis.cpp
 * @brief Simple example demonstrating x3dna library usage
 *
 * This example shows how to:
 * 1. Initialize the library
 * 2. Load a PDB structure
 * 3. Find base pairs
 * 4. Calculate step parameters
 */

#include <x3dna/all.hpp>
#include <iostream>
#include <iomanip>

int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file>\n";
        std::cerr << "Example: " << argv[0] << " 1ehz.pdb\n";
        return 1;
    }

    // Initialize the library (auto-detects resources)
    if (!x3dna::init()) {
        std::cerr << "Failed to initialize x3dna library.\n";
        std::cerr << "Make sure resources directory is accessible.\n";
        return 1;
    }

    std::cout << "x3dna library version: " << x3dna::version() << "\n\n";

    try {
        // Load the structure
        std::cout << "Loading structure: " << argv[1] << "\n";
        auto structure = x3dna::load_structure(argv[1]);

        std::cout << "Structure loaded:\n";
        std::cout << "  Chains: " << structure.chains().size() << "\n";
        std::cout << "  Residues: " << structure.num_residues() << "\n";
        std::cout << "  Atoms: " << structure.num_atoms() << "\n\n";

        // Find base pairs
        std::cout << "Finding base pairs...\n";
        auto pairs = x3dna::find_base_pairs(structure);

        std::cout << "Found " << pairs.size() << " base pairs:\n";
        for (size_t i = 0; i < std::min(pairs.size(), size_t(10)); ++i) {
            const auto& bp = pairs[i];
            std::cout << "  " << (i + 1) << ". "
                      << "residue " << bp.residue_idx1()
                      << " - residue " << bp.residue_idx2()
                      << " (" << bp.bp_type() << ")\n";
        }
        if (pairs.size() > 10) {
            std::cout << "  ... and " << (pairs.size() - 10) << " more\n";
        }
        std::cout << "\n";

        // Calculate step parameters
        if (pairs.size() >= 2) {
            std::cout << "Calculating step parameters...\n";
            auto step_params = x3dna::calculate_step_parameters(pairs);

            std::cout << "Step parameters (first 5 steps):\n";
            std::cout << std::fixed << std::setprecision(2);
            std::cout << "  Step  Shift  Slide   Rise   Tilt   Roll  Twist\n";
            std::cout << "  ----  -----  -----  -----  -----  -----  -----\n";

            for (size_t i = 0; i < std::min(step_params.size(), size_t(5)); ++i) {
                const auto& p = step_params[i];
                std::cout << "  " << std::setw(4) << (i + 1)
                          << std::setw(7) << p.shift
                          << std::setw(7) << p.slide
                          << std::setw(7) << p.rise
                          << std::setw(7) << p.tilt
                          << std::setw(7) << p.roll
                          << std::setw(7) << p.twist << "\n";
            }
        }

        std::cout << "\nAnalysis complete!\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

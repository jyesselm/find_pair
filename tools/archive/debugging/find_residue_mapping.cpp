/**
 * @file find_residue_mapping.cpp
 * @brief Find mapping between legacy and modern residue indices
 *
 * Legacy and modern count residues differently. This tool helps find
 * which modern residue index corresponds to a legacy residue index.
 */

#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/residue.hpp>
#include <iostream>
#include <map>
#include <iomanip>

using namespace x3dna::core;
using namespace x3dna::io;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> [legacy_residue_idx]\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/3G8T.pdb 946\n";
        return 1;
    }

    std::filesystem::path pdb_file = argv[1];
    int target_legacy_idx = -1;
    if (argc >= 3) {
        target_legacy_idx = std::stoi(argv[2]);
    }

    // Parse PDB
    PdbParser parser;
    auto structure = parser.parse_file(pdb_file);

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Residue Index Mapping\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "PDB: " << pdb_file << "\n\n";

    // Count residues modern way (how structure.chains() and chain.residues() work)
    size_t modern_idx = 1;
    std::map<std::tuple<std::string, char, int, char>, size_t> modern_residue_map;

    std::cout << "Modern Residue Counting (how modern code counts):\n";
    std::cout << std::string(60, '-') << "\n";

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            std::tuple<std::string, char, int, char> key = std::make_tuple(
                residue.name(), residue.chain_id(), residue.seq_num(), residue.insertion());
            modern_residue_map[key] = modern_idx;

            if (target_legacy_idx > 0 && modern_idx == static_cast<size_t>(target_legacy_idx)) {
                std::cout << ">>> ";
            }
            std::cout << std::setw(4) << modern_idx << ". " << std::setw(3) << residue.name()
                      << " (chain " << residue.chain_id() << ", seq " << std::setw(4)
                      << residue.seq_num();
            if (residue.insertion() != ' ') {
                std::cout << ", ins='" << residue.insertion() << "'";
            }
            std::cout << ")\n";
            modern_idx++;
        }
    }

    std::cout << "\nTotal modern residues: " << (modern_idx - 1) << "\n";

    // Now show what legacy would count (grouping by ResName, ChainID, ResSeq, insertion)
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Legacy Residue Counting (how legacy code counts)\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "Legacy groups atoms by: (ResName, ChainID, ResSeq, insertion_code)\n";
    std::cout << "This is different from modern counting!\n";
    std::cout << "\nTo find legacy residue index, need to:\n";
    std::cout << "  1. Group atoms by (ResName, ChainID, ResSeq, insertion)\n";
    std::cout << "  2. Count unique groups in order\n";
    std::cout << "  3. Map to modern indices\n";

    if (target_legacy_idx > 0) {
        std::cout << "\n" << std::string(60, '=') << "\n";
        std::cout << "Target: Legacy index " << target_legacy_idx << "\n";
        std::cout << std::string(60, '=') << "\n";
        std::cout << "Modern index " << target_legacy_idx << " = ";

        modern_idx = 1;
        for (const auto& chain : structure.chains()) {
            for (const auto& residue : chain.residues()) {
                if (modern_idx == static_cast<size_t>(target_legacy_idx)) {
                    std::cout << residue.name() << " (chain " << residue.chain_id() << ", seq "
                              << residue.seq_num() << ")\n";
                    std::cout << "\nBut legacy might count this differently!\n";
                    std::cout
                        << "Need to use legacy's residue_idx() function to get correct mapping.\n";
                    break;
                }
                modern_idx++;
            }
            if (modern_idx > static_cast<size_t>(target_legacy_idx))
                break;
        }
    }

    return 0;
}

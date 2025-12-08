/**
 * @file compare_residue_identification.cpp
 * @brief Compare residue identification (nucleotide recognition) between modern and legacy
 *
 * This tool focuses specifically on which residues are identified as nucleotides
 * and compares with legacy to find differences in residue type recognition.
 */

#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <map>
#include <set>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;

struct ResidueInfo {
    char chain_id;
    int seq_num;
    char insertion;
    std::string residue_name;
    bool is_nucleotide;
    ResidueType residue_type;
    size_t num_atoms;
    bool has_frame;

    std::string key() const {
        std::string key_str = std::string(1, chain_id) + ":" + std::to_string(seq_num);
        if (insertion != ' ') {
            key_str += insertion;
        }
        return key_str;
    }
};

// Helper function to check if residue is nucleotide (matches BasePairFinder::is_nucleotide logic)
bool is_nucleotide_with_ring_check(const Residue& residue) {
    ResidueType type = residue.residue_type();

    // Check standard nucleotide types
    if (type == ResidueType::ADENINE || type == ResidueType::CYTOSINE || type == ResidueType::GUANINE ||
        type == ResidueType::THYMINE || type == ResidueType::URACIL) {
        return true;
    }

    // Check for modified nucleotides (like XGR, XCR, XTR, XAR in 3KNC)
    // These have ResidueType::UNKNOWN but have ring atoms
    // Legacy code treats residues with ring atoms as nucleotides (RY >= 0)
    if (type == ResidueType::UNKNOWN) {
        // Check for ring atoms (similar to legacy residue_ident and generate_modern_json)
        static const std::vector<std::string> common_ring_atoms = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "};
        int ring_atom_count = 0;
        for (const auto& atom_name : common_ring_atoms) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == atom_name) {
                    ring_atom_count++;
                    break;
                }
            }
        }
        // If has >= 3 ring atoms, treat as nucleotide (matches legacy and generate_modern_json)
        if (ring_atom_count >= 3) {
            return true;
        }
    }

    return false;
}

std::vector<ResidueInfo> extract_modern_residues(const Structure& structure) {
    std::vector<ResidueInfo> residues;

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            ResidueInfo info;
            info.chain_id = residue.chain_id();
            info.seq_num = residue.seq_num();
            info.insertion = residue.insertion();
            info.residue_name = residue.name();
            // Use ring atom detection logic (matches BasePairFinder::is_nucleotide)
            info.is_nucleotide = is_nucleotide_with_ring_check(residue);
            info.residue_type = residue.residue_type();
            info.num_atoms = residue.num_atoms();
            info.has_frame = residue.reference_frame().has_value();

            residues.push_back(info);
        }
    }

    return residues;
}

std::vector<ResidueInfo> extract_legacy_residues(const std::filesystem::path& json_file) {
    std::vector<ResidueInfo> residues;
    std::map<std::string, bool> has_frame_map; // key -> has_frame

    if (!std::filesystem::exists(json_file)) {
        return residues;
    }

    std::ifstream file(json_file);
    if (!file.is_open()) {
        return residues;
    }

    nlohmann::json legacy_json;
    file >> legacy_json;

    // First pass: collect all residues from base_frame_calc and frame_calc
    std::map<std::string, ResidueInfo> residue_map;

    if (legacy_json.contains("calculations")) {
        for (const auto& calc : legacy_json["calculations"]) {
            if (calc.contains("type")) {
                std::string type = calc["type"];

                if (type == "base_frame_calc" || type == "frame_calc" || type == "ref_frame") {
                    if (calc.contains("chain_id") && calc.contains("residue_seq") && calc.contains("residue_name")) {
                        char chain_id = calc["chain_id"].get<std::string>()[0];
                        int residue_seq = calc["residue_seq"];
                        char insertion = calc.contains("insertion") ? calc["insertion"].get<std::string>()[0] : ' ';
                        std::string residue_name = calc["residue_name"];

                        ResidueInfo info;
                        info.chain_id = chain_id;
                        info.seq_num = residue_seq;
                        info.insertion = insertion;
                        info.residue_name = residue_name;
                        info.is_nucleotide = (type == "base_frame_calc" || type == "frame_calc");
                        info.has_frame = true;
                        info.num_atoms = 0; // Not available in legacy JSON

                        // Determine residue type from name
                        std::string name_upper = residue_name;
                        std::transform(name_upper.begin(), name_upper.end(), name_upper.begin(), ::toupper);
                        if (name_upper.find("A") != std::string::npos || name_upper == "  A" || name_upper == " DA") {
                            info.residue_type = ResidueType::ADENINE;
                        } else if (name_upper.find("C") != std::string::npos || name_upper == "  C" ||
                                   name_upper == " DC") {
                            info.residue_type = ResidueType::CYTOSINE;
                        } else if (name_upper.find("G") != std::string::npos || name_upper == "  G" ||
                                   name_upper == " DG") {
                            info.residue_type = ResidueType::GUANINE;
                        } else if (name_upper.find("T") != std::string::npos || name_upper == "  T" ||
                                   name_upper == " DT") {
                            info.residue_type = ResidueType::THYMINE;
                        } else if (name_upper.find("U") != std::string::npos || name_upper == "  U" ||
                                   name_upper == " DU") {
                            info.residue_type = ResidueType::URACIL;
                        } else {
                            info.residue_type = ResidueType::UNKNOWN;
                        }

                        std::string key = info.key();
                        if (residue_map.find(key) == residue_map.end()) {
                            residue_map[key] = info;
                        }
                    }
                }
            }
        }
    }

    // Convert map to vector
    for (const auto& [key, info] : residue_map) {
        residues.push_back(info);
    }

    return residues;
}

void print_residue_comparison(const std::vector<ResidueInfo>& modern, const std::vector<ResidueInfo>& legacy,
                              const std::string& pdb_id) {
    std::cout << "\n========================================\n";
    std::cout << "Residue Identification Comparison: " << pdb_id << "\n";
    std::cout << "========================================\n\n";

    // Create maps for easy lookup
    std::map<std::string, const ResidueInfo*> modern_map;
    std::map<std::string, const ResidueInfo*> legacy_map;

    for (const auto& res : modern) {
        modern_map[res.key()] = &res;
    }

    for (const auto& res : legacy) {
        legacy_map[res.key()] = &res;
    }

    // Count nucleotides
    size_t modern_nucleotides = 0;
    size_t legacy_nucleotides = 0;

    for (const auto& res : modern) {
        if (res.is_nucleotide)
            modern_nucleotides++;
    }

    for (const auto& res : legacy) {
        if (res.is_nucleotide)
            legacy_nucleotides++;
    }

    std::cout << "Nucleotide Recognition:\n";
    std::cout << "  Modern: " << modern_nucleotides << " / " << modern.size() << " residues\n";
    std::cout << "  Legacy: " << legacy_nucleotides << " / " << legacy.size() << " residues\n";
    if (modern_nucleotides != legacy_nucleotides) {
        std::cout << "  ⚠️  DIFFERENCE: " << std::abs((long long)modern_nucleotides - (long long)legacy_nucleotides)
                  << " nucleotides\n";
    } else {
        std::cout << "  ✅ Match\n";
    }
    std::cout << "\n";

    // Find residues that differ in nucleotide recognition
    std::vector<std::string> modern_only_nuc, legacy_only_nuc, different_recognition;

    for (const auto& [key, modern_res] : modern_map) {
        if (legacy_map.find(key) != legacy_map.end()) {
            const auto& legacy_res = *legacy_map[key];
            if (modern_res->is_nucleotide != legacy_res.is_nucleotide) {
                different_recognition.push_back(key);
            }
        } else if (modern_res->is_nucleotide) {
            modern_only_nuc.push_back(key);
        }
    }

    for (const auto& [key, legacy_res] : legacy_map) {
        if (modern_map.find(key) == modern_map.end() && legacy_res->is_nucleotide) {
            legacy_only_nuc.push_back(key);
        }
    }

    if (!different_recognition.empty()) {
        std::cout << "Residues with Different Nucleotide Recognition (" << different_recognition.size() << "):\n";
        for (const auto& key : different_recognition) {
            const auto& modern_res = *modern_map[key];
            const auto& legacy_res = *legacy_map[key];
            std::cout << "  " << key << " (" << modern_res.residue_name << "):\n";
            std::cout << "    Modern: " << (modern_res.is_nucleotide ? "nucleotide" : "not nucleotide")
                      << " (type=" << static_cast<int>(modern_res.residue_type) << ", atoms=" << modern_res.num_atoms
                      << ")\n";
            std::cout << "    Legacy: " << (legacy_res.is_nucleotide ? "nucleotide" : "not nucleotide")
                      << " (type=" << static_cast<int>(legacy_res.residue_type) << ")\n";
        }
        std::cout << "\n";
    }

    if (!modern_only_nuc.empty()) {
        std::cout << "Modern-only Nucleotides (" << modern_only_nuc.size() << "):\n";
        for (const auto& key : modern_only_nuc) {
            const auto& res = *modern_map[key];
            std::cout << "  " << key << " (" << res.residue_name << ", type=" << static_cast<int>(res.residue_type)
                      << ", atoms=" << res.num_atoms << ")\n";
        }
        std::cout << "\n";
    }

    if (!legacy_only_nuc.empty()) {
        std::cout << "Legacy-only Nucleotides (" << legacy_only_nuc.size() << "):\n";
        for (const auto& key : legacy_only_nuc) {
            const auto* res = legacy_map[key];
            std::cout << "  " << key << " (" << res->residue_name << ", type=" << static_cast<int>(res->residue_type)
                      << ")\n";
        }
        std::cout << "\n";
    }

    // Find residues missing frames
    std::vector<std::string> modern_no_frame;
    for (const auto& res : modern) {
        if (res.is_nucleotide && !res.has_frame) {
            modern_no_frame.push_back(res.key());
        }
    }

    if (!modern_no_frame.empty()) {
        std::cout << "Nucleotides without Frames in Modern (" << modern_no_frame.size() << "):\n";
        for (const auto& key : modern_no_frame) {
            const auto& res = *modern_map[key];
            std::cout << "  " << key << " (" << res.residue_name << ", atoms=" << res.num_atoms << ")\n";
        }
        std::cout << "\n";
    }

    // Summary
    std::cout << "Summary:\n";
    std::cout << "  Total residues (modern): " << modern.size() << "\n";
    std::cout << "  Total residues (legacy): " << legacy.size() << "\n";
    std::cout << "  Nucleotides (modern): " << modern_nucleotides << "\n";
    std::cout << "  Nucleotides (legacy): " << legacy_nucleotides << "\n";
    std::cout << "  Recognition differences: " << different_recognition.size() << "\n";
    std::cout << "  Modern-only nucleotides: " << modern_only_nuc.size() << "\n";
    std::cout << "  Legacy-only nucleotides: " << legacy_only_nuc.size() << "\n";
    std::cout << "  Modern nucleotides without frames: " << modern_no_frame.size() << "\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> [legacy_json_file]\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/3KNC.pdb data/json_legacy/3KNC.json\n";
        std::cerr << "\n";
        std::cerr << "If legacy_json_file is not provided, will look for it in data/json_legacy/\n";
        return 1;
    }

    std::filesystem::path pdb_file = argv[1];
    std::filesystem::path legacy_json_file;

    if (argc == 3) {
        legacy_json_file = argv[2];
    } else {
        // Auto-detect legacy JSON file
        std::filesystem::path legacy_dir = "data/json_legacy";
        legacy_json_file = legacy_dir / pdb_file.filename().replace_extension(".json");
    }

    if (!std::filesystem::exists(pdb_file)) {
        std::cerr << "Error: PDB file not found: " << pdb_file << "\n";
        return 1;
    }

    // Extract PDB ID from filename
    std::string pdb_id = pdb_file.stem().string();

    try {
        // Parse with modern parser
        std::cout << "Parsing PDB file with modern parser: " << pdb_file << "\n";
        PdbParser parser;
        parser.set_include_hetatm(true);
        parser.set_include_waters(true);

        Structure structure = parser.parse_file(pdb_file);

        // Calculate frames for all residues (needed for accurate comparison)
        BaseFrameCalculator calculator("data/templates");
        calculator.calculate_all_frames(structure);

        std::vector<ResidueInfo> modern_residues = extract_modern_residues(structure);

        // Extract legacy residues
        std::cout << "Extracting residues from legacy JSON: " << legacy_json_file << "\n";
        std::vector<ResidueInfo> legacy_residues = extract_legacy_residues(legacy_json_file);

        // Print comparison
        print_residue_comparison(modern_residues, legacy_residues, pdb_id);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

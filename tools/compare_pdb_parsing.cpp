/**
 * @file compare_pdb_parsing.cpp
 * @brief Compare PDB parsing between modern and legacy code
 * 
 * This tool parses a PDB file with the modern parser and compares the results
 * with legacy JSON output to identify parsing differences.
 */

#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <map>
#include <set>
#include <iomanip>

using namespace x3dna::core;
using namespace x3dna::io;

struct ParsingStats {
    size_t total_atoms = 0;
    size_t total_residues = 0;
    size_t nucleotide_residues = 0;
    std::map<std::string, size_t> residue_types; // residue name -> count
    std::map<char, size_t> chain_counts; // chain_id -> residue count
    std::set<std::tuple<char, int, char>> residues; // (chain_id, seq_num, insertion)
};

ParsingStats analyze_structure(const Structure& structure) {
    ParsingStats stats;
    
    for (const auto& chain : structure.chains()) {
        stats.chain_counts[chain.chain_id()] = chain.num_residues();
        
        for (const auto& residue : chain.residues()) {
            stats.total_residues++;
            stats.total_atoms += residue.num_atoms();
            
            auto key = std::make_tuple(residue.chain_id(), residue.seq_num(), residue.insertion());
            stats.residues.insert(key);
            
            stats.residue_types[residue.name()]++;
            
            if (residue.is_nucleotide()) {
                stats.nucleotide_residues++;
            }
        }
    }
    
    return stats;
}

ParsingStats analyze_legacy_json(const std::filesystem::path& json_file) {
    ParsingStats stats;
    
    if (!std::filesystem::exists(json_file)) {
        std::cerr << "Warning: Legacy JSON file not found: " << json_file << std::endl;
        return stats;
    }
    
    std::ifstream file(json_file);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open legacy JSON file: " << json_file << std::endl;
        return stats;
    }
    
    nlohmann::json legacy_json;
    file >> legacy_json;
    
    // Extract atom count from pdb_atoms record
    if (legacy_json.contains("calculations")) {
        for (const auto& calc : legacy_json["calculations"]) {
            if (calc.contains("type") && calc["type"] == "pdb_atoms") {
                if (calc.contains("atoms") && calc["atoms"].is_array()) {
                    stats.total_atoms = calc["atoms"].size();
                }
                break;
            }
        }
    }
    
    // Extract residue information from base_frame_calc and frame_calc records
    std::set<std::tuple<char, int, char>> seen_residues;
    
    if (legacy_json.contains("calculations")) {
        for (const auto& calc : legacy_json["calculations"]) {
            if (calc.contains("type")) {
                std::string type = calc["type"];
                
                if (type == "base_frame_calc" || type == "frame_calc" || type == "ref_frame") {
                    if (calc.contains("chain_id") && calc.contains("residue_seq") && 
                        calc.contains("residue_name")) {
                        char chain_id = calc["chain_id"].get<std::string>()[0];
                        int residue_seq = calc["residue_seq"];
                        char insertion = calc.contains("insertion") ? 
                                        calc["insertion"].get<std::string>()[0] : ' ';
                        std::string residue_name = calc["residue_name"];
                        
                        auto key = std::make_tuple(chain_id, residue_seq, insertion);
                        if (seen_residues.find(key) == seen_residues.end()) {
                            seen_residues.insert(key);
                            stats.total_residues++;
                            stats.residue_types[residue_name]++;
                            
                            // Check if it's a nucleotide (has base_frame_calc)
                            if (type == "base_frame_calc" || type == "frame_calc") {
                                stats.nucleotide_residues++;
                            }
                        }
                    }
                }
            }
        }
    }
    
    stats.residues = seen_residues;
    
    return stats;
}

void print_comparison(const ParsingStats& modern, const ParsingStats& legacy, 
                     const std::string& pdb_id) {
    std::cout << "\n========================================\n";
    std::cout << "PDB Parsing Comparison: " << pdb_id << "\n";
    std::cout << "========================================\n\n";
    
    // Atom count comparison
    std::cout << "Atom Counts:\n";
    std::cout << "  Modern: " << modern.total_atoms << "\n";
    std::cout << "  Legacy: " << legacy.total_atoms << "\n";
    if (modern.total_atoms != legacy.total_atoms) {
        std::cout << "  ⚠️  DIFFERENCE: " << std::abs((long long)modern.total_atoms - (long long)legacy.total_atoms) 
                  << " atoms\n";
    } else {
        std::cout << "  ✅ Match\n";
    }
    std::cout << "\n";
    
    // Residue count comparison
    std::cout << "Residue Counts:\n";
    std::cout << "  Modern: " << modern.total_residues << "\n";
    std::cout << "  Legacy: " << legacy.total_residues << "\n";
    if (modern.total_residues != legacy.total_residues) {
        std::cout << "  ⚠️  DIFFERENCE: " << std::abs((long long)modern.total_residues - (long long)legacy.total_residues) 
                  << " residues\n";
    } else {
        std::cout << "  ✅ Match\n";
    }
    std::cout << "\n";
    
    // Nucleotide count comparison
    std::cout << "Nucleotide Residues:\n";
    std::cout << "  Modern: " << modern.nucleotide_residues << "\n";
    std::cout << "  Legacy: " << legacy.nucleotide_residues << "\n";
    if (modern.nucleotide_residues != legacy.nucleotide_residues) {
        std::cout << "  ⚠️  DIFFERENCE: " << std::abs((long long)modern.nucleotide_residues - (long long)legacy.nucleotide_residues) 
                  << " nucleotides\n";
    } else {
        std::cout << "  ✅ Match\n";
    }
    std::cout << "\n";
    
    // Residue differences
    std::set<std::tuple<char, int, char>> modern_only, legacy_only;
    
    for (const auto& res : modern.residues) {
        if (legacy.residues.find(res) == legacy.residues.end()) {
            modern_only.insert(res);
        }
    }
    
    for (const auto& res : legacy.residues) {
        if (modern.residues.find(res) == modern.residues.end()) {
            legacy_only.insert(res);
        }
    }
    
    if (!modern_only.empty() || !legacy_only.empty()) {
        std::cout << "Residue Differences:\n";
        
        if (!modern_only.empty()) {
            std::cout << "  Modern-only residues (" << modern_only.size() << "):\n";
            for (const auto& [chain_id, seq_num, insertion] : modern_only) {
                std::cout << "    (" << chain_id << ", " << seq_num;
                if (insertion != ' ') std::cout << insertion;
                std::cout << ")\n";
            }
        }
        
        if (!legacy_only.empty()) {
            std::cout << "  Legacy-only residues (" << legacy_only.size() << "):\n";
            for (const auto& [chain_id, seq_num, insertion] : legacy_only) {
                std::cout << "    (" << chain_id << ", " << seq_num;
                if (insertion != ' ') std::cout << insertion;
                std::cout << ")\n";
            }
        }
    } else {
        std::cout << "Residue Sets: ✅ Match (all residues found in both)\n";
    }
    std::cout << "\n";
    
    // Residue type distribution
    std::cout << "Residue Type Distribution (Top 10):\n";
    std::vector<std::pair<std::string, size_t>> modern_types(modern.residue_types.begin(), 
                                                             modern.residue_types.end());
    std::sort(modern_types.begin(), modern_types.end(), 
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    for (size_t i = 0; i < std::min(modern_types.size(), size_t(10)); i++) {
        std::cout << "  " << std::setw(6) << modern_types[i].first << ": " 
                  << std::setw(4) << modern_types[i].second;
        if (legacy.residue_types.count(modern_types[i].first)) {
            std::cout << " (legacy: " << legacy.residue_types.at(modern_types[i].first) << ")";
        } else {
            std::cout << " (legacy: 0) ⚠️";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
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
        
        ParsingStats modern_stats = analyze_structure(structure);
        
        // Analyze legacy JSON
        std::cout << "Analyzing legacy JSON: " << legacy_json_file << "\n";
        ParsingStats legacy_stats = analyze_legacy_json(legacy_json_file);
        
        // Print comparison
        print_comparison(modern_stats, legacy_stats, pdb_id);
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}


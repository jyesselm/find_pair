/**
 * @file list_all_hbonds.cpp
 * @brief List all H-bonds detected for a specific pair
 * 
 * This tool lists ALL H-bonds found between two residues, showing:
 * - Atom names (donor/acceptor)
 * - Distance
 * - Type (standard '-' or non-standard '*')
 * - Whether it's "good" (distance in [2.5, 3.5] and type='-')
 * - Detailed information for debugging
 */

#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <iomanip>
#include <vector>
#include <algorithm>

using json = nlohmann::json;

struct HydrogenBondDetail {
    std::string donor_atom;
    std::string acceptor_atom;
    double distance;
    char type;  // '-' for standard, '*' for non-standard
    bool is_good;  // distance in [2.5, 3.5] and type='-'
    int hbond_idx;
};

void print_hbond_details(const std::vector<HydrogenBondDetail>& hbonds, const std::string& source) {
    std::cout << "\n" << source << " H-bonds (" << hbonds.size() << " total):\n";
    std::cout << "========================================\n";
    
    if (hbonds.empty()) {
        std::cout << "  (none)\n";
        return;
    }
    
    int good_count = 0;
    for (size_t i = 0; i < hbonds.size(); i++) {
        const auto& hb = hbonds[i];
        std::cout << "  " << std::setw(3) << (i+1) << ". ";
        std::cout << std::setw(6) << hb.donor_atom << " -> " << std::setw(6) << hb.acceptor_atom;
        std::cout << "  dist=" << std::fixed << std::setprecision(6) << std::setw(10) << hb.distance;
        std::cout << "  type=" << hb.type;
        std::cout << "  good=" << (hb.is_good ? "YES" : "NO ");
        if (hb.hbond_idx > 0) {
            std::cout << "  idx=" << hb.hbond_idx;
        }
        std::cout << "\n";
        
        if (hb.is_good) {
            good_count++;
        }
    }
    
    std::cout << "\n  Summary:\n";
    std::cout << "    Total H-bonds: " << hbonds.size() << "\n";
    std::cout << "    Good H-bonds (type='-' and dist in [2.5, 3.5]): " << good_count << "\n";
    std::cout << "    adjust_pairQuality: ";
    if (good_count >= 2) {
        std::cout << "-3.0 (2+ good H-bonds)\n";
    } else {
        std::cout << "-" << good_count << ".0 (" << good_count << " good H-bond";
        if (good_count != 1) std::cout << "s";
        std::cout << ")\n";
    }
}

std::vector<HydrogenBondDetail> extract_hbonds_from_json(const json& hbond_record) {
    std::vector<HydrogenBondDetail> hbonds;
    
    if (hbond_record.contains("hbonds") && hbond_record["hbonds"].is_array()) {
        for (const auto& hb : hbond_record["hbonds"]) {
            HydrogenBondDetail detail;
            
            // Handle both field name formats
            if (hb.contains("donor_atom") && hb.contains("acceptor_atom")) {
                detail.donor_atom = hb.value("donor_atom", "");
                detail.acceptor_atom = hb.value("acceptor_atom", "");
            } else {
                detail.donor_atom = hb.value("atom1_name", "");
                detail.acceptor_atom = hb.value("atom2_name", "");
            }
            
            detail.distance = hb.value("distance", 0.0);
            std::string type_str = hb.value("type", "-");
            detail.type = type_str.empty() ? '-' : type_str[0];
            detail.hbond_idx = hb.value("hbond_idx", 0);
            
            // Good H-bond: type='-' AND distance in [2.5, 3.5]
            detail.is_good = (detail.type == '-' && detail.distance >= 2.5 && detail.distance <= 3.5);
            
            hbonds.push_back(detail);
        }
    }
    
    return hbonds;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_id> <residue1_idx> <residue2_idx> [modern|legacy|both]\n";
        std::cerr << "Example: " << argv[0] << " 3G8T 92 160\n";
        std::cerr << "Example: " << argv[0] << " 3G8T 92 160 modern\n";
        std::cerr << "Example: " << argv[0] << " 6CAQ 75 78 both\n";
        std::cerr << "\n";
        std::cerr << "This tool lists ALL H-bonds detected for a specific pair.\n";
        std::cerr << "Default: shows modern only (if legacy file exists, shows both)\n";
        return 1;
    }
    
    std::string pdb_id = argv[1];
    int idx1 = std::stoi(argv[2]);
    int idx2 = std::stoi(argv[3]);
    std::string mode = (argc >= 5) ? argv[4] : "auto";
    
    std::filesystem::path modern_file = "data/json/" + pdb_id + "_hbond_list.json";
    std::filesystem::path legacy_file = "data/json_legacy/" + pdb_id + "_hbond_list.json";
    
    bool show_modern = (mode == "modern" || mode == "both" || mode == "auto");
    bool show_legacy = (mode == "legacy" || mode == "both");
    
    if (mode == "auto") {
        // Auto mode: show modern, and legacy if file exists
        show_legacy = std::filesystem::exists(legacy_file);
    }
    
    std::cout << "========================================\n";
    std::cout << "All H-bonds Detected\n";
    std::cout << "========================================\n";
    std::cout << "PDB: " << pdb_id << "\n";
    std::cout << "Pair: (" << idx1 << ", " << idx2 << ")\n";
    
    try {
        // Load modern H-bond records
        if (show_modern) {
            if (!std::filesystem::exists(modern_file)) {
                std::cerr << "Error: Modern H-bond JSON not found: " << modern_file << "\n";
                return 1;
            }
            
            std::ifstream modern_stream(modern_file);
            json modern_data = json::parse(modern_stream);
            
            // Find matching record
            std::vector<HydrogenBondDetail> modern_hbonds;
            for (const auto& record : modern_data) {
                int r1 = -1, r2 = -1;
                if (record.contains("base_i") && record.contains("base_j")) {
                    r1 = record["base_i"];
                    r2 = record["base_j"];
                } else if (record.contains("residue1_idx") && record.contains("residue2_idx")) {
                    r1 = record["residue1_idx"];
                    r2 = record["residue2_idx"];
                }
                
                if ((r1 == idx1 && r2 == idx2) || (r1 == idx2 && r2 == idx1)) {
                    modern_hbonds = extract_hbonds_from_json(record);
                    break;
                }
            }
            
            print_hbond_details(modern_hbonds, "Modern");
        }
        
        // Load legacy H-bond records
        if (show_legacy) {
            if (!std::filesystem::exists(legacy_file)) {
                std::cerr << "Warning: Legacy H-bond JSON not found: " << legacy_file << "\n";
            } else {
                std::ifstream legacy_stream(legacy_file);
                json legacy_data = json::parse(legacy_stream);
                
                // Find matching record
                std::vector<HydrogenBondDetail> legacy_hbonds;
                for (const auto& record : legacy_data) {
                    int r1 = -1, r2 = -1;
                    if (record.contains("base_i") && record.contains("base_j")) {
                        r1 = record["base_i"];
                        r2 = record["base_j"];
                    } else if (record.contains("residue1_idx") && record.contains("residue2_idx")) {
                        r1 = record["residue1_idx"];
                        r2 = record["residue2_idx"];
                    }
                    
                    if ((r1 == idx1 && r2 == idx2) || (r1 == idx2 && r2 == idx1)) {
                        legacy_hbonds = extract_hbonds_from_json(record);
                        break;
                    }
                }
                
                print_hbond_details(legacy_hbonds, "Legacy");
            }
        }
        
        std::cout << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}


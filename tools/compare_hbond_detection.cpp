/**
 * @file compare_hbond_detection.cpp
 * @brief Compare H-bond detection between legacy and modern for specific pairs
 * 
 * This tool focuses ONLY on H-bond detection, separate from base pair validation.
 * It compares:
 * 1. H-bonds found by legacy vs modern
 * 2. H-bond distances and types
 * 3. How many "good" H-bonds (distance in [2.5, 3.5]) are found
 * 4. The resulting adjust_pairQuality value
 * 
 * NOTE: This is separate from polygon overlap adjustment, which is a different step.
 */

#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>

using json = nlohmann::json;

struct HydrogenBondInfo {
    std::string atom1_name;
    std::string atom2_name;
    double distance;
    char type;  // '-' for standard, '*' for non-standard
    bool is_good;  // distance in [2.5, 3.5]
};

struct HBondComparison {
    int residue1_idx;
    int residue2_idx;
    std::vector<HydrogenBondInfo> modern_hbonds;
    std::vector<HydrogenBondInfo> legacy_hbonds;
    int modern_good_count;
    int legacy_good_count;
    double modern_adjust_pairQuality;
    double legacy_adjust_pairQuality;
};

std::vector<HydrogenBondInfo> extract_hbonds_from_json(const json& hbond_record) {
    std::vector<HydrogenBondInfo> hbonds;
    
    if (hbond_record.contains("hbonds") && hbond_record["hbonds"].is_array()) {
        for (const auto& hb : hbond_record["hbonds"]) {
            HydrogenBondInfo info;
            // Handle both field name formats
            if (hb.contains("donor_atom") && hb.contains("acceptor_atom")) {
                info.atom1_name = hb.value("donor_atom", "");
                info.atom2_name = hb.value("acceptor_atom", "");
            } else {
                info.atom1_name = hb.value("atom1_name", "");
                info.atom2_name = hb.value("atom2_name", "");
            }
            info.distance = hb.value("distance", 0.0);
            std::string type_str = hb.value("type", "-");
            info.type = type_str.empty() ? '-' : type_str[0];
            info.is_good = (info.distance >= 2.5 && info.distance <= 3.5 && info.type == '-');
            hbonds.push_back(info);
        }
    }
    
    return hbonds;
}

HBondComparison extract_hbond_comparison(const json& modern_hbond_records,
                                         const json& legacy_hbond_records,
                                         int idx1, int idx2) {
    HBondComparison comp;
    comp.residue1_idx = idx1;
    comp.residue2_idx = idx2;
    comp.modern_good_count = 0;
    comp.legacy_good_count = 0;
    comp.modern_adjust_pairQuality = 0.0;
    comp.legacy_adjust_pairQuality = 0.0;
    
    // Find modern H-bonds
    for (const auto& record : modern_hbond_records) {
        int r1 = -1, r2 = -1;
        if (record.contains("base_i") && record.contains("base_j")) {
            r1 = record["base_i"];
            r2 = record["base_j"];
        } else if (record.contains("residue1_idx") && record.contains("residue2_idx")) {
            r1 = record["residue1_idx"];
            r2 = record["residue2_idx"];
        }
        
        if ((r1 == idx1 && r2 == idx2) || (r1 == idx2 && r2 == idx1)) {
            comp.modern_hbonds = extract_hbonds_from_json(record);
            for (const auto& hb : comp.modern_hbonds) {
                if (hb.is_good) {
                    comp.modern_good_count++;
                }
            }
            // Calculate adjust_pairQuality: -3.0 if >= 2 good, else -count
            if (comp.modern_good_count >= 2) {
                comp.modern_adjust_pairQuality = -3.0;
            } else {
                comp.modern_adjust_pairQuality = -static_cast<double>(comp.modern_good_count);
            }
            break;
        }
    }
    
    // Find legacy H-bonds
    for (const auto& record : legacy_hbond_records) {
        int r1 = -1, r2 = -1;
        if (record.contains("base_i") && record.contains("base_j")) {
            r1 = record["base_i"];
            r2 = record["base_j"];
        } else if (record.contains("residue1_idx") && record.contains("residue2_idx")) {
            r1 = record["residue1_idx"];
            r2 = record["residue2_idx"];
        }
        
        if ((r1 == idx1 && r2 == idx2) || (r1 == idx2 && r2 == idx1)) {
            comp.legacy_hbonds = extract_hbonds_from_json(record);
            for (const auto& hb : comp.legacy_hbonds) {
                if (hb.is_good) {
                    comp.legacy_good_count++;
                }
            }
            // Calculate adjust_pairQuality: -3.0 if >= 2 good, else -count
            if (comp.legacy_good_count >= 2) {
                comp.legacy_adjust_pairQuality = -3.0;
            } else {
                comp.legacy_adjust_pairQuality = -static_cast<double>(comp.legacy_good_count);
            }
            break;
        }
    }
    
    return comp;
}

void print_hbond_comparison(const HBondComparison& comp) {
    std::cout << "\n========================================\n";
    std::cout << "H-bond Detection Comparison\n";
    std::cout << "========================================\n";
    std::cout << "Pair: (" << comp.residue1_idx << ", " << comp.residue2_idx << ")\n\n";
    
    std::cout << "Modern H-bonds found: " << comp.modern_hbonds.size() << "\n";
    std::cout << "Legacy H-bonds found: " << comp.legacy_hbonds.size() << "\n";
    std::cout << "Difference: " << (static_cast<int>(comp.modern_hbonds.size()) - static_cast<int>(comp.legacy_hbonds.size())) << "\n\n";
    
    std::cout << "Good H-bonds (distance in [2.5, 3.5]):\n";
    std::cout << "  Modern: " << comp.modern_good_count << "\n";
    std::cout << "  Legacy: " << comp.legacy_good_count << "\n";
    std::cout << "  Difference: " << (comp.modern_good_count - comp.legacy_good_count) << "\n\n";
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "adjust_pairQuality:\n";
    std::cout << "  Modern: " << comp.modern_adjust_pairQuality << "\n";
    std::cout << "  Legacy: " << comp.legacy_adjust_pairQuality << "\n";
    std::cout << "  Difference: " << (comp.modern_adjust_pairQuality - comp.legacy_adjust_pairQuality) << "\n\n";
    
    // Print detailed H-bond lists
    std::cout << "Modern H-bonds:\n";
    if (comp.modern_hbonds.empty()) {
        std::cout << "  (none)\n";
    } else {
        for (size_t i = 0; i < comp.modern_hbonds.size(); i++) {
            const auto& hb = comp.modern_hbonds[i];
            std::cout << "  " << (i+1) << ". " << hb.atom1_name << " - " << hb.atom2_name
                      << " (dist=" << std::setprecision(3) << hb.distance
                      << ", type=" << hb.type << ", good=" << (hb.is_good ? "yes" : "no") << ")\n";
        }
    }
    
    std::cout << "\nLegacy H-bonds:\n";
    if (comp.legacy_hbonds.empty()) {
        std::cout << "  (none)\n";
    } else {
        for (size_t i = 0; i < comp.legacy_hbonds.size(); i++) {
            const auto& hb = comp.legacy_hbonds[i];
            std::cout << "  " << (i+1) << ". " << hb.atom1_name << " - " << hb.atom2_name
                      << " (dist=" << std::setprecision(3) << hb.distance
                      << ", type=" << hb.type << ", good=" << (hb.is_good ? "yes" : "no") << ")\n";
        }
    }
    
    // Find missing/extra H-bonds
    std::cout << "\nH-bond Differences:\n";
    std::map<std::string, HydrogenBondInfo> modern_map;
    for (const auto& hb : comp.modern_hbonds) {
        std::string key = hb.atom1_name + "-" + hb.atom2_name;
        modern_map[key] = hb;
    }
    
    std::map<std::string, HydrogenBondInfo> legacy_map;
    for (const auto& hb : comp.legacy_hbonds) {
        std::string key = hb.atom1_name + "-" + hb.atom2_name;
        legacy_map[key] = hb;
    }
    
    std::vector<std::string> missing_in_modern;
    for (const auto& [key, hb] : legacy_map) {
        if (modern_map.find(key) == modern_map.end()) {
            missing_in_modern.push_back(key);
        }
    }
    
    std::vector<std::string> extra_in_modern;
    for (const auto& [key, hb] : modern_map) {
        if (legacy_map.find(key) == legacy_map.end()) {
            extra_in_modern.push_back(key);
        }
    }
    
    if (missing_in_modern.empty() && extra_in_modern.empty()) {
        std::cout << "  ✓ All H-bonds match\n";
    } else {
        if (!missing_in_modern.empty()) {
            std::cout << "  Missing in modern (" << missing_in_modern.size() << "):\n";
            for (const auto& key : missing_in_modern) {
                const auto& hb = legacy_map[key];
                std::cout << "    - " << key << " (dist=" << std::setprecision(3) << hb.distance
                          << ", type=" << hb.type << ", good=" << (hb.is_good ? "yes" : "no") << ")\n";
            }
        }
        if (!extra_in_modern.empty()) {
            std::cout << "  Extra in modern (" << extra_in_modern.size() << "):\n";
            for (const auto& key : extra_in_modern) {
                const auto& hb = modern_map[key];
                std::cout << "    + " << key << " (dist=" << std::setprecision(3) << hb.distance
                          << ", type=" << hb.type << ", good=" << (hb.is_good ? "yes" : "no") << ")\n";
            }
        }
    }
    
    if (comp.modern_adjust_pairQuality != comp.legacy_adjust_pairQuality) {
        std::cout << "\n⚠️  adjust_pairQuality MISMATCH!\n";
        std::cout << "   This will cause quality score differences\n";
        std::cout << "   Quality score difference: " << (comp.modern_adjust_pairQuality - comp.legacy_adjust_pairQuality) << "\n";
    } else {
        std::cout << "\n✓ adjust_pairQuality matches\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_id> <residue1_idx> <residue2_idx>\n";
        std::cerr << "Example: " << argv[0] << " 3G8T 92 160\n";
        std::cerr << "Example: " << argv[0] << " 6CAQ 75 78\n";
        std::cerr << "\n";
        std::cerr << "This tool compares H-bond detection ONLY (separate from polygon overlap).\n";
        return 1;
    }
    
    std::string pdb_id = argv[1];
    int idx1 = std::stoi(argv[2]);
    int idx2 = std::stoi(argv[3]);
    
    std::filesystem::path modern_file = "data/json/" + pdb_id + "_hbond_list.json";
    std::filesystem::path legacy_file = "data/json_legacy/" + pdb_id + "_hbond_list.json";
    
    if (!std::filesystem::exists(modern_file)) {
        std::cerr << "Error: Modern H-bond JSON not found: " << modern_file << "\n";
        return 1;
    }
    
    json legacy_data = json::array(); // Empty if legacy file doesn't exist
    if (std::filesystem::exists(legacy_file)) {
        std::ifstream legacy_stream(legacy_file);
        legacy_data = json::parse(legacy_stream);
    } else {
        std::cerr << "Warning: Legacy H-bond JSON not found: " << legacy_file << "\n";
        std::cerr << "         Will only show modern H-bond detection.\n";
    }
    
    try {
        // Load modern H-bond records
        std::ifstream modern_stream(modern_file);
        json modern_data = json::parse(modern_stream);
        
        // Extract H-bond comparison
        HBondComparison comp = extract_hbond_comparison(modern_data, legacy_data, idx1, idx2);
        
        // Print comparison
        print_hbond_comparison(comp);
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

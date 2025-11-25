/**
 * @file compare_hbond_stages.cpp
 * @brief Compare H-bond detection at ALL stages between legacy and modern
 *
 * This tool provides comprehensive comparison at multiple stages:
 * 1. Initial detection (before conflict resolution)
 * 2. After conflict resolution
 * 3. After validation
 * 4. Atom-by-atom comparison
 */

#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/algorithms/hydrogen_bond_finder.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <iomanip>
#include <map>
#include <set>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using json = nlohmann::json;

struct HBondInfo {
    std::string donor_atom;
    std::string acceptor_atom;
    double distance;
    char type;
    int linkage_type;

    // Normalize atom name (trim spaces, handle PDB format)
    static std::string normalize_atom_name(const std::string& name) {
        std::string normalized = name;
        // Remove leading/trailing spaces
        while (!normalized.empty() && normalized[0] == ' ') {
            normalized.erase(0, 1);
        }
        while (!normalized.empty() && normalized.back() == ' ') {
            normalized.pop_back();
        }
        return normalized;
    }

    bool operator==(const HBondInfo& other) const {
        std::string norm_donor1 = normalize_atom_name(donor_atom);
        std::string norm_acceptor1 = normalize_atom_name(acceptor_atom);
        std::string norm_donor2 = normalize_atom_name(other.donor_atom);
        std::string norm_acceptor2 = normalize_atom_name(other.acceptor_atom);

        return norm_donor1 == norm_donor2 && norm_acceptor1 == norm_acceptor2 &&
               std::abs(distance - other.distance) < 0.01;
    }
};

void print_stage_comparison(const std::string& stage_name,
                            const std::vector<HBondInfo>& modern_hbonds,
                            const std::vector<HBondInfo>& legacy_hbonds) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << stage_name << "\n";
    std::cout << std::string(60, '=') << "\n";

    std::cout << "\nModern: " << modern_hbonds.size() << " H-bonds\n";
    for (size_t i = 0; i < modern_hbonds.size(); i++) {
        const auto& hb = modern_hbonds[i];
        std::cout << "  " << (i + 1) << ". " << std::setw(6) << hb.donor_atom << " -> "
                  << std::setw(6) << hb.acceptor_atom << ", dist=" << std::fixed
                  << std::setprecision(3) << hb.distance;
        if (hb.type != ' ') {
            std::cout << ", type=" << hb.type;
        }
        if (hb.linkage_type != 0) {
            std::cout << ", lkg=" << hb.linkage_type;
        }
        std::cout << " [repr: donor='" << hb.donor_atom << "' acceptor='" << hb.acceptor_atom
                  << "']";
        std::cout << "\n";
    }

    std::cout << "\nLegacy: " << legacy_hbonds.size() << " H-bonds\n";
    for (size_t i = 0; i < legacy_hbonds.size(); i++) {
        const auto& hb = legacy_hbonds[i];
        std::cout << "  " << (i + 1) << ". " << std::setw(6) << hb.donor_atom << " -> "
                  << std::setw(6) << hb.acceptor_atom << ", dist=" << std::fixed
                  << std::setprecision(3) << hb.distance;
        if (hb.type != ' ') {
            std::cout << ", type=" << hb.type;
        }
        if (hb.linkage_type != 0) {
            std::cout << ", lkg=" << hb.linkage_type;
        }
        std::cout << " [repr: donor='" << hb.donor_atom << "' acceptor='" << hb.acceptor_atom
                  << "']";
        std::cout << "\n";
    }

    // Find matches
    std::vector<bool> modern_matched(modern_hbonds.size(), false);
    std::vector<bool> legacy_matched(legacy_hbonds.size(), false);
    int matches = 0;

    for (size_t i = 0; i < modern_hbonds.size(); i++) {
        for (size_t j = 0; j < legacy_hbonds.size(); j++) {
            if (!legacy_matched[j] && modern_hbonds[i] == legacy_hbonds[j]) {
                modern_matched[i] = true;
                legacy_matched[j] = true;
                matches++;
                break;
            }
        }
    }

    std::cout << "\nMatches: " << matches << " / "
              << std::max(modern_hbonds.size(), legacy_hbonds.size()) << "\n";

    if (static_cast<size_t>(matches) < modern_hbonds.size() ||
        static_cast<size_t>(matches) < legacy_hbonds.size()) {
        std::cout << "\nMissing in modern:\n";
        for (size_t i = 0; i < legacy_hbonds.size(); i++) {
            if (!legacy_matched[i]) {
                const auto& hb = legacy_hbonds[i];
                std::cout << "  - " << hb.donor_atom << " -> " << hb.acceptor_atom
                          << " (dist=" << std::fixed << std::setprecision(3) << hb.distance
                          << ")\n";
            }
        }

        std::cout << "\nExtra in modern:\n";
        for (size_t i = 0; i < modern_hbonds.size(); i++) {
            if (!modern_matched[i]) {
                const auto& hb = modern_hbonds[i];
                std::cout << "  + " << hb.donor_atom << " -> " << hb.acceptor_atom
                          << " (dist=" << std::fixed << std::setprecision(3) << hb.distance
                          << ")\n";
            }
        }
    }
}

std::vector<HBondInfo>
extract_modern_initial(const x3dna::algorithms::DetailedHBondResult& result) {
    std::vector<HBondInfo> hbonds;
    for (const auto& hb : result.initial_hbonds) {
        HBondInfo info;
        info.donor_atom = hb.donor_atom;
        info.acceptor_atom = hb.acceptor_atom;
        info.distance = std::abs(hb.distance);
        info.type = '-'; // Initial H-bonds don't have types yet
        info.linkage_type = 0;
        hbonds.push_back(info);
    }
    return hbonds;
}

std::vector<HBondInfo>
extract_modern_after_conflict(const x3dna::algorithms::DetailedHBondResult& result) {
    std::vector<HBondInfo> hbonds;
    for (const auto& hb : result.after_conflict_resolution) {
        HBondInfo info;
        info.donor_atom = hb.donor_atom;
        info.acceptor_atom = hb.acceptor_atom;
        info.distance = std::abs(hb.distance);
        info.type = '-'; // Types assigned in validation
        info.linkage_type = hb.linkage_type;
        hbonds.push_back(info);
    }
    return hbonds;
}

std::vector<HBondInfo>
extract_modern_after_validation(const x3dna::algorithms::DetailedHBondResult& result) {
    std::vector<HBondInfo> hbonds;
    for (const auto& hb : result.after_validation) {
        HBondInfo info;
        info.donor_atom = hb.donor_atom;
        info.acceptor_atom = hb.acceptor_atom;
        info.distance = std::abs(hb.distance);
        info.type = hb.type;
        info.linkage_type = hb.linkage_type;
        hbonds.push_back(info);
    }
    return hbonds;
}

std::vector<HBondInfo> extract_legacy_from_json(const json& legacy_record) {
    std::vector<HBondInfo> hbonds;

    if (legacy_record.contains("hbonds") && legacy_record["hbonds"].is_array()) {
        for (const auto& hb : legacy_record["hbonds"]) {
            HBondInfo info;
            info.donor_atom = hb.value("donor_atom", "");
            info.acceptor_atom = hb.value("acceptor_atom", "");
            info.distance = std::abs(hb.value("distance", 0.0));
            std::string type_str = hb.value("type", " ");
            info.type = type_str.empty() ? ' ' : type_str[0];
            info.linkage_type = hb.value("linkage_type", 0);
            hbonds.push_back(info);
        }
    }

    return hbonds;
}

json find_legacy_pair(const std::string& content, int residue1_idx, int residue2_idx) {
    std::string search1 = "\"base_i\": " + std::to_string(residue1_idx);
    std::string search2 = "\"base_j\": " + std::to_string(residue2_idx);

    size_t pair_pos = content.find(search1);
    if (pair_pos != std::string::npos) {
        size_t check_pos2 = content.find(search2, pair_pos);
        if (check_pos2 != std::string::npos && check_pos2 < pair_pos + 200) {
            size_t obj_start = content.rfind('{', pair_pos);
            if (obj_start != std::string::npos) {
                size_t type_pos = content.find("\"type\": \"hbond_list\"", obj_start);
                if (type_pos == std::string::npos) {
                    type_pos = content.find("\"type\":\"hbond_list\"", obj_start);
                }
                if (type_pos != std::string::npos && type_pos < pair_pos + 500) {
                    int brace_count = 0;
                    bool in_string = false;
                    bool escape_next = false;
                    size_t obj_end = obj_start;

                    for (size_t i = obj_start; i < content.length(); i++) {
                        char c = content[i];
                        if (escape_next) {
                            escape_next = false;
                            continue;
                        }
                        if (c == '\\') {
                            escape_next = true;
                            continue;
                        }
                        if (c == '"') {
                            in_string = !in_string;
                            continue;
                        }
                        if (!in_string) {
                            if (c == '{')
                                brace_count++;
                            else if (c == '}') {
                                brace_count--;
                                if (brace_count == 0) {
                                    obj_end = i + 1;
                                    break;
                                }
                            }
                        }
                    }

                    std::string obj_str = content.substr(obj_start, obj_end - obj_start);
                    try {
                        json obj = json::parse(obj_str);
                        if (obj.value("base_i", -1) == residue1_idx &&
                            obj.value("base_j", -1) == residue2_idx) {
                            return obj;
                        }
                    } catch (const json::parse_error&) {
                        // Continue
                    }
                }
            }
        }
    }

    // Try reverse order
    search1 = "\"base_i\": " + std::to_string(residue2_idx);
    search2 = "\"base_j\": " + std::to_string(residue1_idx);

    pair_pos = content.find(search1);
    if (pair_pos != std::string::npos) {
        size_t check_pos2 = content.find(search2, pair_pos);
        if (check_pos2 != std::string::npos && check_pos2 < pair_pos + 200) {
            size_t obj_start = content.rfind('{', pair_pos);
            if (obj_start != std::string::npos) {
                size_t type_pos = content.find("\"type\": \"hbond_list\"", obj_start);
                if (type_pos == std::string::npos) {
                    type_pos = content.find("\"type\":\"hbond_list\"", obj_start);
                }
                if (type_pos != std::string::npos && type_pos < pair_pos + 500) {
                    int brace_count = 0;
                    bool in_string = false;
                    bool escape_next = false;
                    size_t obj_end = obj_start;

                    for (size_t i = obj_start; i < content.length(); i++) {
                        char c = content[i];
                        if (escape_next) {
                            escape_next = false;
                            continue;
                        }
                        if (c == '\\') {
                            escape_next = true;
                            continue;
                        }
                        if (c == '"') {
                            in_string = !in_string;
                            continue;
                        }
                        if (!in_string) {
                            if (c == '{')
                                brace_count++;
                            else if (c == '}') {
                                brace_count--;
                                if (brace_count == 0) {
                                    obj_end = i + 1;
                                    break;
                                }
                            }
                        }
                    }

                    std::string obj_str = content.substr(obj_start, obj_end - obj_start);
                    try {
                        json obj = json::parse(obj_str);
                        if (obj.value("base_i", -1) == residue2_idx &&
                            obj.value("base_j", -1) == residue1_idx) {
                            return obj;
                        }
                    } catch (const json::parse_error&) {
                        // Continue
                    }
                }
            }
        }
    }

    return json();
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <pdb_file> <residue1_idx> <residue2_idx> [legacy_json]\n";
        return 1;
    }

    std::filesystem::path pdb_file = argv[1];
    int residue1_idx = std::stoi(argv[2]);
    int residue2_idx = std::stoi(argv[3]);

    // Parse PDB - include all residues to match legacy (HETATMs, waters, etc.)
    PdbParser parser;
    parser.set_include_hetatm(true); // Include HETATM records (legacy includes all)
    parser.set_include_waters(true); // Include water molecules (legacy includes all)
    auto structure = parser.parse_file(pdb_file);

    // Get residues using legacy residue index (PDB file order)
    // Use Structure's built-in method to get residues in legacy order
    const Residue* res1 = structure.get_residue_by_legacy_idx(residue1_idx);
    const Residue* res2 = structure.get_residue_by_legacy_idx(residue2_idx);

    if (!res1 || !res2) {
        std::cerr << "Error: Could not find residues " << residue1_idx << " and " << residue2_idx
                  << "\n";
        return 1;
    }

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Comprehensive H-bond Comparison\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "Pair: (" << residue1_idx << ", " << residue2_idx << ")\n";
    std::cout << "Residue 1: " << res1->name() << " (chain " << res1->chain_id() << ", seq "
              << res1->seq_num() << ")\n";
    std::cout << "Residue 2: " << res2->name() << " (chain " << res2->chain_id() << ", seq "
              << res2->seq_num() << ")\n";

    // Calculate frames
    BaseFrameCalculator calculator;
    calculator.calculate_frame_const(*res1);
    calculator.calculate_frame_const(*res2);

    // Get modern H-bonds at all stages
    ValidationParameters params = ValidationParameters::defaults();
    auto detailed = HydrogenBondFinder::find_hydrogen_bonds_detailed(
        *res1, *res2, params.hb_lower, params.hb_dist1,
        params.hb_dist1); // hb_dist2 = hb_dist1 for now

    // Extract at each stage
    auto modern_initial = extract_modern_initial(detailed);
    auto modern_after_conflict = extract_modern_after_conflict(detailed);
    auto modern_after_validation = extract_modern_after_validation(detailed);

    // Get legacy H-bonds (only have after validation from JSON)
    std::vector<HBondInfo> legacy_after_validation;
    if (argc >= 5) {
        std::filesystem::path legacy_json = argv[4];
        if (std::filesystem::exists(legacy_json)) {
            std::ifstream f(legacy_json);
            std::string content((std::istreambuf_iterator<char>(f)),
                                std::istreambuf_iterator<char>());

            json legacy_obj = find_legacy_pair(content, residue1_idx, residue2_idx);
            if (!legacy_obj.empty()) {
                legacy_after_validation = extract_legacy_from_json(legacy_obj);
            }
        }
    }

    // Print comparisons
    print_stage_comparison("Stage 1: Initial Detection (before conflict resolution)",
                           modern_initial, {}); // Legacy initial not available from JSON

    print_stage_comparison("Stage 2: After Conflict Resolution", modern_after_conflict,
                           {}); // Legacy after conflict not available from JSON

    print_stage_comparison("Stage 3: After Validation (final)", modern_after_validation,
                           legacy_after_validation);

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Note: Legacy initial and after-conflict stages not available from JSON.\n";
    std::cout << "      Need to add debug output to legacy code to compare those stages.\n";
    std::cout << std::string(60, '=') << "\n";

    return 0;
}

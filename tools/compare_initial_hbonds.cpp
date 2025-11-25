/**
 * @file compare_initial_hbonds.cpp
 * @brief Compare H-bond detection at the SAME step between legacy and modern
 *
 * This tool compares H-bonds at the step that legacy records to JSON:
 * - Legacy records ALL H-bonds AFTER hb_atompair (conflict resolution) and validate_hbonds
 * - This includes H-bonds with type=' ' (which legacy records but doesn't use)
 * - Modern equivalent: after_validation (all H-bonds with types assigned)
 *
 * This is the "initial batch" of H-bonds that were found, after conflict resolution
 * and validation, which is what legacy records to JSON.
 */

#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/hydrogen_bond_finder.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>

using json = nlohmann::json;
using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;

struct InitialHBond {
    std::string atom1_name;
    std::string atom2_name;
    double distance;

    bool operator<(const InitialHBond& other) const {
        if (atom1_name != other.atom1_name)
            return atom1_name < other.atom1_name;
        if (atom2_name != other.atom2_name)
            return atom2_name < other.atom2_name;
        return distance < other.distance;
    }

    bool operator==(const InitialHBond& other) const {
        return atom1_name == other.atom1_name && atom2_name == other.atom2_name &&
               std::abs(distance - other.distance) < 0.001;
    }
};

std::vector<InitialHBond> get_modern_initial_hbonds(const Residue& res1, const Residue& res2,
                                                    double hb_lower, double hb_dist1) {
    std::vector<InitialHBond> hbonds;

    // Get H-bonds at the SAME step as legacy records
    // Legacy records ALL H-bonds AFTER validation (with types assigned)
    // So we should compare modern's after_validation (not initial_hbonds)
    DetailedHBondResult detailed =
        HydrogenBondFinder::find_hydrogen_bonds_detailed(res1, res2, hb_lower, hb_dist1, 4.5);

    // Extract after_validation H-bonds (matches what legacy records to JSON)
    // Legacy records ALL num_hbonds H-bonds after hb_atompair and validate_hbonds
    // This includes H-bonds with type=' ' (which legacy records but doesn't use)
    for (const auto& hbond : detailed.after_validation) {
        InitialHBond hb;
        hb.atom1_name = hbond.donor_atom;
        hb.atom2_name = hbond.acceptor_atom;
        // Use absolute distance (legacy uses fabs(hb_dist[i]) in JSON)
        hb.distance = std::abs(hbond.distance);
        hbonds.push_back(hb);
    }

    return hbonds;
}

std::vector<InitialHBond> extract_legacy_initial_hbonds(const json& legacy_record) {
    std::vector<InitialHBond> hbonds;

    // Legacy records all H-bonds in hbond_list, including those with type=' '
    // The initial H-bonds are those found before validation
    // In legacy, initial H-bonds are those that pass distance and good_hbatoms checks
    // We can identify them as all H-bonds in the record (legacy doesn't filter them out)

    if (legacy_record.contains("hbonds") && legacy_record["hbonds"].is_array()) {
        for (const auto& hb : legacy_record["hbonds"]) {
            InitialHBond hbond;

            // Try different field name formats
            if (hb.contains("donor_atom") && hb.contains("acceptor_atom")) {
                hbond.atom1_name = hb.value("donor_atom", "");
                hbond.atom2_name = hb.value("acceptor_atom", "");
            } else if (hb.contains("atom1_name") && hb.contains("atom2_name")) {
                hbond.atom1_name = hb.value("atom1_name", "");
                hbond.atom2_name = hb.value("atom2_name", "");
            } else if (hb.contains("donor") && hb.contains("acceptor")) {
                hbond.atom1_name = hb.value("donor", "");
                hbond.atom2_name = hb.value("acceptor", "");
            } else {
                continue; // Skip if we can't find atom names
            }

            hbond.distance = hb.value("distance", 0.0);
            hbonds.push_back(hbond);
        }
    }

    return hbonds;
}

void print_comparison(const std::vector<InitialHBond>& modern_hbonds,
                      const std::vector<InitialHBond>& legacy_hbonds, int residue1_idx,
                      int residue2_idx) {
    std::cout << "\n========================================\n";
    std::cout << "H-bond Detection Comparison (After Validation)\n";
    std::cout << "========================================\n";
    std::cout << "Pair: (" << residue1_idx << ", " << residue2_idx << ")\n\n";
    std::cout << "Comparing at SAME step: After conflict resolution and validation\n";
    std::cout << "(Legacy records ALL H-bonds to JSON after validate_hbonds)\n\n";

    std::cout << "Modern H-bonds (after_validation): " << modern_hbonds.size() << "\n";
    std::cout << "Legacy H-bonds (from JSON): " << legacy_hbonds.size() << "\n";
    std::cout << "Difference: "
              << (static_cast<int>(modern_hbonds.size()) - static_cast<int>(legacy_hbonds.size()))
              << "\n\n";

    // Create sets for comparison
    std::set<InitialHBond> modern_set(modern_hbonds.begin(), modern_hbonds.end());
    std::set<InitialHBond> legacy_set(legacy_hbonds.begin(), legacy_hbonds.end());

    // Find common, missing, and extra
    std::vector<InitialHBond> common, missing, extra;

    for (const auto& hb : modern_set) {
        if (legacy_set.find(hb) != legacy_set.end()) {
            common.push_back(hb);
        } else {
            extra.push_back(hb);
        }
    }

    for (const auto& hb : legacy_set) {
        if (modern_set.find(hb) == modern_set.end()) {
            missing.push_back(hb);
        }
    }

    std::cout << "Common H-bonds: " << common.size() << "\n";
    std::cout << "Missing in modern: " << missing.size() << "\n";
    std::cout << "Extra in modern: " << extra.size() << "\n\n";

    if (!missing.empty()) {
        std::cout << "Missing in modern (found in legacy but not modern):\n";
        for (const auto& hb : missing) {
            std::cout << "  - " << hb.atom1_name << " -> " << hb.atom2_name
                      << " (dist=" << std::fixed << std::setprecision(3) << hb.distance << ")\n";
        }
        std::cout << "\n";
    }

    if (!extra.empty()) {
        std::cout << "Extra in modern (found in modern but not legacy):\n";
        for (const auto& hb : extra) {
            std::cout << "  + " << hb.atom1_name << " -> " << hb.atom2_name
                      << " (dist=" << std::fixed << std::setprecision(3) << hb.distance << ")\n";
        }
        std::cout << "\n";
    }

    if (common.size() == modern_hbonds.size() && common.size() == legacy_hbonds.size()) {
        std::cout << "✓ H-bond detection matches perfectly!\n";
    } else {
        std::cout << "⚠️  H-bond detection differs\n";
        std::cout
            << "\nNote: This compares H-bonds AFTER validation (what legacy records to JSON).\n";
        std::cout << "If there are differences, check:\n";
        std::cout << "  1. Atom selection (seidx range vs all atoms)\n";
        std::cout << "  2. Distance calculations\n";
        std::cout << "  3. Conflict resolution logic\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <pdb_file> <residue1_idx> <residue2_idx> [legacy_hbond_json]\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/3G8T.pdb 946 947\n";
        std::cerr << "Example: " << argv[0]
                  << " data/pdb/3G8T.pdb 946 947 data/json_legacy/3G8T.json\n";
        return 1;
    }

    std::filesystem::path pdb_file = argv[1];
    int residue1_idx = std::stoi(argv[2]);
    int residue2_idx = std::stoi(argv[3]);

    // Load PDB - include all residues to match legacy (HETATMs, waters, etc.)
    PdbParser parser;
    parser.set_include_hetatm(true); // Include HETATM records (legacy includes all)
    parser.set_include_waters(true); // Include water molecules (legacy includes all)
    Structure structure = parser.parse_file(pdb_file);

    // Get residues using legacy residue index (PDB file order)
    // Use Structure's built-in method to get residues in legacy order
    // Note: get_residue_by_legacy_idx returns const, but we need non-const for frame calculation
    // So we'll find the residue and then get a non-const pointer
    const Residue* res1_const = structure.get_residue_by_legacy_idx(residue1_idx);
    const Residue* res2_const = structure.get_residue_by_legacy_idx(residue2_idx);

    // Get non-const pointers (safe since we're not modifying the structure)
    Residue* res1 = const_cast<Residue*>(res1_const);
    Residue* res2 = const_cast<Residue*>(res2_const);

    if (!res1 || !res2) {
        std::cerr << "Error: Could not find residues " << residue1_idx << " and/or " << residue2_idx
                  << "\n";
        return 1;
    }

    // Calculate frames
    BaseFrameCalculator calculator;
    calculator.calculate_frame(*res1);
    calculator.calculate_frame(*res2);

    if (!res1->reference_frame().has_value() || !res2->reference_frame().has_value()) {
        std::cerr << "Error: Could not calculate frames for residues\n";
        return 1;
    }

    // Get modern initial H-bonds
    ValidationParameters params = ValidationParameters::defaults();
    std::vector<InitialHBond> modern_hbonds =
        get_modern_initial_hbonds(*res1, *res2, params.hb_lower, params.hb_dist1);

    // Get legacy initial H-bonds if JSON provided
    std::vector<InitialHBond> legacy_hbonds;
    if (argc >= 5) {
        std::filesystem::path legacy_json = argv[4];
        if (std::filesystem::exists(legacy_json)) {
            // Legacy JSON has multiple JSON objects - read entire file and search for pair
            std::ifstream f(legacy_json);
            std::string content((std::istreambuf_iterator<char>(f)),
                                std::istreambuf_iterator<char>());

            // Search for the specific pair using string search
            // Legacy JSON format: "base_i": 946 (with space after colon)
            std::string search1 = "\"base_i\": " + std::to_string(residue1_idx);
            std::string search2 = "\"base_j\": " + std::to_string(residue2_idx);
            std::string search3 = "\"base_i\": " + std::to_string(residue2_idx);
            std::string search4 = "\"base_j\": " + std::to_string(residue1_idx);

            // Search for the pair directly
            size_t pair_pos = content.find(search1);
            if (pair_pos != std::string::npos) {
                // Check if base_j matches nearby
                size_t check_pos2 = content.find(search2, pair_pos);
                if (check_pos2 != std::string::npos && check_pos2 < pair_pos + 200) {
                    // Found the pair - find the JSON object containing it
                    size_t obj_start = content.rfind('{', pair_pos);
                    if (obj_start != std::string::npos) {
                        // Find the type field to confirm it's hbond_list
                        // Legacy format: "type": "hbond_list" (with space after colon)
                        size_t type_pos = content.find("\"type\": \"hbond_list\"", obj_start);
                        if (type_pos == std::string::npos) {
                            // Try without space
                            type_pos = content.find("\"type\":\"hbond_list\"", obj_start);
                        }
                        if (type_pos != std::string::npos && type_pos < pair_pos + 500) {
                            // Found our pair - extract the JSON object
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
                                int found_base_i = obj.value("base_i", -1);
                                int found_base_j = obj.value("base_j", -1);
                                if (found_base_i == residue1_idx && found_base_j == residue2_idx) {
                                    legacy_hbonds = extract_legacy_initial_hbonds(obj);
                                }
                            } catch (const json::parse_error& e) {
                                // Continue searching
                            }
                        }
                    }
                }
            }

            // Also try reverse order
            pair_pos = content.find(search3);
            if (pair_pos != std::string::npos) {
                size_t check_pos2 = content.find(search4, pair_pos);
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
                                    legacy_hbonds = extract_legacy_initial_hbonds(obj);
                                }
                            } catch (const json::parse_error&) {
                                // Continue
                            }
                        }
                    }
                }
            }
        }
    }

    // Print comparison
    print_comparison(modern_hbonds, legacy_hbonds, residue1_idx, residue2_idx);

    return 0;
}

/**
 * @file test_single_pdb.cpp
 * @brief Test tool to compare generated JSON with legacy JSON for a single PDB file
 *
 * Usage: ./test_single_pdb <pdb_name>
 * Example: ./test_single_pdb 2GQ4
 */

#include <iostream>
#include <fstream>
#include <filesystem>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <nlohmann/json.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/io/serializers.hpp>
#include <x3dna/core/structure.hpp>

namespace fs = std::filesystem;

void print_atom_diff(const nlohmann::json& gen_atom, const nlohmann::json& leg_atom, const std::string& key_str) {
    std::cout << "  Key: " << key_str << "\n";

    // Helper to safely get string
    auto get_str = [](const nlohmann::json& j, const std::string& key) {
        if (j.is_null() || j.empty())
            return std::string("");
        if (!j.contains(key))
            return std::string("");
        if (j[key].is_string())
            return j[key].get<std::string>();
        if (j[key].is_number())
            return std::to_string(j[key].get<int>());
        return std::string("");
    };

    if (!gen_atom.is_null() && !gen_atom.empty()) {
        std::cout << "    Generated: " << gen_atom.dump() << "\n";
    } else {
        std::cout << "    Generated: (missing)\n";
    }
    if (!leg_atom.is_null() && !leg_atom.empty()) {
        std::cout << "    Legacy:     " << leg_atom.dump() << "\n";
    } else {
        std::cout << "    Legacy:     (missing)\n";
    }

    // Check each field
    std::vector<std::string> fields = {"atom_name", "residue_name", "chain_id", "record_type", "insertion", "alt_loc"};
    for (const auto& field : fields) {
        std::string gen_val = get_str(gen_atom, field);
        std::string leg_val = get_str(leg_atom, field);
        if (gen_val != leg_val) {
            std::cout << "    ✗ " << field << ": gen='" << gen_val << "' leg='" << leg_val << "'\n";
        }
    }

    // Check residue_seq separately
    int gen_seq = (gen_atom.is_null() || gen_atom.empty())
                      ? 0
                      : (gen_atom.contains("residue_seq") && gen_atom["residue_seq"].is_number()
                             ? gen_atom["residue_seq"].get<int>()
                             : 0);
    int leg_seq = (leg_atom.is_null() || leg_atom.empty())
                      ? 0
                      : (leg_atom.contains("residue_seq") && leg_atom["residue_seq"].is_number()
                             ? leg_atom["residue_seq"].get<int>()
                             : 0);
    if (gen_seq != leg_seq) {
        std::cout << "    ✗ residue_seq: gen=" << gen_seq << " leg=" << leg_seq << "\n";
    }

    // Check coordinates
    std::vector<double> gen_xyz, leg_xyz;
    if (!gen_atom.is_null() && !gen_atom.empty() && gen_atom.contains("xyz") && gen_atom["xyz"].is_array()) {
        gen_xyz = gen_atom["xyz"].get<std::vector<double>>();
    }
    if (!leg_atom.is_null() && !leg_atom.empty() && leg_atom.contains("xyz") && leg_atom["xyz"].is_array()) {
        leg_xyz = leg_atom["xyz"].get<std::vector<double>>();
    }

    if (gen_xyz.size() == 3 && leg_xyz.size() == 3) {
        bool coord_match = true;
        for (size_t i = 0; i < 3; ++i) {
            if (std::abs(gen_xyz[i] - leg_xyz[i]) > 0.0001) {
                coord_match = false;
                std::cout << "    ✗ xyz[" << i << "]: gen=" << gen_xyz[i] << " leg=" << leg_xyz[i]
                          << " (diff=" << std::abs(gen_xyz[i] - leg_xyz[i]) << ")\n";
            }
        }
        if (coord_match) {
            std::cout << "    ✓ Coordinates match\n";
        }
    }
    std::cout << "\n";
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <pdb_name>\n";
        std::cerr << "Example: " << argv[0] << " 2GQ4\n";
        return 1;
    }

    std::string pdb_name = argv[1];

    // File paths
    fs::path pdb_file = fs::path("data/pdb") / (pdb_name + ".pdb");
    fs::path gen_json_file = fs::path("data/json") / (pdb_name + ".json");
    fs::path leg_json_file = fs::path("data/json_legacy") / (pdb_name + ".json");

    // Check if files exist
    if (!fs::exists(pdb_file)) {
        std::cerr << "Error: PDB file not found: " << pdb_file << "\n";
        return 1;
    }

    if (!fs::exists(leg_json_file)) {
        std::cerr << "Error: Legacy JSON file not found: " << leg_json_file << "\n";
        return 1;
    }

    std::cout << "=" << std::string(70, '=') << "\n";
    std::cout << "Testing PDB: " << pdb_name << "\n";
    std::cout << "=" << std::string(70, '=') << "\n\n";

    // Parse PDB file
    std::cout << "1. Parsing PDB file: " << pdb_file << "\n";
    x3dna::io::PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(true);

    x3dna::core::Structure structure;
    try {
        structure = parser.parse_file(pdb_file);
        std::cout << "   ✓ Parsed successfully: " << structure.num_atoms() << " atoms\n\n";
    } catch (const std::exception& e) {
        std::cerr << "   ✗ Error parsing PDB: " << e.what() << "\n";
        return 1;
    }

    // Generate JSON
    std::cout << "2. Generating JSON...\n";
    nlohmann::json gen_json;
    gen_json["pdb_file"] = pdb_file.string();
    gen_json["pdb_name"] = pdb_name;
    gen_json["calculations"] = nlohmann::json::array();

    nlohmann::json metadata;
    metadata["version"] = "X3DNA Modernized C++ Library";
    gen_json["metadata"] = metadata;

    nlohmann::json structure_json = x3dna::io::StructureSerializer::to_legacy_json(structure);
    nlohmann::json pdb_atoms_record;
    pdb_atoms_record["type"] = "pdb_atoms";
    pdb_atoms_record["num_atoms"] = structure_json["num_atoms"];
    pdb_atoms_record["atoms"] = structure_json["atoms"];
    gen_json["calculations"].push_back(pdb_atoms_record);

    // Write generated JSON
    std::ofstream gen_out(gen_json_file);
    if (!gen_out.is_open()) {
        std::cerr << "Error: Cannot write to " << gen_json_file << "\n";
        return 1;
    }
    gen_out << gen_json.dump(2);
    gen_out.close();
    std::cout << "   ✓ Generated JSON: " << gen_json_file << "\n\n";

    // Load legacy JSON
    std::cout << "3. Loading legacy JSON: " << leg_json_file << "\n";
    std::ifstream leg_in(leg_json_file);
    if (!leg_in.is_open()) {
        std::cerr << "   ✗ Error: Cannot open " << leg_json_file << "\n";
        return 1;
    }
    nlohmann::json leg_json;
    leg_in >> leg_json;
    leg_in.close();
    std::cout << "   ✓ Loaded successfully\n\n";

    // Compare
    std::cout << "4. Comparing pdb_atoms records...\n";

    auto gen_atoms = std::find_if(gen_json["calculations"].begin(), gen_json["calculations"].end(),
                                  [](const nlohmann::json& j) {
                                      return j.value("type", "") == "pdb_atoms";
                                  });
    auto leg_atoms = std::find_if(leg_json["calculations"].begin(), leg_json["calculations"].end(),
                                  [](const nlohmann::json& j) {
                                      return j.value("type", "") == "pdb_atoms";
                                  });

    if (gen_atoms == gen_json["calculations"].end() || leg_atoms == leg_json["calculations"].end()) {
        std::cerr << "   ✗ Error: pdb_atoms record not found\n";
        return 1;
    }

    auto gen_atom_list = (*gen_atoms)["atoms"];
    auto leg_atom_list = (*leg_atoms)["atoms"];

    std::cout << "   Generated atoms: " << gen_atom_list.size() << "\n";
    std::cout << "   Legacy atoms:     " << leg_atom_list.size() << "\n";
    std::cout << "   Difference:       " << (int(gen_atom_list.size()) - int(leg_atom_list.size())) << "\n\n";

    // Helper function to safely get string value
    auto get_string = [](const nlohmann::json& j, const std::string& key, const std::string& default_val = "") {
        if (!j.contains(key))
            return default_val;
        if (j[key].is_string())
            return j[key].get<std::string>();
        if (j[key].is_number())
            return std::to_string(j[key].get<int>());
        return default_val;
    };

    auto get_int = [](const nlohmann::json& j, const std::string& key, int default_val = 0) {
        if (!j.contains(key))
            return default_val;
        if (j[key].is_number())
            return j[key].get<int>();
        if (j[key].is_string()) {
            try {
                return std::stoi(j[key].get<std::string>());
            } catch (...) {
                return default_val;
            }
        }
        return default_val;
    };

    // Build maps
    std::map<std::tuple<std::string, int, std::string, std::string>, nlohmann::json> gen_map;
    for (const auto& atom : gen_atom_list) {
        std::string chain_id = get_string(atom, "chain_id");
        int residue_seq = get_int(atom, "residue_seq");
        std::string insertion = get_string(atom, "insertion");
        std::string atom_name = get_string(atom, "atom_name");
        auto key = std::make_tuple(chain_id, residue_seq, insertion, atom_name);
        gen_map[key] = atom;
    }

    std::map<std::tuple<std::string, int, std::string, std::string>, nlohmann::json> leg_map;
    for (const auto& atom : leg_atom_list) {
        std::string chain_id = get_string(atom, "chain_id");
        int residue_seq = get_int(atom, "residue_seq");
        std::string insertion = get_string(atom, "insertion");
        std::string atom_name = get_string(atom, "atom_name");
        auto key = std::make_tuple(chain_id, residue_seq, insertion, atom_name);
        leg_map[key] = atom;
    }

    std::set<std::tuple<std::string, int, std::string, std::string>> gen_keys;
    std::set<std::tuple<std::string, int, std::string, std::string>> leg_keys;
    for (const auto& pair : gen_map)
        gen_keys.insert(pair.first);
    for (const auto& pair : leg_map)
        leg_keys.insert(pair.first);

    auto missing = leg_keys;
    for (const auto& k : gen_keys)
        missing.erase(k);

    auto extra = gen_keys;
    for (const auto& k : leg_keys)
        extra.erase(k);

    auto common = gen_keys;
    for (const auto& k : leg_keys) {
        if (gen_keys.find(k) == gen_keys.end()) {
            common.erase(k);
        }
    }

    std::cout << "   Missing atoms (in legacy but not generated): " << missing.size() << "\n";
    std::cout << "   Extra atoms (in generated but not legacy):   " << extra.size() << "\n";
    std::cout << "   Common atoms:                                 " << common.size() << "\n\n";

    // Check for mismatches in common atoms
    size_t mismatches = 0;
    for (const auto& key : common) {
        const auto& gen_atom = gen_map[key];
        const auto& leg_atom = leg_map[key];

        bool match = true;
        for (const auto& field : {"atom_name", "residue_name", "chain_id", "record_type", "insertion", "alt_loc"}) {
            std::string gen_val = get_string(gen_atom, field);
            std::string leg_val = get_string(leg_atom, field);
            if (gen_val != leg_val) {
                match = false;
                break;
            }
        }
        // Check residue_seq separately (it's a number)
        if (match && get_int(gen_atom, "residue_seq") != get_int(leg_atom, "residue_seq")) {
            match = false;
        }

        auto gen_xyz = gen_atom.value("xyz", std::vector<double>());
        auto leg_xyz = leg_atom.value("xyz", std::vector<double>());
        if (gen_xyz.size() == 3 && leg_xyz.size() == 3) {
            for (size_t i = 0; i < 3; ++i) {
                if (std::abs(gen_xyz[i] - leg_xyz[i]) > 0.0001) {
                    match = false;
                    break;
                }
            }
        }

        if (!match) {
            mismatches++;
        }
    }

    std::cout << "   Field/coordinate mismatches in common atoms: " << mismatches << "\n\n";

    // Show details
    if (!missing.empty() || !extra.empty() || mismatches > 0) {
        std::cout << "5. Detailed differences:\n";
        std::cout << "=" << std::string(70, '=') << "\n\n";

        if (!missing.empty()) {
            std::cout << "Missing atoms (first 10):\n";
            size_t count = 0;
            for (const auto& key : missing) {
                if (count++ >= 10)
                    break;
                std::string key_str = std::get<0>(key) + ":" + std::to_string(std::get<1>(key)) + ":" +
                                      std::get<2>(key) + ":" + std::get<3>(key);
                print_atom_diff(nlohmann::json::object(), leg_map[key], key_str);
            }
            if (missing.size() > 10) {
                std::cout << "  ... and " << (missing.size() - 10) << " more\n\n";
            }
        }

        if (!extra.empty()) {
            std::cout << "Extra atoms (first 10):\n";
            size_t count = 0;
            for (const auto& key : extra) {
                if (count++ >= 10)
                    break;
                std::string key_str = std::get<0>(key) + ":" + std::to_string(std::get<1>(key)) + ":" +
                                      std::get<2>(key) + ":" + std::get<3>(key);
                print_atom_diff(gen_map[key], nlohmann::json::object(), key_str);
            }
            if (extra.size() > 10) {
                std::cout << "  ... and " << (extra.size() - 10) << " more\n\n";
            }
        }

        if (mismatches > 0) {
            std::cout << "Mismatched atoms (first 10):\n";
            size_t count = 0;
            for (const auto& key : common) {
                const auto& gen_atom = gen_map[key];
                const auto& leg_atom = leg_map[key];

                bool match = true;
                for (const auto& field :
                     {"atom_name", "residue_name", "chain_id", "record_type", "insertion", "alt_loc"}) {
                    std::string gen_val = get_string(gen_atom, field);
                    std::string leg_val = get_string(leg_atom, field);
                    if (gen_val != leg_val) {
                        match = false;
                        break;
                    }
                }
                // Check residue_seq separately (it's a number)
                if (match && get_int(gen_atom, "residue_seq") != get_int(leg_atom, "residue_seq")) {
                    match = false;
                }

                auto gen_xyz = gen_atom.value("xyz", std::vector<double>());
                auto leg_xyz = leg_atom.value("xyz", std::vector<double>());
                if (gen_xyz.size() == 3 && leg_xyz.size() == 3) {
                    for (size_t i = 0; i < 3; ++i) {
                        if (std::abs(gen_xyz[i] - leg_xyz[i]) > 0.0001) {
                            match = false;
                            break;
                        }
                    }
                }

                if (!match) {
                    if (count++ >= 10)
                        break;
                    std::string key_str = std::get<0>(key) + ":" + std::to_string(std::get<1>(key)) + ":" +
                                          std::get<2>(key) + ":" + std::get<3>(key);
                    print_atom_diff(gen_atom, leg_atom, key_str);
                }
            }
            if (mismatches > 10) {
                std::cout << "  ... and " << (mismatches - 10) << " more\n\n";
            }
        }
    } else {
        std::cout << "5. Result: ✓ PERFECT MATCH!\n";
    }

    std::cout << "=" << std::string(70, '=') << "\n";

    return (missing.empty() && extra.empty() && mismatches == 0) ? 0 : 1;
}

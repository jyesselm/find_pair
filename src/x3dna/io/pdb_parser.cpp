/**
 * @file pdb_parser.cpp
 * @brief Implementation of PDB file parser using GEMMI library
 */

#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/constants.hpp>
#include <x3dna/core/modified_nucleotide_registry.hpp>
#include <x3dna/core/structure/residue.hpp>  // Polymorphic types
#include <gemmi/pdb.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/gz.hpp>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <map>
#include <vector>

namespace x3dna {
namespace io {

// ParseError implementation
PdbParser::ParseError::ParseError(const std::string& message, size_t line_number)
    : std::runtime_error(line_number > 0 ? message + " (line " + std::to_string(line_number) + ")" : message),
      line_number_(line_number) {}

// Public methods
core::Structure PdbParser::parse_file(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) {
        throw ParseError("PDB file does not exist: " + path.string());
    }

    try {
        // Use GEMMI to read PDB file (handles .pdb and .pdb.gz)
        gemmi::Structure gemmi_struct = gemmi::read_structure(gemmi::MaybeGzipped(path.string()));

        std::string pdb_id = path.stem().string();
        if (!gemmi_struct.name.empty()) {
            pdb_id = gemmi_struct.name;
        }

        return convert_gemmi_structure(gemmi_struct, pdb_id);

    } catch (const ParseError&) {
        throw;
    } catch (const std::exception& e) {
        throw ParseError("Error parsing PDB file " + path.string() + ": " + e.what());
    }
}

core::Structure PdbParser::parse_stream(std::istream& stream) {
    if (!stream.good()) {
        throw ParseError("Input stream is not valid");
    }

    std::stringstream buffer;
    buffer << stream.rdbuf();
    return parse_string(buffer.str());
}

core::Structure PdbParser::parse_string(const std::string& content) {
    try {
        if (content.empty()) {
            throw ParseError("Empty PDB content");
        }

        // Use GEMMI to parse PDB string
        gemmi::Structure gemmi_struct = gemmi::read_pdb_string(content, "input");

        std::string pdb_id = gemmi_struct.name;
        if (pdb_id.empty()) {
            pdb_id = "unknown";
        }

        return convert_gemmi_structure(gemmi_struct, pdb_id);

    } catch (const ParseError&) {
        throw;
    } catch (const std::exception& e) {
        throw ParseError("Error parsing PDB content: " + std::string(e.what()));
    }
}

// Convert GEMMI Structure to our Structure
core::Structure PdbParser::convert_gemmi_structure(const gemmi::Structure& gemmi_struct, const std::string& pdb_id) {
    std::map<ResidueKey, std::vector<core::Atom>> residue_atoms;
    std::vector<std::string> chain_order; // Track chain order as first encountered

    // Legacy indices: assign sequentially as atoms are encountered (1-based)
    int legacy_atom_idx = 1;
    std::map<ResidueKey, int> legacy_residue_idx_map;
    int legacy_residue_idx = 1;

    // Process only the first model (consistent with legacy behavior)
    if (gemmi_struct.models.empty()) {
        return core::Structure(pdb_id);
    }
    const gemmi::Model& model = gemmi_struct.models[0];
    int model_number = 1;

    for (const gemmi::Chain& gemmi_chain : model.chains) {
        // Get chain ID (use full string for CIF compatibility)
        std::string chain_id = gemmi_chain.name;

        // Track chain order on first encounter
        if (std::find(chain_order.begin(), chain_order.end(), chain_id) == chain_order.end()) {
            chain_order.push_back(chain_id);
        }

        for (const gemmi::Residue& gemmi_residue : gemmi_chain.residues) {
            // Get residue properties
            std::string original_residue_name = gemmi_residue.name;
            std::string residue_name = normalize_residue_name_from_gemmi(original_residue_name);

            // Get sequence number
            int residue_seq = gemmi_residue.seqid.num.value;

            // Get insertion code (convert char to string)
            std::string insertion;
            if (gemmi_residue.seqid.icode != ' ' && gemmi_residue.seqid.icode != '\0') {
                insertion = std::string(1, gemmi_residue.seqid.icode);
            }

            // Determine if this is a HETATM residue
            bool is_hetatm = (gemmi_residue.het_flag == 'H');

            // Check if we should process this residue
            if (is_hetatm) {
                bool is_modified_nuc = is_modified_nucleotide_name(residue_name);
                if (!include_hetatm_ && !is_modified_nuc) {
                    continue;
                }
                if (!include_waters_ && is_water(residue_name)) {
                    continue;
                }
            }

            for (const gemmi::Atom& gemmi_atom : gemmi_residue.atoms) {
                // Get alternate location
                char alt_loc = gemmi_atom.altloc;
                if (alt_loc == '\0') {
                    alt_loc = ' ';
                }

                // Check alt_loc filter
                if (!check_alt_loc_filter(alt_loc)) {
                    continue;
                }

                // Normalize atom name to PDB 4-character format
                std::string original_atom_name = gemmi_atom.name;
                std::string atom_name = normalize_atom_name_from_gemmi(original_atom_name);

                // Create atom using Builder pattern
                auto builder = core::Atom::create(atom_name, geometry::Vector3D(gemmi_atom.pos.x, gemmi_atom.pos.y,
                                                                                gemmi_atom.pos.z))
                                   .alt_loc(alt_loc)
                                   .occupancy(gemmi_atom.occ)
                                   .b_factor(gemmi_atom.b_iso)
                                   .atom_serial(gemmi_atom.serial)
                                   .model_number(model_number);

                // Set element if available
                if (gemmi_atom.element != gemmi::El::X) {
                    builder.element(gemmi_atom.element.name());
                }

                core::Atom atom = builder.build();

                // Assign legacy atom index
                atom.set_legacy_atom_idx(legacy_atom_idx++);

                // Track legacy residue index (will be set on Residue later)
                // Convert het_flag to record_type: 'H' for HETATM, 'A' for ATOM
                char record_type = is_hetatm ? 'H' : 'A';
                ResidueKey key{residue_name, chain_id, residue_seq, insertion, record_type};
                if (legacy_residue_idx_map.find(key) == legacy_residue_idx_map.end()) {
                    legacy_residue_idx_map[key] = legacy_residue_idx++;
                }

                // Add to residue group
                residue_atoms[key].push_back(atom);
            }
        }
    }

    return build_structure_from_residues(pdb_id, residue_atoms, legacy_residue_idx_map, chain_order);
}

// Normalize atom name from GEMMI format to legacy PDB 4-character format
std::string PdbParser::normalize_atom_name_from_gemmi(const std::string& name) const {
    if (name.empty()) {
        return "    ";
    }

    // GEMMI atom names are trimmed, convert to PDB 4-char format
    std::string normalized;
    size_t len = name.length();

    if (len == 1) {
        normalized = " " + name + "  ";
    } else if (len == 2) {
        normalized = " " + name + " ";
    } else if (len == 3) {
        normalized = " " + name;
    } else {
        normalized = name.substr(0, 4);
    }

    // Handle character replacement: '*' becomes '\''
    for (char& c : normalized) {
        if (c == '*') {
            c = '\'';
        }
    }

    // Apply exact matches for phosphate atoms, etc.
    normalized = apply_atom_name_exact_matches(normalized);

    return ensure_atom_name_length(normalized);
}

// Normalize residue name from GEMMI
std::string PdbParser::normalize_residue_name_from_gemmi(const std::string& name) const {
    if (name.empty()) {
        return "";
    }

    // Trim whitespace
    size_t start = name.find_first_not_of(" \t");
    if (start == std::string::npos) {
        return "";
    }
    size_t end = name.find_last_not_of(" \t");
    return name.substr(start, end - start + 1);
}

// Legacy helper methods (still needed for normalization)
bool PdbParser::is_water(const std::string& residue_name) const {
    if (residue_name.length() != 3) {
        return false;
    }

    char c0 = residue_name[0] & ~0x20;
    char c1 = residue_name[1] & ~0x20;
    char c2 = residue_name[2] & ~0x20;

    if (c0 == 'H' && c1 == 'O' && c2 == 'H')
        return true;
    if (c0 == 'W' && c1 == 'A' && c2 == 'T')
        return true;
    if (c0 == 'H' && c1 == '2' && c2 == 'O')
        return true;
    if (c0 == 'O' && c1 == 'H' && c2 == '2')
        return true;
    if (c0 == 'S' && c1 == 'O' && c2 == 'L')
        return true;

    return false;
}

bool PdbParser::is_modified_nucleotide_name(const std::string& residue_name) const {
    // Use centralized registry instead of hardcoded list
    return core::ModifiedNucleotideRegistry::contains(residue_name);
}

bool PdbParser::check_alt_loc_filter(char alt_loc) const {
    return (alt_loc == ' ' || alt_loc == 'A' || alt_loc == '1');
}

std::string PdbParser::apply_atom_name_exact_matches(const std::string& name) const {
    if (name == " O1'")
        return " O4'";
    if (name == " OL ")
        return " O1P";
    if (name == " OP1")
        return " O1P";
    if (name == " OR ")
        return " O2P";
    if (name == " OP2")
        return " O2P";
    if (name == " OP3")
        return " O3P";
    if (name == " C5A")
        return " C5M";
    if (name == " O5T")
        return " O5'";
    if (name == " O3T")
        return " O3'";
    if (name == "   P" || name == "P   ")
        return " P  ";

    std::string trimmed = name;
    trimmed.erase(0, trimmed.find_first_not_of(" \t"));
    trimmed.erase(trimmed.find_last_not_of(" \t") + 1);

    if (trimmed == "OP1")
        return " O1P";
    if (trimmed == "OP2")
        return " O2P";
    if (trimmed == "OP3")
        return " O3P";
    if (trimmed == "P")
        return " P  ";

    return name;
}

std::string PdbParser::ensure_atom_name_length(const std::string& name) const {
    if (name.length() == 4) {
        return name;
    }
    if (name.length() < 4) {
        return name + std::string(4 - name.length(), ' ');
    }
    return name.substr(0, 4);
}

std::string PdbParser::normalize_atom_name(const std::string& name) const {
    return normalize_atom_name_from_gemmi(name);
}

std::string PdbParser::normalize_residue_name(const std::string& name) const {
    return normalize_residue_name_from_gemmi(name);
}

core::Structure PdbParser::build_structure_from_residues(
    const std::string& pdb_id, const std::map<ResidueKey, std::vector<core::Atom>>& residue_atoms,
    const std::map<ResidueKey, int>& legacy_idx_map, const std::vector<std::string>& chain_order) const {

    core::Structure structure(pdb_id);
    std::map<std::string, core::Chain> chains;

    for (const auto& [key, atoms] : residue_atoms) {
        if (atoms.empty()) {
            continue;
        }

        core::Residue residue = core::Residue::create_from_atoms(key.residue_name, key.residue_seq, key.chain_id,
                                                                 key.insertion_code, atoms);

        // Set legacy_residue_idx from passed map (preserves PDB file order)
        residue.set_legacy_residue_idx(legacy_idx_map.at(key));

        // Store record_type in Structure map instead of on Residue
        structure.set_residue_record_type(key.chain_id, key.residue_seq, key.insertion_code, key.record_type);

        auto [it, inserted] = chains.try_emplace(key.chain_id, key.chain_id);
        it->second.add_residue(residue);
    }

    // Add chains in file encounter order (passed from parsing)
    for (const auto& chain_id : chain_order) {
        auto it = chains.find(chain_id);
        if (it != chains.end()) {
            structure.add_chain(it->second);
        }
    }

    return structure;
}

// === Polymorphic Structure parsing ===

core::structure::Structure PdbParser::parse_file_poly(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) {
        throw ParseError("PDB file does not exist: " + path.string());
    }

    try {
        gemmi::Structure gemmi_struct = gemmi::read_structure(gemmi::MaybeGzipped(path.string()));

        std::string pdb_id = path.stem().string();
        if (!gemmi_struct.name.empty()) {
            pdb_id = gemmi_struct.name;
        }

        // Reuse the atom collection logic from convert_gemmi_structure
        std::map<ResidueKey, std::vector<core::Atom>> residue_atoms;
        std::vector<std::string> chain_order;
        int legacy_atom_idx = 1;
        std::map<ResidueKey, int> legacy_residue_idx_map;
        int legacy_residue_idx = 1;

        if (gemmi_struct.models.empty()) {
            return core::structure::Structure(pdb_id);
        }

        const gemmi::Model& model = gemmi_struct.models[0];
        int model_number = 1;

        for (const gemmi::Chain& gemmi_chain : model.chains) {
            std::string chain_id = gemmi_chain.name;

            if (std::find(chain_order.begin(), chain_order.end(), chain_id) == chain_order.end()) {
                chain_order.push_back(chain_id);
            }

            for (const gemmi::Residue& gemmi_residue : gemmi_chain.residues) {
                std::string original_residue_name = gemmi_residue.name;
                std::string residue_name = normalize_residue_name_from_gemmi(original_residue_name);
                int residue_seq = gemmi_residue.seqid.num.value;

                std::string insertion;
                if (gemmi_residue.seqid.icode != ' ' && gemmi_residue.seqid.icode != '\0') {
                    insertion = std::string(1, gemmi_residue.seqid.icode);
                }

                bool is_hetatm = (gemmi_residue.het_flag == 'H');

                if (is_hetatm) {
                    bool is_modified_nuc = is_modified_nucleotide_name(residue_name);
                    if (!include_hetatm_ && !is_modified_nuc) {
                        continue;
                    }
                    if (!include_waters_ && is_water(residue_name)) {
                        continue;
                    }
                }

                for (const gemmi::Atom& gemmi_atom : gemmi_residue.atoms) {
                    char alt_loc = gemmi_atom.altloc;
                    if (alt_loc == '\0') {
                        alt_loc = ' ';
                    }

                    if (!check_alt_loc_filter(alt_loc)) {
                        continue;
                    }

                    std::string original_atom_name = gemmi_atom.name;
                    std::string atom_name = normalize_atom_name_from_gemmi(original_atom_name);

                    auto builder = core::Atom::create(atom_name, geometry::Vector3D(gemmi_atom.pos.x, gemmi_atom.pos.y,
                                                                                    gemmi_atom.pos.z))
                                       .alt_loc(alt_loc)
                                       .occupancy(gemmi_atom.occ)
                                       .b_factor(gemmi_atom.b_iso)
                                       .atom_serial(gemmi_atom.serial)
                                       .model_number(model_number);

                    if (gemmi_atom.element != gemmi::El::X) {
                        builder.element(gemmi_atom.element.name());
                    }

                    core::Atom atom = builder.build();
                    atom.set_legacy_atom_idx(legacy_atom_idx++);

                    char record_type = is_hetatm ? 'H' : 'A';
                    ResidueKey key{residue_name, chain_id, residue_seq, insertion, record_type};
                    if (legacy_residue_idx_map.find(key) == legacy_residue_idx_map.end()) {
                        legacy_residue_idx_map[key] = legacy_residue_idx++;
                    }

                    residue_atoms[key].push_back(atom);
                }
            }
        }

        return build_poly_structure_from_residues(pdb_id, residue_atoms, legacy_residue_idx_map, chain_order);

    } catch (const ParseError&) {
        throw;
    } catch (const std::exception& e) {
        throw ParseError("Error parsing PDB file " + path.string() + ": " + e.what());
    }
}

core::structure::Structure PdbParser::parse_stream_poly(std::istream& stream) {
    if (!stream.good()) {
        throw ParseError("Input stream is not valid");
    }

    std::stringstream buffer;
    buffer << stream.rdbuf();
    return parse_string_poly(buffer.str());
}

core::structure::Structure PdbParser::parse_string_poly(const std::string& content) {
    try {
        if (content.empty()) {
            throw ParseError("Empty PDB content");
        }

        gemmi::Structure gemmi_struct = gemmi::read_pdb_string(content, "input");

        std::string pdb_id = gemmi_struct.name;
        if (pdb_id.empty()) {
            pdb_id = "unknown";
        }

        // Reuse the atom collection logic
        std::map<ResidueKey, std::vector<core::Atom>> residue_atoms;
        std::vector<std::string> chain_order;
        int legacy_atom_idx = 1;
        std::map<ResidueKey, int> legacy_residue_idx_map;
        int legacy_residue_idx = 1;

        if (gemmi_struct.models.empty()) {
            return core::structure::Structure(pdb_id);
        }

        const gemmi::Model& model = gemmi_struct.models[0];
        int model_number = 1;

        for (const gemmi::Chain& gemmi_chain : model.chains) {
            std::string chain_id = gemmi_chain.name;

            if (std::find(chain_order.begin(), chain_order.end(), chain_id) == chain_order.end()) {
                chain_order.push_back(chain_id);
            }

            for (const gemmi::Residue& gemmi_residue : gemmi_chain.residues) {
                std::string original_residue_name = gemmi_residue.name;
                std::string residue_name = normalize_residue_name_from_gemmi(original_residue_name);
                int residue_seq = gemmi_residue.seqid.num.value;

                std::string insertion;
                if (gemmi_residue.seqid.icode != ' ' && gemmi_residue.seqid.icode != '\0') {
                    insertion = std::string(1, gemmi_residue.seqid.icode);
                }

                bool is_hetatm = (gemmi_residue.het_flag == 'H');

                if (is_hetatm) {
                    bool is_modified_nuc = is_modified_nucleotide_name(residue_name);
                    if (!include_hetatm_ && !is_modified_nuc) {
                        continue;
                    }
                    if (!include_waters_ && is_water(residue_name)) {
                        continue;
                    }
                }

                for (const gemmi::Atom& gemmi_atom : gemmi_residue.atoms) {
                    char alt_loc = gemmi_atom.altloc;
                    if (alt_loc == '\0') {
                        alt_loc = ' ';
                    }

                    if (!check_alt_loc_filter(alt_loc)) {
                        continue;
                    }

                    std::string original_atom_name = gemmi_atom.name;
                    std::string atom_name = normalize_atom_name_from_gemmi(original_atom_name);

                    auto builder = core::Atom::create(atom_name, geometry::Vector3D(gemmi_atom.pos.x, gemmi_atom.pos.y,
                                                                                    gemmi_atom.pos.z))
                                       .alt_loc(alt_loc)
                                       .occupancy(gemmi_atom.occ)
                                       .b_factor(gemmi_atom.b_iso)
                                       .atom_serial(gemmi_atom.serial)
                                       .model_number(model_number);

                    if (gemmi_atom.element != gemmi::El::X) {
                        builder.element(gemmi_atom.element.name());
                    }

                    core::Atom atom = builder.build();
                    atom.set_legacy_atom_idx(legacy_atom_idx++);

                    char record_type = is_hetatm ? 'H' : 'A';
                    ResidueKey key{residue_name, chain_id, residue_seq, insertion, record_type};
                    if (legacy_residue_idx_map.find(key) == legacy_residue_idx_map.end()) {
                        legacy_residue_idx_map[key] = legacy_residue_idx++;
                    }

                    residue_atoms[key].push_back(atom);
                }
            }
        }

        return build_poly_structure_from_residues(pdb_id, residue_atoms, legacy_residue_idx_map, chain_order);

    } catch (const ParseError&) {
        throw;
    } catch (const std::exception& e) {
        throw ParseError("Error parsing PDB content: " + std::string(e.what()));
    }
}

core::structure::Structure PdbParser::build_poly_structure_from_residues(
    const std::string& pdb_id,
    const std::map<ResidueKey, std::vector<core::Atom>>& residue_atoms,
    const std::map<ResidueKey, int>& legacy_idx_map,
    const std::vector<std::string>& chain_order) const {

    core::structure::Structure structure(pdb_id);
    std::map<std::string, core::structure::Chain> chains;

    for (const auto& [key, atoms] : residue_atoms) {
        if (atoms.empty()) {
            continue;
        }

        // Use ResidueFactory to create the appropriate polymorphic type
        auto residue = core::structure::ResidueFactory::create(
            key.residue_name, key.residue_seq, key.chain_id, key.insertion_code, atoms);

        // Set legacy_residue_idx from passed map
        residue->set_legacy_residue_idx(legacy_idx_map.at(key));

        // Store record_type in Structure map
        structure.set_residue_record_type(key.chain_id, key.residue_seq, key.insertion_code, key.record_type);

        auto [it, inserted] = chains.try_emplace(key.chain_id, key.chain_id);
        it->second.add_residue(std::move(residue));
    }

    // Add chains in file encounter order
    for (const auto& chain_id : chain_order) {
        auto it = chains.find(chain_id);
        if (it != chains.end()) {
            structure.add_chain(std::move(it->second));
        }
    }

    return structure;
}

} // namespace io
} // namespace x3dna

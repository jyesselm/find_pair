/**
 * @file cif_parser.cpp
 * @brief Implementation of CIF/mmCIF file parser using GEMMI library
 */

#include <x3dna/io/cif_parser.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/typing/type_registry.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/gz.hpp>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <fstream>
#include <sstream>

namespace x3dna {
namespace io {

// ParseError implementation
CifParser::ParseError::ParseError(const std::string& message) : std::runtime_error(message) {}

core::Structure CifParser::parse_file(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) {
        throw ParseError("CIF file does not exist: " + path.string());
    }

    try {
        // Use GEMMI to read CIF file (handles .cif and .cif.gz)
        gemmi::Structure gemmi_struct = gemmi::read_structure(gemmi::MaybeGzipped(path.string()));

        std::string pdb_id = path.stem().string();
        if (!gemmi_struct.name.empty()) {
            pdb_id = gemmi_struct.name;
        }
        // Remove .cif extension if present (for .cif.gz files)
        if (pdb_id.size() > 4 && pdb_id.substr(pdb_id.size() - 4) == ".cif") {
            pdb_id = pdb_id.substr(0, pdb_id.size() - 4);
        }

        return convert_gemmi_structure(gemmi_struct, pdb_id);

    } catch (const ParseError&) {
        throw;
    } catch (const std::exception& e) {
        throw ParseError("Error parsing CIF file " + path.string() + ": " + e.what());
    }
}

core::Structure CifParser::parse_string(const std::string& content) {
    try {
        if (content.empty()) {
            throw ParseError("Empty CIF content");
        }

        // Use GEMMI to parse CIF string
        gemmi::cif::Document doc = gemmi::cif::read_string(content);
        if (doc.blocks.empty()) {
            throw ParseError("No data blocks found in CIF content");
        }

        gemmi::Structure gemmi_struct = gemmi::make_structure_from_block(doc.blocks[0]);

        std::string pdb_id = gemmi_struct.name;
        if (pdb_id.empty()) {
            pdb_id = "unknown";
        }

        return convert_gemmi_structure(gemmi_struct, pdb_id);

    } catch (const ParseError&) {
        throw;
    } catch (const std::exception& e) {
        throw ParseError("Error parsing CIF content: " + std::string(e.what()));
    }
}

// Convert GEMMI Structure to our Structure
core::Structure CifParser::convert_gemmi_structure(const gemmi::Structure& gemmi_struct, const std::string& pdb_id) {
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
            std::string residue_name = normalize_residue_name(original_residue_name);

            // Get sequence number (use auth_seq_id for PDB compatibility)
            int residue_seq = gemmi_residue.seqid.num.value;

            // Get insertion code (convert char to string)
            std::string insertion;
            if (gemmi_residue.seqid.icode != ' ' && gemmi_residue.seqid.icode != '\0') {
                insertion = std::string(1, gemmi_residue.seqid.icode);
            }

            // Determine if this is a HETATM residue
            bool is_hetatm = (gemmi_residue.het_flag == 'H');

            // Check if we should process this residue
            if (!should_keep_atom(is_hetatm, ' ', residue_name)) {
                // Check at residue level (alt_loc checked per-atom)
                if (is_hetatm && !is_modified_nucleotide_name(residue_name)) {
                    if (!include_hetatm_)
                        continue;
                    if (!include_waters_ && is_water(residue_name))
                        continue;
                }
            }

            for (const gemmi::Atom& gemmi_atom : gemmi_residue.atoms) {
                // Get alternate location
                char alt_loc = gemmi_atom.altloc;
                if (alt_loc == '\0') {
                    alt_loc = ' ';
                }

                // Check if we should keep this atom
                if (!should_keep_atom(is_hetatm, alt_loc, residue_name)) {
                    continue;
                }

                // Normalize atom name to PDB 4-character format
                std::string original_atom_name = gemmi_atom.name;
                std::string atom_name = normalize_atom_name(original_atom_name);

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

bool CifParser::should_keep_atom(bool is_hetatm, char alt_loc, const std::string& residue_name) const {
    // Check alt_loc filter first
    if (!check_alt_loc_filter(alt_loc)) {
        return false;
    }

    // For HETATM records
    if (is_hetatm) {
        // Auto-include modified nucleotides even without include_hetatm_ flag
        if (is_modified_nucleotide_name(residue_name)) {
            return true;
        }

        // Skip if HETATM not enabled
        if (!include_hetatm_) {
            return false;
        }

        // Skip waters if not enabled
        if (!include_waters_ && is_water(residue_name)) {
            return false;
        }
    }

    return true;
}

bool CifParser::check_alt_loc_filter(char alt_loc) const {
    // Keep atoms with alt_loc = ' ', 'A', or '1' (same as PdbParser)
    return (alt_loc == ' ' || alt_loc == '\0' || alt_loc == 'A' || alt_loc == '1');
}

bool CifParser::is_water(const std::string& residue_name) const {
    std::string name = residue_name;
    name.erase(0, name.find_first_not_of(" \t"));
    name.erase(name.find_last_not_of(" \t") + 1);

    if (name.length() != 3) {
        return false;
    }

    std::string upper;
    for (char c : name) {
        upper += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }

    return (upper == "HOH" || upper == "WAT" || upper == "H2O" || upper == "OH2" || upper == "SOL");
}

bool CifParser::is_modified_nucleotide_name(const std::string& residue_name) const {
    // Use centralized registry instead of hardcoded list
    return core::TypeRegistry::instance().is_nucleotide(residue_name);
}

std::string CifParser::normalize_atom_name(const std::string& name) const {
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

std::string CifParser::apply_atom_name_formatting_rules(const std::string& name) const {
    // Not needed with GEMMI-based parsing, but kept for API compatibility
    return name;
}

// Apply exact matches for atom name normalization during CIF parsing.
// NOTE: This function uses PDB 4-character format with position-specific padding.
// This is intentional - the padding encodes position information from the PDB file format
// (columns 13-16). The output is later trimmed when stored in the Atom class.
// Examples:
//   " OP1" (CIF format) -> " O1P" (legacy format)
//   "   P" or "P   " (malformed) -> " P  " (standard format)
std::string CifParser::apply_atom_name_exact_matches(const std::string& name) const {
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
    size_t start = trimmed.find_first_not_of(" \t");
    size_t end = trimmed.find_last_not_of(" \t");
    if (start != std::string::npos) {
        trimmed = trimmed.substr(start, end - start + 1);
    }

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

std::string CifParser::ensure_atom_name_length(const std::string& name) const {
    if (name.length() == 4) {
        return name;
    }
    if (name.length() < 4) {
        return name + std::string(4 - name.length(), ' ');
    }
    return name.substr(0, 4);
}

std::string CifParser::normalize_residue_name(const std::string& name) const {
    if (name.empty()) {
        return "";
    }
    size_t start = name.find_first_not_of(" \t");
    if (start == std::string::npos) {
        return "";
    }
    size_t end = name.find_last_not_of(" \t");
    return name.substr(start, end - start + 1);
}

core::Structure CifParser::build_structure_from_residues(
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

} // namespace io
} // namespace x3dna

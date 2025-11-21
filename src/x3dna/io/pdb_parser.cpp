/**
 * @file pdb_parser.cpp
 * @brief Implementation of PDB file parser
 */

#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/chain.hpp>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace x3dna {
namespace io {

// ParseError implementation
PdbParser::ParseError::ParseError(const std::string& message, size_t line_number)
    : std::runtime_error(line_number > 0 
                         ? message + " (line " + std::to_string(line_number) + ")"
                         : message), 
      line_number_(line_number) {
}

// Public methods
core::Structure PdbParser::parse_file(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) {
        throw ParseError("PDB file does not exist: " + path.string());
    }

    std::ifstream file(path);
    if (!file.is_open()) {
        throw ParseError("Cannot open PDB file: " + path.string());
    }

    // Extract PDB ID from filename (without extension)
    std::string pdb_id = path.stem().string();
    
    return parse_impl(file, pdb_id);
}

core::Structure PdbParser::parse_stream(std::istream& stream) {
    if (!stream.good()) {
        throw ParseError("Input stream is not valid");
    }
    return parse_impl(stream);
}

core::Structure PdbParser::parse_string(const std::string& content) {
    std::istringstream stream(content);
    return parse_impl(stream);
}

// Private helper methods
std::string PdbParser::parse_atom_name(const std::string& line) {
    // PDB format: columns 13-16 (1-indexed)
    // Atom name is exactly 4 characters (may have leading/trailing spaces)
    // We normalize to match legacy JSON format
    if (line.length() < 16) {
        return "    ";  // Return 4 spaces if line too short
    }
    
    // Extract exactly 4 characters
    std::string name = line.substr(12, 4);  // 0-indexed: [12, 16)
    
    // Normalize to match legacy format (handles OP1->O1P, OP2->O2P, etc.)
    return normalize_atom_name(name);
}

std::string PdbParser::parse_residue_name(const std::string& line) {
    // PDB format: columns 18-20 (1-indexed)
    if (line.length() < 20) {
        return "";
    }
    
    std::string name = line.substr(17, 3);  // 0-indexed: [17, 20)
    
    // Trim whitespace
    name.erase(0, name.find_first_not_of(" \t"));
    name.erase(name.find_last_not_of(" \t") + 1);
    
    return normalize_residue_name(name);
}

char PdbParser::parse_chain_id(const std::string& line) {
    // PDB format: column 22 (1-indexed)
    if (line.length() < 22) {
        return ' ';  // Default to space if not specified
    }
    
    char chain_id = line[21];  // 0-indexed: position 21
    // Normalize '0' to space to match legacy behavior
    if (chain_id == '0') {
        return ' ';
    }
    return (chain_id == ' ') ? ' ' : chain_id;
}

char PdbParser::parse_alt_loc(const std::string& line) {
    // PDB format: column 17 (1-indexed)
    if (line.length() < 17) {
        return ' ';  // Default to space if not specified
    }
    
    char alt_loc = line[16];  // 0-indexed: position 16
    return (alt_loc == ' ') ? ' ' : alt_loc;
}

char PdbParser::parse_insertion(const std::string& line) {
    // PDB format: column 27 (1-indexed)
    if (line.length() < 27) {
        return ' ';  // Default to space if not specified
    }
    
    char insertion = line[26];  // 0-indexed: position 26
    return (insertion == ' ') ? ' ' : insertion;
}

double PdbParser::parse_occupancy(const std::string& line) {
    // PDB format: columns 55-60 (1-indexed)
    if (line.length() < 60) {
        return 1.0;  // Default to 1.0 if not specified
    }
    
    std::string occ_str = line.substr(54, 6);  // 0-indexed: [54, 60)
    
    // Trim whitespace
    occ_str.erase(0, occ_str.find_first_not_of(" \t"));
    occ_str.erase(occ_str.find_last_not_of(" \t") + 1);
    
    if (occ_str.empty()) {
        return 1.0;  // Default to 1.0 if empty
    }
    
    try {
        return std::stod(occ_str);
    } catch (const std::exception& e) {
        return 1.0;  // Default to 1.0 if parsing fails
    }
}

int PdbParser::parse_residue_seq(const std::string& line) {
    // PDB format: columns 23-26 (1-indexed)
    if (line.length() < 26) {
        throw ParseError("Line too short to contain residue sequence number");
    }
    
    std::string seq_str = line.substr(22, 4);  // 0-indexed: [22, 26)
    
    // Trim whitespace
    seq_str.erase(0, seq_str.find_first_not_of(" \t"));
    seq_str.erase(seq_str.find_last_not_of(" \t") + 1);
    
    if (seq_str.empty()) {
        throw ParseError("Residue sequence number is empty");
    }
    
    try {
        return std::stoi(seq_str);
    } catch (const std::exception& e) {
        throw ParseError("Cannot parse residue sequence number: " + seq_str);
    }
}

geometry::Vector3D PdbParser::parse_coordinates(const std::string& line) {
    // PDB format: columns 31-38 (x), 39-46 (y), 47-54 (z) (1-indexed)
    if (line.length() < 54) {
        throw ParseError("Line too short to contain coordinates");
    }
    
    std::string x_str = line.substr(30, 8);  // 0-indexed: [30, 38)
    std::string y_str = line.substr(38, 8);  // 0-indexed: [38, 46)
    std::string z_str = line.substr(46, 8);  // 0-indexed: [46, 54)
    
    // Trim whitespace
    x_str.erase(0, x_str.find_first_not_of(" \t"));
    x_str.erase(x_str.find_last_not_of(" \t") + 1);
    y_str.erase(0, y_str.find_first_not_of(" \t"));
    y_str.erase(y_str.find_last_not_of(" \t") + 1);
    z_str.erase(0, z_str.find_first_not_of(" \t"));
    z_str.erase(z_str.find_last_not_of(" \t") + 1);
    
    try {
        double x = std::stod(x_str);
        double y = std::stod(y_str);
        double z = std::stod(z_str);
        return geometry::Vector3D(x, y, z);
    } catch (const std::exception& e) {
        throw ParseError("Cannot parse coordinates: x=" + x_str + ", y=" + y_str + ", z=" + z_str);
    }
}

core::Atom PdbParser::parse_atom_line(const std::string& line, size_t line_number) {
    try {
        std::string atom_name = parse_atom_name(line);
        std::string residue_name = parse_residue_name(line);
        char chain_id = parse_chain_id(line);
        int residue_seq = parse_residue_seq(line);
        char alt_loc = parse_alt_loc(line);
        char insertion = parse_insertion(line);
        double occupancy = parse_occupancy(line);
        geometry::Vector3D coords = parse_coordinates(line);
        
        core::Atom atom(atom_name, coords, residue_name, chain_id, residue_seq, 'A');
        atom.set_alt_loc(alt_loc);
        atom.set_insertion(insertion);
        atom.set_occupancy(occupancy);
        return atom;
    } catch (const std::exception& e) {
        throw ParseError("Error parsing ATOM line: " + std::string(e.what()), line_number);
    }
}

core::Atom PdbParser::parse_hetatm_line(const std::string& line, size_t line_number) {
    try {
        std::string atom_name = parse_atom_name(line);
        std::string residue_name = parse_residue_name(line);
        char chain_id = parse_chain_id(line);
        int residue_seq = parse_residue_seq(line);
        char alt_loc = parse_alt_loc(line);
        char insertion = parse_insertion(line);
        double occupancy = parse_occupancy(line);
        geometry::Vector3D coords = parse_coordinates(line);
        
        core::Atom atom(atom_name, coords, residue_name, chain_id, residue_seq, 'H');
        atom.set_alt_loc(alt_loc);
        atom.set_insertion(insertion);
        atom.set_occupancy(occupancy);
        return atom;
    } catch (const std::exception& e) {
        throw ParseError("Error parsing HETATM line: " + std::string(e.what()), line_number);
    }
}

bool PdbParser::is_water(const std::string& residue_name) const {
    // Common water residue names
    std::string upper_name = residue_name;
    std::transform(upper_name.begin(), upper_name.end(), upper_name.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    
    return (upper_name == "HOH" || upper_name == "WAT" || upper_name == "H2O" ||
            upper_name == "OH2" || upper_name == "SOL");
}

std::string PdbParser::normalize_atom_name(const std::string& name) const {
    // PDB atom names are typically 4 characters with leading space if needed
    // Examples: " C1'", " N3 ", " P  ", " O5'"
    if (name.length() == 0) {
        return "    ";  // 4 spaces
    }
    
    std::string normalized = name;
    
    // Normalize phosphate atom names to match legacy format
    // Legacy uses " O1P" and " O2P", but PDB files may use " OP1" and " OP2"
    std::string trimmed = normalized;
    trimmed.erase(0, trimmed.find_first_not_of(" \t"));
    trimmed.erase(trimmed.find_last_not_of(" \t") + 1);
    
    if (trimmed == "OP1") {
        normalized = " O1P";
    } else if (trimmed == "OP2") {
        normalized = " O2P";
    } else if (trimmed == "OP3") {
        normalized = " O3P";
    }
    
    if (normalized.length() == 4) {
        return normalized;  // Already 4 characters
    }
    
    if (normalized.length() < 4) {
        // Pad to 4 characters with leading space
        return std::string(4 - normalized.length(), ' ') + normalized;
    }
    
    // If longer than 4, truncate (shouldn't happen in standard PDB)
    return normalized.substr(0, 4);
}

std::string PdbParser::normalize_residue_name(const std::string& name) const {
    // PDB residue names are typically 3 characters with leading space if needed
    // Examples: "  A", "  C", "  G", "  T", "  U", "HOH"
    if (name.length() == 0) {
        return "   ";  // 3 spaces
    }
    
    if (name.length() == 3) {
        return name;  // Already 3 characters
    }
    
    if (name.length() < 3) {
        // Pad to 3 characters with leading space
        return std::string(3 - name.length(), ' ') + name;
    }
    
    // If longer than 3, truncate (shouldn't happen in standard PDB)
    return name.substr(0, 3);
}

core::Structure PdbParser::parse_impl(std::istream& stream, const std::string& pdb_id) {
    core::Structure structure(pdb_id);
    
    // Map to group atoms by chain and residue
    // Key: (chain_id, residue_seq), Value: vector of atoms
    std::map<std::pair<char, int>, std::vector<core::Atom>> residue_atoms;
    
    // Track atoms by position to filter alternate conformations
    // Key: (chain_id, residue_seq, insertion, atom_name), Value: (atom, occupancy)
    // We'll keep the atom with highest occupancy for each position
    std::map<std::tuple<char, int, char, std::string>, std::pair<core::Atom, double>> atom_candidates;
    
    size_t line_number = 0;
    std::string line;
    
    while (std::getline(stream, line)) {
        line_number++;
        
        // Skip empty lines
        if (line.empty() || line.find_first_not_of(" \t") == std::string::npos) {
            continue;
        }
        
        // Check for ATOM record
        if (line.length() >= 4 && line.substr(0, 4) == "ATOM") {
            try {
                core::Atom atom = parse_atom_line(line, line_number);
                
                // Use (chain_id, residue_seq, insertion, atom_name) as key to handle insertion codes
                // For alternate conformations with same key, keep the one with highest occupancy
                std::tuple<char, int, char, std::string> atom_key = 
                    std::make_tuple(atom.chain_id(), atom.residue_seq(), atom.insertion(), atom.name());
                
                double occupancy = atom.occupancy();
                char alt_loc = atom.alt_loc();
                
                // Filter: legacy code only filters by occupancy if Gvars.OCCUPANCY is TRUE
                // Default is FALSE, so by default it doesn't filter by occupancy
                // Also filter by alt_loc: default ALT_LIST is "A1" (keep 'A', '1', or ' ')
                bool keep_by_occupancy = !filter_by_occupancy_ || occupancy > 0.0;
                
                if (keep_by_occupancy) {
                    // Check alt_loc: default ALT_LIST "A1" means keep 'A', '1', or ' '
                    bool keep_alt_loc = (alt_loc == ' ' || alt_loc == 'A' || alt_loc == '1');
                    
                    if (keep_alt_loc) {
                        // For alternate conformations, keep the one with highest occupancy
                        if (atom_candidates.find(atom_key) != atom_candidates.end()) {
                            // Already have a candidate, keep the one with higher occupancy
                            if (occupancy > atom_candidates[atom_key].second) {
                                atom_candidates[atom_key] = {atom, occupancy};
                            }
                        } else {
                            // First candidate for this position
                            atom_candidates[atom_key] = {atom, occupancy};
                        }
                    }
                }
                // Skip atoms with occupancy <= 0 or alt_loc not in "A1"
            } catch (const ParseError& e) {
                // Skip malformed lines or rethrow if critical
                // For now, we'll skip and continue
                continue;
            }
        }
        // Check for HETATM record
        else if (line.length() >= 6 && line.substr(0, 6) == "HETATM") {
            if (!include_hetatm_) {
                continue;  // Skip HETATM if not included
            }
            
            try {
                core::Atom atom = parse_hetatm_line(line, line_number);
                
                // Check if water and should be skipped
                if (!include_waters_ && is_water(atom.residue_name())) {
                    continue;
                }
                
                // Use (chain_id, residue_seq, insertion, atom_name) as key to handle insertion codes
                // For alternate conformations with same key, keep the one with highest occupancy
                std::tuple<char, int, char, std::string> atom_key = 
                    std::make_tuple(atom.chain_id(), atom.residue_seq(), atom.insertion(), atom.name());
                
                double occupancy = atom.occupancy();
                char alt_loc = atom.alt_loc();
                
                // Filter: legacy code only filters by occupancy if Gvars.OCCUPANCY is TRUE
                // Default is FALSE, so by default it doesn't filter by occupancy
                // Also filter by alt_loc: default ALT_LIST is "A1" (keep 'A', '1', or ' ')
                bool keep_by_occupancy = !filter_by_occupancy_ || occupancy > 0.0;
                
                if (keep_by_occupancy) {
                    // Check alt_loc: default ALT_LIST "A1" means keep 'A', '1', or ' '
                    bool keep_alt_loc = (alt_loc == ' ' || alt_loc == 'A' || alt_loc == '1');
                    
                    if (keep_alt_loc) {
                        // For alternate conformations, keep the one with highest occupancy
                        if (atom_candidates.find(atom_key) != atom_candidates.end()) {
                            // Already have a candidate, keep the one with higher occupancy
                            if (occupancy > atom_candidates[atom_key].second) {
                                atom_candidates[atom_key] = {atom, occupancy};
                            }
                        } else {
                            // First candidate for this position
                            atom_candidates[atom_key] = {atom, occupancy};
                        }
                    }
                }
                // Skip atoms with occupancy <= 0 or alt_loc not in "A1"
            } catch (const ParseError& e) {
                // Skip malformed lines
                continue;
            }
        }
        // Check for HEADER record to extract PDB ID if not provided
        else if (line.length() >= 6 && line.substr(0, 6) == "HEADER" && pdb_id.empty()) {
            // PDB ID is in columns 63-66 (1-indexed)
            if (line.length() >= 66) {
                std::string header_id = line.substr(62, 4);
                // Trim whitespace
                header_id.erase(0, header_id.find_first_not_of(" \t"));
                header_id.erase(header_id.find_last_not_of(" \t") + 1);
                if (!header_id.empty()) {
                    // Update structure PDB ID (we'll need to add a setter or reconstruct)
                    // For now, we'll use the extracted ID
                }
            }
        }
    }
    
    // Process atom candidates: add all candidates to residue_atoms
    // (We've selected highest occupancy for alternate conformations)
    for (const auto& [atom_key, atom_occ_pair] : atom_candidates) {
        const auto& atom = atom_occ_pair.first;
        char chain_id = atom.chain_id();
        int residue_seq = atom.residue_seq();
        residue_atoms[{chain_id, residue_seq}].push_back(atom);
    }
    
    // Create chains and residues from grouped atoms
    std::map<char, core::Chain> chains;
    
    for (const auto& [key, atoms] : residue_atoms) {
        if (atoms.empty()) {
            continue;
        }
        
        char chain_id = key.first;
        int residue_seq = key.second;
        
        // Get residue name from first atom
        std::string residue_name = atoms[0].residue_name();
        
        // Create residue
        core::Residue residue(residue_name, residue_seq, chain_id);
        
        // Add all atoms to residue
        for (const auto& atom : atoms) {
            residue.add_atom(atom);
        }
        
        // Add residue to chain
        if (chains.find(chain_id) == chains.end()) {
            chains[chain_id] = core::Chain(chain_id);
        }
        chains[chain_id].add_residue(residue);
    }
    
    // Add chains to structure
    for (auto& [chain_id, chain] : chains) {
        structure.add_chain(chain);
    }
    
    return structure;
}

} // namespace io
} // namespace x3dna


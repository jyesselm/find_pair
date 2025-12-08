/**
 * @file pdb_parser.cpp
 * @brief Implementation of PDB file parser
 */

#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/residue_factory.hpp>
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

    std::ifstream file(path);
    if (!file.is_open()) {
        throw ParseError("Cannot open PDB file: " + path.string());
    }

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

// Private helper methods for parsing individual fields
std::string PdbParser::parse_atom_name(const std::string& line) {
    if (line.length() < 16) {
        return "    ";
    }
    std::string name = line.substr(12, 4);
    return normalize_atom_name(name);
}

std::string PdbParser::parse_residue_name(const std::string& line) {
    if (line.length() < 20) {
        return "";
    }
    std::string name = line.substr(17, 3);

    // Legacy normalization: shift left if last char is space (up to 2 times)
    for (int i = 0; i < 2; i++) {
        if (name.length() == 3 && name[2] == ' ') {
            name[2] = name[1];
            name[1] = name[0];
            name[0] = ' ';
        }
    }

    return normalize_residue_name(name);
}

char PdbParser::parse_chain_id(const std::string& line) {
    if (line.length() < 22) {
        return ' ';
    }
    return line[21];
}

char PdbParser::parse_alt_loc(const std::string& line) {
    if (line.length() < 17) {
        return ' ';
    }
    return line[16];
}

char PdbParser::parse_insertion(const std::string& line) {
    if (line.length() < 27) {
        return ' ';
    }
    return line[26];
}

double PdbParser::parse_occupancy(const std::string& line) {
    if (line.length() < 60) {
        return 1.0;
    }

    // Optimized: parse directly from line without substring copy
    // Occupancy is in columns 55-60 (0-indexed: 54-60)
    size_t start = 54;
    size_t end = 60;

    // Find first non-whitespace
    while (start < end && (line[start] == ' ' || line[start] == '\t')) {
        start++;
    }
    // Find last non-whitespace
    while (end > start && (line[end - 1] == ' ' || line[end - 1] == '\t')) {
        end--;
    }

    if (start >= end) {
        return 1.0;
    }

    // Parse directly from substring view
    std::string occ_str = line.substr(start, end - start);
    try {
        return std::stod(occ_str);
    } catch (const std::exception&) {
        return 1.0;
    }
}

int PdbParser::parse_residue_seq(const std::string& line) {
    if (line.length() < 26) {
        throw ParseError("Line too short to contain residue sequence number");
    }

    // Optimized: parse directly without multiple string operations
    // Residue sequence is in columns 23-26 (0-indexed: 22-26)
    size_t start = 22;
    size_t end = 26;

    // Skip leading whitespace
    while (start < end && (line[start] == ' ' || line[start] == '\t')) {
        start++;
    }
    // Skip trailing whitespace
    while (end > start && (line[end - 1] == ' ' || line[end - 1] == '\t')) {
        end--;
    }

    if (start >= end) {
        throw ParseError("Residue sequence number is empty");
    }

    std::string seq_str = line.substr(start, end - start);
    try {
        return std::stoi(seq_str);
    } catch (const std::exception&) {
        throw ParseError("Cannot parse residue sequence number: " + seq_str);
    }
}

geometry::Vector3D PdbParser::parse_coordinates(const std::string& line) {
    if (line.length() < 54) {
        throw ParseError("Line too short to contain coordinates");
    }

    // Optimized: parse coordinates directly without substring copies
    // X: columns 31-38 (0-indexed: 30-38)
    // Y: columns 39-46 (0-indexed: 38-46)
    // Z: columns 47-54 (0-indexed: 46-54)
    auto parse_coord = [&line](size_t start_col, size_t len) -> double {
        size_t start = start_col;
        size_t end = start_col + len;

        // Skip leading whitespace
        while (start < end && (line[start] == ' ' || line[start] == '\t')) {
            start++;
        }
        // Skip trailing whitespace
        while (end > start && (line[end - 1] == ' ' || line[end - 1] == '\t')) {
            end--;
        }

        if (start >= end) {
            return 0.0;
        }

        // Parse directly from line substring
        std::string coord_str = line.substr(start, end - start);
        try {
            return std::stod(coord_str);
        } catch (const std::exception&) {
            return 0.0;
        }
    };

    try {
        double x = parse_coord(30, 8);
        double y = parse_coord(38, 8);
        double z = parse_coord(46, 8);
        return geometry::Vector3D(x, y, z);
    } catch (const std::exception&) {
        throw ParseError("Cannot parse coordinates");
    }
}

// Helper function to parse atom metadata (serial, b-factor, element)
PdbParser::atom_metadata PdbParser::parse_atom_metadata(const std::string& line, size_t line_number) {
    atom_metadata metadata;
    metadata.atom_serial = 0;
    metadata.b_factor = 0.0;
    metadata.element = "";
    metadata.original_atom_name = line.length() >= 16 ? line.substr(12, 4) : "    ";
    metadata.original_residue_name = line.length() >= 20 ? line.substr(17, 3) : "   ";

    // Parse atom serial number (columns 7-11)
    if (line.length() >= 11) {
        std::string serial_str = line.substr(6, 5);
        try {
            metadata.atom_serial = std::stoi(serial_str);
        } catch (...) {
            metadata.atom_serial = static_cast<int>(line_number);
        }
    }

    // Parse B-factor (columns 61-66)
    if (line.length() >= 66) {
        std::string bfactor_str = line.substr(60, 6);
        try {
            metadata.b_factor = std::stod(bfactor_str);
        } catch (...) {
            metadata.b_factor = 0.0;
        }
    }

    // Parse element symbol (columns 77-78)
    if (line.length() >= 78) {
        metadata.element = line.substr(76, 2);
        metadata.element.erase(0, metadata.element.find_first_not_of(" \t"));
        metadata.element.erase(metadata.element.find_last_not_of(" \t") + 1);
    }

    return metadata;
}

// Build atom from parsed data
core::Atom PdbParser::build_atom_from_parsed_data(const std::string& line, size_t line_number, char atom_type) {
    std::string atom_name = parse_atom_name(line);
    std::string residue_name = parse_residue_name(line);
    char chain_id = parse_chain_id(line);
    int residue_seq = parse_residue_seq(line);
    char alt_loc = parse_alt_loc(line);
    char insertion = parse_insertion(line);
    double occupancy = parse_occupancy(line);
    geometry::Vector3D coords = parse_coordinates(line);

    atom_metadata metadata = parse_atom_metadata(line, line_number);

    core::Atom atom(atom_name, coords, residue_name, chain_id, residue_seq, atom_type);
    atom.set_alt_loc(alt_loc);
    atom.set_insertion(insertion);
    atom.set_occupancy(occupancy);
    atom.set_atom_serial(metadata.atom_serial);
    atom.set_line_number(line_number);
    atom.set_b_factor(metadata.b_factor);
    if (!metadata.element.empty()) {
        atom.set_element(metadata.element);
    }
    atom.set_original_atom_name(metadata.original_atom_name);
    atom.set_original_residue_name(metadata.original_residue_name);

    return atom;
}

core::Atom PdbParser::parse_atom_line(const std::string& line, size_t line_number) {
    try {
        return build_atom_from_parsed_data(line, line_number, 'A');
    } catch (const std::exception& e) {
        throw ParseError("Error parsing ATOM line: " + std::string(e.what()), line_number);
    }
}

core::Atom PdbParser::parse_hetatm_line(const std::string& line, size_t line_number) {
    try {
        return build_atom_from_parsed_data(line, line_number, 'H');
    } catch (const std::exception& e) {
        throw ParseError("Error parsing HETATM line: " + std::string(e.what()), line_number);
    }
}

bool PdbParser::is_water(const std::string& residue_name) const {
    // Optimized: Fast case-insensitive comparison for common water names
    // Avoid string copy and transformation for better performance
    if (residue_name.length() != 3) {
        return false;
    }

    // Convert to uppercase using bit manipulation (faster than std::toupper)
    char c0 = residue_name[0] & ~0x20; // Convert to uppercase
    char c1 = residue_name[1] & ~0x20;
    char c2 = residue_name[2] & ~0x20;

    // Check common water residue names
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
    // Trim the residue name
    std::string trimmed = residue_name;
    while (!trimmed.empty() && trimmed[0] == ' ') {
        trimmed.erase(0, 1);
    }
    while (!trimmed.empty() && trimmed.back() == ' ') {
        trimmed.pop_back();
    }

    // Common modified nucleotides that appear as HETATM
    // These should be auto-included even without -T flag
    static const std::vector<std::string> modified_nucleotides = {
        // Pseudouridine and variants
        "PSU", "P5P", "PU",
        // Modified Adenine
        "A2M", "1MA", "2MA", "6MA", "OMA", "MIA", "I6A", "T6A", "M6A",
        "A23", // 2-aminoadenine
        "DA",  // Deoxyadenosine
        // Modified Cytosine
        "5MC", "OMC", "S4C", "5IC", "5FC", "CBR",
        "CVC", // Cytidine derivative
        "CM0", // Modified cytosine
        // Modified Guanine
        "OMG", "1MG", "2MG", "7MG", "M2G", "YYG", "YG", "QUO",
        // Modified Uracil/Thymine
        "5MU", "H2U", "DHU", "OMU", "4SU", "S4U", "5BU", "2MU", "UR3", "RT", "70U",
        "2YR", // Dihydrouridine derivative
        // Other common modified bases
        "I", "DI", // Inosine
        "EPE",     // Ethylpseudouridine or similar
        "J48",     // Modified base
        "KIR",     // Modified base
        "NCA",     // N-carboxyaminoadenine or similar
        "NF2",     // Modified base
        "NMN",     // Nicotinamide mononucleotide or modified base
        "NNR",     // Modified base
        "WVQ"      // Modified base
    };

    for (const auto& mod : modified_nucleotides) {
        if (trimmed == mod) {
            return true;
        }
    }

    return false;
}

// Atom name normalization helpers
std::string PdbParser::apply_atom_name_formatting_rules(const std::string& name) const {
    if (name.length() != 4) {
        return name;
    }

    // Rule 1: Move first 3 chars to " X" format if applicable
    if (name[0] != ' ' && !std::isdigit(static_cast<unsigned char>(name[0])) &&
        (name[1] == ' ' || std::isdigit(static_cast<unsigned char>(name[1]))) && name[3] == ' ') {
        return " " + name.substr(0, 3);
    }

    // Rule 2: Move chars 2-3 to " X " format if applicable
    if (name[0] == ' ' && name[1] == ' ' && std::isdigit(static_cast<unsigned char>(name[3]))) {
        return " " + name.substr(2, 2) + " ";
    }

    // Rule 3: Normalize P formatting
    if (name == "   P" || name == "P   ") {
        return " P  ";
    }

    // Rule 4: Handle "OP1 " and "OP2 " formats
    if (name == "OP1 ") {
        return " O1P";
    }
    if (name == "OP2 ") {
        return " O2P";
    }

    return name;
}

std::string PdbParser::apply_atom_name_exact_matches(const std::string& name) const {
    if (name == " O1'") {
        return " O4'";
    }
    if (name == " OL ") {
        return " O1P";
    }
    if (name == " OP1") {
        return " O1P";
    }
    if (name == " OR ") {
        return " O2P";
    }
    if (name == " OP2") {
        return " O2P";
    }
    if (name == " OP3") {
        return " O3P";
    }
    if (name == " C5A") {
        return " C5M";
    }
    if (name == " O5T") {
        return " O5'";
    }
    if (name == " O3T") {
        return " O3'";
    }
    if (name == "   P" || name == "P   ") {
        return " P  ";
    }

    // Check trimmed version for phosphate names
    std::string trimmed = name;
    trimmed.erase(0, trimmed.find_first_not_of(" \t"));
    trimmed.erase(trimmed.find_last_not_of(" \t") + 1);

    if (trimmed == "OP1") {
        return " O1P";
    }
    if (trimmed == "OP2") {
        return " O2P";
    }
    if (trimmed == "OP3") {
        return " O3P";
    }
    if (trimmed == "OP1 ") {
        return " O1P";
    }
    if (trimmed == "OP2 ") {
        return " O2P";
    }

    return name;
}

std::string PdbParser::ensure_atom_name_length(const std::string& name) const {
    if (name.length() == 4) {
        return name;
    }
    if (name.length() < 4) {
        return std::string(4 - name.length(), ' ') + name;
    }
    return name.substr(0, 4);
}

std::string PdbParser::normalize_atom_name(const std::string& name) const {
    if (name.empty()) {
        return "    ";
    }

    // Ensure 4 characters first
    std::string normalized = ensure_atom_name_length(name);

    // Apply formatting rules
    normalized = apply_atom_name_formatting_rules(normalized);

    // Handle character replacement: '*' becomes '\''
    if (normalized.length() == 4 && normalized[3] == '*') {
        normalized[3] = '\'';
    }

    // Apply exact matches
    normalized = apply_atom_name_exact_matches(normalized);

    // Ensure final length
    return ensure_atom_name_length(normalized);
}

std::string PdbParser::normalize_residue_name(const std::string& name) const {
    if (name.empty()) {
        return "   ";
    }
    if (name.length() == 3) {
        return name;
    }
    if (name.length() < 3) {
        return std::string(3 - name.length(), ' ') + name;
    }
    return name.substr(0, 3);
}

// Filtering helpers
bool PdbParser::check_alt_loc_filter(char alt_loc) const {
    // Default ALT_LIST is "A1" (with space added) = " A1"
    // Keep atoms with alt_loc = ' ', 'A', or '1'
    return (alt_loc == ' ' || alt_loc == 'A' || alt_loc == '1');
}

bool PdbParser::should_keep_atom(const core::Atom& atom) const {
    // Check occupancy filter
    bool keep_by_occupancy = !filter_by_occupancy_ || atom.occupancy() > 0.0;
    if (!keep_by_occupancy) {
        return false;
    }

    // Check alt_loc filter
    return check_alt_loc_filter(atom.alt_loc());
}

// Record processing helpers
void PdbParser::process_atom_record(
    const std::string& line, size_t line_number, int model_number,
    std::map<std::tuple<std::string, char, int, char>, std::vector<core::Atom>>& residue_atoms, int& legacy_atom_idx,
    int& legacy_residue_idx, std::map<std::tuple<std::string, char, int, char>, int>& legacy_residue_idx_map) {
    try {
        core::Atom atom = parse_atom_line(line, line_number);
        atom.set_model_number(model_number);

        // Only assign legacy indices to atoms that pass filtering
        // This matches the legacy code's behavior: atom_idx is only assigned to kept atoms
        if (should_keep_atom(atom)) {
            // Assign legacy atom index sequentially (1-based) as kept atoms are encountered
            // This matches the legacy code's behavior of assigning indices in PDB file order
            atom.set_legacy_atom_idx(legacy_atom_idx++);

            // Assign legacy residue index (1-based) - assign new index when residue changes
            // CRITICAL: Legacy groups by (ResName, ChainID, ResSeq, insertion) - must include
            // ResName! Legacy code: sprintf(bidx, "%3s%c%4ld%c", ResName, ChainID, ResSeq, iCode)
            auto residue_key = std::make_tuple(atom.residue_name(), atom.chain_id(), atom.residue_seq(),
                                               atom.insertion());
            auto residue_it = legacy_residue_idx_map.find(residue_key);
            if (residue_it == legacy_residue_idx_map.end()) {
                // New residue - assign next residue index
                legacy_residue_idx_map[residue_key] = legacy_residue_idx++;
            }
            atom.set_legacy_residue_idx(legacy_residue_idx_map[residue_key]);

            // Use (residue_name, chain_id, residue_seq, insertion_code) as key to match legacy
            // behavior Legacy uses ResName + ChainID + ResSeq + iCode to identify residues
            residue_atoms[{atom.residue_name(), atom.chain_id(), atom.residue_seq(), atom.insertion()}].push_back(atom);
        }
    } catch (const ParseError&) {
        // Skip malformed lines
    }
}

void PdbParser::process_hetatm_record(
    const std::string& line, size_t line_number, int model_number,
    std::map<std::tuple<std::string, char, int, char>, std::vector<core::Atom>>& residue_atoms, int& legacy_atom_idx,
    int& legacy_residue_idx, std::map<std::tuple<std::string, char, int, char>, int>& legacy_residue_idx_map) {

    try {
        core::Atom atom = parse_hetatm_line(line, line_number);

        // Auto-include modified nucleotides even if include_hetatm_ is false
        // Check if residue name is a known modified nucleotide
        bool is_modified_nucleotide = is_modified_nucleotide_name(atom.residue_name());

        if (!include_hetatm_ && !is_modified_nucleotide) {
            return;
        }
        atom.set_model_number(model_number);

        // Check if water and should be skipped
        if (!include_waters_ && is_water(atom.residue_name())) {
            return;
        }

        // Only assign legacy indices to atoms that pass filtering
        // This matches the legacy code's behavior: atom_idx is only assigned to kept atoms
        if (should_keep_atom(atom)) {
            // Assign legacy atom index sequentially (1-based) as kept atoms are encountered
            // This matches the legacy code's behavior of assigning indices in PDB file order
            atom.set_legacy_atom_idx(legacy_atom_idx++);

            // Assign legacy residue index (1-based) - assign new index when residue changes
            // CRITICAL: Legacy groups by (ResName, ChainID, ResSeq, insertion) - must include
            // ResName! Legacy code: sprintf(bidx, "%3s%c%4ld%c", ResName, ChainID, ResSeq, iCode)
            auto residue_key = std::make_tuple(atom.residue_name(), atom.chain_id(), atom.residue_seq(),
                                               atom.insertion());
            auto residue_it = legacy_residue_idx_map.find(residue_key);
            if (residue_it == legacy_residue_idx_map.end()) {
                // New residue - assign next residue index
                legacy_residue_idx_map[residue_key] = legacy_residue_idx++;
            }
            atom.set_legacy_residue_idx(legacy_residue_idx_map[residue_key]);

            // Use (residue_name, chain_id, residue_seq, insertion_code) as key to match legacy
            // behavior Legacy uses ResName + ChainID + ResSeq + iCode to identify residues
            residue_atoms[{atom.residue_name(), atom.chain_id(), atom.residue_seq(), atom.insertion()}].push_back(atom);
        }
    } catch (const ParseError&) {
        // Skip malformed lines
    }
}

// Model/END record handling
void PdbParser::handle_model_record(const std::string& line, int& current_model_number, bool& /* all_models */) {
    if (line.length() > 10) {
        try {
            std::string model_str = line.substr(6, 4);
            current_model_number = std::stoi(model_str);
        } catch (...) {
            current_model_number++;
        }
    } else {
        current_model_number++;
    }
    // Note: all_models is not currently configurable, stays false
}

bool PdbParser::handle_end_record(const std::string& line, bool all_models) {
    // Optimized: Use direct character comparison instead of substr
    // If it's ENDMDL and we're processing all models, continue to next model
    if (all_models && line.length() >= 6 && line[0] == 'E' && line[1] == 'N' && line[2] == 'D' && line[3] == 'M' &&
        line[4] == 'D' && line[5] == 'L') {
        return false; // Don't stop, continue
    }
    // Otherwise, stop at first END
    return true; // Stop processing
}

// Structure building - optimized for speed
core::Structure PdbParser::build_structure_from_residues(
    const std::string& pdb_id,
    const std::map<std::tuple<std::string, char, int, char>, std::vector<core::Atom>>& residue_atoms) const {

    core::Structure structure(pdb_id);
    std::map<char, core::Chain> chains;

    for (const auto& [key, atoms] : residue_atoms) {
        if (atoms.empty()) {
            continue;
        }

        auto [residue_name, chain_id, residue_seq, insertion_code] = key;

        // Use ResidueFactory to create residue with all properties initialized
        core::Residue residue = core::ResidueFactory::create(residue_name, residue_seq, chain_id, insertion_code,
                                                             atoms);

        // Use try_emplace for efficiency (C++17) - avoids unnecessary chain copy
        auto [it, inserted] = chains.try_emplace(chain_id, chain_id);
        it->second.add_residue(residue);
    }

    // Add all chains to structure
    for (auto& [chain_id, chain] : chains) {
        structure.add_chain(chain);
    }

    return structure;
}

// Main parsing implementation - optimized for speed
core::Structure PdbParser::parse_impl(std::istream& stream, const std::string& pdb_id) {
    // Key: (residue_name, chain_id, residue_seq, insertion_code) - matches legacy residue
    // identification Legacy groups by (ResName, ChainID, ResSeq, insertion) - sprintf(bidx,
    // "%3s%c%4ld%c", ...)
    std::map<std::tuple<std::string, char, int, char>, std::vector<core::Atom>> residue_atoms;
    size_t line_number = 0;
    std::string line;
    bool all_models = false;

    // Legacy indices: assign sequentially as atoms/residues are encountered (1-based)
    int legacy_atom_idx = 1;
    // CRITICAL: Legacy groups by (ResName, ChainID, ResSeq, insertion) - must include ResName!
    // Legacy code: sprintf(bidx, "%3s%c%4ld%c", ResName, ChainID, ResSeq, iCode)
    std::map<std::tuple<std::string, char, int, char>, int> legacy_residue_idx_map;
    int legacy_residue_idx = 1;
    int current_model_number = 0;

    // Reserve space for line to avoid reallocations
    line.reserve(120); // Typical PDB line length (80 chars + buffer)

    while (std::getline(stream, line)) {
        line_number++;

        // Fast empty line check - skip immediately if empty
        if (line.empty()) {
            continue;
        }

        // Quick whitespace-only check - only check first char
        if (line[0] == ' ' || line[0] == '\t' || line[0] == '\0') {
            // If starts with whitespace, check if entire line is whitespace
            if (line.find_first_not_of(" \t") == std::string::npos) {
                continue;
            }
        }

        // Fast prefix checks using direct character comparison instead of substr
        if (line.length() >= 4) {
            // ATOM records - most common case first
            if (line[0] == 'A' && line[1] == 'T' && line[2] == 'O' && line[3] == 'M' && line.length() >= 52) {
                process_atom_record(line, line_number, current_model_number, residue_atoms, legacy_atom_idx,
                                    legacy_residue_idx, legacy_residue_idx_map);
                continue;
            }

            // END records
            if (line.length() >= 3 && line[0] == 'E' && line[1] == 'N' && line[2] == 'D') {
                if (handle_end_record(line, all_models)) {
                    break;
                }
                continue;
            }

            // HETATM records
            if (line.length() >= 6 && line[0] == 'H' && line[1] == 'E' && line[2] == 'T' && line[3] == 'A' &&
                line[4] == 'T' && line[5] == 'M' && line.length() >= 52) {
                process_hetatm_record(line, line_number, current_model_number, residue_atoms, legacy_atom_idx,
                                      legacy_residue_idx, legacy_residue_idx_map);
                continue;
            }

            // MODEL records
            if (line.length() >= 6 && line[0] == 'M' && line[1] == 'O' && line[2] == 'D' && line[3] == 'E' &&
                line[4] == 'L' && line[5] == ' ') {
                handle_model_record(line, current_model_number, all_models);
                continue;
            }

            // HEADER records
            if (line.length() >= 6 && line[0] == 'H' && line[1] == 'E' && line[2] == 'A' && line[3] == 'D' &&
                line[4] == 'E' && line[5] == 'R' && pdb_id.empty() && line.length() >= 66) {
                std::string header_id = line.substr(62, 4);
                header_id.erase(0, header_id.find_first_not_of(" \t"));
                header_id.erase(header_id.find_last_not_of(" \t") + 1);
                // Note: PDB ID extraction from HEADER not fully implemented
            }
        }
    }

    return build_structure_from_residues(pdb_id, residue_atoms);
}

} // namespace io
} // namespace x3dna

/**
 * @file input_file_parser.cpp
 * @brief Implementation of input file parser
 */

#include <x3dna/io/input_file_parser.hpp>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

namespace x3dna {
namespace io {

InputData InputFileParser::parse(const std::filesystem::path& input_file) {
    std::ifstream file(input_file);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open input file: " + input_file.string());
    }
    return parse_stream(file);
}

InputData InputFileParser::parse_stream(std::istream& stream) {
    InputData data;
    std::string line;
    size_t line_number = 0;
    
    // Line 1: PDB file path
    if (!std::getline(stream, line)) {
        throw std::runtime_error("Input file is empty");
    }
    line_number++;
    // Trim whitespace
    line.erase(0, line.find_first_not_of(" \t"));
    line.erase(line.find_last_not_of(" \t") + 1);
    data.pdb_file = line;
    
    // Line 2: Output file name
    if (!std::getline(stream, line)) {
        throw std::runtime_error("Input file missing output file name");
    }
    line_number++;
    line.erase(0, line.find_first_not_of(" \t"));
    line.erase(line.find_last_not_of(" \t") + 1);
    data.output_file = line;
    
    // Line 3: Duplex number
    if (!std::getline(stream, line)) {
        throw std::runtime_error("Input file missing duplex number");
    }
    line_number++;
    std::istringstream iss(line);
    iss >> data.duplex_number;
    
    // Line 4: Number of base pairs
    if (!std::getline(stream, line)) {
        throw std::runtime_error("Input file missing number of base pairs");
    }
    line_number++;
    iss.clear();
    iss.str(line);
    iss >> data.num_base_pairs;
    
    // Line 5: Flags
    if (!std::getline(stream, line)) {
        throw std::runtime_error("Input file missing flags");
    }
    line_number++;
    iss.clear();
    iss.str(line);
    iss >> data.flags;
    
    // Remaining lines: Base pair data
    // Format: bp_num res1 res2 flag # comment
    // Note: res1 and res2 are 1-based residue indices in the original
    // We convert to 0-based for our internal representation
    data.base_pairs.reserve(data.num_base_pairs);
    
    while (std::getline(stream, line)) {
        line_number++;
        
        // Skip empty lines
        if (line.empty() || line.find_first_not_of(" \t") == std::string::npos) {
            continue;
        }
        
        // Check for comment lines (criteria, helix info, etc.)
        if (line.find("#####") == 0) {
            if (line.find("Base-pair criteria") != std::string::npos) {
                data.criteria_line = line;
            } else if (line.find("Helix #") != std::string::npos) {
                data.helix_info.push_back(line);
            }
            continue;
        }
        
        // Parse base pair line
        // Format: bp_num res1 res2 flag # comment
        try {
            auto [res1, res2] = parse_base_pair_line(line, line_number);
            
            // Create base pair (type will be determined later during analysis)
            core::BasePair bp(res1, res2, core::BasePairType::UNKNOWN);
            data.base_pairs.push_back(bp);
        } catch (const std::exception& e) {
            // Skip malformed lines, but log warning
            continue;
        }
    }
    
    return data;
}

std::pair<size_t, size_t> InputFileParser::parse_base_pair_line(const std::string& line, size_t line_number) {
    std::istringstream iss(line);
    
    // Skip base pair number
    int bp_num;
    if (!(iss >> bp_num)) {
        throw std::runtime_error("Cannot parse base pair number at line " + std::to_string(line_number));
    }
    
    // Read residue indices (1-based in file, convert to 0-based)
    int res1, res2;
    if (!(iss >> res1 >> res2)) {
        throw std::runtime_error("Cannot parse residue indices at line " + std::to_string(line_number));
    }
    
    // Convert from 1-based to 0-based
    size_t res1_0based = static_cast<size_t>(res1 - 1);
    size_t res2_0based = static_cast<size_t>(res2 - 1);
    
    return {res1_0based, res2_0based};
}

} // namespace io
} // namespace x3dna


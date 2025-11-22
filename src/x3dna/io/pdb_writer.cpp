/**
 * @file pdb_writer.cpp
 * @brief Implementation of PDB writer
 */

#include <x3dna/io/pdb_writer.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/atom.hpp>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace x3dna {
namespace io {

void PdbWriter::write_file(const core::Structure& structure, const std::filesystem::path& path) {
    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + path.string());
    }
    write_stream(structure, file);
}

void PdbWriter::write_stream(const core::Structure& structure, std::ostream& stream) {
    int atom_serial = 1; // 1-based serial numbers
    
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            for (const auto& atom : residue.atoms()) {
                std::string line = format_atom_line(atom, atom_serial++);
                stream << line << "\n";
            }
        }
    }
    
    // Write END record
    stream << "END\n";
}

std::string PdbWriter::to_string(const core::Structure& structure) {
    std::ostringstream stream;
    write_stream(structure, stream);
    return stream.str();
}

std::string PdbWriter::format_atom_line(const core::Atom& atom, int atom_serial) const {
    std::ostringstream line;
    
    // Determine record type
    std::string record_type = (atom.record_type() == 'H') ? "HETATM" : "ATOM  ";
    
    // Format: RECORD_TYPE SERIAL NAME RESNAME CHAIN RESSEQ X Y Z OCCUPANCY B_FACTOR
    // Columns: 1-6, 7-11, 13-16, 18-20, 22, 23-26, 31-38, 39-46, 47-54, 55-60, 61-66
    
    line << std::left << std::setw(6) << record_type;
    line << std::right << std::setw(5) << atom_serial;
    line << " ";
    line << std::left << std::setw(4) << atom.name();
    line << std::left << std::setw(3) << atom.residue_name();
    line << " ";
    line << atom.chain_id();
    line << std::right << std::setw(4) << atom.residue_seq();
    line << "    ";
    
    // Coordinates (8.3 format)
    const auto& pos = atom.position();
    line << format_coordinate(pos.x());
    line << format_coordinate(pos.y());
    line << format_coordinate(pos.z());
    
    // Occupancy (6.2 format)
    line << std::fixed << std::setprecision(2) << std::setw(6) << atom.occupancy();
    
    // B-factor (6.2 format)
    line << std::fixed << std::setprecision(2) << std::setw(6) << atom.b_factor();
    
    // Element (optional, columns 77-78)
    if (!atom.element().empty()) {
        line << "          " << std::setw(2) << atom.element();
    }
    
    return line.str();
}

std::string PdbWriter::format_coordinate(double coord) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << std::setw(8) << coord;
    return oss.str();
}

} // namespace io
} // namespace x3dna


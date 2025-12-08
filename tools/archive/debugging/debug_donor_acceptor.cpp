/**
 * @file debug_donor_acceptor.cpp
 * @brief Debug tool to test donor_acceptor function for specific atom pairs
 *
 * This tool helps debug why H-bond types are classified as '*' vs '-'
 */

#include <iostream>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <iomanip>

using namespace x3dna::algorithms;

void test_donor_acceptor(char base1, char base2, const std::string& atom1, const std::string& atom2) {
    char type = BasePairValidator::donor_acceptor(base1, base2, atom1, atom2);

    std::cout << "  Base pair: " << base1 << "-" << base2 << "\n";
    std::cout << "  Atoms: " << atom1 << " -> " << atom2 << "\n";
    std::cout << "  Type: " << type << " (" << (type == '-' ? "standard" : "non-standard") << ")\n";
    std::cout << "\n";
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <base1> <base2> <atom1> <atom2>\n";
        std::cerr << "Example: " << argv[0] << " C G \" N1 \" \" O2'\"\n";
        std::cerr << "Example: " << argv[0] << " C G \" N3 \" \" N2 \"\n";
        return 1;
    }

    char base1 = argv[1][0];
    char base2 = argv[2][0];
    std::string atom1 = argv[3];
    std::string atom2 = argv[4];

    std::cout << "========================================\n";
    std::cout << "Donor-Acceptor Type Test\n";
    std::cout << "========================================\n";

    test_donor_acceptor(base1, base2, atom1, atom2);

    return 0;
}

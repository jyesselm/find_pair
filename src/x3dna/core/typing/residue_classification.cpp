/**
 * @file residue_classification.cpp
 * @brief Implementation of ResidueClassification legacy compatibility methods
 */

#include "x3dna/core/typing/residue_classification.hpp"
#include "x3dna/core/residue_type.hpp"

namespace x3dna {
namespace core {
namespace typing {

ResidueType ResidueClassification::to_legacy_type() const {
    if (molecule_type == MoleculeType::WATER) {
        return ResidueType::WATER;
    }
    if (molecule_type == MoleculeType::ION) {
        return ResidueType::ION;
    }
    if (molecule_type == MoleculeType::LIGAND) {
        return ResidueType::LIGAND;
    }
    if (molecule_type == MoleculeType::PROTEIN) {
        return ResidueType::AMINO_ACID;
    }
    if (molecule_type != MoleculeType::NUCLEIC_ACID) {
        return ResidueType::UNKNOWN;
    }

    // For nucleic acids, map base type to legacy enum
    switch (base_type) {
        case BaseType::ADENINE: return ResidueType::ADENINE;
        case BaseType::GUANINE: return ResidueType::GUANINE;
        case BaseType::CYTOSINE: return ResidueType::CYTOSINE;
        case BaseType::THYMINE: return ResidueType::THYMINE;
        case BaseType::URACIL: return ResidueType::URACIL;
        case BaseType::PSEUDOURIDINE: return ResidueType::PSEUDOURIDINE;
        case BaseType::INOSINE: return ResidueType::INOSINE;
        default:
            return is_modified_nucleotide ? ResidueType::NONCANONICAL_RNA : ResidueType::NUCLEOTIDE;
    }
}

ResidueClassification ResidueClassification::from_legacy(
    ResidueType type,
    const std::string& residue_name,
    bool is_purine_hint) {

    ResidueClassification result;
    result.residue_name = residue_name;

    switch (type) {
        case ResidueType::WATER:
            result.molecule_type = MoleculeType::WATER;
            result.solvent_type = SolventType::WATER;
            return result;

        case ResidueType::ION:
            result.molecule_type = MoleculeType::ION;
            return result;

        case ResidueType::LIGAND:
            result.molecule_type = MoleculeType::LIGAND;
            return result;

        case ResidueType::AMINO_ACID:
            result.molecule_type = MoleculeType::PROTEIN;
            return result;

        case ResidueType::UNKNOWN:
            return result;

        default:
            break;
    }

    // Nucleotide types
    result.molecule_type = MoleculeType::NUCLEIC_ACID;

    // Detect DNA vs RNA from residue name
    if (residue_name.size() >= 2 && residue_name[0] == 'D') {
        result.nucleic_acid_type = NucleicAcidType::DNA;
    } else {
        result.nucleic_acid_type = NucleicAcidType::RNA;
    }

    // Map legacy type to base type
    switch (type) {
        case ResidueType::ADENINE:
            result.base_type = BaseType::ADENINE;
            result.canonical_code = 'A';
            result.one_letter_code = 'A';
            break;
        case ResidueType::GUANINE:
            result.base_type = BaseType::GUANINE;
            result.canonical_code = 'G';
            result.one_letter_code = 'G';
            break;
        case ResidueType::CYTOSINE:
            result.base_type = BaseType::CYTOSINE;
            result.canonical_code = 'C';
            result.one_letter_code = 'C';
            break;
        case ResidueType::THYMINE:
            result.base_type = BaseType::THYMINE;
            result.canonical_code = 'T';
            result.one_letter_code = 'T';
            result.nucleic_acid_type = NucleicAcidType::DNA;
            break;
        case ResidueType::URACIL:
            result.base_type = BaseType::URACIL;
            result.canonical_code = 'U';
            result.one_letter_code = 'U';
            result.nucleic_acid_type = NucleicAcidType::RNA;
            break;
        case ResidueType::PSEUDOURIDINE:
            result.base_type = BaseType::PSEUDOURIDINE;
            result.canonical_code = 'U';
            result.one_letter_code = 'P';
            result.is_modified_nucleotide = true;
            break;
        case ResidueType::INOSINE:
            result.base_type = BaseType::INOSINE;
            result.canonical_code = 'I';
            result.one_letter_code = 'I';
            result.is_modified_nucleotide = true;
            break;
        case ResidueType::NONCANONICAL_RNA:
            result.is_modified_nucleotide = true;
            result.base_category = is_purine_hint ? BaseCategory::PURINE : BaseCategory::PYRIMIDINE;
            break;
        default:
            break;
    }

    // Set base category from base type
    result.base_category = get_base_category(result.base_type);

    return result;
}

} // namespace typing
} // namespace core
} // namespace x3dna

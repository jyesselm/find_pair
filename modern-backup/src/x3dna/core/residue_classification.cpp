/**
 * @file residue_classification.cpp
 * @brief Implementation of hierarchical residue classification system
 */

#include "x3dna/core/residue_classification.hpp"

namespace x3dna {
namespace core {

bool ResidueClassification::is_nucleotide() const {
    return molecule_type == MoleculeType::NUCLEIC_ACID;
}

bool ResidueClassification::is_rna() const {
    return molecule_type == MoleculeType::NUCLEIC_ACID && nucleic_acid_type == NucleicAcidType::RNA;
}

bool ResidueClassification::is_dna() const {
    return molecule_type == MoleculeType::NUCLEIC_ACID && nucleic_acid_type == NucleicAcidType::DNA;
}

bool ResidueClassification::is_purine() const {
    return base_category == BaseCategory::PURINE;
}

bool ResidueClassification::is_pyrimidine() const {
    return base_category == BaseCategory::PYRIMIDINE;
}

bool ResidueClassification::is_canonical() const {
    return is_nucleotide() && !is_modified;
}

bool ResidueClassification::is_protein() const {
    return molecule_type == MoleculeType::PROTEIN;
}

bool ResidueClassification::is_water() const {
    return molecule_type == MoleculeType::WATER;
}

bool ResidueClassification::is_ion() const {
    return molecule_type == MoleculeType::ION;
}

bool ResidueClassification::is_ligand() const {
    return molecule_type == MoleculeType::LIGAND;
}

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
    if (base_type == BaseType::ADENINE) {
        return ResidueType::ADENINE;
    }
    if (base_type == BaseType::GUANINE) {
        return ResidueType::GUANINE;
    }
    if (base_type == BaseType::CYTOSINE) {
        return ResidueType::CYTOSINE;
    }
    if (base_type == BaseType::THYMINE) {
        return ResidueType::THYMINE;
    }
    if (base_type == BaseType::URACIL) {
        return ResidueType::URACIL;
    }
    if (base_type == BaseType::PSEUDOURIDINE) {
        return ResidueType::PSEUDOURIDINE;
    }
    if (base_type == BaseType::INOSINE) {
        return ResidueType::INOSINE;
    }

    return is_modified ? ResidueType::NONCANONICAL_RNA : ResidueType::NUCLEOTIDE;
}

ResidueClassification ResidueClassification::from_legacy(ResidueType type, const std::string& residue_name,
                                                         bool is_purine_hint) {

    ResidueClassification result;
    result.residue_name = residue_name;

    if (type == ResidueType::WATER) {
        result.molecule_type = MoleculeType::WATER;
        return result;
    }

    if (type == ResidueType::ION) {
        result.molecule_type = MoleculeType::ION;
        return result;
    }

    if (type == ResidueType::LIGAND) {
        result.molecule_type = MoleculeType::LIGAND;
        return result;
    }

    if (type == ResidueType::AMINO_ACID) {
        result.molecule_type = MoleculeType::PROTEIN;
        return result;
    }

    if (type == ResidueType::UNKNOWN) {
        return result;
    }

    // Nucleotide types
    result.molecule_type = MoleculeType::NUCLEIC_ACID;

    // Detect DNA vs RNA from residue name
    // DNA residues typically start with 'D' (DA, DC, DG, DT)
    if (residue_name.size() >= 2 && residue_name[0] == 'D') {
        result.nucleic_acid_type = NucleicAcidType::DNA;
    } else {
        result.nucleic_acid_type = NucleicAcidType::RNA;
    }

    // Map legacy type to base type
    if (type == ResidueType::ADENINE) {
        result.base_type = BaseType::ADENINE;
        result.canonical_code = 'A';
    } else if (type == ResidueType::GUANINE) {
        result.base_type = BaseType::GUANINE;
        result.canonical_code = 'G';
    } else if (type == ResidueType::CYTOSINE) {
        result.base_type = BaseType::CYTOSINE;
        result.canonical_code = 'C';
    } else if (type == ResidueType::THYMINE) {
        result.base_type = BaseType::THYMINE;
        result.canonical_code = 'T';
        result.nucleic_acid_type = NucleicAcidType::DNA;
    } else if (type == ResidueType::URACIL) {
        result.base_type = BaseType::URACIL;
        result.canonical_code = 'U';
        result.nucleic_acid_type = NucleicAcidType::RNA;
    } else if (type == ResidueType::PSEUDOURIDINE) {
        result.base_type = BaseType::PSEUDOURIDINE;
        result.canonical_code = 'U';
        result.is_modified = true;
    } else if (type == ResidueType::INOSINE) {
        result.base_type = BaseType::INOSINE;
        result.canonical_code = 'I';
        result.is_modified = true;
    } else if (type == ResidueType::NONCANONICAL_RNA) {
        result.is_modified = true;
        // Use purine hint to set category
        result.base_category = is_purine_hint ? BaseCategory::PURINE : BaseCategory::PYRIMIDINE;
    }

    // Set base category from base type
    result.base_category = get_base_category(result.base_type);

    return result;
}

} // namespace core
} // namespace x3dna

/**
 * @file residue_factory.hpp
 * @brief Factory for creating polymorphic residue types
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <x3dna/core/atom.hpp>
#include <x3dna/core/residue/iresidue.hpp>
#include <x3dna/core/residue/rna.hpp>
#include <x3dna/core/residue/dna.hpp>
#include <x3dna/core/residue/protein.hpp>
#include <x3dna/core/residue/ligand.hpp>
#include <x3dna/core/typing/type_registry.hpp>
#include <x3dna/core/modified_nucleotide_registry.hpp>
#include <x3dna/core/string_utils.hpp>

namespace x3dna {
namespace core {
namespace poly {


/**
 * @class ResidueFactory
 * @brief Creates polymorphic residue objects based on residue name classification
 */
class ResidueFactory {
public:
    /**
     * @brief Create a residue from parsed PDB/CIF data
     *
     * Uses TypeRegistry to classify the residue and creates the appropriate
     * concrete type (RNA, DNA, Protein, or Ligand).
     *
     * @param name Residue name (e.g., "A", "DA", "ALA", "HOH")
     * @param seq_num Sequence number
     * @param chain_id Chain identifier
     * @param insertion Insertion code (default "")
     * @param atoms Atoms belonging to this residue
     * @return unique_ptr to the appropriate residue type
     */
    [[nodiscard]] static std::unique_ptr<IResidue> create(
        const std::string& name,
        int seq_num,
        const std::string& chain_id,
        const std::string& insertion,
        const std::vector<Atom>& atoms) {

        std::string trimmed_name = trim(name);
        auto classification = typing::TypeRegistry::instance().classify_residue(trimmed_name);

        std::unique_ptr<IResidue> residue;

        if (classification.is_rna()) {
            auto rna = std::make_unique<RNA>(trimmed_name, seq_num, chain_id, insertion);
            rna->set_classification(classification);
            rna->set_one_letter_code(ModifiedNucleotideRegistry::get_one_letter_code(trimmed_name));
            for (const auto& atom : atoms) {
                rna->add_atom(atom);
            }
            residue = std::move(rna);
        } else if (classification.is_dna()) {
            auto dna = std::make_unique<DNA>(trimmed_name, seq_num, chain_id, insertion);
            dna->set_classification(classification);
            dna->set_one_letter_code(ModifiedNucleotideRegistry::get_one_letter_code(trimmed_name));
            for (const auto& atom : atoms) {
                dna->add_atom(atom);
            }
            residue = std::move(dna);
        } else if (classification.is_protein()) {
            auto protein = std::make_unique<Protein>(trimmed_name, seq_num, chain_id, insertion);
            protein->set_classification(classification);
            protein->set_one_letter_code(classification.one_letter_code);
            for (const auto& atom : atoms) {
                protein->add_atom(atom);
            }
            residue = std::move(protein);
        } else {
            // Water, ions, ligands, unknown
            auto ligand = std::make_unique<Ligand>(trimmed_name, seq_num, chain_id, insertion);
            ligand->set_classification(classification);
            for (const auto& atom : atoms) {
                ligand->add_atom(atom);
            }
            residue = std::move(ligand);
        }

        return residue;
    }

    /**
     * @brief Create a residue without atoms (atoms added later)
     */
    [[nodiscard]] static std::unique_ptr<IResidue> create(
        const std::string& name,
        int seq_num,
        const std::string& chain_id,
        const std::string& insertion = "") {

        return create(name, seq_num, chain_id, insertion, {});
    }

    /**
     * @brief Create an RNA nucleotide directly
     */
    [[nodiscard]] static std::unique_ptr<RNA> create_rna(
        const std::string& name,
        int seq_num,
        const std::string& chain_id,
        const std::string& insertion = "") {

        std::string trimmed_name = trim(name);
        auto classification = typing::TypeRegistry::instance().classify_residue(trimmed_name);
        auto rna = std::make_unique<RNA>(trimmed_name, seq_num, chain_id, insertion);
        rna->set_classification(classification);
        rna->set_one_letter_code(ModifiedNucleotideRegistry::get_one_letter_code(trimmed_name));
        return rna;
    }

    /**
     * @brief Create a DNA nucleotide directly
     */
    [[nodiscard]] static std::unique_ptr<DNA> create_dna(
        const std::string& name,
        int seq_num,
        const std::string& chain_id,
        const std::string& insertion = "") {

        std::string trimmed_name = trim(name);
        auto classification = typing::TypeRegistry::instance().classify_residue(trimmed_name);
        auto dna = std::make_unique<DNA>(trimmed_name, seq_num, chain_id, insertion);
        dna->set_classification(classification);
        dna->set_one_letter_code(ModifiedNucleotideRegistry::get_one_letter_code(trimmed_name));
        return dna;
    }

    /**
     * @brief Create a nucleotide, automatically determining RNA vs DNA
     */
    [[nodiscard]] static std::unique_ptr<INucleotide> create_nucleotide(
        const std::string& name,
        int seq_num,
        const std::string& chain_id,
        const std::string& insertion = "") {

        std::string trimmed_name = trim(name);
        auto classification = typing::TypeRegistry::instance().classify_residue(trimmed_name);

        if (classification.is_dna()) {
            return create_dna(name, seq_num, chain_id, insertion);
        }
        // Default to RNA for nucleotides
        return create_rna(name, seq_num, chain_id, insertion);
    }
};

} // namespace poly
} // namespace core
} // namespace x3dna

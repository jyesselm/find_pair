/**
 * @file iresidue.hpp
 * @brief Pure interface for all residue types
 */

#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <x3dna/core/atom.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/core/typing/residue_classification.hpp>
#include <x3dna/core/typing/nucleotide_type.hpp>

namespace x3dna {
namespace core {
namespace structure {

/**
 * @class IResidue
 * @brief Pure interface for all residue types (RNA, DNA, Protein, Ligand)
 */
class IResidue {
public:
    virtual ~IResidue() = default;

    // === Identity ===
    [[nodiscard]] virtual const std::string& name() const = 0;
    [[nodiscard]] virtual int seq_num() const = 0;
    [[nodiscard]] virtual const std::string& chain_id() const = 0;
    [[nodiscard]] virtual const std::string& insertion() const = 0;

    // === Atoms ===
    [[nodiscard]] virtual const std::vector<Atom>& atoms() const = 0;
    [[nodiscard]] virtual std::vector<Atom>& atoms() = 0;
    [[nodiscard]] virtual size_t num_atoms() const = 0;
    [[nodiscard]] virtual std::optional<Atom> find_atom(const std::string& atom_name) const = 0;
    virtual void add_atom(const Atom& atom) = 0;

    // === Type queries ===
    [[nodiscard]] virtual bool is_nucleotide() const = 0;
    [[nodiscard]] virtual bool is_rna() const = 0;
    [[nodiscard]] virtual bool is_dna() const = 0;
    [[nodiscard]] virtual bool is_protein() const = 0;
    [[nodiscard]] virtual bool is_ligand() const = 0;

    // === Legacy support ===
    [[nodiscard]] virtual int legacy_residue_idx() const = 0;
    virtual void set_legacy_residue_idx(int idx) = 0;
    [[nodiscard]] virtual std::pair<int, int> atom_range() const = 0;

    // === Classification ===
    [[nodiscard]] virtual const typing::ResidueClassification& classification() const = 0;

    // === Clone for copying ===
    [[nodiscard]] virtual std::unique_ptr<IResidue> clone() const = 0;
};

/**
 * @class INucleotide
 * @brief Interface for nucleotide residues (RNA and DNA)
 */
class INucleotide : public virtual IResidue {
public:
    // === Nucleotide-specific ===
    [[nodiscard]] virtual char one_letter_code() const = 0;
    [[nodiscard]] virtual bool is_purine() const = 0;
    [[nodiscard]] virtual bool is_pyrimidine() const = 0;
    [[nodiscard]] virtual typing::BaseType base_type() const = 0;
    [[nodiscard]] virtual int ry_classification() const = 0;

    // === Reference frame ===
    [[nodiscard]] virtual std::optional<ReferenceFrame> reference_frame() const = 0;
    virtual void set_reference_frame(const ReferenceFrame& frame) = 0;

    // === Ring atoms ===
    [[nodiscard]] virtual std::vector<Atom> ring_atoms() const = 0;

    // === IResidue defaults for nucleotides ===
    [[nodiscard]] bool is_nucleotide() const override { return true; }
    [[nodiscard]] bool is_protein() const override { return false; }
    [[nodiscard]] bool is_ligand() const override { return false; }
};

} // namespace structure
} // namespace core
} // namespace x3dna

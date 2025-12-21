/**
 * @file protein.hpp
 * @brief Protein (amino acid) residue class
 */

#pragma once

#include <limits>
#include <x3dna/core/residue/iresidue.hpp>
#include <x3dna/core/string_utils.hpp>

namespace x3dna {
namespace core {
namespace poly {


/**
 * @class Protein
 * @brief Represents a protein residue (amino acid)
 */
class Protein final : public IResidue {
public:
    Protein() = default;

    Protein(const std::string& name, int seq_num, const std::string& chain_id,
            const std::string& insertion = "")
        : name_(trim(name)), seq_num_(seq_num), chain_id_(chain_id), insertion_(insertion) {}

    // === Identity (IResidue) ===
    [[nodiscard]] const std::string& name() const override { return name_; }
    [[nodiscard]] int seq_num() const override { return seq_num_; }
    [[nodiscard]] const std::string& chain_id() const override { return chain_id_; }
    [[nodiscard]] const std::string& insertion() const override { return insertion_; }

    // === Atoms (IResidue) ===
    [[nodiscard]] const std::vector<Atom>& atoms() const override { return atoms_; }
    [[nodiscard]] std::vector<Atom>& atoms() override { return atoms_; }
    [[nodiscard]] size_t num_atoms() const override { return atoms_.size(); }

    void add_atom(const Atom& atom) override { atoms_.push_back(atom); }

    [[nodiscard]] std::optional<Atom> find_atom(const std::string& atom_name) const override {
        for (const auto& atom : atoms_) {
            if (atom.name() == atom_name) {
                return atom;
            }
            std::string trimmed = trim(atom_name);
            std::string atom_trimmed = trim(atom.name());
            if (atom_trimmed == trimmed) {
                return atom;
            }
        }
        return std::nullopt;
    }

    // === Type queries (IResidue) ===
    [[nodiscard]] bool is_nucleotide() const override { return false; }
    [[nodiscard]] bool is_rna() const override { return false; }
    [[nodiscard]] bool is_dna() const override { return false; }
    [[nodiscard]] bool is_protein() const override { return true; }
    [[nodiscard]] bool is_ligand() const override { return false; }

    // === Legacy support (IResidue) ===
    [[nodiscard]] int legacy_residue_idx() const override { return legacy_residue_idx_; }
    void set_legacy_residue_idx(int idx) override { legacy_residue_idx_ = idx; }

    [[nodiscard]] std::pair<int, int> atom_range() const override {
        if (atoms_.empty()) return {0, 0};
        int min_idx = std::numeric_limits<int>::max();
        int max_idx = 0;
        for (const auto& atom : atoms_) {
            int legacy_idx = atom.legacy_atom_idx();
            if (legacy_idx > 0) {
                min_idx = std::min(min_idx, legacy_idx);
                max_idx = std::max(max_idx, legacy_idx);
            }
        }
        if (min_idx == std::numeric_limits<int>::max()) return {0, 0};
        return {min_idx, max_idx};
    }

    // === Classification (IResidue) ===
    [[nodiscard]] const typing::ResidueClassification& classification() const override {
        return classification_;
    }
    void set_classification(const typing::ResidueClassification& c) { classification_ = c; }

    // === Clone (IResidue) ===
    [[nodiscard]] std::unique_ptr<IResidue> clone() const override {
        return std::unique_ptr<IResidue>(new Protein(*this));
    }

    // === Protein-specific ===
    [[nodiscard]] char one_letter_code() const { return one_letter_code_; }
    void set_one_letter_code(char code) { one_letter_code_ = code; }

private:
    std::string name_;
    int seq_num_ = 0;
    std::string chain_id_;
    std::string insertion_;
    std::vector<Atom> atoms_;
    typing::ResidueClassification classification_;
    int legacy_residue_idx_ = 0;
    char one_letter_code_ = '?';
};

} // namespace poly
} // namespace core
} // namespace x3dna

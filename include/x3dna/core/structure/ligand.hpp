/**
 * @file ligand.hpp
 * @brief Ligand residue class (water, ions, small molecules)
 */

#pragma once

#include <limits>
#include <x3dna/core/structure/iresidue.hpp>
#include <x3dna/core/string_utils.hpp>

namespace x3dna {
namespace core {
namespace structure {


/**
 * @class Ligand
 * @brief Represents a ligand residue (water, ion, or small molecule)
 */
class Ligand final : public IResidue {
public:
    Ligand() = default;

    Ligand(const std::string& name, int seq_num, const std::string& chain_id,
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
    [[nodiscard]] bool is_protein() const override { return false; }
    [[nodiscard]] bool is_ligand() const override { return true; }

    // === More specific ligand type queries ===
    [[nodiscard]] bool is_water() const { return classification_.is_water(); }
    [[nodiscard]] bool is_ion() const { return classification_.is_ion(); }

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
        return std::unique_ptr<IResidue>(new Ligand(*this));
    }

private:
    std::string name_;
    int seq_num_ = 0;
    std::string chain_id_;
    std::string insertion_;
    std::vector<Atom> atoms_;
    typing::ResidueClassification classification_;
    int legacy_residue_idx_ = 0;
};

} // namespace structure
} // namespace core
} // namespace x3dna

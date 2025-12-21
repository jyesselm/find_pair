/**
 * @file chain.hpp
 * @brief Chain class using polymorphic residue types
 */

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <x3dna/core/residue/iresidue.hpp>

namespace x3dna {
namespace core {
namespace poly {


/**
 * @class Chain
 * @brief Represents a chain of polymorphic residues
 *
 * This Chain stores unique_ptr<IResidue> for polymorphic residue handling.
 * Use dynamic_cast to access type-specific methods.
 */
class Chain {
public:
    Chain() = default;

    explicit Chain(const std::string& id) : chain_id_(id) {}

    // Move-only (no copying due to unique_ptr)
    Chain(const Chain&) = delete;
    Chain& operator=(const Chain&) = delete;
    Chain(Chain&&) = default;
    Chain& operator=(Chain&&) = default;

    // Deep copy via clone()
    [[nodiscard]] Chain clone() const {
        Chain copy(chain_id_);
        copy.residues_.reserve(residues_.size());
        for (const auto& res : residues_) {
            copy.residues_.push_back(res->clone());
        }
        return copy;
    }

    // === Identity ===
    [[nodiscard]] const std::string& chain_id() const { return chain_id_; }
    void set_chain_id(const std::string& id) { chain_id_ = id; }

    // === Residue access ===
    [[nodiscard]] size_t num_residues() const { return residues_.size(); }
    [[nodiscard]] size_t size() const { return residues_.size(); }
    [[nodiscard]] bool empty() const { return residues_.empty(); }

    [[nodiscard]] const IResidue& operator[](size_t idx) const { return *residues_[idx]; }
    [[nodiscard]] IResidue& operator[](size_t idx) { return *residues_[idx]; }

    [[nodiscard]] const IResidue& at(size_t idx) const { return *residues_.at(idx); }
    [[nodiscard]] IResidue& at(size_t idx) { return *residues_.at(idx); }

    // === Residue ownership ===
    void add_residue(std::unique_ptr<IResidue> residue) {
        residues_.push_back(std::move(residue));
    }

    // Get raw pointer (non-owning)
    [[nodiscard]] IResidue* get_residue(size_t idx) {
        return idx < residues_.size() ? residues_[idx].get() : nullptr;
    }

    [[nodiscard]] const IResidue* get_residue(size_t idx) const {
        return idx < residues_.size() ? residues_[idx].get() : nullptr;
    }

    // === Iteration (via raw pointers for convenience) ===
    class Iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = IResidue;
        using difference_type = std::ptrdiff_t;
        using pointer = IResidue*;
        using reference = IResidue&;

        Iterator(std::vector<std::unique_ptr<IResidue>>::iterator it) : it_(it) {}

        reference operator*() const { return **it_; }
        pointer operator->() const { return it_->get(); }

        Iterator& operator++() { ++it_; return *this; }
        Iterator operator++(int) { Iterator tmp = *this; ++it_; return tmp; }

        bool operator==(const Iterator& other) const { return it_ == other.it_; }
        bool operator!=(const Iterator& other) const { return it_ != other.it_; }

    private:
        std::vector<std::unique_ptr<IResidue>>::iterator it_;
    };

    class ConstIterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = const IResidue;
        using difference_type = std::ptrdiff_t;
        using pointer = const IResidue*;
        using reference = const IResidue&;

        ConstIterator(std::vector<std::unique_ptr<IResidue>>::const_iterator it) : it_(it) {}

        reference operator*() const { return **it_; }
        pointer operator->() const { return it_->get(); }

        ConstIterator& operator++() { ++it_; return *this; }
        ConstIterator operator++(int) { ConstIterator tmp = *this; ++it_; return tmp; }

        bool operator==(const ConstIterator& other) const { return it_ == other.it_; }
        bool operator!=(const ConstIterator& other) const { return it_ != other.it_; }

    private:
        std::vector<std::unique_ptr<IResidue>>::const_iterator it_;
    };

    [[nodiscard]] Iterator begin() { return Iterator(residues_.begin()); }
    [[nodiscard]] Iterator end() { return Iterator(residues_.end()); }
    [[nodiscard]] ConstIterator begin() const { return ConstIterator(residues_.begin()); }
    [[nodiscard]] ConstIterator end() const { return ConstIterator(residues_.end()); }

    // === Atom count ===
    [[nodiscard]] size_t num_atoms() const {
        size_t count = 0;
        for (const auto& res : residues_) {
            count += res->num_atoms();
        }
        return count;
    }

    // === Sequence ===
    [[nodiscard]] std::string sequence() const {
        std::string seq;
        for (const auto& res : residues_) {
            auto* nuc = dynamic_cast<const INucleotide*>(res.get());
            if (nuc) {
                char code = nuc->one_letter_code();
                if (code != '?') {
                    seq += code;
                }
            }
        }
        return seq;
    }

    // === Nucleotide access ===
    [[nodiscard]] std::vector<INucleotide*> nucleotides() {
        std::vector<INucleotide*> nts;
        for (auto& res : residues_) {
            if (auto* nuc = dynamic_cast<INucleotide*>(res.get())) {
                nts.push_back(nuc);
            }
        }
        return nts;
    }

    [[nodiscard]] std::vector<const INucleotide*> nucleotides() const {
        std::vector<const INucleotide*> nts;
        for (const auto& res : residues_) {
            if (auto* nuc = dynamic_cast<const INucleotide*>(res.get())) {
                nts.push_back(nuc);
            }
        }
        return nts;
    }

    // === Find residue ===
    [[nodiscard]] IResidue* find_residue(int seq_num) {
        for (auto& res : residues_) {
            if (res->seq_num() == seq_num) {
                return res.get();
            }
        }
        return nullptr;
    }

    [[nodiscard]] const IResidue* find_residue(int seq_num) const {
        for (const auto& res : residues_) {
            if (res->seq_num() == seq_num) {
                return res.get();
            }
        }
        return nullptr;
    }

private:
    std::string chain_id_;
    std::vector<std::unique_ptr<IResidue>> residues_;
};

} // namespace poly
} // namespace core
} // namespace x3dna

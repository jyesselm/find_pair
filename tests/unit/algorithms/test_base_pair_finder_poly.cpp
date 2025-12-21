/**
 * @file test_base_pair_finder_poly.cpp
 * @brief Unit tests for polymorphic BasePairFinder methods
 */

#include <gtest/gtest.h>
#include <x3dna/algorithms/base_pair_finder.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/core/structure/residue.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>
#include <filesystem>

using namespace x3dna::algorithms;
using namespace x3dna::geometry;

// Use structure namespace types explicitly to avoid collisions with legacy types
using x3dna::core::Atom;
using x3dna::core::ReferenceFrame;
using x3dna::core::structure::IResidue;
using x3dna::core::structure::INucleotide;
using x3dna::core::structure::RNA;
using x3dna::core::structure::DNA;
using x3dna::core::structure::Protein;
using x3dna::core::structure::Ligand;
using x3dna::core::structure::Chain;
using x3dna::core::structure::Structure;

namespace {

// Helper to create a standard adenine with ring atoms
std::unique_ptr<IResidue> create_adenine(const std::string& chain_id, int seq_num,
                                          const std::string& name = "A") {
    // RNA constructor: name, seq_num, chain_id, insertion
    auto rna = std::make_unique<RNA>(name, seq_num, chain_id);

    // Add ring atoms for adenine (purine)
    rna->add_atom(Atom("N9", Vector3D(0.0, 0.0, 0.0)));
    rna->add_atom(Atom("C8", Vector3D(1.2, 0.3, 0.0)));
    rna->add_atom(Atom("N7", Vector3D(2.1, -0.5, 0.0)));
    rna->add_atom(Atom("C5", Vector3D(1.5, -1.7, 0.0)));
    rna->add_atom(Atom("C6", Vector3D(2.0, -3.0, 0.0)));
    rna->add_atom(Atom("N6", Vector3D(3.3, -3.3, 0.0)));
    rna->add_atom(Atom("N1", Vector3D(1.1, -4.0, 0.0)));
    rna->add_atom(Atom("C2", Vector3D(-0.2, -3.7, 0.0)));
    rna->add_atom(Atom("N3", Vector3D(-0.7, -2.5, 0.0)));
    rna->add_atom(Atom("C4", Vector3D(0.2, -1.5, 0.0)));
    rna->add_atom(Atom("C1'", Vector3D(-1.0, 1.0, 0.0)));

    return rna;
}

// Helper to create a standard uracil with ring atoms
std::unique_ptr<IResidue> create_uracil(const std::string& chain_id, int seq_num,
                                         const std::string& name = "U") {
    // RNA constructor: name, seq_num, chain_id, insertion
    auto rna = std::make_unique<RNA>(name, seq_num, chain_id);

    // Add ring atoms for uracil (pyrimidine)
    rna->add_atom(Atom("N1", Vector3D(0.0, 0.0, 0.0)));
    rna->add_atom(Atom("C2", Vector3D(1.2, 0.3, 0.0)));
    rna->add_atom(Atom("O2", Vector3D(1.5, 1.5, 0.0)));
    rna->add_atom(Atom("N3", Vector3D(2.1, -0.7, 0.0)));
    rna->add_atom(Atom("C4", Vector3D(1.8, -2.0, 0.0)));
    rna->add_atom(Atom("O4", Vector3D(2.6, -2.9, 0.0)));
    rna->add_atom(Atom("C5", Vector3D(0.4, -2.3, 0.0)));
    rna->add_atom(Atom("C6", Vector3D(-0.5, -1.3, 0.0)));
    rna->add_atom(Atom("C1'", Vector3D(-1.0, 1.0, 0.0)));

    return rna;
}

// Create a simple reference frame
ReferenceFrame create_simple_frame(double x, double y, double z) {
    return ReferenceFrame(Matrix3D::identity(), Vector3D(x, y, z));
}

} // namespace

class BasePairFinderPolyTest : public ::testing::Test {
protected:
    void SetUp() override {
        finder_ = std::make_unique<BasePairFinder>();
    }

    std::unique_ptr<BasePairFinder> finder_;
};

TEST_F(BasePairFinderPolyTest, IsNucleotideForRNA) {
    auto adenine = create_adenine("A", 1);
    EXPECT_TRUE(BasePairFinder::is_nucleotide(*adenine));
}

TEST_F(BasePairFinderPolyTest, IsNucleotideForDNA) {
    // Constructor: name, seq_num, chain_id
    DNA dna("DT", 1, "A");
    EXPECT_TRUE(BasePairFinder::is_nucleotide(dna));
}

TEST_F(BasePairFinderPolyTest, IsNucleotideForProtein) {
    // Constructor: name, seq_num, chain_id
    Protein protein("ALA", 1, "A");
    EXPECT_FALSE(BasePairFinder::is_nucleotide(protein));
}

TEST_F(BasePairFinderPolyTest, IsNucleotideForLigand) {
    // Constructor: name, seq_num, chain_id
    Ligand ligand("HOH", 1, "A");
    EXPECT_FALSE(BasePairFinder::is_nucleotide(ligand));
}

TEST_F(BasePairFinderPolyTest, FindPairsEmptyStructure) {
    Structure structure("TEST");

    auto pairs = finder_->find_pairs(structure);
    EXPECT_TRUE(pairs.empty());
}

TEST_F(BasePairFinderPolyTest, FindPairsSingleChain) {
    Structure structure("TEST");
    Chain chain("A");

    auto adenine = create_adenine("A", 1);
    adenine->set_legacy_residue_idx(1);
    chain.add_residue(std::move(adenine));

    auto uracil = create_uracil("A", 2);
    uracil->set_legacy_residue_idx(2);
    chain.add_residue(std::move(uracil));

    structure.add_chain(std::move(chain));

    // Without frames, no pairs should be found
    auto pairs = finder_->find_pairs(structure);
    EXPECT_TRUE(pairs.empty());
}

TEST_F(BasePairFinderPolyTest, FindPairsWithFrames) {
    Structure structure("TEST");
    Chain chain("A");

    // Create adenine with frame
    auto adenine = create_adenine("A", 1);
    adenine->set_legacy_residue_idx(1);
    auto* nuc1 = dynamic_cast<INucleotide*>(adenine.get());
    if (nuc1) {
        nuc1->set_reference_frame(create_simple_frame(0, 0, 0));
    }
    chain.add_residue(std::move(adenine));

    // Create uracil with frame close to adenine
    auto uracil = create_uracil("A", 2);
    uracil->set_legacy_residue_idx(2);
    auto* nuc2 = dynamic_cast<INucleotide*>(uracil.get());
    if (nuc2) {
        nuc2->set_reference_frame(create_simple_frame(5, 0, 0));
    }
    chain.add_residue(std::move(uracil));

    structure.add_chain(std::move(chain));

    // With frames but no proper pair geometry, might not find pairs
    // This test verifies the method runs without crashing
    auto pairs = finder_->find_pairs(structure);
    // Pairs may or may not be found depending on geometry validation
}

TEST_F(BasePairFinderPolyTest, FindPairsAllPairsStrategy) {
    finder_->set_strategy(PairFindingStrategy::ALL_PAIRS);
    EXPECT_EQ(finder_->strategy(), PairFindingStrategy::ALL_PAIRS);

    Structure structure("TEST");
    Chain chain("A");

    auto adenine = create_adenine("A", 1);
    adenine->set_legacy_residue_idx(1);
    auto* nuc1 = dynamic_cast<INucleotide*>(adenine.get());
    if (nuc1) {
        nuc1->set_reference_frame(create_simple_frame(0, 0, 0));
    }
    chain.add_residue(std::move(adenine));

    auto uracil = create_uracil("A", 2);
    uracil->set_legacy_residue_idx(2);
    auto* nuc2 = dynamic_cast<INucleotide*>(uracil.get());
    if (nuc2) {
        nuc2->set_reference_frame(create_simple_frame(5, 0, 0));
    }
    chain.add_residue(std::move(uracil));

    structure.add_chain(std::move(chain));

    auto pairs = finder_->find_pairs(structure);
    // Verifies the ALL_PAIRS strategy runs
}

TEST_F(BasePairFinderPolyTest, FindPairsConstStructure) {
    Structure structure("TEST");
    Chain chain("A");

    auto adenine = create_adenine("A", 1);
    adenine->set_legacy_residue_idx(1);
    chain.add_residue(std::move(adenine));

    structure.add_chain(std::move(chain));

    // Test const overload
    const Structure& const_structure = structure;
    auto pairs = finder_->find_pairs(const_structure);
    EXPECT_TRUE(pairs.empty());
}

TEST_F(BasePairFinderPolyTest, FindPairsMultipleChains) {
    Structure structure("TEST");

    // Chain A with adenine
    Chain chain_a("A");
    auto adenine = create_adenine("A", 1);
    adenine->set_legacy_residue_idx(1);
    auto* nuc1 = dynamic_cast<INucleotide*>(adenine.get());
    if (nuc1) {
        nuc1->set_reference_frame(create_simple_frame(0, 0, 0));
    }
    chain_a.add_residue(std::move(adenine));
    structure.add_chain(std::move(chain_a));

    // Chain B with uracil
    Chain chain_b("B");
    auto uracil = create_uracil("B", 1);
    uracil->set_legacy_residue_idx(2);
    auto* nuc2 = dynamic_cast<INucleotide*>(uracil.get());
    if (nuc2) {
        nuc2->set_reference_frame(create_simple_frame(5, 0, 0));
    }
    chain_b.add_residue(std::move(uracil));
    structure.add_chain(std::move(chain_b));

    // Should search across chains
    auto pairs = finder_->find_pairs(structure);
    // Result depends on geometry validation
}

TEST_F(BasePairFinderPolyTest, ParameterAccess) {
    // Test that parameters can be set and retrieved
    auto params = ValidationParameters::defaults();
    finder_->set_parameters(params);

    const auto& retrieved_params = finder_->parameters();
    EXPECT_DOUBLE_EQ(retrieved_params.max_dorg, params.max_dorg);
}

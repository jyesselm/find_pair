/**
 * @file test_chain_detector.cpp
 * @brief Unit tests for ChainDetector class
 */

#include <gtest/gtest.h>
#include <x3dna/algorithms/chain_detector.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>

using namespace x3dna::algorithms;
using namespace x3dna::core;
using namespace x3dna::geometry;

class ChainDetectorTest : public ::testing::Test {
protected:
    void SetUp() override {
        detector_ = std::make_unique<ChainDetector>();
    }

    std::unique_ptr<ChainDetector> detector_;
};

// ============================================================================
// RNA Connectivity Tests
// ============================================================================

TEST_F(ChainDetectorTest, RNAResiduesConnected_Forward) {
    // Create two connected RNA residues: res1.O3' -> res2.P (5' to 3')
    std::vector<Atom> atoms1 = {Atom("O3'", Vector3D(0.0, 0.0, 0.0)), Atom("C1'", Vector3D(1.0, 0.0, 0.0))};
    auto res1 = Residue::create_from_atoms("G", 1, "A", "", atoms1);

    std::vector<Atom> atoms2 = {Atom("P", Vector3D(0.0, 0.0, 2.0)), // 2.0 Å from res1.O3'
                                Atom("C1'", Vector3D(2.0, 0.0, 2.0))};
    auto res2 = Residue::create_from_atoms("C", 2, "A", "", atoms2);

    int connection = detector_->are_rna_residues_connected(res1, res2);
    EXPECT_EQ(connection, 1); // Forward connection (5' to 3')
}

TEST_F(ChainDetectorTest, RNAResiduesConnected_Reverse) {
    // Create two connected RNA residues: res2.O3' -> res1.P (3' to 5')
    std::vector<Atom> atoms1 = {Atom("P", Vector3D(0.0, 0.0, 2.0)), Atom("C1'", Vector3D(1.0, 0.0, 0.0))};
    auto res1 = Residue::create_from_atoms("G", 1, "A", "", atoms1);

    std::vector<Atom> atoms2 = {Atom("O3'", Vector3D(0.0, 0.0, 0.0)), // 2.0 Å from res1.P
                                Atom("C1'", Vector3D(2.0, 0.0, 2.0))};
    auto res2 = Residue::create_from_atoms("C", 2, "A", "", atoms2);

    int connection = detector_->are_rna_residues_connected(res1, res2);
    EXPECT_EQ(connection, -1); // Reverse connection (3' to 5')
}

TEST_F(ChainDetectorTest, RNAResiduesConnected_NotConnected) {
    // Create two disconnected RNA residues (distance > cutoff)
    std::vector<Atom> atoms1 = {Atom("O3'", Vector3D(0.0, 0.0, 0.0)), Atom("C1'", Vector3D(1.0, 0.0, 0.0))};
    auto res1 = Residue::create_from_atoms("G", 1, "A", "", atoms1);

    std::vector<Atom> atoms2 = {Atom("P", Vector3D(0.0, 0.0, 10.0)), // 10.0 Å - too far
                                Atom("C1'", Vector3D(2.0, 0.0, 10.0))};
    auto res2 = Residue::create_from_atoms("C", 2, "A", "", atoms2);

    int connection = detector_->are_rna_residues_connected(res1, res2);
    EXPECT_EQ(connection, 0); // Not connected
}

TEST_F(ChainDetectorTest, RNAResiduesConnected_Triphosphate) {
    // Test triphosphate support (PA instead of P)
    std::vector<Atom> atoms1 = {Atom("O3'", Vector3D(0.0, 0.0, 0.0)), Atom("C1'", Vector3D(1.0, 0.0, 0.0))};
    auto res1 = Residue::create_from_atoms("G", 1, "A", "", atoms1);

    std::vector<Atom> atoms2 = {Atom("PA", Vector3D(0.0, 0.0, 2.0)), // Triphosphate
                                Atom("C1'", Vector3D(2.0, 0.0, 2.0))};
    auto res2 = Residue::create_from_atoms("A", 2, "A", "", atoms2);

    int connection = detector_->are_rna_residues_connected(res1, res2);
    EXPECT_EQ(connection, 1); // Forward connection using PA
}

// ============================================================================
// Protein Connectivity Tests
// ============================================================================

TEST_F(ChainDetectorTest, ProteinResiduesConnected_Forward) {
    // Create two connected protein residues: res1.C -> res2.N (N-term to C-term)
    std::vector<Atom> atoms1 = {Atom("C", Vector3D(0.0, 0.0, 0.0)), Atom("CA", Vector3D(1.0, 0.0, 0.0))};
    auto res1 = Residue::create_from_atoms("ALA", 1, "A", "", atoms1);

    std::vector<Atom> atoms2 = {Atom("N", Vector3D(0.0, 0.0, 1.5)), // 1.5 Å from res1.C (peptide bond)
                                Atom("CA", Vector3D(2.0, 0.0, 1.5))};
    auto res2 = Residue::create_from_atoms("GLY", 2, "A", "", atoms2);

    int connection = detector_->are_protein_residues_connected(res1, res2);
    EXPECT_EQ(connection, 1); // Forward connection (N to C)
}

TEST_F(ChainDetectorTest, ProteinResiduesConnected_Reverse) {
    // Create two connected protein residues: res2.C -> res1.N (C-term to N-term)
    std::vector<Atom> atoms1 = {Atom("N", Vector3D(0.0, 0.0, 1.5)), Atom("CA", Vector3D(1.0, 0.0, 0.0))};
    auto res1 = Residue::create_from_atoms("ALA", 1, "A", "", atoms1);

    std::vector<Atom> atoms2 = {Atom("C", Vector3D(0.0, 0.0, 0.0)), // 1.5 Å from res1.N
                                Atom("CA", Vector3D(2.0, 0.0, 1.5))};
    auto res2 = Residue::create_from_atoms("GLY", 2, "A", "", atoms2);

    int connection = detector_->are_protein_residues_connected(res1, res2);
    EXPECT_EQ(connection, -1); // Reverse connection (C to N)
}

TEST_F(ChainDetectorTest, ProteinResiduesConnected_NotConnected) {
    // Create two disconnected protein residues (distance > cutoff)
    std::vector<Atom> atoms1 = {Atom("C", Vector3D(0.0, 0.0, 0.0)), Atom("CA", Vector3D(1.0, 0.0, 0.0))};
    auto res1 = Residue::create_from_atoms("ALA", 1, "A", "", atoms1);

    std::vector<Atom> atoms2 = {Atom("N", Vector3D(0.0, 0.0, 5.0)), // 5.0 Å - too far
                                Atom("CA", Vector3D(2.0, 0.0, 5.0))};
    auto res2 = Residue::create_from_atoms("GLY", 2, "A", "", atoms2);

    int connection = detector_->are_protein_residues_connected(res1, res2);
    EXPECT_EQ(connection, 0); // Not connected
}

// ============================================================================
// Chain Detection Tests
// ============================================================================

TEST_F(ChainDetectorTest, DetectRNAChains_SimpleChain) {
    // Create a simple RNA chain: G1 - C2 - A3
    Structure structure;
    Chain chain("A");

    // Residue 1: G (has O3')
    std::vector<Atom> atoms1 = {Atom("P", Vector3D(0.0, 0.0, 0.0)), Atom("O3'", Vector3D(0.0, 0.0, 2.5)),
                                Atom("C1'", Vector3D(1.0, 0.0, 0.0))};
    auto res1 = Residue::create_from_atoms("G", 1, "A", "", atoms1);

    // Residue 2: C (has O3' and P connected to res1)
    std::vector<Atom> atoms2 = {Atom("P", Vector3D(0.0, 0.0, 2.5)), Atom("O3'", Vector3D(0.0, 0.0, 5.0)),
                                Atom("C1'", Vector3D(2.0, 0.0, 2.5))};
    auto res2 = Residue::create_from_atoms("C", 2, "A", "", atoms2);

    // Residue 3: A (has P connected to res2)
    std::vector<Atom> atoms3 = {Atom("P", Vector3D(0.0, 0.0, 5.0)), Atom("O3'", Vector3D(0.0, 0.0, 7.5)),
                                Atom("C1'", Vector3D(3.0, 0.0, 5.0))};
    auto res3 = Residue::create_from_atoms("A", 3, "A", "", atoms3);

    chain.add_residue(res1);
    chain.add_residue(res2);
    chain.add_residue(res3);
    structure.add_chain(chain);

    auto chains = detector_->detect_rna_chains(structure);

    EXPECT_EQ(chains.size(), 1);
    EXPECT_EQ(chains[0].residues.size(), 3);
    EXPECT_TRUE(chains[0].is_rna);
    EXPECT_FALSE(chains[0].is_protein);
    EXPECT_EQ(chains[0].chain_id, "A");
}

TEST_F(ChainDetectorTest, DetectRNAChains_MultipleChains) {
    // Create two separate RNA chains in same PDB chain
    Structure structure;
    Chain chain("A");

    // Chain 1: G1 - C2 (connected)
    std::vector<Atom> atoms1 = {Atom("O3'", Vector3D(0.0, 0.0, 0.0)), Atom("C1'", Vector3D(1.0, 0.0, 0.0))};
    auto res1 = Residue::create_from_atoms("G", 1, "A", "", atoms1);

    std::vector<Atom> atoms2 = {Atom("P", Vector3D(0.0, 0.0, 0.0)), Atom("O3'", Vector3D(0.0, 0.0, 2.5)),
                                Atom("C1'", Vector3D(2.0, 0.0, 0.0))};
    auto res2 = Residue::create_from_atoms("C", 2, "A", "", atoms2);

    // Chain 2: A10 - U11 (connected, but separate from chain 1)
    std::vector<Atom> atoms10 = {Atom("O3'", Vector3D(10.0, 0.0, 0.0)), Atom("C1'", Vector3D(11.0, 0.0, 0.0))};
    auto res10 = Residue::create_from_atoms("A", 10, "A", "", atoms10);

    std::vector<Atom> atoms11 = {Atom("P", Vector3D(10.0, 0.0, 0.0)), Atom("C1'", Vector3D(12.0, 0.0, 0.0))};
    auto res11 = Residue::create_from_atoms("U", 11, "A", "", atoms11);

    chain.add_residue(res1);
    chain.add_residue(res2);
    chain.add_residue(res10);
    chain.add_residue(res11);
    structure.add_chain(chain);

    auto chains = detector_->detect_rna_chains(structure);

    EXPECT_EQ(chains.size(), 2);
    EXPECT_EQ(chains[0].residues.size(), 2);
    EXPECT_EQ(chains[1].residues.size(), 2);
}

TEST_F(ChainDetectorTest, DetectProteinChains_SimpleChain) {
    // Create a simple protein chain: ALA1 - GLY2
    Structure structure;
    Chain chain("A");

    std::vector<Atom> atoms1 = {Atom("N", Vector3D(0.0, 0.0, 0.0)), Atom("CA", Vector3D(1.0, 0.0, 0.0)),
                                Atom("C", Vector3D(2.0, 0.0, 0.0))};
    auto res1 = Residue::create_from_atoms("ALA", 1, "A", "", atoms1);

    std::vector<Atom> atoms2 = {Atom("N", Vector3D(2.0, 0.0, 0.0)), // Connected to res1.C
                                Atom("CA", Vector3D(3.0, 0.0, 0.0)), Atom("C", Vector3D(4.0, 0.0, 0.0))};
    auto res2 = Residue::create_from_atoms("GLY", 2, "A", "", atoms2);

    chain.add_residue(res1);
    chain.add_residue(res2);
    structure.add_chain(chain);

    auto chains = detector_->detect_protein_chains(structure);

    EXPECT_EQ(chains.size(), 1);
    EXPECT_EQ(chains[0].residues.size(), 2);
    EXPECT_FALSE(chains[0].is_rna);
    EXPECT_TRUE(chains[0].is_protein);
    EXPECT_EQ(chains[0].chain_id, "A");
}

TEST_F(ChainDetectorTest, DetectAllChains_Mixed) {
    // Create structure with both RNA and protein
    Structure structure;

    // RNA chain
    Chain rna_chain("R");
    std::vector<Atom> rna_atoms = {Atom("P", Vector3D(0.0, 0.0, 0.0)), Atom("C1'", Vector3D(1.0, 0.0, 0.0))};
    auto rna_res = Residue::create_from_atoms("G", 1, "R", "", rna_atoms);
    rna_chain.add_residue(rna_res);
    structure.add_chain(rna_chain);

    // Protein chain
    Chain protein_chain("P");
    std::vector<Atom> protein_atoms = {Atom("N", Vector3D(0.0, 0.0, 0.0)), Atom("CA", Vector3D(1.0, 0.0, 0.0)),
                                       Atom("C", Vector3D(2.0, 0.0, 0.0))};
    auto protein_res = Residue::create_from_atoms("ALA", 1, "P", "", protein_atoms);
    protein_chain.add_residue(protein_res);
    structure.add_chain(protein_chain);

    auto chains = detector_->detect_all_chains(structure);

    EXPECT_EQ(chains.size(), 2);

    // Find RNA and protein chains
    bool found_rna = false;
    bool found_protein = false;
    for (const auto& chain : chains) {
        if (chain.is_rna)
            found_rna = true;
        if (chain.is_protein)
            found_protein = true;
    }

    EXPECT_TRUE(found_rna);
    EXPECT_TRUE(found_protein);
}

// ============================================================================
// Configuration Tests
// ============================================================================

TEST_F(ChainDetectorTest, CustomConfiguration) {
    ChainDetector::Config config;
    config.rna_connectivity_cutoff = 3.0; // Relaxed cutoff
    config.protein_connectivity_cutoff = 2.5;
    config.merge_adjacent_chains = false;

    ChainDetector custom_detector(config);

    // Test with distance that would fail with default cutoff but passes with custom
    std::vector<Atom> atoms1 = {Atom("O3'", Vector3D(0.0, 0.0, 0.0)), Atom("C1'", Vector3D(1.0, 0.0, 0.0))};
    auto res1 = Residue::create_from_atoms("G", 1, "A", "", atoms1);

    std::vector<Atom> atoms2 = {Atom("P", Vector3D(0.0, 0.0, 2.9)), // 2.9 Å - fails default (2.75), passes custom (3.0)
                                Atom("C1'", Vector3D(2.0, 0.0, 2.9))};
    auto res2 = Residue::create_from_atoms("C", 2, "A", "", atoms2);

    int connection_default = detector_->are_rna_residues_connected(res1, res2);
    int connection_custom = custom_detector.are_rna_residues_connected(res1, res2);

    EXPECT_EQ(connection_default, 0); // Not connected with default cutoff
    EXPECT_EQ(connection_custom, 1);  // Connected with custom cutoff
}

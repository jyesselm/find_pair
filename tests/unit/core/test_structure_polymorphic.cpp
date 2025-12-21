/**
 * @file test_structure_polymorphic.cpp
 * @brief Unit tests for the polymorphic Structure class
 */

#include <gtest/gtest.h>
#include <x3dna/core/structure/residue.hpp>
#include <x3dna/geometry/vector3d.hpp>

using namespace x3dna::core::structure;
using namespace x3dna::core;  // For Atom
using namespace x3dna::geometry;

class StructurePolymorphicTest : public ::testing::Test {
protected:
    void SetUp() override {}

    // Helper to create a simple chain with residues
    Chain make_chain(const std::string& id, const std::vector<std::string>& residue_names) {
        Chain chain(id);
        int seq = 1;
        for (const auto& name : residue_names) {
            chain.add_residue(ResidueFactory::create(name, seq++, id, ""));
        }
        return chain;
    }
};

// === Basic Construction ===

TEST_F(StructurePolymorphicTest, DefaultConstruction) {
    Structure structure;
    EXPECT_TRUE(structure.empty());
    EXPECT_EQ(structure.size(), 0);
    EXPECT_EQ(structure.pdb_id(), "");
}

TEST_F(StructurePolymorphicTest, ConstructWithPdbId) {
    Structure structure("1ABC");
    EXPECT_EQ(structure.pdb_id(), "1ABC");
    EXPECT_TRUE(structure.empty());
}

// === Adding Chains ===

TEST_F(StructurePolymorphicTest, AddChain) {
    Structure structure("1ABC");

    structure.add_chain(make_chain("A", {"A", "G", "C", "U"}));

    EXPECT_EQ(structure.size(), 1);
    EXPECT_EQ(structure[0].chain_id(), "A");
    EXPECT_EQ(structure[0].size(), 4);
}

TEST_F(StructurePolymorphicTest, AddMultipleChains) {
    Structure structure("1ABC");

    structure.add_chain(make_chain("A", {"A", "G"}));
    structure.add_chain(make_chain("B", {"C", "U"}));

    EXPECT_EQ(structure.size(), 2);
    EXPECT_EQ(structure[0].chain_id(), "A");
    EXPECT_EQ(structure[1].chain_id(), "B");
}

// === Counts ===

TEST_F(StructurePolymorphicTest, NumResidues) {
    Structure structure("1ABC");
    structure.add_chain(make_chain("A", {"A", "G", "C"}));
    structure.add_chain(make_chain("B", {"U", "A"}));

    EXPECT_EQ(structure.num_residues(), 5);
}

TEST_F(StructurePolymorphicTest, NumAtoms) {
    Structure structure("1ABC");

    Chain chain("A");
    auto res = ResidueFactory::create("A", 1, "A", "");
    res->add_atom(Atom("N9", Vector3D(0, 0, 0)));
    res->add_atom(Atom("C8", Vector3D(1, 0, 0)));
    chain.add_residue(std::move(res));
    structure.add_chain(std::move(chain));

    EXPECT_EQ(structure.num_atoms(), 2);
}

// === Residue Access ===

TEST_F(StructurePolymorphicTest, AllResidues) {
    Structure structure("1ABC");
    structure.add_chain(make_chain("A", {"A", "G"}));
    structure.add_chain(make_chain("B", {"C"}));

    auto residues = structure.all_residues();
    EXPECT_EQ(residues.size(), 3);
    EXPECT_EQ(residues[0]->name(), "A");
    EXPECT_EQ(residues[1]->name(), "G");
    EXPECT_EQ(residues[2]->name(), "C");
}

TEST_F(StructurePolymorphicTest, Nucleotides) {
    Structure structure("1ABC");

    Chain chain("A");
    chain.add_residue(ResidueFactory::create("A", 1, "A", ""));     // RNA
    chain.add_residue(ResidueFactory::create("ALA", 2, "A", ""));   // Protein
    chain.add_residue(ResidueFactory::create("G", 3, "A", ""));     // RNA
    structure.add_chain(std::move(chain));

    auto nucs = structure.nucleotides();
    EXPECT_EQ(nucs.size(), 2);
    EXPECT_EQ(nucs[0]->one_letter_code(), 'A');
    EXPECT_EQ(nucs[1]->one_letter_code(), 'G');
}

// === Find Chain ===

TEST_F(StructurePolymorphicTest, FindChain) {
    Structure structure("1ABC");
    structure.add_chain(make_chain("A", {"A", "G"}));
    structure.add_chain(make_chain("B", {"C", "U"}));

    auto* found = structure.find_chain("B");
    ASSERT_NE(found, nullptr);
    EXPECT_EQ(found->chain_id(), "B");
    EXPECT_EQ(found->size(), 2);

    auto* not_found = structure.find_chain("Z");
    EXPECT_EQ(not_found, nullptr);
}

// === Clone ===

TEST_F(StructurePolymorphicTest, CloneStructure) {
    Structure original("1ABC");
    original.add_chain(make_chain("A", {"A", "G"}));

    Structure cloned = original.clone();

    EXPECT_EQ(cloned.pdb_id(), "1ABC");
    EXPECT_EQ(cloned.size(), 1);
    EXPECT_EQ(cloned[0].size(), 2);

    // Verify deep copy
    EXPECT_NE(&cloned[0], &original[0]);
}

// === Move Semantics ===

TEST_F(StructurePolymorphicTest, MoveConstruction) {
    Structure original("1ABC");
    original.add_chain(make_chain("A", {"A", "G"}));

    Structure moved(std::move(original));

    EXPECT_EQ(moved.pdb_id(), "1ABC");
    EXPECT_EQ(moved.size(), 1);
}

// === Iteration ===

TEST_F(StructurePolymorphicTest, IterateChains) {
    Structure structure("1ABC");
    structure.add_chain(make_chain("A", {"A"}));
    structure.add_chain(make_chain("B", {"G"}));

    int count = 0;
    for (const auto& chain : structure) {
        EXPECT_FALSE(chain.empty());
        count++;
    }
    EXPECT_EQ(count, 2);
}

// === Legacy Index Support ===

TEST_F(StructurePolymorphicTest, SetLegacyIndices) {
    Structure structure("1ABC");

    Chain chain("A");
    auto res = ResidueFactory::create("A", 1, "A", "");
    res->add_atom(Atom("N9", Vector3D(0, 0, 0)));
    chain.add_residue(std::move(res));
    structure.add_chain(std::move(chain));

    // Set legacy indices
    std::map<std::tuple<std::string, int, std::string>, int> residue_idx_map;
    residue_idx_map[std::make_tuple("A", 1, "")] = 42;

    std::map<std::tuple<std::string, int, std::string, std::string>, int> atom_idx_map;
    atom_idx_map[std::make_tuple("A", 1, "", "N9")] = 100;

    structure.set_legacy_indices(atom_idx_map, residue_idx_map);

    EXPECT_EQ(structure[0][0].legacy_residue_idx(), 42);
    EXPECT_EQ(structure[0][0].atoms()[0].legacy_atom_idx(), 100);
}

TEST_F(StructurePolymorphicTest, GetResidueByLegacyIdx) {
    Structure structure("1ABC");

    Chain chain("A");
    auto res1 = ResidueFactory::create("A", 1, "A", "");
    auto res2 = ResidueFactory::create("G", 2, "A", "");
    chain.add_residue(std::move(res1));
    chain.add_residue(std::move(res2));
    structure.add_chain(std::move(chain));

    // Set legacy indices
    std::map<std::tuple<std::string, int, std::string>, int> residue_idx_map;
    residue_idx_map[std::make_tuple("A", 1, "")] = 10;
    residue_idx_map[std::make_tuple("A", 2, "")] = 20;
    std::map<std::tuple<std::string, int, std::string, std::string>, int> atom_idx_map;

    structure.set_legacy_indices(atom_idx_map, residue_idx_map);

    auto* found = structure.get_residue_by_legacy_idx(20);
    ASSERT_NE(found, nullptr);
    EXPECT_EQ(found->name(), "G");

    auto* not_found = structure.get_residue_by_legacy_idx(999);
    EXPECT_EQ(not_found, nullptr);
}

// === Record Type Support ===

TEST_F(StructurePolymorphicTest, RecordType) {
    Structure structure("1ABC");

    structure.set_residue_record_type("A", 1, "", 'H');

    EXPECT_EQ(structure.get_residue_record_type("A", 1, ""), 'H');
    EXPECT_EQ(structure.get_residue_record_type("A", 2, ""), 'A'); // Default
}

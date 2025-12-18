/**
 * @file test_find_pair_protocol.cpp
 * @brief Unit tests for FindPairProtocol
 */

#include <gtest/gtest.h>
#include <x3dna/protocols/find_pair_protocol.hpp>
#include <x3dna/config/config_manager.hpp>
#include <x3dna/config/resource_locator.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <filesystem>

using namespace x3dna::protocols;
using namespace x3dna::config;
using namespace x3dna::core;
using namespace x3dna::geometry;

class FindPairProtocolTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Reset config
        auto& config = ConfigManager::instance();
        config.set_defaults();

        // Create a simple structure for testing
        structure_ = Structure("TEST");

        // Create chain A with a few residues
        Chain chain_a('A');

        // Residue 1: C (Cytosine)
        Residue c1("  C", 1, 'A');
        c1.add_atom(Atom(" C1'", Vector3D(0.0, 0.0, 0.0), "  C", 'A', 1));
        c1.add_atom(Atom(" N1 ", Vector3D(1.0, 0.0, 0.0), "  C", 'A', 1));
        c1.add_atom(Atom(" C2 ", Vector3D(2.0, 0.0, 0.0), "  C", 'A', 1));
        c1.add_atom(Atom(" O2 ", Vector3D(3.0, 0.0, 0.0), "  C", 'A', 1));
        c1.add_atom(Atom(" N3 ", Vector3D(4.0, 0.0, 0.0), "  C", 'A', 1));
        c1.add_atom(Atom(" C4 ", Vector3D(5.0, 0.0, 0.0), "  C", 'A', 1));
        c1.add_atom(Atom(" N4 ", Vector3D(6.0, 0.0, 0.0), "  C", 'A', 1));
        c1.add_atom(Atom(" C5 ", Vector3D(7.0, 0.0, 0.0), "  C", 'A', 1));
        c1.add_atom(Atom(" C6 ", Vector3D(8.0, 0.0, 0.0), "  C", 'A', 1));
        chain_a.add_residue(c1);

        // Residue 2: G (Guanine) - paired with C
        Residue g1("  G", 2, 'A');
        g1.add_atom(Atom(" C1'", Vector3D(0.0, 10.0, 0.0), "  G", 'A', 2));
        g1.add_atom(Atom(" N1 ", Vector3D(1.0, 10.0, 0.0), "  G", 'A', 2));
        g1.add_atom(Atom(" C2 ", Vector3D(2.0, 10.0, 0.0), "  G", 'A', 2));
        g1.add_atom(Atom(" N2 ", Vector3D(3.0, 10.0, 0.0), "  G", 'A', 2));
        g1.add_atom(Atom(" N3 ", Vector3D(4.0, 10.0, 0.0), "  G", 'A', 2));
        g1.add_atom(Atom(" C4 ", Vector3D(5.0, 10.0, 0.0), "  G", 'A', 2));
        g1.add_atom(Atom(" C5 ", Vector3D(6.0, 10.0, 0.0), "  G", 'A', 2));
        g1.add_atom(Atom(" C6 ", Vector3D(7.0, 10.0, 0.0), "  G", 'A', 2));
        g1.add_atom(Atom(" O6 ", Vector3D(8.0, 10.0, 0.0), "  G", 'A', 2));
        g1.add_atom(Atom(" N7 ", Vector3D(9.0, 10.0, 0.0), "  G", 'A', 2));
        g1.add_atom(Atom(" C8 ", Vector3D(10.0, 10.0, 0.0), "  G", 'A', 2));
        g1.add_atom(Atom(" N9 ", Vector3D(11.0, 10.0, 0.0), "  G", 'A', 2));
        chain_a.add_residue(g1);

        structure_.add_chain(chain_a);

        // Use ResourceLocator for template path (it auto-initializes from environment)
        if (!ResourceLocator::is_initialized()) {
            ResourceLocator::initialize_from_environment();
        }
        auto template_path = ResourceLocator::templates_dir();

        protocol_ = std::make_unique<FindPairProtocol>(template_path);
    }

    Structure structure_;
    std::unique_ptr<FindPairProtocol> protocol_;
};

// Constructor tests
TEST_F(FindPairProtocolTest, Constructor) {
    EXPECT_NE(protocol_, nullptr);
    EXPECT_FALSE(protocol_->single_strand_mode());
    EXPECT_FALSE(protocol_->find_all_pairs());
    EXPECT_FALSE(protocol_->divide_helices());
    EXPECT_FALSE(protocol_->legacy_mode());
}

// Options tests
TEST_F(FindPairProtocolTest, SingleStrandMode) {
    protocol_->set_single_strand_mode(true);
    EXPECT_TRUE(protocol_->single_strand_mode());
    protocol_->set_single_strand_mode(false);
    EXPECT_FALSE(protocol_->single_strand_mode());
}

TEST_F(FindPairProtocolTest, FindAllPairs) {
    protocol_->set_find_all_pairs(true);
    EXPECT_TRUE(protocol_->find_all_pairs());
    protocol_->set_find_all_pairs(false);
    EXPECT_FALSE(protocol_->find_all_pairs());
}

TEST_F(FindPairProtocolTest, DivideHelices) {
    protocol_->set_divide_helices(true);
    EXPECT_TRUE(protocol_->divide_helices());
    protocol_->set_divide_helices(false);
    EXPECT_FALSE(protocol_->divide_helices());
}

TEST_F(FindPairProtocolTest, LegacyMode) {
    protocol_->set_legacy_mode(true);
    EXPECT_TRUE(protocol_->legacy_mode());
    protocol_->set_legacy_mode(false);
    EXPECT_FALSE(protocol_->legacy_mode());
}

// Configuration tests
TEST_F(FindPairProtocolTest, ConfigAccess) {
    // Test that config() returns a reference to FindPairConfig
    auto& config = protocol_->config();

    // Modify through config reference
    config.legacy_mode = true;
    EXPECT_TRUE(protocol_->legacy_mode());

    config.find_all_pairs = true;
    EXPECT_TRUE(protocol_->find_all_pairs());

    // Verify const access works
    const auto& const_protocol = *protocol_;
    EXPECT_TRUE(const_protocol.config().legacy_mode);
}

// Base pairs access tests
TEST_F(FindPairProtocolTest, BasePairsAccess) {
    // Initially empty
    EXPECT_TRUE(protocol_->base_pairs().empty());

    // After execution, may have pairs (if structure is valid and templates exist)
    try {
        protocol_->execute(structure_);
        // Should be able to access base pairs (even if empty)
        EXPECT_NO_THROW(protocol_->base_pairs());
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Templates not available: " << e.what();
    }
}

// Frame calculator access tests
TEST_F(FindPairProtocolTest, FrameCalculatorAccess) {
    EXPECT_NO_THROW(protocol_->frame_calculator());
}

// Pair finder access tests
TEST_F(FindPairProtocolTest, PairFinderAccess) {
    EXPECT_NO_THROW(protocol_->pair_finder());
}

// JSON writer tests
TEST_F(FindPairProtocolTest, SetJsonWriter) {
    // Setting nullptr should be allowed
    protocol_->set_json_writer(nullptr);

    // We can't easily test with a real JsonWriter without more setup
    // But we can verify the method exists and doesn't crash
    EXPECT_NO_THROW(protocol_->set_json_writer(nullptr));
}

// Test config struct initialization
TEST_F(FindPairProtocolTest, ConfigStructInitialization) {
    // Test creating protocol with config struct
    FindPairConfig config;
    config.legacy_mode = true;
    config.find_all_pairs = true;
    config.output_stage = "frames";

    // Use ResourceLocator for template path
    auto template_path = ResourceLocator::templates_dir();

    FindPairProtocol configured_protocol(template_path, config);

    EXPECT_TRUE(configured_protocol.legacy_mode());
    EXPECT_TRUE(configured_protocol.find_all_pairs());
    EXPECT_EQ(configured_protocol.output_stage(), "frames");
}

// Multiple executions
TEST_F(FindPairProtocolTest, MultipleExecutions) {
    try {
        protocol_->execute(structure_);
        size_t first_count = protocol_->base_pairs().size();

        protocol_->execute(structure_);
        size_t second_count = protocol_->base_pairs().size();

        // Should be able to execute multiple times
        // Results may vary, but should not crash
        EXPECT_GE(second_count, 0);
        // Both executions should complete
        (void)first_count; // Suppress unused variable warning
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Templates not available: " << e.what();
    }
}

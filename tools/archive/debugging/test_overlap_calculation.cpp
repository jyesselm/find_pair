/**
 * Test executable to test overlap calculation with various inputs
 * This helps verify the modern implementation matches legacy behavior
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/geometry/vector3d.hpp>

using namespace x3dna::core;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;

// Test with simple rectangular polygons
void test_simple_overlap() {
    std::cout << "=== Test 1: Simple Rectangular Overlap ===\n";

    // Create two residues with rectangular base atoms
    Residue res1("  A", 1, 'A');
    Residue res2("  A", 2, 'A');

    // Residue 1: rectangle from (0,0,0) to (10,10,0)
    res1.add_atom(Atom(" C4 ", Vector3D(0, 0, 0)));
    res1.add_atom(Atom(" N3 ", Vector3D(10, 0, 0)));
    res1.add_atom(Atom(" C2 ", Vector3D(10, 10, 0)));
    res1.add_atom(Atom(" N1 ", Vector3D(0, 10, 0)));

    // Residue 2: rectangle from (5,0,0) to (15,10,0) - overlaps by 5x10
    res2.add_atom(Atom(" C4 ", Vector3D(5, 0, 0)));
    res2.add_atom(Atom(" N3 ", Vector3D(15, 0, 0)));
    res2.add_atom(Atom(" C2 ", Vector3D(15, 10, 0)));
    res2.add_atom(Atom(" N1 ", Vector3D(5, 10, 0)));

    Vector3D oave(0, 0, 0);
    Vector3D zave(0, 0, 1); // z-axis

    BasePairValidator validator;
    double overlap = validator.calculate_overlap_area(res1, res2, oave, zave);

    std::cout << "  Expected overlap: ~50.0 (5x10 rectangle)\n";
    std::cout << "  Calculated overlap: " << std::fixed << std::setprecision(6) << overlap << "\n\n";
}

// Test with non-overlapping polygons
void test_no_overlap() {
    std::cout << "=== Test 2: Non-Overlapping Polygons ===\n";

    Residue res1("  A", 1, 'A');
    Residue res2("  A", 2, 'A');

    // Residue 1: rectangle at (0,0)
    res1.add_atom(Atom(" C4 ", Vector3D(0, 0, 0)));
    res1.add_atom(Atom(" N3 ", Vector3D(5, 0, 0)));
    res1.add_atom(Atom(" C2 ", Vector3D(5, 5, 0)));
    res1.add_atom(Atom(" N1 ", Vector3D(0, 5, 0)));

    // Residue 2: rectangle at (10,10) - no overlap
    res2.add_atom(Atom(" C4 ", Vector3D(10, 10, 0)));
    res2.add_atom(Atom(" N3 ", Vector3D(15, 10, 0)));
    res2.add_atom(Atom(" C2 ", Vector3D(15, 15, 0)));
    res2.add_atom(Atom(" N1 ", Vector3D(10, 15, 0)));

    Vector3D oave(0, 0, 0);
    Vector3D zave(0, 0, 1);

    BasePairValidator validator;
    double overlap = validator.calculate_overlap_area(res1, res2, oave, zave);

    std::cout << "  Expected overlap: 0.0\n";
    std::cout << "  Calculated overlap: " << std::fixed << std::setprecision(6) << overlap << "\n\n";
}

// Test with one polygon inside another
void test_inside_overlap() {
    std::cout << "=== Test 3: One Polygon Inside Another ===\n";

    Residue res1("  A", 1, 'A');
    Residue res2("  A", 2, 'A');

    // Residue 1: large rectangle (0,0) to (10,10)
    res1.add_atom(Atom(" C4 ", Vector3D(0, 0, 0)));
    res1.add_atom(Atom(" N3 ", Vector3D(10, 0, 0)));
    res1.add_atom(Atom(" C2 ", Vector3D(10, 10, 0)));
    res1.add_atom(Atom(" N1 ", Vector3D(0, 10, 0)));

    // Residue 2: small rectangle (2,2) to (8,8) - completely inside
    res2.add_atom(Atom(" C4 ", Vector3D(2, 2, 0)));
    res2.add_atom(Atom(" N3 ", Vector3D(8, 2, 0)));
    res2.add_atom(Atom(" C2 ", Vector3D(8, 8, 0)));
    res2.add_atom(Atom(" N1 ", Vector3D(2, 8, 0)));

    Vector3D oave(0, 0, 0);
    Vector3D zave(0, 0, 1);

    BasePairValidator validator;
    double overlap = validator.calculate_overlap_area(res1, res2, oave, zave);

    std::cout << "  Expected overlap: ~36.0 (6x6 rectangle)\n";
    std::cout << "  Calculated overlap: " << std::fixed << std::setprecision(6) << overlap << "\n\n";
}

// Test with complex polygon shapes
void test_complex_shapes() {
    std::cout << "=== Test 4: Complex Polygon Shapes ===\n";

    Residue res1("  A", 1, 'A');
    Residue res2("  A", 2, 'A');

    // Test 4a: Irregular pentagon vs triangle
    std::cout << "Test 4a: Irregular pentagon vs triangle\n";
    res1 = Residue("  A", 1, 'A');
    res2 = Residue("  A", 2, 'A');

    // Pentagon: (0,0), (5,0), (6,3), (3,5), (0,3)
    res1.add_atom(Atom(" C4 ", Vector3D(0, 0, 0)));
    res1.add_atom(Atom(" N3 ", Vector3D(5, 0, 0)));
    res1.add_atom(Atom(" C2 ", Vector3D(6, 3, 0)));
    res1.add_atom(Atom(" N1 ", Vector3D(3, 5, 0)));
    res1.add_atom(Atom(" C6 ", Vector3D(0, 3, 0)));

    // Triangle: (2,1), (4,1), (3,3) - overlaps with pentagon
    res2.add_atom(Atom(" C4 ", Vector3D(2, 1, 0)));
    res2.add_atom(Atom(" N3 ", Vector3D(4, 1, 0)));
    res2.add_atom(Atom(" C2 ", Vector3D(3, 3, 0)));

    Vector3D oave(0, 0, 0);
    Vector3D zave(0, 0, 1);

    BasePairValidator validator;
    double overlap = validator.calculate_overlap_area(res1, res2, oave, zave);

    std::cout << "  Expected overlap: > 0.0 (triangle inside pentagon)\n";
    std::cout << "  Calculated overlap: " << std::fixed << std::setprecision(6) << overlap << "\n\n";

    // Test 4b: Star-like shape vs rectangle
    std::cout << "Test 4b: Star-like shape vs rectangle\n";
    res1 = Residue("  A", 1, 'A');
    res2 = Residue("  A", 2, 'A');

    // Star-like shape: (5,5), (6,2), (9,5), (6,8), (5,5) - center at (6,5)
    res1.add_atom(Atom(" C4 ", Vector3D(5, 5, 0)));
    res1.add_atom(Atom(" N3 ", Vector3D(6, 2, 0)));
    res1.add_atom(Atom(" C2 ", Vector3D(9, 5, 0)));
    res1.add_atom(Atom(" N1 ", Vector3D(6, 8, 0)));
    res1.add_atom(Atom(" C6 ", Vector3D(5, 5, 0))); // Close polygon

    // Rectangle: (4,4) to (8,6) - overlaps with star
    res2.add_atom(Atom(" C4 ", Vector3D(4, 4, 0)));
    res2.add_atom(Atom(" N3 ", Vector3D(8, 4, 0)));
    res2.add_atom(Atom(" C2 ", Vector3D(8, 6, 0)));
    res2.add_atom(Atom(" N1 ", Vector3D(4, 6, 0)));

    overlap = validator.calculate_overlap_area(res1, res2, oave, zave);

    std::cout << "  Expected overlap: > 0.0 (partial overlap)\n";
    std::cout << "  Calculated overlap: " << std::fixed << std::setprecision(6) << overlap << "\n\n";

    // Test 4c: Concave polygon vs convex polygon
    std::cout << "Test 4c: Concave (L-shape) vs convex polygon\n";
    res1 = Residue("  A", 1, 'A');
    res2 = Residue("  A", 2, 'A');

    // L-shape: (0,0), (4,0), (4,2), (2,2), (2,4), (0,4)
    res1.add_atom(Atom(" C4 ", Vector3D(0, 0, 0)));
    res1.add_atom(Atom(" N3 ", Vector3D(4, 0, 0)));
    res1.add_atom(Atom(" C2 ", Vector3D(4, 2, 0)));
    res1.add_atom(Atom(" N1 ", Vector3D(2, 2, 0)));
    res1.add_atom(Atom(" C6 ", Vector3D(2, 4, 0)));
    res1.add_atom(Atom(" C5 ", Vector3D(0, 4, 0)));

    // Triangle: (1,1), (3,1), (2,3) - overlaps with L-shape
    res2.add_atom(Atom(" C4 ", Vector3D(1, 1, 0)));
    res2.add_atom(Atom(" N3 ", Vector3D(3, 1, 0)));
    res2.add_atom(Atom(" C2 ", Vector3D(2, 3, 0)));

    overlap = validator.calculate_overlap_area(res1, res2, oave, zave);

    std::cout << "  Expected overlap: > 0.0 (triangle overlaps L-shape)\n";
    std::cout << "  Calculated overlap: " << std::fixed << std::setprecision(6) << overlap << "\n\n";

    // Test 4d: Rotated shapes (test projection)
    std::cout << "Test 4d: Rotated shapes (test z-axis projection)\n";
    res1 = Residue("  A", 1, 'A');
    res2 = Residue("  A", 2, 'A');

    // Square in XY plane: (0,0,0), (2,0,0), (2,2,0), (0,2,0)
    res1.add_atom(Atom(" C4 ", Vector3D(0, 0, 0)));
    res1.add_atom(Atom(" N3 ", Vector3D(2, 0, 0)));
    res1.add_atom(Atom(" C2 ", Vector3D(2, 2, 0)));
    res1.add_atom(Atom(" N1 ", Vector3D(0, 2, 0)));

    // Square rotated 45 degrees in XY plane: (1,0,0), (2,1,0), (1,2,0), (0,1,0)
    res2.add_atom(Atom(" C4 ", Vector3D(1, 0, 0)));
    res2.add_atom(Atom(" N3 ", Vector3D(2, 1, 0)));
    res2.add_atom(Atom(" C2 ", Vector3D(1, 2, 0)));
    res2.add_atom(Atom(" N1 ", Vector3D(0, 1, 0)));

    // Use z-axis pointing in different direction to test projection
    zave = Vector3D(0, 0, 1); // Normal z-axis
    overlap = validator.calculate_overlap_area(res1, res2, oave, zave);

    std::cout << "  Expected overlap: > 0.0 (rotated squares overlap)\n";
    std::cout << "  Calculated overlap: " << std::fixed << std::setprecision(6) << overlap << "\n\n";

    // Test 4e: Very small overlap (near threshold)
    std::cout << "Test 4e: Very small overlap (near 0.01 threshold)\n";
    res1 = Residue("  A", 1, 'A');
    res2 = Residue("  A", 2, 'A');

    // Small square: (0,0), (0.1,0), (0.1,0.1), (0,0.1)
    res1.add_atom(Atom(" C4 ", Vector3D(0, 0, 0)));
    res1.add_atom(Atom(" N3 ", Vector3D(0.1, 0, 0)));
    res1.add_atom(Atom(" C2 ", Vector3D(0.1, 0.1, 0)));
    res1.add_atom(Atom(" N1 ", Vector3D(0, 0.1, 0)));

    // Overlapping small square: (0.05,0.05), (0.15,0.05), (0.15,0.15), (0.05,0.15)
    res2.add_atom(Atom(" C4 ", Vector3D(0.05, 0.05, 0)));
    res2.add_atom(Atom(" N3 ", Vector3D(0.15, 0.05, 0)));
    res2.add_atom(Atom(" C2 ", Vector3D(0.15, 0.15, 0)));
    res2.add_atom(Atom(" N1 ", Vector3D(0.05, 0.15, 0)));

    overlap = validator.calculate_overlap_area(res1, res2, oave, zave);

    std::cout << "  Expected overlap: ~0.0025 (0.05 x 0.05)\n";
    std::cout << "  Calculated overlap: " << std::fixed << std::setprecision(8) << overlap << "\n";
    std::cout << "  Below threshold (0.01): " << (overlap < 0.01 ? "YES" : "NO") << "\n\n";
}

// Test with legacy verify_oarea cases (matches org/src/ana_fncs.c lines 3360-3404)
void test_legacy_verify_oarea() {
    std::cout << "=== Test 5: Legacy verify_oarea Test Cases ===\n";

    // Test case from legacy verify_oarea function
    // These match the exact test cases used in the legacy code

    Residue res1("  A", 1, 'A');
    Residue res2("  A", 2, 'A');

    // Case 1: Square a: (0,0), (10,0), (10,10), (0,10)
    std::cout << "Test 5a: Square self-intersection (should equal area)\n";
    res1 = Residue("  A", 1, 'A');
    res1.add_atom(Atom(" C4 ", Vector3D(0, 0, 0)));
    res1.add_atom(Atom(" N3 ", Vector3D(10, 0, 0)));
    res1.add_atom(Atom(" C2 ", Vector3D(10, 10, 0)));
    res1.add_atom(Atom(" N1 ", Vector3D(0, 10, 0)));

    res2 = res1; // Same polygon

    Vector3D oave(0, 0, 0);
    Vector3D zave(0, 0, 1);

    BasePairValidator validator;
    double overlap = validator.calculate_overlap_area(res1, res2, oave, zave);

    std::cout << "  Expected: ~100.0 (10x10 square area)\n";
    std::cout << "  Calculated: " << std::fixed << std::setprecision(6) << overlap << "\n\n";

    // Case 2: Square b: (5,0), (10,0), (10,5), (5,5)
    std::cout << "Test 5b: Smaller square self-intersection\n";
    res1 = Residue("  A", 1, 'A');
    res1.add_atom(Atom(" C4 ", Vector3D(5, 0, 0)));
    res1.add_atom(Atom(" N3 ", Vector3D(10, 0, 0)));
    res1.add_atom(Atom(" C2 ", Vector3D(10, 5, 0)));
    res1.add_atom(Atom(" N1 ", Vector3D(5, 5, 0)));

    res2 = res1;

    overlap = validator.calculate_overlap_area(res1, res2, oave, zave);

    std::cout << "  Expected: ~25.0 (5x5 square area)\n";
    std::cout << "  Calculated: " << std::fixed << std::setprecision(6) << overlap << "\n\n";

    // Case 3: Overlap between a and b
    std::cout << "Test 5c: Overlap between squares a and b (legacy test)\n";
    res1 = Residue("  A", 1, 'A');
    res1.add_atom(Atom(" C4 ", Vector3D(0, 0, 0)));
    res1.add_atom(Atom(" N3 ", Vector3D(10, 0, 0)));
    res1.add_atom(Atom(" C2 ", Vector3D(10, 10, 0)));
    res1.add_atom(Atom(" N1 ", Vector3D(0, 10, 0)));

    res2 = Residue("  A", 2, 'A');
    res2.add_atom(Atom(" C4 ", Vector3D(5, 0, 0)));
    res2.add_atom(Atom(" N3 ", Vector3D(10, 0, 0)));
    res2.add_atom(Atom(" C2 ", Vector3D(10, 5, 0)));
    res2.add_atom(Atom(" N1 ", Vector3D(5, 5, 0)));

    overlap = validator.calculate_overlap_area(res1, res2, oave, zave);

    std::cout << "  Expected: ~25.0 (5x5 overlap region)\n";
    std::cout << "  Calculated: " << std::fixed << std::setprecision(6) << overlap << "\n\n";

    // Case 4: Complex polygon c: (-3,-2), (-1,4), (6,1), (3,10), (-4,9)
    std::cout << "Test 5d: Complex polygon self-intersection\n";
    res1 = Residue("  A", 1, 'A');
    res1.add_atom(Atom(" C4 ", Vector3D(-3, -2, 0)));
    res1.add_atom(Atom(" N3 ", Vector3D(-1, 4, 0)));
    res1.add_atom(Atom(" C2 ", Vector3D(6, 1, 0)));
    res1.add_atom(Atom(" N1 ", Vector3D(3, 10, 0)));
    res1.add_atom(Atom(" C6 ", Vector3D(-4, 9, 0)));

    res2 = res1;

    overlap = validator.calculate_overlap_area(res1, res2, oave, zave);

    std::cout << "  Expected: > 0.0 (polygon area)\n";
    std::cout << "  Calculated: " << std::fixed << std::setprecision(6) << overlap << "\n\n";

    // Case 5: Overlap between a and c
    std::cout << "Test 5e: Overlap between square a and complex polygon c\n";
    res1 = Residue("  A", 1, 'A');
    res1.add_atom(Atom(" C4 ", Vector3D(0, 0, 0)));
    res1.add_atom(Atom(" N3 ", Vector3D(10, 0, 0)));
    res1.add_atom(Atom(" C2 ", Vector3D(10, 10, 0)));
    res1.add_atom(Atom(" N1 ", Vector3D(0, 10, 0)));

    res2 = Residue("  A", 2, 'A');
    res2.add_atom(Atom(" C4 ", Vector3D(-3, -2, 0)));
    res2.add_atom(Atom(" N3 ", Vector3D(-1, 4, 0)));
    res2.add_atom(Atom(" C2 ", Vector3D(6, 1, 0)));
    res2.add_atom(Atom(" N1 ", Vector3D(3, 10, 0)));
    res2.add_atom(Atom(" C6 ", Vector3D(-4, 9, 0)));

    overlap = validator.calculate_overlap_area(res1, res2, oave, zave);

    std::cout << "  Expected: > 0.0 (partial overlap)\n";
    std::cout << "  Calculated: " << std::fixed << std::setprecision(6) << overlap << "\n\n";
}

// Test with real PDB data
void test_pdb_data(const std::string& pdb_file, size_t res1_idx, size_t res2_idx) {
    std::cout << "=== Test 5: Real PDB Data ===\n";
    std::cout << "  PDB: " << pdb_file << "\n";
    std::cout << "  Residues: " << res1_idx << ", " << res2_idx << "\n";

    // TODO: Load PDB and extract residues
    // For now, just show the structure
    std::cout << "  (Not yet implemented - need PDB parser)\n\n";
}

int main(int argc, char* argv[]) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Overlap Calculation Test Suite\n";
    std::cout << "==============================\n\n";

    test_simple_overlap();
    test_no_overlap();
    test_inside_overlap();
    test_complex_shapes();
    test_legacy_verify_oarea();

    if (argc >= 4) {
        test_pdb_data(argv[1], std::stoul(argv[2]), std::stoul(argv[3]));
    }

    return 0;
}

# Stage 1: Core Geometry Classes

## Objectives

Implement the fundamental geometry classes (Vector3D, Matrix3D, Quaternion) and geometric utilities that form the foundation for all subsequent work.

## Duration

**1 week**

## Dependencies

- ✅ Stage 0: Project Setup & Infrastructure

## Tasks

### Task 1.1: Implement Vector3D
- [x] Create `include/x3dna/geometry/vector3d.hpp`
- [x] Implement basic constructors
- [x] Implement getters/setters
- [x] Implement arithmetic operations (+, -, *, /)
- [x] Implement dot product
- [x] Implement cross product
- [x] Implement length/normalize
- [x] Implement distance calculations
- [x] Add JSON serialization
- [x] Write comprehensive unit tests

**Key Methods**:
```cpp
class Vector3D {
    Vector3D(double x, double y, double z);
    double x(), y(), z() const;
    double length() const;
    Vector3D normalized() const;
    double dot(const Vector3D& other) const;
    Vector3D cross(const Vector3D& other) const;
    double distance_to(const Vector3D& other) const;
    nlohmann::json to_json() const;
};
```

**Deliverable**: Fully tested Vector3D class

### Task 1.2: Implement Matrix3D
- [x] Create `include/x3dna/geometry/matrix3d.hpp`
- [x] Implement constructors (identity, from array, from rows/columns)
- [x] Implement matrix-vector multiplication
- [x] Implement matrix-matrix multiplication
- [x] Implement transpose
- [x] Implement inverse
- [x] Implement determinant
- [x] Implement row/column access
- [x] Add JSON serialization (including legacy format)
- [x] Write comprehensive unit tests

**Key Methods**:
```cpp
class Matrix3D {
    Matrix3D();  // Identity
    Matrix3D(const std::array<double, 9>& values);
    Vector3D operator*(const Vector3D& vec) const;
    Matrix3D operator*(const Matrix3D& other) const;
    Matrix3D transpose() const;
    Matrix3D inverse() const;
    double determinant() const;
    std::array<double, 9> as_array() const;  // Row-major
    nlohmann::json to_json_legacy() const;  // 3x3 nested array
};
```

**Deliverable**: Fully tested Matrix3D class

### Task 1.3: Implement Quaternion (Optional but Recommended)
- [ ] Create `include/x3dna/geometry/Quaternion.hpp`
- [ ] Implement quaternion operations
- [ ] Implement conversion to/from rotation matrix
- [ ] Implement slerp (spherical interpolation)
- [ ] Add unit tests

**Note**: Useful for smooth rotations and some algorithms

**Deliverable**: Quaternion class (if implemented)

### Task 1.4: Implement LeastSquaresFitter
- [x] Create `include/x3dna/geometry/least_squares_fitter.hpp`
- [x] Implement quaternion-based algorithm (Jacobi eigenvalue decomposition)
- [x] Implement point set alignment
- [x] Return rotation matrix, translation, and RMS
- [x] Add JSON serialization (matches legacy `ls_fitting` format)
- [x] Write comprehensive unit tests

**Key Methods**:
```cpp
class LeastSquaresFitter {
    struct FitResult {
        Matrix3D rotation;
        Vector3D translation;
        double rms;
    };
    
    FitResult fit(
        const std::vector<Vector3D>& points1,
        const std::vector<Vector3D>& points2
    );
};
```

**Critical**: This must match the original `ls_fitting()` algorithm exactly!

**Deliverable**: LeastSquaresFitter with tests

### Task 1.5: Implement Geometry Utilities
- [ ] Create `include/x3dna/geometry/GeometryUtils.hpp`
- [ ] Implement angle calculations
- [ ] Implement plane fitting
- [ ] Implement distance calculations
- [ ] Implement coordinate transformations
- [ ] Add helper functions for common operations
- [ ] Write unit tests

**Deliverable**: Geometry utility functions

### Task 1.6: Integration Testing
- [x] Test Vector3D with Matrix3D operations
- [x] Test LeastSquaresFitter with known point sets
- [x] Verify numerical stability
- [x] Test edge cases (zero vectors, singular matrices, etc.)
- [x] Regression tests with real-world data from legacy JSON files

**Deliverable**: Integration tests passing

## Testing Plan

### Unit Tests

#### Vector3D Tests
- [x] Test constructors
- [x] Test getters/setters
- [x] Test arithmetic operations
- [x] Test dot product (known values)
- [x] Test cross product (orthogonality)
- [x] Test length calculation
- [x] Test normalization
- [x] Test distance calculations
- [x] Test JSON serialization/deserialization
- [x] Test edge cases (zero vector, etc.)

#### Matrix3D Tests
- [x] Test identity matrix
- [x] Test constructors
- [x] Test matrix-vector multiplication
- [x] Test matrix-matrix multiplication
- [x] Test transpose (verify property: (AB)^T = B^T A^T)
- [x] Test inverse (verify: A * A^-1 = I)
- [x] Test determinant calculation
- [x] Test row/column access
- [x] Test JSON serialization (both formats)
- [x] Test edge cases (singular matrix, etc.)

#### LeastSquaresFitter Tests
- [x] Test with known point sets (simple translation)
- [x] Test with known point sets (simple rotation)
- [x] Test with known point sets (combined transform)
- [x] Test RMS calculation accuracy
- [x] Test with noisy data
- [x] Test with degenerate cases (collinear points)
- [x] Compare with reference implementation (legacy JSON regression tests)
- [x] Test JSON serialization (legacy format)

### Integration Tests
- [x] Test complete transformation pipeline
- [x] Test LeastSquaresFitter with real PDB atom coordinates
- [x] Verify numerical precision (compare with original code)

### Regression Tests
- [x] Load legacy JSON with `ls_fitting` records
- [x] Recreate point sets from JSON
- [x] Run LeastSquaresFitter
- [x] Compare rotation matrices and translations
- [x] Verify RMS values match within tolerance (0.001)

## Test Data

### Vector3D Test Cases
```cpp
// Known dot product
Vector3D v1(1, 2, 3);
Vector3D v2(4, 5, 6);
EXPECT_DOUBLE_EQ(v1.dot(v2), 32.0);

// Known cross product
Vector3D v3(1, 0, 0);
Vector3D v4(0, 1, 0);
EXPECT_EQ(v3.cross(v4), Vector3D(0, 0, 1));
```

### Matrix3D Test Cases
```cpp
// Identity matrix
Matrix3D I;
Vector3D v(1, 2, 3);
EXPECT_EQ(I * v, v);

// Inverse property
Matrix3D M = create_test_matrix();
Matrix3D M_inv = M.inverse();
Matrix3D I_result = M * M_inv;
EXPECT_TRUE(I_result.approximately_equals(Matrix3D(), 1e-6));
```

### LeastSquaresFitter Test Cases
```cpp
// Simple translation
std::vector<Vector3D> points1 = {{0,0,0}, {1,0,0}, {0,1,0}};
std::vector<Vector3D> points2 = {{1,1,1}, {2,1,1}, {1,2,1}};
auto result = fitter.fit(points1, points2);
EXPECT_NEAR(result.translation.x(), 1.0, 1e-6);
EXPECT_NEAR(result.rms, 0.0, 1e-6);
```

## Success Criteria

- [x] All unit tests pass (100% for geometry classes)
- [x] All integration tests pass
- [ ] Code coverage > 90% for geometry classes (pending measurement)
- [x] LeastSquaresFitter produces results matching legacy JSON within 0.001 tolerance
- [ ] No memory leaks (valgrind clean - pending)
- [x] Performance is acceptable (benchmark if needed)
- [x] JSON serialization works for all classes
- [x] Documentation generated and reviewed

## Code Quality Checks

- [ ] No compiler warnings
- [ ] Passes clang-tidy checks
- [ ] Follows code style (clang-format)
- [ ] All public methods documented
- [ ] Const-correctness throughout

## Deliverables

1. ✅ `Vector3D` class fully implemented and tested
2. ✅ `Matrix3D` class fully implemented and tested
3. ✅ `LeastSquaresFitter` class fully implemented and tested
4. ✅ `GeometryUtils` utility functions
5. ✅ Comprehensive unit test suite
6. ✅ Integration tests
7. ✅ Regression tests comparing with legacy JSON
8. ✅ Documentation (Doxygen)

## Files Created

```
include/x3dna/geometry/
├── vector3d.hpp ✅
├── matrix3d.hpp ✅
├── least_squares_fitter.hpp ✅
└── (Quaternion and GeometryUtils - optional, not yet implemented)

tests/unit/geometry/
├── test_vector3d.cpp ✅
├── test_matrix3d.cpp ✅
└── test_least_squares_fitter.cpp ✅
```

## Risks & Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| LeastSquaresFitter doesn't match original | High | Compare with legacy JSON early, use reference implementation |
| Numerical precision issues | Medium | Use appropriate tolerance, test with known values |
| Matrix inverse fails for singular cases | Medium | Add proper error handling, test edge cases |
| Performance issues | Low | Profile if needed, optimize hot paths |

## Validation Checklist

Before moving to Stage 2:
- [x] All geometry classes compile and link
- [x] All unit tests pass
- [x] Integration tests pass
- [x] Regression tests match legacy JSON
- [x] Code reviewed
- [x] Documentation complete
- [ ] No memory leaks (pending valgrind)
- [x] Performance acceptable

## Next Stage

After completing Stage 1, proceed to **Stage 2: Core Domain Objects** (`STAGE_02_CORE_OBJECTS.md`)

---

*Estimated Completion: Week 2*


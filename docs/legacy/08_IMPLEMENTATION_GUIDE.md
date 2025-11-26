# Legacy to Modern Implementation Guide

**Date**: 2025-01-XX  
**Purpose**: Step-by-step guide for matching legacy behavior exactly in modern C++ implementation  
**Status**: Comprehensive conversion reference

---

## Table of Contents

1. [Index Conversion](#index-conversion)
2. [Matching Base Frames](#matching-base-frames)
3. [Matching H-Bond Detection](#matching-h-bond-detection)
4. [Matching Validation Logic](#matching-validation-logic)
5. [Matching Greedy Matching](#matching-greedy-matching)
6. [Common Conversion Pitfalls](#common-conversion-pitfalls)
7. [Verification Checklist](#verification-checklist)
8. [Debugging Strategies](#debugging-strategies)

---

## Index Conversion

### The Critical Challenge

**Legacy**: Everything uses 1-based indexing  
**Modern**: Standard C++ uses 0-based indexing  
**Solution**: Convert at boundaries only

### Conversion Strategy

**Golden Rule**: Store and compute in 0-based internally, convert only when interfacing with legacy or JSON output.

#### Residue Index Conversion

```cpp
// Legacy code uses:
for (i = 1; i <= num_residue; i++) {
    process_residue(i);  // i ∈ [1, num_residue]
}

// Modern equivalent:
for (size_t i = 0; i < residues.size(); i++) {
    process_residue(i);  // i ∈ [0, residues.size()-1]
    // Legacy index: legacy_idx = i + 1
}
```

**When to Convert**:
- ✅ **Store legacy_residue_idx on atoms** during PDB parsing
- ✅ **Convert to legacy (1-based) when writing JSON** for comparison
- ✅ **Use legacy indices when matching residues** from legacy output
- ❌ **Don't convert inside algorithms** - use 0-based consistently

#### Array Access Conversion

```cpp
// Legacy: orien[i][j] where i=1..n, j=1..9
double zx_legacy = orien[residue_idx][7];  // z-axis x-component

// Modern: Use std::vector or Eigen::Matrix, 0-based
double zx_modern = frames[residue_idx - 1][6];  // Convert legacy→modern
// or store 0-based: frames[residue_idx][6]
```

**Frame Matrix Mapping**:
```cpp
// Legacy stores flattened 3×3 as [1..9]
// Modern should also flatten, but 0-indexed
struct Frame {
    double data[9];  // [0..8], not [1..9]
};

// Conversion:
// Legacy [1..9] → Modern [0..8]
// Access pattern: modern[i] = legacy[i+1]
// Or: legacy[i] = modern[i-1] when i ∈ [1,9]
```

#### Frame Element Access

```cpp
// Legacy z-axis access (1-based):
double zx = orien[i][7];  // Actually indices [6..8] in 0-based terms
double zy = orien[i][8];
double zz = orien[i][9];

// Modern equivalent (0-based):
double zx = frame.data[6];  // Row 3, Col 1 (z-axis x-component)
double zy = frame.data[7];  // Row 3, Col 2
double zz = frame.data[8];  // Row 3, Col 3

// Mapping:
// Legacy orien[i][j] where j ∈ [1,9]
// Modern frame.data[k] where k = j - 1 ∈ [0,8]
```

---

## Matching Base Frames

### Critical Requirements

1. **Template Loading**: Exact same template files
2. **Ring Atom Matching**: Same matching order and criteria
3. **Least-Squares Fitting**: Exact algorithm
4. **Matrix Storage**: Same flattening pattern

### Step-by-Step Matching

#### Step 1: Template File Path

```cpp
// Legacy constructs path:
// {BDIR}/Atomic_{base}.pdb
std::string template_path = base_dir + "/Atomic_" + base_letter + ".pdb";

// Must match exactly - same files, same path construction
```

#### Step 2: Ring Atom Matching

```cpp
// Legacy order (RA_LIST):
static const char* RING_ATOMS[] = {
    " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
};

// Must match in this exact order
// For purines: all 9 atoms
// For pyrimidines: first 6 only
int ring_count = (is_purine) ? 9 : 6;

// Match atoms in this order
std::vector<Atom*> matched_standard;
std::vector<Atom*> matched_experimental;
for (int i = 0; i < ring_count; i++) {
    auto std_atom = find_atom(standard_residue, RING_ATOMS[i]);
    auto exp_atom = find_atom(experimental_residue, RING_ATOMS[i]);
    if (std_atom && exp_atom) {
        matched_standard.push_back(std_atom);
        matched_experimental.push_back(exp_atom);
    }
}

// Require minimum 3 matches
if (matched_standard.size() < 3) {
    // Skip this residue
}
```

#### Step 3: Least-Squares Fitting

```cpp
// Use exact same algorithm from ls_fitting()
// See 04_ALGORITHMS.md for complete implementation

// Key requirements:
// 1. Calculate covariance matrix
// 2. Build 4×4 quaternion matrix
// 3. Eigenvalue decomposition (largest eigenvalue)
// 4. Convert quaternion to rotation matrix
// 5. Calculate translation
// 6. Return RMS

// Must match exactly - use same numerical recipes algorithm
```

#### Step 4: Matrix Storage

```cpp
// Legacy stores as flattened array [1..9]
// Modern: store as [0..8] but convert for JSON

struct RotationMatrix {
    double data[9];  // 0-indexed internally
    
    // Convert to legacy format for JSON
    std::vector<double> to_legacy_format() const {
        std::vector<double> result(10);  // [0] unused, [1..9] used
        for (int i = 0; i < 9; i++) {
            result[i + 1] = data[i];
        }
        return result;
    }
};

// Mapping:
// Legacy orien[residue][j] where j ∈ [1,9]
// Modern frame.data[k] where k = j - 1 ∈ [0,8]
```

### Verification Points

```cpp
// 1. Same RMS values (within floating-point precision)
assert(std::abs(modern_rms - legacy_rms) < 1e-10);

// 2. Rotation matrices match (element-by-element)
for (int i = 0; i < 9; i++) {
    assert(std::abs(modern_matrix[i] - legacy_matrix[i+1]) < 1e-10);
}

// 3. Origins match
for (int i = 0; i < 3; i++) {
    assert(std::abs(modern_origin[i] - legacy_origin[i+1]) < 1e-10);
}
```

---

## Matching H-Bond Detection

### Critical Requirements

1. **Exact same distance limits**: `hb_lower = 1.8`, `hb_dist1 = 4.0`
2. **Exact same conflict resolution**: Iterative algorithm
3. **Exact same validation**: Donor/acceptor roles
4. **Exact same output format**: H-bond string format

### Step-by-Step Matching

#### Step 1: Initial Detection

```cpp
// Legacy uses nested loops:
std::vector<HBond> initial_hbonds;
for (int m = seidx[i][1]; m <= seidx[i][2]; m++) {
    for (int n = seidx[j][1]; n <= seidx[j][2]; n++) {
        if (good_hbatoms(atom_i, atom_j) && 
            within_limits(atom_i->xyz, atom_j->xyz, 1.8, 4.0)) {
            initial_hbonds.push_back({
                atom_i, atom_j, 
                distance(atom_i->xyz, atom_j->xyz)
            });
        }
    }
}

// Must check atoms in same order
// Use residue start/end indices to iterate
```

#### Step 2: Conflict Resolution

```cpp
// Use exact algorithm from hb_atompair()
// See 04_ALGORITHMS.md

// Critical points:
// 1. Iterative until convergence
// 2. Find best H-bond for each atom (by distance)
// 3. Mark conflicts by negating distances
// 4. Calculate linkage types

// Must match iteration order and conflict marking exactly
```

#### Step 3: Validation

```cpp
// Use exact logic from validate_hbonds()
for (auto& hbond : hbonds) {
    if (hbond.distance <= 0.0) {  // Conflicted
        hbond.distance = std::abs(hbond.distance);
        hbond.type = donor_acceptor(base_i, base_j, 
                                    hbond.atom1, hbond.atom2);
    } else {
        hbond.type = ' ';  // Normal
    }
}

// Filter invalid H-bonds
// Apply same distance and type checks
```

### Verification Points

```cpp
// 1. Same number of H-bonds found
assert(modern_count == legacy_count);

// 2. Same atoms involved
for (int i = 0; i < count; i++) {
    assert(modern_hbonds[i].atom1 == legacy_hbonds[i].atom1);
    assert(modern_hbonds[i].atom2 == legacy_hbonds[i].atom2);
}

// 3. Same distances (exact match)
assert(std::abs(modern_dist - legacy_dist) < 1e-10);

// 4. Same types (' ', '-', '*')
assert(modern_type == legacy_type);
```

---

## Matching Validation Logic

### Critical Parameters

**Must match exactly**:
```cpp
constexpr double MIN_DORG = 0.0;
constexpr double MAX_DORG = 15.0;
constexpr double MIN_DV = 0.0;
constexpr double MAX_DV = 2.5;
constexpr double MAX_PLANE_ANGLE = 65.0;
constexpr double MIN_DNN = 4.5;
constexpr double MAX_DNN = 1e18;  // XBIG
constexpr int MIN_BASE_HB = 1;
constexpr double OVERLAP_THRESHOLD = 0.01;
```

### Step-by-Step Matching

#### Step 1: Geometric Constraints

```cpp
// Calculate metrics exactly as legacy
double dorg = distance(org_i, org_j);
double dv = std::abs(dot(dorg_vector, zave));
double plane_angle = angle_between(z_axis_i, z_axis_j);
double dNN = distance(NC1xyz_i, NC1xyz_j);

// Check constraints in same order
bool passes = 
    (dorg >= MIN_DORG && dorg <= MAX_DORG) &&
    (dv >= MIN_DV && dv <= MAX_DV) &&
    (plane_angle <= MAX_PLANE_ANGLE) &&
    (dNN >= MIN_DNN && dNN <= MAX_DNN);
```

#### Step 2: Overlap Check

```cpp
// Use exact algorithm from get_oarea()
// See 04_ALGORITHMS.md

// Critical: Must use same projection and intersection algorithm
double overlap = calculate_overlap_area(residue_i, residue_j, 
                                        oave, zave);

if (overlap >= OVERLAP_THRESHOLD) {
    return false;  // Reject
}
```

#### Step 3: H-Bond Requirement

```cpp
// Count H-bonds using exact detection algorithm
int hbond_count = detect_and_count_hbonds(residue_i, residue_j);

// Check requirement
if (hbond_count < MIN_BASE_HB) {
    return false;  // Reject
}
```

#### Step 4: Quality Score

```cpp
// Calculate exactly as legacy
double quality_score = 
    dorg + 
    2.0 * dv + 
    plane_angle / 20.0;

// Add H-bond adjustment
quality_score += adjust_pair_quality(residue_i, residue_j, hbond_count);

// Lower score = better pair (must match legacy scoring)
```

### Verification Points

```cpp
// 1. Same validation decision (accept/reject)
assert(modern_valid == legacy_valid);

// 2. Same quality scores
assert(std::abs(modern_score - legacy_score) < 1e-10);

// 3. Same metrics
assert(std::abs(modern_dorg - legacy_dorg) < 1e-10);
assert(std::abs(modern_dv - legacy_dv) < 1e-10);
// ... etc
```

---

## Matching Greedy Matching

### Critical Requirements

1. **Iteration order**: Process residues in same order
2. **Mutual matching**: Both bases must choose each other
3. **Quality scoring**: Lower score = better
4. **Tie-breaking**: Must match legacy behavior

### Step-by-Step Matching

#### Step 1: Iteration Structure

```cpp
// Legacy uses nested iterations
std::vector<bool> matched(residues.size(), false);
int num_pairs = 0;

while (true) {
    size_t old_count = std::count(matched.begin(), matched.end(), true);
    
    // Process all residues
    for (size_t i = 0; i < residues.size(); i++) {
        if (matched[i] || !is_nucleotide(residues[i])) {
            continue;
        }
        
        // Find best partner
        auto best_partner = find_best_partner(i, matched);
        
        if (best_partner.has_value()) {
            size_t j = best_partner.value();
            
            // Check mutual match
            auto mutual_partner = find_best_partner(j, matched);
            if (mutual_partner.has_value() && mutual_partner.value() == i) {
                // Match found!
                matched[i] = true;
                matched[j] = true;
                pairs.push_back({i, j});
                num_pairs++;
            }
        }
    }
    
    size_t new_count = std::count(matched.begin(), matched.end(), true);
    if (new_count == old_count) {
        break;  // No more pairs found
    }
}
```

#### Step 2: Best Partner Selection

```cpp
std::optional<size_t> find_best_partner(size_t i, 
                                        const std::vector<bool>& matched) {
    double best_score = std::numeric_limits<double>::max();
    std::optional<size_t> best_j;
    
    // Check all other residues
    for (size_t j = 0; j < residues.size(); j++) {
        if (i == j || matched[j] || !is_nucleotide(residues[j])) {
            continue;
        }
        
        auto validation_result = check_pair(i, j);
        if (validation_result.valid && 
            validation_result.quality_score < best_score) {
            best_score = validation_result.quality_score;
            best_j = j;
        }
    }
    
    return best_j;
}
```

### Verification Points

```cpp
// 1. Same number of pairs found
assert(modern_pairs.size() == legacy_num_pairs);

// 2. Same pairs (same residue indices)
for (size_t i = 0; i < modern_pairs.size(); i++) {
    assert(modern_pairs[i].first == legacy_pairs[i].residue_i);
    assert(modern_pairs[i].second == legacy_pairs[i].residue_j);
}

// 3. Same order (if order matters)
```

---

## Common Conversion Pitfalls

### Pitfall 1: Residue Index Mismatch

**Problem**: Using simple counter instead of legacy_residue_idx

```cpp
// ❌ WRONG
size_t idx = 0;
for (auto& residue : residues) {
    if (is_nucleotide(residue)) {
        record_frame(idx++);  // Wrong - skips amino acids
    }
}

// ✅ CORRECT
for (auto& residue : residues) {
    if (is_nucleotide(residue)) {
        int legacy_idx = residue.atoms()[0].legacy_residue_idx();
        record_frame(legacy_idx - 1);  // Convert to 0-based
    }
}
```

### Pitfall 2: Frame Matrix Indexing

**Problem**: Wrong element access for z-axis

```cpp
// ❌ WRONG
double zx = frame[6];  // If frame is 0-indexed, this is wrong
double zx = frame[7];  // Confusing legacy [7] with modern [6]

// ✅ CORRECT
// Legacy stores z-axis at indices [6,7,8] (0-based view of [7,8,9])
// Modern should store at [6,7,8] in 0-indexed array
double zx = frame.data[6];  // Row 3, Col 1
double zy = frame.data[7];  // Row 3, Col 2
double zz = frame.data[8];  // Row 3, Col 3
```

### Pitfall 3: H-Bond Distance Sign

**Problem**: Not handling negative distances (conflicts)

```cpp
// ❌ WRONG
for (auto& hbond : hbonds) {
    double dist = hbond.distance;  // May be negative!
    if (dist < threshold) { ... }
}

// ✅ CORRECT
for (auto& hbond : hbonds) {
    bool is_conflict = (hbond.distance <= 0.0);
    double dist = std::abs(hbond.distance);  // Always positive
    if (dist < threshold && !is_conflict) { ... }
}
```

### Pitfall 4: Parameter Value Mismatch

**Problem**: Using different parameter values

```cpp
// ❌ WRONG
constexpr double HB_LOWER = 1.5;  // Different from legacy

// ✅ CORRECT
constexpr double HB_LOWER = 1.8;  // Match legacy exactly
constexpr double HB_DIST1 = 4.0;
constexpr double HB_DIST2 = 0.0;  // Critical: must be 0.0
```

### Pitfall 5: Iteration Order

**Problem**: Different iteration order leads to different results

```cpp
// ❌ WRONG - may give different results
std::set<Residue*> processed;
for (auto& residue : residues) {
    if (processed.count(&residue)) continue;
    // Process...
}

// ✅ CORRECT - match legacy order
for (size_t i = 0; i < residues.size(); i++) {
    // Process residues[i] in order
}
```

---

## Verification Checklist

### Frame Calculation

- [ ] Same RMS values (within 1e-10)
- [ ] Same rotation matrices (element-by-element)
- [ ] Same origins (x, y, z coordinates)
- [ ] Same template files used
- [ ] Same ring atoms matched

### H-Bond Detection

- [ ] Same number of H-bonds found
- [ ] Same atoms involved
- [ ] Same distances (exact match)
- [ ] Same conflict resolution
- [ ] Same H-bond types (' ', '-', '*')

### Validation

- [ ] Same pairs validated/rejected
- [ ] Same quality scores
- [ ] Same geometric metrics (dorg, dv, plane_angle, dNN)
- [ ] Same overlap calculations
- [ ] Same H-bond counts

### Base Pair Finding

- [ ] Same number of pairs found
- [ ] Same pairs (same residue indices)
- [ ] Same iteration count
- [ ] Same mutual matching decisions

### Overall

- [ ] JSON output matches legacy format
- [ ] All indices match (residue, atom, pair)
- [ ] All parameters match (constants, thresholds)
- [ ] All algorithms match (exact same logic)

---

## Debugging Strategies

### Strategy 1: Step-by-Step Comparison

```cpp
// Add detailed logging at each step
void calculate_frame(size_t residue_idx) {
    LOG("Frame calculation for residue {}", residue_idx);
    LOG("  Template: {}", template_path);
    LOG("  Matched atoms: {}", matched_count);
    LOG("  RMS: {}", rms_value);
    LOG("  Rotation matrix: {}", rotation_matrix);
    LOG("  Origin: {}", origin);
    
    // Compare with legacy JSON output
}
```

### Strategy 2: Isolated Function Testing

```cpp
// Test each function independently
TEST(LeastSquaresFitting, MatchesLegacy) {
    // Load test data from legacy JSON
    auto standard_points = load_from_json("legacy_standard.json");
    auto experimental_points = load_from_json("legacy_experimental.json");
    
    auto result = ls_fitting(standard_points, experimental_points);
    
    // Compare with legacy result
    assert(result.rms == legacy_rms);
    assert(result.rotation_matrix == legacy_matrix);
}
```

### Strategy 3: Binary Comparison

```cpp
// Compare binary representations
void compare_frames(const Frame& modern, const LegacyFrame& legacy) {
    for (int i = 0; i < 9; i++) {
        double diff = std::abs(modern.data[i] - legacy.data[i+1]);
        if (diff > 1e-10) {
            LOG_ERROR("Frame mismatch at index {}: {} vs {}", 
                      i, modern.data[i], legacy.data[i+1]);
        }
    }
}
```

### Strategy 4: Incremental Validation

```cpp
// Validate after each major step
void find_pairs() {
    // Step 1: Calculate frames
    calculate_all_frames();
    validate_frames();  // Compare with legacy
    
    // Step 2: Detect H-bonds
    detect_all_hbonds();
    validate_hbonds();  // Compare with legacy
    
    // Step 3: Validate pairs
    validate_all_pairs();
    validate_pair_validation();  // Compare with legacy
    
    // Step 4: Find pairs
    run_greedy_matching();
    validate_pairs();  // Compare with legacy
}
```

---

**Next**: [Core Functions](03_CORE_FUNCTIONS.md) for detailed function specifications


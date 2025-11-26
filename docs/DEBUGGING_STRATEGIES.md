# Debugging Strategies for 100% Match

## Quick Reference

### Current Status
- **Match Rate**: 99.5% (1044/1048 pairs)
- **Perfect PDBs**: 94/100 (94%)
- **Issues**: 6 PDBs with differences

### Top Priorities
1. **Residue Inclusion** (3KNC: 0→8 pairs, 5UJ2: 5→6 pairs)
2. **Quality Score Calculation** (3G8T, 6CAQ: 10-12 missing pairs)
3. **H-bond Code Refactoring** (Code organization)
4. **bp_type_id Full Implementation** (Requires Stage 6)

---

## Strategy 1: Break Up PDBs into Smaller Chunks

### Why?
- Large PDBs are hard to debug
- Isolate problematic sections
- Faster iteration cycles

### Implementation

**Tool: `tools/create_pdb_subset.cpp`**
```cpp
// Extract specific residues from PDB
// Usage: create_pdb_subset input.pdb output.pdb --residues 90-95,160-165
```

**Script: `scripts/generate_test_cases.py`**
```python
# Automatically create subsets around problematic pairs
# Example: Extract residues around (92, 160) in 3G8T
# Includes context (neighboring residues)
```

### Use Cases

**3G8T Missing Pairs**:
- Extract residues 85-100, 155-170 (around pair 92, 160)
- Extract residues 940-950 (around pair 946, 947)

**6CAQ Missing Pairs**:
- Extract residues 70-85 (around pair 75, 78)
- Extract residues 960-1030 (around pair 968, 1024)

**3KNC No Pairs**:
- Extract first 20 residues to see why only 16/66 are recognized

### Benefits
- Faster debugging cycles
- Easier to understand issues
- Can test fixes on smaller datasets

---

## Strategy 2: Create Legacy Test Executables

### Why?
- Isolate each algorithm for testing
- Compare step-by-step with modern
- Understand exactly how legacy works

### Proposed Executables

#### 1. `org/bin/test_pdb_parsing`
**Purpose**: Test PDB parsing only  
**Input**: PDB file  
**Output**: JSON with atom/residue counts, residue identification  
**Compare**: Modern parsing output

**Key Outputs**:
- Total atoms
- Total residues
- Nucleotide residues (with indices)
- Residue types (A, C, G, T, amino acid, other)
- HETATM handling

#### 2. `org/bin/test_frame_calculation`
**Purpose**: Test base frame calculation only  
**Input**: PDB file + residue index  
**Output**: JSON with frame (orien, org), RMS fit, matched atoms  
**Compare**: Modern frame calculation

**Key Outputs**:
- Template file used
- Matched atoms
- RMS fit value
- Rotation matrix (orien)
- Origin (org)

#### 3. `org/bin/test_hbond_detection`
**Purpose**: Test H-bond detection only  
**Input**: PDB file + two residue indices  
**Output**: JSON with H-bonds found, validation steps  
**Compare**: Modern H-bond detection

**Key Outputs**:
- Initial H-bonds (before conflict resolution)
- After conflict resolution
- After validation
- Final H-bonds (type != ' ')

#### 4. `org/bin/test_pair_validation`
**Purpose**: Test pair validation only  
**Input**: PDB file + two residue indices  
**Output**: JSON with all validation checks, quality_score, bp_type_id  
**Compare**: Modern validation

**Key Outputs**:
- Geometric values (dorg, d_v, plane_angle, dNN)
- Validation checks (distance, d_v, plane_angle, dNN, overlap)
- H-bond counts (num_base_hb, num_o2_hb)
- Quality score (base + adjustments)
- bp_type_id

#### 5. `org/bin/test_quality_score`
**Purpose**: Test quality score calculation only  
**Input**: PDB file + two residue indices  
**Output**: JSON with base score, adjust_pairQuality, bp_type_id, final score  
**Compare**: Modern quality score

**Key Outputs**:
- Base quality_score (dorg + 2.0*d_v + plane_angle/20.0)
- adjust_pairQuality value
- bp_type_id value
- Final adjusted quality_score (rtn_val[5])

### Implementation

**File Structure**:
```
org/src/
├── test_pdb_parsing.c      # New
├── test_frame_calculation.c # New
├── test_hbond_detection.c   # New
├── test_pair_validation.c   # New
└── test_quality_score.c     # New
```

**Key Features**:
- Use existing functions (no algorithm changes)
- Add JSON output for easy comparison
- Focused output (only relevant information)
- Keep existing `find_pair` executable unchanged

### Benefits
- Can test each algorithm independently
- Easy to compare with modern
- Understand legacy behavior exactly
- Faster debugging

---

## Strategy 3: Algorithm-Specific Debugging

### For Each Algorithm, Create Focused Comparison Tools

#### PDB Parsing
**Tool**: `tools/compare_pdb_parsing.cpp`
- Compare atom-by-atom parsing
- Compare residue identification
- Compare HETATM handling
- Output: Differences found

#### Frame Calculation
**Tool**: `tools/compare_frame_calculation.cpp`
- Compare template selection
- Compare atom matching
- Compare RMS fit values
- Compare final frames
- Output: Differences found

#### H-bond Detection
**Tool**: `tools/compare_hbond_detection.cpp`
- Compare initial H-bond finding
- Compare conflict resolution
- Compare validation
- Compare final H-bond list
- Output: Differences found

#### Pair Validation
**Tool**: `tools/compare_pair_validation.cpp`
- Compare geometric calculations
- Compare validation checks
- Compare H-bond counting
- Compare quality scores
- Output: Differences found

#### Pair Selection
**Tool**: `tools/compare_pair_selection.cpp`
- Compare best partner selection
- Compare mutual matching
- Compare iteration order
- Output: Differences found

### Benefits
- Focused debugging
- Easy to identify root causes
- Can fix issues one at a time

---

## Strategy 4: Incremental Testing

### Approach
Test each algorithm in isolation before full pipeline:

1. **Stage 1: PDB Parsing**
   - Compare atom/residue counts
   - Fix any differences
   - Verify 100% match

2. **Stage 2: Frame Calculation**
   - Compare frames for each residue
   - Fix any differences
   - Verify 100% match

3. **Stage 3: H-bond Detection**
   - Compare H-bonds for test pairs
   - Fix any differences
   - Verify 100% match

4. **Stage 4: Pair Validation**
   - Compare validation for test pairs
   - Fix any differences
   - Verify 100% match

5. **Stage 5: Pair Selection**
   - Compare selected pairs
   - Fix any differences
   - Verify 100% match

### Benefits
- Catch issues early
- Fix one stage at a time
- Verify 100% match at each stage
- Easier to debug

---

## Strategy 5: Focused Test Cases

### Create Test Cases for Specific Issues

#### Test Case 1: 3KNC Residue Recognition
**Goal**: Understand why only 16/66 residues are recognized  
**Test**: Extract first 20 residues, compare residue identification  
**Expected**: All 20 should be recognized as nucleotides

#### Test Case 2: 5UJ2 Missing Residue
**Goal**: Understand why residue 2 is missing  
**Test**: Extract residues 1-10, compare parsing  
**Expected**: Residue 2 should be found

#### Test Case 3: 3G8T Missing Pair (92, 160)
**Goal**: Understand why pair (92, 160) is missing  
**Test**: Extract residues 85-100, 155-170, compare validation  
**Expected**: Pair should be found and validated

#### Test Case 4: Quality Score Differences
**Goal**: Understand quality score calculation differences  
**Test**: Compare quality scores for mismatched pairs  
**Expected**: Quality scores should match exactly

### Benefits
- Focused debugging
- Faster iteration
- Easier to understand issues

---

## Recommended Workflow

### For Each Problematic PDB:

1. **Create Subset**: Extract relevant residues
2. **Run Legacy Test**: Use focused test executable
3. **Run Modern Test**: Use modern equivalent
4. **Compare Outputs**: Use comparison tool
5. **Identify Differences**: Find root cause
6. **Fix Issue**: Update modern code
7. **Verify Fix**: Re-run comparison
8. **Test on Full PDB**: Verify fix works on full dataset

### Example: Fixing 3G8T Missing Pair (92, 160)

1. **Create Subset**:
   ```bash
   create_pdb_subset 3G8T.pdb 3G8T_subset.pdb --residues 85-100,155-170
   ```

2. **Run Legacy Test**:
   ```bash
   org/bin/test_pair_validation 3G8T_subset.pdb 92 160 > legacy_validation.json
   ```

3. **Run Modern Test**:
   ```bash
   tools/test_pair_validation 3G8T_subset.pdb 92 160 > modern_validation.json
   ```

4. **Compare**:
   ```bash
   tools/compare_pair_validation legacy_validation.json modern_validation.json
   ```

5. **Identify Differences**: Check output for mismatches

6. **Fix Issue**: Update modern code

7. **Verify**: Re-run comparison, should match

8. **Test on Full PDB**: Run full pipeline on 3G8T, verify pair is found

---

## Tools to Create

### High Priority
1. ✅ `tools/create_pdb_subset.cpp` - Extract residues from PDB
2. ✅ `tools/compare_pdb_parsing.cpp` - Compare parsing
3. ✅ `tools/compare_quality_scores.cpp` - Compare quality scores
4. ✅ `org/bin/test_pair_validation` - Legacy validation test
5. ✅ `org/bin/test_quality_score` - Legacy quality score test

### Medium Priority
1. `tools/compare_frame_calculation.cpp` - Compare frame calculation
2. `tools/compare_hbond_detection.cpp` - Compare H-bond detection
3. `org/bin/test_pdb_parsing` - Legacy parsing test
4. `org/bin/test_frame_calculation` - Legacy frame test
5. `org/bin/test_hbond_detection` - Legacy H-bond test

### Low Priority
1. `scripts/generate_test_cases.py` - Auto-generate test cases
2. `tools/compare_pair_selection.cpp` - Compare pair selection

---

## Success Metrics

### Current
- Match Rate: 99.5%
- Perfect PDBs: 94/100

### Target
- Match Rate: 100%
- Perfect PDBs: 100/100

### Progress Tracking
- Track match rate improvement weekly
- Track number of perfect PDBs
- Track number of issues resolved
- Document any new issues found

---

## Next Steps

1. **Start with H-bond Refactoring** (Week 1)
   - Clean up code organization
   - No functional changes expected

2. **Investigate Residue Inclusion** (Week 2)
   - Create comparison tools
   - Fix 3KNC and 5UJ2

3. **Investigate Quality Scores** (Week 3)
   - Create comparison tools
   - Fix 3G8T and 6CAQ

4. **Create Legacy Test Executables** (Week 4)
   - Enable focused testing
   - Improve debugging capabilities

5. **Iterate Until 100% Match** (Weeks 5-7)
   - Fix remaining issues
   - Verify 100% match


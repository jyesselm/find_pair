# Difference Investigation Plan

## Executive Summary

Based on analysis of **50 PDBs**, we've identified the major patterns in differences between legacy and modern base pair detection.

## Statistics (50 PDBs Analyzed)

| Category | Count | Percentage | Notes |
|----------|-------|------------|-------|
| **Perfect Matches** | 13 | 26% | Pairs match exactly |
| **Modern Finds More** | 15 | 30% | Modern finds additional pairs |
| **Both Different** | 7 | 14% | Both found pairs but sets differ |
| **Legacy Finds None** | 15 | 30% | Modern finds pairs, legacy doesn't |
| **Legacy Finds More** | 0 | 0% | Rare occurrence |

**Key Finding**: Modern code is more comprehensive, finding pairs in 30% of cases where legacy found none, and finding additional pairs in 30% of cases where both found some.

## Investigation Priorities

### Priority 1: Verify Parameters for Perfect Matches ✅

**Goal**: Ensure that when pairs match exactly, all calculated parameters also match.

**Test Cases**:
- 6V9Q (7 pairs)
- 7EH2 (24 pairs)  
- 1A34 (9 pairs)
- 1AV6, 1B2M, 1BMV, 1C9S (0 pairs - no structures)

**Actions**:
1. Run parameter comparison script on all perfect matches
2. Compare step parameters (Shift, Slide, Rise, Tilt, Roll, Twist)
3. Compare helical parameters (X-disp, Y-disp, h-Rise, Inclination, Tip, h-Twist)
4. Verify numerical accuracy (tolerance: 0.01)

**Status**: Script created, ready to run

### Priority 2: Investigate Why Modern Finds More Pairs

**Goal**: Understand why modern code finds additional pairs that legacy misses.

**Top Priority Cases** (all have +10 additional pairs):
- **157D**: Legacy=2, Modern=12 (+10)
- **1BNA**: Legacy=2, Modern=12 (+10)
- **1CSL**: Legacy=0, Modern=10 (+10)
- **1DFU**: Legacy=0, Modern=10 (+10)
- **1DUQ**: Legacy=0, Modern=10 (+10)

**Average Additional Pairs**: 7.1

**Investigation Steps**:

1. **Validation Thresholds Comparison**
   - Compare `ValidationParameters` between legacy and modern
   - Check distance thresholds (dorg, dNN, hb_dist, etc.)
   - Check angle thresholds (plane_angle, etc.)
   - Check overlap thresholds

2. **Hydrogen Bond Detection**
   - Compare hydrogen bond detection algorithms
   - Check minimum hydrogen bond requirements
   - Verify hydrogen bond distance calculations

3. **Frame Calculation Differences**
   - Compare reference frame calculations
   - Check template matching
   - Verify RMS fitting thresholds

4. **PDB Parsing Differences**
   - Compare how residues are identified
   - Check heteroatom handling
   - Verify chain handling

**Expected Outcome**: Identify specific threshold or algorithm differences that cause modern to find more pairs.

### Priority 3: Investigate Cases Where Legacy Found None

**Goal**: Understand why modern finds pairs when legacy finds none.

**Cases** (15 PDBs):
- 100D, 161D, 165D, 168D (10 pairs each)
- 1D87, 1D88, 1D96, 1D9H, 1DI2, 1DNO (various counts)

**Investigation Steps**:

1. **Check if pairs are valid**
   - Visual inspection in PyMOL
   - Verify hydrogen bond distances
   - Check geometry

2. **Compare PDB parsing**
   - Check if legacy failed to parse certain residues
   - Verify chain identification
   - Check for parsing errors

3. **Compare initial filtering**
   - Check if legacy filters out structures early
   - Compare residue type detection
   - Verify nucleotide identification

**Expected Outcome**: Determine if modern's additional pairs are valid or if legacy's filtering is too strict.

### Priority 4: Investigate Both Found, Different Sets

**Goal**: Understand why both codes find pairs but the sets differ.

**Cases** (7 PDBs):
- 1ASY: Legacy=49, Modern=48 (-1 legacy-only)
- 1ASZ: Legacy=47, Modern=48 (+1 modern-only)
- 1B23: Legacy=20, Modern=26 (+6 modern-only)
- 1C0A: Legacy=19, Modern=25 (+6 modern-only)
- 1D4R: Various differences

**Investigation Steps**:

1. **Identify specific differing pairs**
   - Extract pairs only in legacy
   - Extract pairs only in modern
   - Compare their validation scores

2. **Check edge cases**
   - Pairs near validation thresholds
   - Pairs with unusual geometry
   - Pairs with modified nucleotides

3. **Compare conflict resolution**
   - Check hydrogen bond conflict resolution
   - Verify overlap calculations
   - Compare tie-breaking algorithms

**Expected Outcome**: Identify specific pairs and reasons for differences.

## Action Items

### Immediate (Next Session)

1. ✅ **Run parameter comparison on perfect matches**
   - Verify all parameters match for 6V9Q, 7EH2, 1A34
   - Document any numerical differences

2. **Deep dive on top 5 "modern finds more" cases**
   - Run detailed validation comparison
   - Check specific pairs modern found
   - Verify if pairs are valid

3. **Analyze "legacy found none" cases**
   - Check PDB parsing differences
   - Verify if modern's pairs are valid

### Medium Term

4. **Compare validation thresholds systematically**
   - Extract all thresholds from legacy code
   - Compare with modern defaults
   - Document differences

5. **Create detailed comparison for each category**
   - Write scripts to extract validation details
   - Compare algorithms step-by-step
   - Document findings

### Long Term

6. **Make thresholds configurable**
   - Allow matching legacy behavior exactly
   - Add flag for "legacy mode" validation
   - Test with legacy thresholds

7. **Validate modern improvements**
   - Determine if additional pairs found by modern are valid
   - Decide if we want to match legacy exactly or keep improvements

## Tools Created

1. **`compare_parameters_for_matching_pairs.py`**
   - Verifies parameters match for perfect pair matches
   - Extracts and compares step and helical parameters

2. **`analyze_differences_plan.py`**
   - Categorizes differences
   - Identifies priority cases
   - Generates investigation plan

3. **`batch_compare.py`**
   - Batch comparison across multiple PDBs
   - Generates statistics

4. **`compare_side_by_side.py`**
   - Detailed side-by-side comparison
   - Shows parameters

## Next Steps

1. Run parameter verification: `python3 scripts/compare_parameters_for_matching_pairs.py 6V9Q 7EH2 1A34`
2. Deep dive on 157D: Why does modern find 10 more pairs?
3. Analyze 1ASY: Why did legacy find 1 more pair (49 vs 48)?
4. Compare validation thresholds systematically


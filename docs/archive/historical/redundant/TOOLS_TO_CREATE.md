# Tools to Create for 100% Match

## Existing Tools ✅

These tools already exist and are built:
- `compare_quality_scores` - Compare quality scores between legacy and modern
- `compare_pdb_parsing` - Compare PDB parsing (atom/residue counts)
- `compare_residue_identification` - Compare nucleotide recognition
- `compare_hbond_detection` - Compare H-bond detection
- `debug_donor_acceptor` - Debug donor/acceptor type determination
- `detect_hbonds_standalone` - Detect H-bonds independently
- `list_all_hbonds` - List all H-bonds from JSON files
- `debug_frame_calculation` - Debug frame calculation

---

## High Priority Tools to Create

### 1. `create_pdb_subset.cpp` ⭐ HIGHEST PRIORITY

**Purpose**: Extract specific residues/chains from a PDB to create smaller test cases for debugging.

**Usage**:
```bash
build/create_pdb_subset input.pdb output.pdb --residues A:90-95,A:155-165
build/create_pdb_subset input.pdb output.pdb --chain A --start 90 --end 95
build/create_pdb_subset input.pdb output.pdb --around-pair A:92 A:160 --context 10
```

**Features**:
- Extract specific residues by chain:seq_num
- Extract around a specific pair with context (neighboring residues)
- Preserve atom numbering and chain information
- Useful for isolating problematic pairs

**Use Cases**:
- Extract residues around 3G8T pair (92, 160) for focused debugging
- Extract first 20 residues from 3KNC to debug residue recognition
- Create minimal test cases for specific issues

**Estimated Time**: 2-3 hours

---

### 2. `analyze_mismatched_pairs.cpp` ⭐ HIGH PRIORITY

**Purpose**: Automatically analyze all mismatched pairs and generate a comprehensive report.

**Usage**:
```bash
build/analyze_mismatched_pairs --pdb 3G8T --legacy-json data/json_legacy/3G8T.json --modern-json data/json/3G8T.json
build/analyze_mismatched_pairs --test-set 100 --output report.md
```

**Features**:
- Identify all missing pairs (in legacy but not modern)
- Identify all extra pairs (in modern but not legacy)
- For each mismatched pair:
  - Check if it exists in validation records
  - Compare quality scores (base, adjust_pairQuality, bp_type_id, final)
  - Compare H-bond counts and types
  - Compare geometric values (dorg, d_v, plane_angle, dNN)
  - Check if it's a tie-breaking issue
- Generate markdown report with prioritized issues

**Output Format**:
```markdown
# Mismatched Pairs Analysis: 3G8T

## Missing Pairs (2)

### Pair (92, 160)
- **Status**: ✅ EXISTS in validation (is_valid=1)
- **Quality Score**: Legacy=7.64, Modern=10.64 (diff=3.0)
- **adjust_pairQuality**: Legacy=-3.0, Modern=0.0 ⚠️ DIFFERENCE
- **H-bonds**: Legacy=2 good, Modern=0 good ⚠️ DIFFERENCE
- **Root Cause**: H-bond detection difference (likely fixed after H-bond fix)
- **Priority**: HIGH - Verify after regeneration

### Pair (946, 947)
- **Status**: ✅ EXISTS in validation (is_valid=1)
- **Quality Score**: Legacy=12.34, Modern=12.35 (diff=0.01)
- **Root Cause**: Tie-breaking or iteration order
- **Priority**: MEDIUM
```

**Estimated Time**: 4-6 hours

---

### 3. `verify_hbond_fix.cpp` ⭐ IMMEDIATE PRIORITY

**Purpose**: Quick verification tool to check if H-bond fix worked for specific pairs.

**Usage**:
```bash
build/verify_hbond_fix data/pdb/3G8T.pdb 92 160
build/verify_hbond_fix data/pdb/3G8T.pdb --json data/json/3G8T.json
```

**Features**:
- Check specific pair for H-bond types
- Verify `adjust_pairQuality` value
- Compare with legacy if JSON provided
- Quick pass/fail output

**Output**:
```
Verifying H-bond fix for 3G8T pair (92, 160):
  Residue 92: A2M (modified nucleotide)
  Residue 160: G (standard)
  
  H-bonds found:
    N1 -> O2': type='*', dist=2.675 ✅
    N3 -> N2: type='-', dist=3.258 ✅ FIXED (was '*' before)
  
  adjust_pairQuality: -3.0 ✅ (was 0.0 before)
  Status: ✅ FIX VERIFIED
```

**Estimated Time**: 1-2 hours

---

### 4. `compare_pair_validation.cpp` ⭐ HIGH PRIORITY

**Purpose**: Compare full pair validation results between legacy and modern for specific pairs.

**Usage**:
```bash
build/compare_pair_validation data/pdb/3G8T.pdb 92 160 --legacy-json data/json_legacy/3G8T.json
build/compare_pair_validation data/pdb/3G8T.pdb --pairs 92:160,946:947
```

**Features**:
- Compare all validation checks (cdns, overlap, hbond_check, etc.)
- Compare geometric values (dorg, d_v, plane_angle, dNN)
- Compare H-bond counts and types
- Compare quality scores (all components)
- Side-by-side comparison output

**Output Format**:
```
Pair Validation Comparison: 3G8T (92, 160)
===========================================

Geometric Values:
  dorg:        Legacy=8.234, Modern=8.234 ✅ MATCH
  d_v:         Legacy=0.123, Modern=0.123 ✅ MATCH
  plane_angle: Legacy=12.5,  Modern=12.5  ✅ MATCH
  dNN:         Legacy=10.2,   Modern=10.2  ✅ MATCH

Validation Checks:
  cdns:        Legacy=1, Modern=1 ✅ MATCH
  overlap:     Legacy=1, Modern=1 ✅ MATCH
  hbond_check: Legacy=1, Modern=1 ✅ MATCH
  is_valid:    Legacy=1, Modern=1 ✅ MATCH

H-bonds:
  num_base_hb: Legacy=2, Modern=2 ✅ MATCH
  num_o2_hb:   Legacy=0, Modern=0 ✅ MATCH
  good_hb:     Legacy=2, Modern=2 ✅ MATCH

Quality Scores:
  base_score:        Legacy=10.64, Modern=10.64 ✅ MATCH
  adjust_pairQuality: Legacy=-3.0, Modern=-3.0  ✅ MATCH
  bp_type_id:        Legacy=0,    Modern=0     ✅ MATCH
  final_score:       Legacy=7.64, Modern=7.64  ✅ MATCH

Selection:
  Legacy: ✅ SELECTED
  Modern: ❌ NOT SELECTED
  Issue:  Quality score difference (likely tie-breaking)
```

**Estimated Time**: 3-4 hours

---

### 5. `batch_verify_fix.cpp` ⭐ MEDIUM PRIORITY

**Purpose**: Batch verify H-bond fix across multiple PDBs.

**Usage**:
```bash
build/batch_verify_fix --test-set 100 --output fix_verification.json
build/batch_verify_fix --pdbs 3G8T 6CAQ --check-modified-nucleotides
```

**Features**:
- Check all PDBs for modified nucleotides
- Verify H-bond types for pairs with modified nucleotides
- Compare `adjust_pairQuality` before/after fix
- Generate summary report

**Estimated Time**: 2-3 hours

---

## Medium Priority Tools

### 6. `find_quality_score_differences.cpp`

**Purpose**: Find all pairs with quality score differences (automated version of compare_quality_scores).

**Usage**:
```bash
build/find_quality_score_differences --pdb 3G8T --threshold 0.01
```

**Features**:
- Compare all pairs between legacy and modern
- Identify pairs with quality score differences above threshold
- Sort by difference magnitude
- Generate prioritized list

**Estimated Time**: 2-3 hours

---

### 7. `debug_tie_breaking.cpp`

**Purpose**: Debug tie-breaking issues when quality scores are equal.

**Usage**:
```bash
build/debug_tie_breaking data/pdb/3G8T.pdb --legacy-json data/json_legacy/3G8T.json
```

**Features**:
- Find all pairs with equal quality scores
- Compare selection order between legacy and modern
- Check iteration order differences
- Identify tie-breaking discrepancies

**Estimated Time**: 2-3 hours

---

## Scripts to Create

### 8. `scripts/generate_test_cases.py` ⭐ HIGH PRIORITY

**Purpose**: Automatically generate test cases around problematic pairs.

**Usage**:
```python
python3 scripts/generate_test_cases.py --pdb 3G8T --pairs 92:160,946:947 --context 10
python3 scripts/generate_test_cases.py --from-report comparison_report.md
```

**Features**:
- Extract residues around problematic pairs
- Generate both legacy and modern JSON for subsets
- Create focused test cases for debugging
- Integrate with `create_pdb_subset`

**Estimated Time**: 2-3 hours

---

### 9. `scripts/quick_fix_verification.py` ⭐ IMMEDIATE PRIORITY

**Purpose**: Quick script to verify H-bond fix worked.

**Usage**:
```bash
python3 scripts/quick_fix_verification.py 3G8T 92 160
python3 scripts/quick_fix_verification.py --regenerate-and-check 3G8T
```

**Features**:
- Regenerate modern JSON if needed
- Check specific pair for H-bond fix
- Compare with legacy
- Quick pass/fail output

**Estimated Time**: 1 hour

---

## Priority Order

1. **IMMEDIATE** (for Step 1: Verify H-bond Fix):
   - `verify_hbond_fix.cpp` or `scripts/quick_fix_verification.py`
   
2. **HIGH** (for Step 2: Quality Score Investigation):
   - `analyze_mismatched_pairs.cpp`
   - `compare_pair_validation.cpp`
   
3. **HIGH** (for Step 3: Residue Inclusion):
   - `create_pdb_subset.cpp`
   - `scripts/generate_test_cases.py`
   
4. **MEDIUM** (for ongoing debugging):
   - `find_quality_score_differences.cpp`
   - `debug_tie_breaking.cpp`
   - `batch_verify_fix.cpp`

---

## Implementation Notes

- All tools should output JSON for programmatic use
- All tools should have human-readable output (markdown or formatted text)
- Tools should be fast and focused (single purpose)
- Tools should integrate with existing comparison infrastructure
- Consider adding `--verbose` and `--quiet` flags to all tools


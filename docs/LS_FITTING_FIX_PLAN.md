# LS_FITTING Fix Plan: Path to 100% Validation

**Date**: December 3, 2025  
**Current Status**: 97% Success Rate (3493/3602 PDBs pass or FP-only)  
**Goal**: Fix remaining 107 count mismatches (3%)

---

## Current Results Summary

| Category | Count | Percentage | Status |
|----------|-------|------------|--------|
| Perfect match | 869 | 24% | ✅ Perfect |
| FP differences only | 2624 | 73% | ✅ Acceptable (1e-6 precision) |
| **Count mismatches** | **107** | **3%** | ⚠️ Need fixing |
| Errors/Skips | 2 | <0.1% | ⚠️ Minor |
| **Total Success** | **3493** | **97%** | ✅ Very good! |

---

## Count Mismatch Analysis

### Category 1: Modern Missing Residues (47 PDBs - Legacy > Modern)

**Pattern**: Modern code is **rejecting** certain modified nucleotides that legacy accepts.

**Examples**:
- `1FIR`: Missing **70U** (7-deoxyuridine, 2-thio modification, has S2 instead of O2)
- `1OB2`: Missing **KIR** (unknown modified base)
- `2G92`, `2Q1O`: Missing **NF2** (appears in multiple PDBs)
- `2XD0` and others: Missing **A23** (2-aminoadenine)

**Root Cause**: 
- These modified nucleotides have non-standard atoms (S2, modified bases)
- Modern RMSD check (threshold 0.2618) rejects them as "not nucleotides"
- Legacy is more permissive or uses different criteria

**Statistics**:
- Total: 47 PDBs
- Avg missing: 1.9 residues per PDB
- Max missing: 8 residues
- Total missing residues: ~89

**Modified Nucleotides Found**:
- `70U` - 7-deoxyuridine (2-thio-uridine with S2)
- `KIR` - Unknown
- `NF2` - Unknown (appears frequently)
- `A23` - 2-aminoadenine (modified adenine)
- Others to be cataloged

### Category 2: Modern Has Extra Residues (60 PDBs - Modern > Legacy)

**Pattern**: Modern code is **including** residues that legacy doesn't.

**Statistics**:
- Total: 60 PDBs  
- Avg extra: 8.8 residues per PDB
- Max extra: 39 residues
- Total extra residues: ~528

**Possible Causes**:
1. Legacy has stricter filtering criteria
2. Legacy might skip certain chain types
3. Deduplication issue in comparison
4. Legacy JSON files incomplete

**Needs Investigation**: Check a few examples to understand pattern.

---

## Fix Strategy

### Phase 1: Investigate Missing Modified Nucleotides ⏳ NEXT

**Goal**: Understand why modern rejects 70U, NF2, A23, KIR, etc.

**Actions**:
1. ✅ Already added `70U` to modified_nucleotides list (line 344 in pdb_parser.cpp)
2. ⏳ Add other modified bases: `NF2`, `KIR`, `A23`
3. ⏳ Check if RMSD threshold is too strict for these
4. ⏳ Consider relaxing RMSD for known modified nucleotides

**Investigation Commands**:
```bash
# Check what atoms each modified base has
for pdb in 1FIR 1OB2 2G92 2Q1O 2XD0; do
    echo "=== $pdb ===" 
    grep -E "70U|KIR|NF2|A23" data/pdb/${pdb}.pdb | head -5
done

# Test with added modified bases
# After adding to modified_nucleotides list, rebuild and test
```

### Phase 2: Investigate Extra Residues in Modern ⏳ NEXT

**Goal**: Understand why modern has MORE residues than legacy in 60 cases.

**Actions**:
1. ⏳ Pick 3-5 examples from Category 2
2. ⏳ Compare residue lists to see what's extra in modern
3. ⏳ Check if legacy is filtering certain residues (amino acids, waters, etc.)
4. ⏳ Verify deduplication is working correctly

**Investigation Commands**:
```bash
# Analyze a specific case
python << 'EOF'
import json
pdb_id = "1H3E"  # 76 legacy vs 82 modern
legacy = json.load(open(f"data/json_legacy/ls_fitting/{pdb_id}.json"))
modern = json.load(open(f"tmp/validation/ls_fitting/{pdb_id}.json"))

# Compare residue lists
# ... (script to find differences)
EOF
```

### Phase 3: Add Missing Modified Nucleotides 🔧 TO DO

**Based on findings, add modified bases to the list:**

File: `src/x3dna/io/pdb_parser.cpp` line ~344

```cpp
// Modified Uracil/Thymine
"5MU", "H2U", "DHU", "OMU", "4SU", "S4U", "5BU", "2MU", "UR3", "RT", 
"70U",  // ✅ ADDED
"NF2",  // TO ADD
"KIR",  // TO ADD
// Modified Adenine  
"A2M", "1MA", "2MA", "6MA", "OMA", "MIA", "I6A", "T6A", "M6A",
"A23",  // TO ADD
```

### Phase 4: Relax RMSD Threshold for Modified Bases (If Needed) 🔧 OPTIONAL

**If modified bases still fail RMSD check:**

File: `src/x3dna/algorithms/base_frame_calculator.cpp` line ~298

**Option A**: Relax threshold for specific modified bases
```cpp
// Check if this is a known modified base with expected differences
bool is_special_modified = (res_name == "70U" || res_name == "NF2" || ...);
double threshold = is_special_modified ? 0.5 : 0.2618;  // More lenient for modified

if (!rmsd_result.has_value() || *rmsd_result > threshold) {
    // Reject
}
```

**Option B**: Skip RMSD check for recognized modified nucleotides
```cpp
if (is_modified_nucleotide && has_ring_atoms && has_nitrogen) {
    // Accept without RMSD check if it's in our modified list
    has_valid_rmsd = true;
}
```

### Phase 5: Fix Category 2 Cases (Modern > Legacy) 🔧 TO DO

**After investigating, apply appropriate fix:**

Possible fixes:
- If legacy is incorrectly filtering: Document as legacy bug, modern is correct
- If modern is incorrectly including: Add filtering to match legacy
- If deduplication issue: Fix deduplication logic

### Phase 6: Rebuild and Revalidate ✅ FINAL

**After all fixes:**

```bash
# Rebuild
cd build && ninja

# Revalidate
python scripts/run_ls_fitting_validation.py

# Should see:
# - Count mismatches: 0
# - Perfect match + FP only: 100%
```

---

## Detailed Investigation Plan

### Step 1: Catalog All Modified Nucleotides in Category 1

```bash
# Extract all missing modified bases
python << 'EOF'
import json
from pathlib import Path

missing_bases = {}

with open("data/validation_results/legacy_more.txt") as f:
    for line in f:
        pdb_id = line.split('\t')[0]
        
        legacy_json = Path(f"data/json_legacy/ls_fitting/{pdb_id}.json")
        modern_json = Path(f"tmp/validation/ls_fitting/{pdb_id}.json")
        
        if legacy_json.exists() and modern_json.exists():
            legacy = json.load(open(legacy_json))
            modern = json.load(open(modern_json))
            
            # Deduplicate
            seen = set()
            unique_legacy = []
            for rec in legacy:
                key = (rec.get('chain_id'), rec.get('residue_seq'), rec.get('residue_name', '').strip())
                if key not in seen:
                    seen.add(key)
                    unique_legacy.append(rec)
            
            # Find missing
            legacy_res = {r['residue_name'].strip() for r in unique_legacy}
            modern_res = {r['residue_name'].strip() for r in modern}
            
            missing = legacy_res - modern_res
            for base in missing:
                if base not in missing_bases:
                    missing_bases[base] = []
                missing_bases[base].append(pdb_id)

print("Modified bases missing in modern:")
for base in sorted(missing_bases.keys()):
    pdbs = missing_bases[base]
    print(f"  {base}: {len(pdbs)} PDBs - {', '.join(pdbs[:5])}")
EOF
```

### Step 2: Check PDB Files for These Bases

```bash
# For each missing base, check its structure
for base in 70U KIR NF2 A23; do
    echo "=== $base ==="
    grep "HETATM.*$base" data/pdb/*.pdb | head -20 | grep -E "N1|N3|N9"
done
```

### Step 3: Add to Modified Nucleotides List

Update `src/x3dna/io/pdb_parser.cpp` to include all found modified bases.

### Step 4: Test Each Fix

```bash
# After each addition, rebuild and test
cd build && ninja

# Test specific PDB
python tests_python/integration/test_ls_fitting.py <PDB_ID>
```

### Step 5: Investigate Category 2 (Modern > Legacy)

```python
# Script to find what modern includes that legacy doesn't
# Check if these are valid nucleotides or amino acids being incorrectly included
```

---

## Expected Timeline

| Phase | Time Estimate | Status |
|-------|--------------|--------|
| 1. Catalog modified bases | 15 min | ⏳ Ready to start |
| 2. Add to modified list | 5 min | ⏳ Waiting |
| 3. Test Category 1 fixes | 10 min | ⏳ Waiting |
| 4. Investigate Category 2 | 30 min | ⏳ Waiting |
| 5. Apply Category 2 fixes | 15 min | ⏳ Waiting |
| 6. Final validation | 15 min | ⏳ Waiting |
| **Total** | **~90 min** | |

---

## Success Criteria

✅ **100% validation**:
- All 3602 PDBs either perfect match or FP differences only
- Count mismatches: 0
- All modified nucleotides correctly recognized

---

## Modified Nucleotides to Add (Preliminary List)

Based on Category 1 analysis:

```cpp
// In src/x3dna/io/pdb_parser.cpp
static const std::vector<std::string> modified_nucleotides = {
    // ... existing ...
    
    // Additional modified bases found in validation:
    "70U",  // 7-deoxyuridine (2-thio modification)
    "NF2",  // TODO: Identify
    "KIR",  // TODO: Identify  
    "A23",  // 2-aminoadenine
    // ... add more as found ...
};
```

---

## Next Immediate Actions

1. ✅ Run complete validation (DONE - 97% success)
2. ✅ Categorize mismatches (DONE - 47 missing, 60 extra)
3. ⏳ Catalog all missing modified bases
4. ⏳ Add them to modified_nucleotides list
5. ⏳ Rebuild and retest
6. ⏳ Investigate Category 2 cases
7. ⏳ Apply final fixes
8. ⏳ Achieve 100%!

---

*The path to 100% is clear: systematically add the ~5-10 missing modified nucleotides and investigate why modern has extra residues in some cases. Most of the work is done - just need to catalog and add the modified bases!*


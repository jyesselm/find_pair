# Investigation: 7EH2 Index Mismatch

**Date**: December 2, 2025  
**Status**: 🔴 **REQUIRES INVESTIGATION**

---

## Problem

During comprehensive validation, PDB **7EH2** failed with a count mismatch:
- **Modern code**: 88 nucleotides
- **Legacy code**: 48 nucleotides  
- **Difference**: 40 nucleotides (45% missing in legacy!)

---

## Detailed Breakdown

| Residue Type | Modern | Legacy | Difference | Notes |
|--------------|--------|--------|------------|-------|
| DA (DNA Adenine) | 20 | 8 | +12 | Legacy missing 60% |
| DC (DNA Cytosine) | 18 | 16 | +2 | Legacy missing 11% |
| DG (DNA Guanine) | 26 | 12 | +14 | Legacy missing 54% |
| DT (DNA Thymine) | 20 | 8 | +12 | Legacy missing 60% |
| G (RNA Guanine) | 4 | 4 | ✅ 0 | Perfect match |
| **TOTAL** | **88** | **48** | **+40** | **45% missing** |

### Pattern

- **RNA residues**: Perfect match (4 G residues)
- **DNA residues**: Legacy finds about half of what modern finds
- DNA purines (DA, DG): ~60% missing in legacy
- DNA pyrimidines (DC, DT): ~11-60% missing in legacy

---

## Matched Residues

The 48 residues that DO match have extremely high legacy indices:
- G17 DA: legacy_idx=**6995**
- G21 DA: legacy_idx=**6999**
- G16 DC: legacy_idx=**6994**

This suggests the PDB has ~7000 total residues, with most being non-nucleotides (proteins, waters, etc.).

---

## Files for Investigation

- **PDB file**: `data/pdb/7EH2.pdb`
- **Modern JSON**: `data/json/base_frame_calc/7EH2.json`
- **Legacy JSON**: `data/json_legacy/base_frame_calc/7EH2.json`
- **Index mapping**: `data/index_mapping/7EH2.json`

---

## Possible Causes

### 1. Modern Code is Correct (Legacy has bug)

**Evidence**:
- Modern consistently finds DNA residues
- Pattern is systematic (not random)
- RNA residues match perfectly

**Possible legacy bug**:
- Legacy might be filtering DNA differently than RNA
- Chain ID filtering issue?
- Insertion code handling?
- Multiple models in PDB?

### 2. Legacy Code is Correct (Modern has bug)

**Evidence**:
- Legacy has been validated on thousands of structures
- Matched indices are consistent

**Possible modern bug**:
- Modern might be accepting invalid nucleotides
- Template matching too permissive?
- Reference frame calculation accepting partial atoms?

### 3. Structural Issue with 7EH2

**Possible PDB issues**:
- Multiple models (NMR structure?)
- Alternate conformations
- Modified nucleotides
- Hybrid structure (RNA/DNA)
- Non-standard residue names

---

## Investigation Steps

### Step 1: Examine the PDB File

```bash
# Check basic stats
grep "^ATOM" data/pdb/7EH2.pdb | wc -l
grep "^HETATM" data/pdb/7EH2.pdb | wc -l

# Check for models
grep "^MODEL" data/pdb/7EH2.pdb

# Count DNA residues manually
grep "^ATOM" data/pdb/7EH2.pdb | awk '{print $4}' | grep "^D[ACGT]$" | sort | uniq -c

# Check chains
grep "^ATOM" data/pdb/7EH2.pdb | awk '{print $5}' | sort | uniq -c
```

### Step 2: Compare Specific Residues

Check which DNA residues modern finds but legacy doesn't:
```bash
python3 scripts/analyze_index_mismatches.py --pdb 7EH2
```

Look at the mapping file to see which specific residues are missing:
```bash
cat data/index_mapping/7EH2.json | jq '.[] | select(.legacy_index == -1)'
```

### Step 3: Check Modern Processing

Run modern code with debug output:
```bash
# Check what modern code sees
./build/generate_modern_json data/pdb/7EH2.pdb data/json 2>&1 | grep -i "nucleotide\|residue"
```

### Step 4: Check Legacy Processing

Run legacy code:
```bash
cd org
./build/bin/find_pair_original ../data/pdb/7EH2.pdb test_output.inp
```

### Step 5: Manual Verification

Pick one of the "modern only" residues and check:
1. Does it exist in the PDB?
2. Does it have all required atoms?
3. Why might legacy skip it?

For example, check G2 DA (modern_idx=0, legacy_idx=-1):
```bash
grep "^ATOM.*\s2\s.*DA" data/pdb/7EH2.pdb
```

---

## Decision Tree

```
Is 7EH2 a multi-model PDB?
├─ YES → Modern needs to filter to single model
└─ NO → Continue

Do "modern only" residues have all required atoms?
├─ NO → Modern template matching too permissive
└─ YES → Continue

Are "modern only" residues in specific chains?
├─ YES → Check chain filtering logic
└─ NO → Continue

Are residues actually DNA/RNA nucleotides?
├─ NO → Modern misidentifying residues
└─ YES → Legacy has bug, needs fix
```

---

## Expected Outcome

One of:
1. **Modern bug found** → Fix modern code → Re-run validation from 7EH2
2. **Legacy bug found** → Document as known difference → Continue validation
3. **PDB issue** → Add to exclusion list → Continue validation

---

## Current Status

**Validation stopped at**: Batch 4, PDB 344 of 1737  
**Progress**: 2,729 / 4,123 PDBs (66%)  
**Results so far**:
- ✅ PASS: 2,464
- ⏭️ SKIP: 264
- ❌ FAIL: 1 (7EH2)

**Resume command** (after fix):
```bash
python3 scripts/validate_all_indices.py --start-from 7EH2
```

---

## References

- Index mapping: `data/index_mapping/7EH2.json`
- Analysis script: `scripts/analyze_index_mismatches.py`
- Validation status: `data/index_validation_status.csv`


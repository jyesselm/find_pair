# 3KNC Debug Report

**Last Updated**: 2025-11-28  
**Status**: ⏳ Pending Investigation  
**Overall Match**: Unknown (residue recognition issues)

---

## Summary

| Metric | Value | Target | Notes |
|--------|-------|--------|-------|
| Residues recognized | 16/66 | 66/66 | Only 24% recognized |
| Nucleotide frames | ? | All | Need investigation |
| Pairs matched | ? | 100% | Blocked by residue issue |

---

## Primary Issue

### Residue Recognition Problem

**Problem**: Only 16 out of 66 residues are recognized as nucleotides.

**Impact**: 
- Missing ~50 residues from analysis
- Likely causing many missing pairs
- Prevents accurate comparison with legacy

**Possible Causes**:
1. **HETATM vs ATOM**: Some nucleotides may be labeled as HETATM
2. **Modified nucleotides**: Non-standard residue names not recognized
3. **PDB parsing issue**: Parser may be skipping some residues
4. **Chain handling**: Multi-chain handling may have issues

---

## Investigation Plan

### Step 1: Examine PDB File

```bash
# Check residue types in PDB
grep -E "^(ATOM|HETATM)" data/pdb/3KNC.pdb | awk '{print $4}' | sort | uniq -c

# Check chains
grep -E "^(ATOM|HETATM)" data/pdb/3KNC.pdb | awk '{print $5}' | sort | uniq -c

# Check legacy residue count
cat data/json_legacy/base_frame_calc/3KNC.json | python3 -c "import json, sys; d=json.load(sys.stdin); print(len(d.get('base_frame_calc', [])))"
```

### Step 2: Compare Residue Lists

```bash
# Generate modern residue list
./build/check_residue_indices data/pdb/3KNC.pdb

# Compare with legacy
python3 scripts/compare_json.py compare 3KNC --verbose
```

### Step 3: Identify Missing Residues

- Get list of residues in legacy JSON
- Get list of residues in modern JSON
- Find the difference
- Analyze why they're missing

---

## Hypotheses

### H1: HETATM Nucleotides

Some nucleotides may be recorded as HETATM instead of ATOM.
- Legacy may handle these differently
- Modern may be filtering them out

**Test**: Check PDB for HETATM nucleotides

### H2: Modified Nucleotides

3KNC may have modified nucleotides with non-standard names.
- Legacy may recognize them
- Modern may not have them in residue type list

**Test**: Check for non-standard residue names in PDB

### H3: Chain ID Issues

Multi-chain PDB may have parsing issues.
- Residues may be in different chains
- Chain ID handling may differ

**Test**: Check chain IDs for all 66 residues

---

## Data to Collect

- [ ] Total residues in PDB
- [ ] Residue types (A, G, C, U, modified)
- [ ] Chain IDs present
- [ ] ATOM vs HETATM counts
- [ ] Legacy residue count
- [ ] Modern residue count
- [ ] Specific missing residues

---

## Potential Fixes

Depending on root cause:

1. **HETATM handling**: Update PdbParser to include HETATM nucleotides
2. **Modified nucleotides**: Add non-standard names to residue type list
3. **Chain handling**: Fix chain ID parsing logic
4. **Residue filtering**: Adjust filtering criteria

---

## Related Documentation

- [IS_NUCLEOTIDE_BUG_ANALYSIS.md](../IS_NUCLEOTIDE_BUG_ANALYSIS.md) - Nucleotide detection logic
- [RESIDUE_INDEXING_ISSUE.md](../RESIDUE_INDEXING_ISSUE.md) - Residue indexing details
- [PDB_PROPERTIES_MATCHING_APPROACH.md](../PDB_PROPERTIES_MATCHING_APPROACH.md) - PDB matching

---

## Investigation Log

### 2025-11-28
- Created this debug report
- Issue known but not yet investigated
- Blocked behind 6CAQ investigation

---

## Next Steps

1. **Wait for 6CAQ completion** (higher priority)
2. **Run PDB analysis** using commands above
3. **Identify root cause** of missing residues
4. **Implement fix**
5. **Re-run comparison**

---

*Debug report for 3KNC - pending investigation of residue recognition issue.*


# LS_FITTING Known Issues

**Date**: December 3, 2025

---

## Issue #1: Modified Nucleotide 70U Not Included in Modern Output

**PDB**: 1FIR  
**Residue**: A34 (70U - 5-(O-METHYLACETO)-2-THIO-2-DEOXY-URIDINE)  
**Status**: ⚠️ INVESTIGATION NEEDED

### Symptom

- Legacy (fixed): 76 ls_fitting records (includes 70U)
- Modern: 75 ls_fitting records (excludes 70U)

### Details

**70U Structure:**
- Has all pyrimidine ring atoms: N1, C2, N3, C4, C5, C6 ✓
- Has C1' atom ✓
- Has nitrogen atoms (N1, N3) ✓
- **Has S2 (sulfur) instead of O2 (oxygen)** - 2-thio modification

**Modern Code Behavior:**
- 70U is correctly parsed (25 atoms including all ring atoms)
- 70U is added to modified_nucleotides list
- BUT: Frame calculation appears to fail (not included in ls_fitting output)

**Legacy Code Behavior:**
- Includes 70U in ls_fitting
- RMS fit: 0.016104 (very good fit)
- Num points: 6 (matched ring atoms)

### Possible Causes

1. **RMSD threshold check failing**: Modern code checks RMSD > 0.2618
   - 70U's S2 modification might cause RMSD to exceed threshold
   - Legacy might use different calculation or threshold

2. **Template matching**: No specific template for 2-thio-uridine
   - Both use standard U template
   - S2 vs O2 difference might affect fit differently

3. **Ring atom matching**: Modern code might be stricter about which atoms to include

### Impact

**Frequency**: Unknown - need to complete batch test to see how many PDBs have this issue

**Workaround Options**:
1. Add specific template for 2-thio-uridine modifications
2. Relax RMSD threshold for known modified nucleotides
3. Special handling for S2-containing uridines
4. Accept that modern is more strict (document difference)

### Next Steps

1. ✅ Document issue (this file)
2. ⏳ Complete batch test to quantify frequency
3. ⏳ Check if other modified nucleotides have similar issues
4. ⏳ Decide on fix strategy based on frequency and impact

---

## Related Files

- `src/x3dna/io/pdb_parser.cpp` - Modified nucleotides list (line 344)
- `src/x3dna/algorithms/base_frame_calculator.cpp` - RMSD threshold check (line 298)
- `src/x3dna/algorithms/base_pair_finder.cpp` - RMSD calculation (line 895)

---

*Status: Documented, awaiting batch test completion to determine frequency and severity*


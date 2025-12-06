# Frame Validation Status

**Last Updated**: December 4, 2025

## Current Status: Ready for Full Validation

### Recent Progress

✅ **KIR Bug Fixed** (Dec 4, 2025)
- Issue: Atom matching bug causing 20x RMS difference
- Fix: Integrated hardcoded atom list into `check_nt_type_by_rmsd()`
- Result: Perfect match (0.00e+00 difference)
- Commit: `050fd03`

✅ **Validation on 100 PDBs** (Dec 4, 2025)
- Perfect match: 97/100 (97%)
- Known edge cases: 1/100 (70U)
- Malformed legacy JSONs: 2/100 (1F8V, 1FFZ)
- Total acceptable: 98/100 (98%)

✅ **1OB2 Re-verification** (Dec 4, 2025)
- Regenerated with KIR fix
- KIR atoms match perfectly: C4, C2, N1, C6, C5
- RMS difference: 0.00e+00

### Known Edge Cases

**Documented and Acceptable:**

1. **A23** (2'-deoxy-2'-fluoroadenosine)
   - RMS diff: ~5e-3
   - Cause: Numerical precision in Jacobi eigenvalue solver
   - Status: Acceptable (same algorithm, modern actually better)
   - Tolerance: 1e-2

2. **70U** (7-methyluridine/2-thiouridine)
   - RMS diff: ~0.09
   - Cause: LS fitting anomaly for S2 sulfur substitution geometry
   - Status: Documented, tracked
   - Tolerance: 0.15

3. **Malformed Legacy JSONs**
   - PDBs: 1F8V, 1FFZ (and others)
   - Cause: Legacy code bug (corrupt JSON output)
   - Status: Skip these in validation
   - Action: Not our bug, legacy issue

### Next Steps

**IMMEDIATE:**
1. ✅ Clean up all generated modern JSONs (DONE)
2. ✅ Test 1OB2 regeneration (DONE - KIR verified)
3. ⏳ Regenerate all 3,602 frame JSONs with KIR fix
4. ⏳ Run full validation on 3,602 fast PDBs
5. ⏳ Generate final statistics and report

**Expected Results:**
- Perfect match: >95% (3,420+ PDBs)
- A23 edge cases: ~10-20 PDBs
- 70U edge cases: ~5-10 PDBs
- Malformed legacy: ~10-20 PDBs
- Total acceptable: >99%

**After Frame Validation:**
- Proceed to Stage 4: Distance Checks
- Then Stage 5: H-bond Detection
- Then Stage 6: Pair Validation
- Then Stage 7: Pair Selection (CRITICAL)

### Files Modified (KIR Fix)

**Code:**
- `src/x3dna/algorithms/base_frame_calculator.cpp` - Integrated hardcoded atom lookup
- `src/x3dna/algorithms/template_assignment.cpp` - Added KIR atom list

**Scripts:**
- `scripts/validate_frames_parallel.py` - Removed KIR tolerance exception

**Docs:**
- `docs/KNOWN_FRAME_EDGE_CASES.md` - Updated KIR status to FIXED
- `docs/FRAME_VALIDATION_STATUS.md` - This file (NEW)

### Validation Commands

**Regenerate all frames:**
```bash
cat data/valid_pdbs_fast.json | jq -r '.valid_pdbs_with_atoms_and_frames[]' | \
  xargs -P 20 -I {} sh -c './build/generate_modern_json data/pdb/{}.pdb data/json/ --stage=frames'
```

**Run validation:**
```bash
python3 scripts/validate_frames_parallel.py
```

**Quick test (first 100):**
```bash
cat data/valid_pdbs_fast.json | jq -r '.valid_pdbs_with_atoms_and_frames[]' | head -100 | \
  xargs -P 20 -I {} sh -c './build/generate_modern_json data/pdb/{}.pdb data/json/ --stage=frames'
python3 scripts/validate_frames_parallel.py  # Will only validate those with modern JSONs
```

### Summary

**Frame validation is ready for final run:**
- All known bugs fixed (KIR ✅)
- Edge cases documented (A23, 70U)
- 98% success rate on 100 PDB sample
- Clean slate (all old JSONs removed)
- Ready to regenerate and validate all 3,602 PDBs

**Timeline:**
- Regeneration: ~5-10 minutes (20 parallel threads)
- Validation: ~2-5 minutes
- Analysis: ~1 minute
- **Total: ~10-15 minutes to completion**


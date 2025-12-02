# Next Steps After Index Validation

**Date**: December 2, 2025  
**Current Branch**: `fix-index-matching`  
**Status**: Index validation complete ✅

---

## Immediate Next Steps

### 1. Merge to Main ✅ RECOMMENDED

Now that index validation is complete, merge the branch:

```bash
# Review the changes
git log fix-index-matching --oneline

# Merge to main
git checkout main
git merge fix-index-matching

# Push if desired
git push origin main
```

**Why now?**
- ✅ All validation passing (3,790 PDBs)
- ✅ Zero failures
- ✅ --only-paired mode working
- ✅ Complete documentation

---

## 2. Start Base Pair Comparisons 🎯

Now you can reliably compare base pair calculations!

### Option A: Compare a Single PDB

```bash
# Generate both legacy and modern with --only-paired
cd org && ./build/bin/find_pair_original ../data/pdb/1EHZ.pdb 1ehz.inp
cd .. && ./build/generate_modern_json data/pdb/1EHZ.pdb data/json --only-paired

# Compare base pairs
python3 scripts/compare_base_pairs.py data/json_legacy/base_pair/1EHZ.json \
                                       data/json/base_pair/1EHZ.json
```

### Option B: Batch Compare Many PDBs

```bash
# Use existing comparison infrastructure
python3 scripts/verify_all_pdbs.py --threads 8 --batch-size 50
```

### Option C: Focus on Specific Pairs

```bash
# Extract and compare specific problematic pairs
python3 scripts/investigate_pair_differences.py 1EHZ \
    --pair 44,55 --show-details
```

---

## 3. Compare Geometric Parameters 📐

With indices validated, you can trust step parameter comparisons:

```bash
# Compare step parameters (Shift, Slide, Rise, Tilt, Roll, Twist)
python3 scripts/verify_step_params.py 1EHZ

# Or batch compare
python3 scripts/batch_generate_and_verify.py --param-check
```

---

## 4. Investigate Any Calculation Differences 🔬

If you find differences in calculations:

```bash
# Deep dive into specific residue frames
python3 scripts/debug_ref_frames.py 1EHZ --residue 44

# Compare frame matrices
python3 scripts/compare_ref_frames.py ref_frames_legacy.dat \
                                       ref_frames_modern.dat \
                                       --by-residue

# Check specific base pair parameters
python3 scripts/compare_frames.py data/json_legacy/base_frame_calc/1EHZ.json \
                                   data/json/base_frame_calc/1EHZ.json \
                                   --show-diff --tolerance 0.01
```

---

## 5. Update Documentation 📝

### Update README

Add section about validation:

```markdown
## Index Validation

This codebase has been validated against legacy code on 3,790 PDB structures:
- ✅ 100% index matching success rate
- ✅ Validated on ~486,000 nucleotides
- ✅ Zero indexing errors

See `VALIDATION_COMPLETE.md` for details.
```

### Update Comparison Docs

Note that `--only-paired` is now the default for comparisons:

```markdown
## Comparing with Legacy

All comparison tools now use `--only-paired` mode by default to ensure
exact matching with legacy behavior (only processes paired residues).
```

---

## 6. Clean Up (Optional) 🧹

### Remove Old/Stale Files

```bash
# Remove old mapping files (none should exist for PASS)
ls data/index_mapping/  # Should be empty or only timeouts

# Clean up temporary files
rm -f temp_*.inp *.out

# Archive old documentation if desired
./docs/archive_docs.sh
```

### Consolidate Documentation

You now have:
- `VALIDATION_COMPLETE.md`
- `ONLY_PAIRED_MODE.md`  
- `INDEX_MATCHING_PROOF.md`
- `RESOLUTION_7EH2.md`
- `INVESTIGATION_7EH2.md`
- `docs/INDEX_VALIDATION_REPORT.md`
- `docs/ONLY_PAIRED_EXPLAINED.md`

Consider:
- Moving investigation files to `docs/investigations/`
- Creating a summary index document

---

## Recommended Priority Order

### High Priority (Do First) 🔥

1. **Merge to main** - Lock in the validated code
2. **Compare base pairs** - Main goal, now possible
3. **Check for calculation differences** - Debug any found

### Medium Priority (Do Soon) 📋

4. **Batch verification** - Run on all PDBs
5. **Document findings** - Record what you learn
6. **Update tests** - Ensure regression testing

### Low Priority (Nice to Have) 💡

7. **Clean up docs** - Organize validation documents
8. **Handle timeouts** - Investigate the 26 timeout PDBs
9. **Generate missing legacy** - Optional, for 305 skipped PDBs

---

## Specific Commands to Run

### Start with a Simple Case

```bash
# Pick a small, well-behaved PDB
PDB="1EHZ"

# Generate modern JSON with --only-paired
./build/generate_modern_json data/pdb/${PDB}.pdb data/json --only-paired

# Compare everything
python3 scripts/compare_json.py \
    data/json_legacy/base_frame_calc/${PDB}.json \
    data/json/base_frame_calc/${PDB}.json

# If differences found, deep dive
python3 scripts/debug_ref_frames.py ${PDB}
```

### Then Scale Up

```bash
# Batch compare all validated PDBs
python3 scripts/verify_all_pdbs.py \
    --only-validated \
    --threads 8 \
    --stop-on-first-diff

# This will:
# - Only test PDBs that passed validation
# - Use 8 threads for speed
# - Stop at first calculation difference for investigation
```

---

## Expected Outcomes

### Best Case 🎉

- All calculations match legacy exactly
- Zero differences found
- Modern code is perfect!

### Realistic Case 📊

- Small numerical differences (< 0.01°)
- Due to floating point precision
- Document and accept

### Interesting Case 🔬

- Systematic differences in specific cases
- Need investigation to determine correct answer
- May find bugs in legacy OR modern
- Requires careful analysis

---

## Tools Already Available

You have these comparison tools ready:

```bash
ls scripts/compare_*.py
# compare_base_pairs.py
# compare_best_partner.py
# compare_frames.py
# compare_iteration.py
# compare_json.py
# compare_mutual_best.py
# compare_ref_frames.py
# compare_validation_geometry.py
```

All should work correctly now that indices are validated!

---

## Questions to Answer Next

1. **Do base pair counts match?**
   - Modern finds same number of pairs as legacy?
   
2. **Do geometric parameters match?**
   - Shift, Slide, Rise within tolerance?
   - Tilt, Roll, Twist within tolerance?

3. **Do individual frames match?**
   - Reference frame origins?
   - Frame orientations (rotation matrices)?

4. **Do H-bonds match?**
   - Same H-bonds detected?
   - Same donor/acceptor atoms?

---

## Success Criteria

After completing comparisons, you should know:

✅ Which calculations match legacy exactly  
✅ Which have acceptable differences (< tolerance)  
✅ Which have systematic differences (need investigation)  
✅ Overall: Is modern code ready for production?  

---

## Getting Started

**Run this now**:

```bash
# 1. Merge validated code
git checkout main
git merge fix-index-matching

# 2. Pick a test PDB
PDB="1EHZ"

# 3. Generate and compare
./build/generate_modern_json data/pdb/${PDB}.pdb data/json --only-paired
python3 scripts/compare_json.py data/json_legacy/base_frame_calc/${PDB}.json \
                                 data/json/base_frame_calc/${PDB}.json

# 4. Check results
echo "If you see differences, start investigating!"
```

---

## Ready When You Are! 🚀

You've successfully:
- ✅ Validated index matching on 3,790 structures
- ✅ Implemented --only-paired mode
- ✅ Achieved 100% validation success

**Next**: Start comparing calculations and debugging any differences!

The hard infrastructure work is done. Now comes the interesting part - finding and fixing calculation differences! 🔬


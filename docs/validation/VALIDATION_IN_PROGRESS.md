# Comprehensive Index Validation: IN PROGRESS

**Started**: December 2, 2025  
**Status**: 🔄 Running...  
**Est. Completion**: ~1 hour

---

## Progress

### Current Status
- **Total PDBs**: 4,123 (from valid_pdbs.json)
- **Tested**: ~300 / 4,123 (~7%)
- **✅ PASS**: ~250
- **❌ FAIL**: 0  
- **⏭️  SKIP**: ~50 (no legacy JSON)

### What's Being Tested

Each PDB validation checks:
1. ✅ Parse PDB file (modern code)
2. ✅ Calculate frames (RMSD filtering applied)
3. ✅ Track nucleotides in ResidueTracker
4. ✅ Load legacy indices from JSON
5. ✅ Validate counts match
6. ✅ Validate each residue matches by PDB properties
7. ✅ Validate atom legacy indices match JSON indices
8. ✅ Export mapping to JSON

### Key Tests Validated So Far

**Filtering Edge Cases**:
- ✅ 1TTT D:16 - Filtered by both (RMSD check failed)
- ✅ Large structures (1FJG, 1HNZ: 1,513 nucleotides each)
- ✅ Small structures (1FXL: 8 nucleotides)
- ✅ Various modified nucleotides (H2U, PSU, 2MG, etc.)

**No Failures Yet**: 0 index mismatches in 250+ PDBs!

---

## Monitoring

**Auto-monitor running**: Checks every 60 seconds, stops on first failure

**Manual check**:
```bash
./check_validation_progress.sh
```

**See latest batch**:
```bash
tail -50 ~/.cursor/projects/.../terminals/4.txt | grep -E "(Batch|FAIL|summary)"
```

---

## If Failure Found

The script will:
1. **Stop immediately**
2. **Report** the failed PDB ID
3. **Save** details to CSV
4. **Create** mapping file for investigation

To investigate:
```bash
# Check the mapping
cat data/index_mapping/{FAILED_PDB}.json

# Check modern vs legacy
cat data/json/base_frame_calc/{FAILED_PDB}.json
cat data/json_legacy/base_frame_calc/{FAILED_PDB}.json

# Resume after fix
python3 scripts/validate_all_indices.py --start-from {FAILED_PDB}
```

---

## Expected Outcome

**Best case**: All 4,123 PDBs pass → **100% validation**  
**Likely case**: Some fails → Investigate and fix root causes  
**Worst case**: Many fails → Review filtering logic

Based on results so far (250+ pass, 0 fail), looking very promising! 🎉

---

## Results Will Show

1. **Which PDBs have matching indices** (can be compared safely)
2. **Which PDBs have mismatches** (need investigation)
3. **Which PDBs are skipped** (no legacy JSON)
4. **Filtering differences** (modern vs legacy)

This creates the **definitive list** of PDBs ready for comparison.

---

## Next Steps After Validation

### If All Pass ✅
→ Proceed to Week 2: Unified comparison framework

### If Some Fail ❌
→ Investigate root causes:
- Different filtering logic?
- Index assignment differences?
- Edge cases in parsing?

→ Fix and re-validate

---

**Estimated completion**: Check back in ~45 minutes or monitor with `check_validation_progress.sh`

